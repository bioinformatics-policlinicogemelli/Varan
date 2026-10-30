#Copyright 2025 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

"""Pre-flight ("dry run") validation for a Varan invocation - see -D/
--dry-run in varan.py.

Runs every check that can be done WITHOUT actually starting the pipeline
(no VEP/vcf2maf, no vendor download, no output folder created, nothing
written to disk) and prints ONLY the problems found, to stdout - no
report file. If nothing is wrong, prints a single "all good" line
instead. The caller (varan.py) turns the returned bool into an exit
code: 0 if there's no FAIL-level issue (WARN-level issues don't block a
real run, they're just worth knowing about upfront), 1 otherwise.

Every check here mirrors a failure this codebase has actually hit at
runtime, deep inside an otherwise-long pipeline run (see git history /
CONVERSATION for the concrete incidents): a BOM-corrupted dict.csv
silently emptying the OncoTree lookup, a blank comb_path producing a
bogus `cp .`, a fold-change of exactly 0.0 crashing math.log2(), a
resumed run appending onto a stale-schema data_cna_hg19.seg.fc.txt, an
OncoKB filter key missing from an old conf.ini, and so on. The point of
-D is to surface all of these BEFORE the expensive VEP/vcf2maf/OncoKB
steps run, not to re-implement the pipeline's own logic.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
from loguru import logger

from config_loader import get_config
from filter_clinvar import check_bool
from vendor_adapters import ADAPTERS
from vendor_adapters.common import SAMPLE_TSV_COLUMNS, load_oncotree_dict
from versioning import get_newest_version, get_version_list

FAIL = "FAIL"
WARN = "WARN"


@dataclass
class Issue:
    level: str  # FAIL or WARN
    message: str


Issues = list


def _fail(issues: Issues, message: str) -> None:
    issues.append(Issue(FAIL, message))


def _warn(issues: Issues, message: str) -> None:
    issues.append(Issue(WARN, message))


def _report(issues: Issues) -> bool:
    """Print only the problems found (or a single all-clear line), to
    stdout only - no file. Returns True if the run could proceed (no
    FAIL-level issue; WARN-level issues are informational).
    """
    if not issues:
        logger.info("Dry run: all checks passed, this run is ready to go.")
        return True

    fails = [i for i in issues if i.level == FAIL]
    warns = [i for i in issues if i.level == WARN]

    if fails:
        logger.critical(
            f"Dry run: {len(fails)} blocking problem(s) found - this run "
            "would fail (or silently produce incomplete output) if started "
            "as-is:")
        for issue in fails:
            logger.critical(f"  [FAIL] {issue.message}")

    if warns:
        logger.warning(
            f"Dry run: {len(warns)} non-blocking issue(s) worth checking:")
        for issue in warns:
            logger.warning(f"  [WARN] {issue.message}")

    return not fails


# --------------------------------------------------------------- conf.ini

def _check_conf_paths(issues: Issues, config, sigma: bool) -> None:
    vcf2maf = config.get("Paths", "VCF2MAF", fallback="").strip()
    if vcf2maf and not Path(vcf2maf).exists():
        _fail(issues, f"conf.ini [Paths] VCF2MAF not found: {vcf2maf}")

    vep_path = config.get("Paths", "VEP_PATH", fallback="").strip()
    if vep_path and not Path(vep_path).exists():
        _fail(issues, f"conf.ini [Paths] VEP_PATH not found: {vep_path}")

    vep_data = config.get("Paths", "VEP_DATA", fallback="").strip()
    if vep_data and not Path(vep_data).exists():
        _fail(issues, f"conf.ini [Paths] VEP_DATA not found: {vep_data}")

    ref_fasta = config.get("Paths", "REF_FASTA", fallback="").strip()
    if ref_fasta:
        if not Path(ref_fasta).exists():
            _fail(issues, f"conf.ini [Paths] REF_FASTA not found: {ref_fasta}")
        elif not Path(f"{ref_fasta}.fai").exists():
            _warn(
                issues,
                f"REF_FASTA has no .fai index next to it: {ref_fasta}.fai")

    clinv = config.get("Paths", "CLINV", fallback="").strip()
    if clinv:
        clinv_path = clinv.split(",")[0].strip()
        if clinv_path and not Path(clinv_path).exists():
            _fail(issues, f"conf.ini [Paths] CLINV not found: {clinv_path}")
        elif clinv_path and clinv_path.endswith((".vcf.gz", ".vcf.bgz")) and (
            not Path(f"{clinv_path}.tbi").exists()):
            _warn(
                issues,
                f"CLINV has no tabix index next to it: {clinv_path}.tbi")

    if sigma:
        rscript = config.get("Paths", "RSCRIPT", fallback="Rscript").strip()
        if not shutil.which(rscript) and not Path(rscript).exists():
            _fail(
                issues,
                f"-g/--sigma was given but RSCRIPT isn't runnable: "
                f"'{rscript}' not found on PATH or as a file")


def _check_oncokb_config(
    issues: Issues, config, oncokb: bool, filters: str,
    analysis_type: Optional[str]) -> None:
    if oncokb and config.get("OncoKB", "ONCOKB", fallback="").strip() == "":
        _fail(
            issues,
            "-k/--oncokb was given but conf.ini's [OncoKB] ONCOKB is empty")

    if not (oncokb and "o" in filters):
        return

    checks = [("Filters", "ONCOKB_FILTER_SNV")]
    if analysis_type not in ["snv", "fus", "tab"]:
        checks.append(("Cna", "ONCOKB_FILTER_CNV"))
    if analysis_type not in ["cnv", "snv", "tab"]:
        checks.append(("FUSION", "ONCOKB_FILTER_FUSION"))

    for section, key in checks:
        if not config.has_option(section, key):
            _fail(
                issues,
                f"-f o/-k was given but conf.ini is missing '{key}' under "
                f"[{section}] (an old conf.ini predating the SNV/CNV/FUSION "
                "OncoKB filter split would trigger this)")


def _check_header_config(issues: Issues, config) -> None:
    for key in ("HEADER_SAMPLE_SHORT", "HEADER_SAMPLE_LONG"):
        value = config.get("ClinicalSample", key, fallback="").strip()
        if not value:
            continue
        declared = len(value.split(","))
        minimum = len(SAMPLE_TSV_COLUMNS) - 3  # snv/cnv/comb_path get dropped
        if declared < minimum:
            _warn(
                issues,
                f"conf.ini [ClinicalSample] {key} declares {declared} "
                f"column name(s), fewer than the {minimum} base "
                "data_clinical_sample.txt columns - extra columns (RUN_ID, "
                "SigMA, vendor passthrough, ...) will likely make the real "
                "count mismatch and abort the run late, right when "
                "data_clinical_sample.txt is written.")


# ------------------------------------------------------------- versioning

def _check_output_folder(issues: Issues, output_folder: str,
                          overwrite_output: bool, resume: bool) -> None:
    if not output_folder:
        _fail(issues, "No output folder given (-o)")
        return

    parent = Path(output_folder).parent
    parent_to_check = parent if str(parent) else Path.cwd()
    if not os.access(parent_to_check, os.W_OK):
        _fail(
            issues,
            f"No write permission on '{parent_to_check}' - the output "
            "folder can't be created there")
        return

    existing_versions = get_version_list(output_folder)
    if not existing_versions:
        logger.info(
            f"Dry run: no existing version of '{output_folder}' - "
            f"'{output_folder}_v1' would be created.")
        return

    latest = Path(output_folder).parent / existing_versions[-1]
    if overwrite_output:
        logger.info(
            f"Dry run: -w was given - '{latest}' would be deleted and "
            "recreated from scratch.")
    elif resume:
        logger.info(f"Dry run: -R was given - '{latest}' would be resumed.")
        segfc = latest / "data_cna_hg19.seg.fc.txt"
        if segfc.exists():
            try:
                with segfc.open() as f:
                    header = f.readline().rstrip("\n").split("\t")
                if len(header) not in (10, 11):
                    _warn(
                        issues,
                        f"{segfc} has an unexpected column count "
                        f"({len(header)}) for a resumed run - appending "
                        "this run's rows onto it risks the exact "
                        "'Expected N fields, saw M' crash seen before; "
                        "consider -w instead of -R for this folder.")
            except OSError as err:
                _warn(issues, f"Could not read {segfc}: {err}")
    else:
        next_version, _ = get_newest_version(output_folder)
        _warn(
            issues,
            f"'{latest}' already exists, and neither -w nor -R was given - "
            f"'{next_version}' would be created fresh instead (existing "
            "data is left untouched). Pass -w to overwrite it or -R to "
            "resume it if that's not what you want.")


# -------------------------------------------------------- vendor adapter

def _check_vendor_input(issues: Issues, pipeline: str,
                         varan_input: Sequence[str]) -> None:
    config = get_config()
    section = f"Vendor.{pipeline}"
    dict_path = (
        dict(config.items(section)).get("dict_path")
        if config.has_section(section) else None
    ) or getattr(ADAPTERS[pipeline], "DICT_PATH", None)

    if dict_path:
        if not Path(dict_path).exists():
            _fail(issues, f"Vendor '{pipeline}' dict_path not found: {dict_path}")
        else:
            name_to_code, _ = load_oncotree_dict(dict_path)
            if not name_to_code:
                _fail(
                    issues,
                    f"Vendor '{pipeline}' dict_path exists but produced zero "
                    f"usable name/code rows: {dict_path} (check its header "
                    "is exactly 'name,code' or 'code,name', and that it "
                    "isn't BOM-prefixed - a common Excel 'CSV UTF-8' export "
                    "artifact)")

    raw_input = varan_input[0]
    if raw_input.startswith("s3://"):
        if not shutil.which("aws"):
            _fail(issues, "Vendor input is an S3 path but the 'aws' CLI isn't installed")
        else:
            try:
                result = subprocess.run(
                    ["aws", "s3", "ls", raw_input.rstrip("/") + "/"],
                    capture_output=True, text=True, timeout=15, check=False)
                if result.returncode != 0:
                    _fail(
                        issues,
                        f"Vendor input S3 path isn't reachable/listable: "
                        f"{raw_input} ({result.stderr.strip() or 'aws s3 ls failed'})")
            except subprocess.TimeoutExpired:
                _warn(
                    issues,
                    f"Timed out checking S3 reachability for {raw_input} - "
                    "could be a slow connection or misconfigured credentials")
    elif not Path(raw_input).exists():
        _fail(issues, f"Vendor input not found: {raw_input}")


# ---------------------------------------------------------------- sample.tsv

def _resolve_sample_tsv(varan_input: Sequence[str]) -> tuple[str, str, str]:
    """Mirror walk.input_extraction_file()/input_extraction_folder() without
    importing walk (which reads conf.ini's [Sample_Type] at import time -
    keeping preflight.py independent of that ordering).
    """
    raw = Path(varan_input[0])
    if raw.is_file():
        sample_tsv = str(raw)
        patient_tsv = varan_input[1].strip() if len(varan_input) > 1 else ""
        fusion_tsv = varan_input[2].strip() if len(varan_input) > 2 else ""
    else:
        sample_tsv = str(raw / "sample.tsv")
        patient_candidate = raw / "patient.tsv"
        patient_tsv = str(patient_candidate) if patient_candidate.exists() else ""
        fusion_candidate = raw / "FUSIONS" / "Fusions.tsv"
        fusion_tsv = str(fusion_candidate) if fusion_candidate.exists() else ""
    return sample_tsv, patient_tsv, fusion_tsv


def _check_sample_tsv(issues: Issues, sample_tsv: str,
                       cnvkit_algorithm: bool) -> Optional[pd.DataFrame]:
    if not Path(sample_tsv).exists():
        _fail(issues, f"sample.tsv not found: {sample_tsv}")
        return None

    try:
        df = pd.read_csv(sample_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception as err:
        _fail(issues, f"sample.tsv couldn't be parsed as TSV: {sample_tsv} ({err})")
        return None

    if df.empty:
        _fail(issues, f"sample.tsv has no rows: {sample_tsv}")
        return None

    required = [
        "SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path",
        "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR"]
    missing_cols = [c for c in required if c not in df.columns]
    if missing_cols:
        _fail(
            issues,
            f"sample.tsv is missing required column(s): {', '.join(missing_cols)}")
        return df

    dup_ids = df["SAMPLE_ID"][df["SAMPLE_ID"].duplicated()].unique().tolist()
    if dup_ids:
        _fail(issues, f"sample.tsv has duplicate SAMPLE_ID(s): {', '.join(dup_ids)}")

    empty_oncotree = df.loc[df["ONCOTREE_CODE"].str.strip() == "", "SAMPLE_ID"].tolist()
    if empty_oncotree:
        _warn(
            issues,
            f"Sample(s) with no ONCOTREE_CODE (a default cancer type will "
            f"be used instead): {', '.join(empty_oncotree)}")

    missing_by_type: dict = {"snv_path": [], "cnv_path": [], "comb_path": []}
    for col in missing_by_type:
        for _, row in df.iterrows():
            value = row[col].strip()
            if value and not Path(value).exists():
                missing_by_type[col].append((row["SAMPLE_ID"], value))

    for col, entries in missing_by_type.items():
        for sample_id, value in entries:
            _fail(
                issues,
                f"Sample '{sample_id}': {col} doesn't exist on disk: {value}")

    if cnvkit_algorithm:
        if "TC" not in df.columns:
            _warn(
                issues,
                "conf.ini's CNVKIT_algorithm is True but sample.tsv has no "
                "TC column - only unadjusted copy number will be computed "
                "for every sample")
        else:
            missing_tc = df.loc[df["TC"].str.strip() == "", "SAMPLE_ID"].tolist()
            if len(missing_tc) == len(df):
                _warn(
                    issues,
                    "conf.ini's CNVKIT_algorithm is True but NO sample has "
                    "a TC value - only unadjusted copy number will be "
                    "computed for all of them")
            elif missing_tc:
                _warn(
                    issues,
                    f"Sample(s) with no TC value (only unadjusted copy "
                    f"number will be computed for them): {', '.join(missing_tc)}")

    return df


def _check_patient_tsv(issues: Issues, patient_tsv: str,
                        sample_df: Optional[pd.DataFrame]) -> None:
    if not patient_tsv:
        return
    if not Path(patient_tsv).exists():
        _fail(issues, f"patient.tsv not found: {patient_tsv}")
        return
    try:
        patient_df = pd.read_csv(patient_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception as err:
        _fail(issues, f"patient.tsv couldn't be parsed as TSV: {patient_tsv} ({err})")
        return
    if "PATIENT_ID" not in patient_df.columns:
        _fail(issues, f"patient.tsv has no PATIENT_ID column: {patient_tsv}")
        return
    if sample_df is None or "PATIENT_ID" not in sample_df.columns:
        return
    missing = sorted(
        set(sample_df["PATIENT_ID"]) - set(patient_df["PATIENT_ID"]))
    if missing:
        _warn(
            issues,
            f"PATIENT_ID(s) referenced in sample.tsv but absent from "
            f"patient.tsv: {', '.join(missing)}")


def _check_fusion_tsv(issues: Issues, fusion_tsv: str,
                       sample_df: Optional[pd.DataFrame]) -> None:
    if not fusion_tsv:
        return
    if not Path(fusion_tsv).exists():
        _fail(issues, f"fusion tsv not found: {fusion_tsv}")
        return
    try:
        fusion_df = pd.read_csv(fusion_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception as err:
        _fail(issues, f"fusion tsv couldn't be parsed as TSV: {fusion_tsv} ({err})")
        return
    id_col = next(
        (c for c in ("Sample_Id", "SAMPLE_ID", "sample_id") if c in fusion_df.columns),
        None)
    if id_col is None:
        _fail(issues, f"fusion tsv has no sample-id column: {fusion_tsv}")
        return
    if sample_df is None or "SAMPLE_ID" not in sample_df.columns:
        return
    unknown = sorted(set(fusion_df[id_col]) - set(sample_df["SAMPLE_ID"]))
    if unknown:
        _warn(
            issues,
            f"fusion tsv references sample(s) not present in sample.tsv: "
            f"{', '.join(unknown)}")


# -------------------------------------------------------------- raw VCFs

def _check_vcf_header(path: Path) -> Optional[str]:
    """Return the `##fileformat` version string, or None if the file
    couldn't be read/has no recognizable header - callers decide severity.
    """
    try:
        with path.open() as f:
            for line in f:
                if line.startswith("##fileformat"):
                    return line.split("=")[-1].strip()
                if not line.startswith("#"):
                    break
    except OSError:
        return None
    return None


def _check_cnv_vcf_rows(issues: Issues, path: Path, sample_id: str) -> None:
    version = _check_vcf_header(path)
    if version is None:
        _fail(issues, f"Sample '{sample_id}': CNV VCF has no ##fileformat line: {path}")
        return
    if version not in ("VCFv4.1", "VCFv4.2"):
        _fail(
            issues,
            f"Sample '{sample_id}': CNV VCF version '{version}' isn't "
            f"supported (only VCFv4.1/VCFv4.2 are): {path}")
        return

    bad_rows = 0
    try:
        with path.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                info = cols[7]
                is_struct = any(
                    tag in info for tag in ("SVTYPE=CNV", "SVTYPE=DUP", "SVTYPE=DEL"))
                if not is_struct:
                    continue
                if "END=" not in info:
                    bad_rows += 1
                    continue
                if version == "VCFv4.2" and len(cols) >= 9 and "SM" not in cols[8].split(":"):
                    bad_rows += 1
    except OSError as err:
        _warn(issues, f"Sample '{sample_id}': couldn't fully read CNV VCF {path}: {err}")
        return

    if bad_rows:
        _fail(
            issues,
            f"Sample '{sample_id}': CNV VCF has {bad_rows} structural "
            f"row(s) missing END/SM - these currently crash CNV "
            f"conversion (StopIteration/ValueError) rather than degrading "
            f"gracefully: {path}")


def _check_referenced_vcfs(issues: Issues, sample_df: pd.DataFrame,
                            analysis_type: Optional[str]) -> None:
    if analysis_type in ("snv", "fus", "tab"):
        return
    for _, row in sample_df.iterrows():
        cnv_path = row.get("cnv_path", "").strip()
        if not cnv_path or not Path(cnv_path).exists():
            continue  # already reported by _check_sample_tsv
        _check_cnv_vcf_rows(issues, Path(cnv_path), row["SAMPLE_ID"])


# ------------------------------------------------- update/extract/remove

def _check_study_folder(issues: Issues, path: str, label: str) -> bool:
    if not path:
        _fail(issues, f"-p/--Path is required for {label}")
        return False
    if not Path(path).is_dir():
        _fail(issues, f"{label}: '{path}' is not a valid folder")
        return False
    clin_sample = Path(path) / "data_clinical_sample.txt"
    if not clin_sample.exists():
        _fail(
            issues,
            f"{label}: '{path}' doesn't look like a Varan output folder - "
            "data_clinical_sample.txt is missing")
        return False
    return True


def _check_sample_list_file(issues: Issues, sample_list: str, path: str,
                             label: str) -> None:
    if not sample_list:
        _fail(issues, f"-s/--SampleList is required for {label}")
        return
    if not Path(sample_list).exists():
        _fail(issues, f"{label}: SampleList not found: {sample_list}")
        return
    try:
        with Path(sample_list).open() as f:
            first_line = f.readline()
            if len(first_line.split("\t")) > 1:
                _fail(
                    issues,
                    f"{label}: SampleList '{sample_list}' has more than one "
                    "column - it should be one sample ID per line")
                return
        with Path(sample_list).open() as f:
            requested = {s.strip() for s in f if s.strip()}
    except OSError as err:
        _fail(issues, f"{label}: couldn't read SampleList '{sample_list}': {err}")
        return

    if not requested:
        _fail(issues, f"{label}: SampleList '{sample_list}' is empty")
        return

    try:
        clin_df = pd.read_csv(
            Path(path) / "data_clinical_sample.txt", sep="\t", skiprows=4,
            dtype=str)
    except Exception as err:
        _warn(
            issues,
            f"{label}: couldn't read data_clinical_sample.txt to "
            f"cross-check SampleList against: {err}")
        return

    known = set(clin_df.get("SAMPLE_ID", []))
    if not requested & known:
        _fail(
            issues,
            f"{label}: none of the sample(s) in '{sample_list}' are present "
            f"in '{path}'s data_clinical_sample.txt")
        return

    missing = sorted(requested - known)
    if missing:
        _warn(
            issues,
            f"{label}: sample(s) in '{sample_list}' not found in "
            f"'{path}': {', '.join(missing)}")


def _run_update_extract_remove_checks(
    issues: Issues, update: bool, extract: bool, remove: bool,
    path: str, new_path: str, sample_list: str, output_folder: str,
    overwrite_output: bool) -> None:
    label = "Update" if update else ("Extract" if extract else "Remove")
    ok = _check_study_folder(issues, path, label)

    if update:
        if not new_path:
            _fail(issues, "-n/--NewPath is required for Update")
        elif not Path(new_path).is_dir():
            _fail(issues, f"Update: '{new_path}' is not a valid folder")
        elif not (Path(new_path) / "data_clinical_sample.txt").exists():
            _warn(
                issues,
                f"Update: '{new_path}' doesn't look like a Varan output "
                "folder - data_clinical_sample.txt is missing")
    else:
        if ok:
            _check_sample_list_file(issues, sample_list, path, label)

    if output_folder:
        _check_output_folder(issues, output_folder, overwrite_output, False)


# -------------------------------------------------------------- entrypoint

def run_preflight_checks(
    *, pipeline: str, varan_input: Optional[Sequence[str]], cancer: Optional[str],
    output_folder: str, oncokb: bool, filters: str,
    analysis_type: Optional[str], overwrite_output: bool, resume: bool,
    sigma: bool, update: bool, extract: bool, remove: bool,
    path: Optional[str], new_path: Optional[str],
    sample_list: Optional[str]) -> bool:
    """Run every pre-flight check applicable to this invocation and print
    the results (problems only, or a single all-clear line) to stdout.

    Returns True if the run could proceed (no FAIL-level issue found).
    """
    issues: Issues = []
    config = get_config()

    if update or extract or remove:
        _run_update_extract_remove_checks(
            issues, update, extract, remove, path, new_path, sample_list,
            output_folder, overwrite_output)
        return _report(issues)

    _check_conf_paths(issues, config, sigma)
    _check_oncokb_config(issues, config, oncokb, filters, analysis_type)
    _check_header_config(issues, config)
    _check_output_folder(issues, output_folder, overwrite_output, resume)

    if not cancer:
        _fail(issues, "-c/--Cancer is required")

    if not varan_input or not varan_input[0].strip():
        _fail(issues, "-i/--varan_input is required")
        return _report(issues)

    if pipeline in ADAPTERS:
        _check_vendor_input(issues, pipeline, varan_input)
        # sample.tsv doesn't exist yet for a vendor pipeline - it's
        # generated by the adapter itself, so the checks below don't apply.
        return _report(issues)

    sample_tsv, patient_tsv, fusion_tsv = _resolve_sample_tsv(varan_input)
    cnvkit_algorithm = check_bool(config.get("Cna", "CNVKIT_algorithm", fallback=""))
    sample_df = _check_sample_tsv(issues, sample_tsv, cnvkit_algorithm)
    _check_patient_tsv(issues, patient_tsv, sample_df)
    _check_fusion_tsv(issues, fusion_tsv, sample_df)
    if sample_df is not None:
        _check_referenced_vcfs(issues, sample_df, analysis_type)

    return _report(issues)
