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

"""Per-sample SigMA (Signature Multivariate Analysis) orchestration.

This is the one function walk.py's per-sample SNV loop calls when the
-g/--sigma CLI flag is set (see varan.py's ANNOTATION BLOCK). It ties
together the pieces already built on this branch:

- sigma_filter.prepare_sigma_maf(): filters one sample's unfiltered,
  vcf2maf-produced MAF down to SigMA-ready rows.
- sigma_cancer_type_map.get_sigma_call_params(): resolves ONCOTREE_CODE ->
  (SigMA tumor_type, do_mva).
- run_sigma.R: the actual R subprocess that loads the SigMA package,
  builds the 96-dim trinucleotide spectrum, and calls SigMA::run().

Every call to run_sigma_for_sample() returns a dict with at least
SAMPLE_ID and SIGMA_STATUS - callers should merge this row into the
SigMA results table regardless of outcome, so a run's skip/failure
reason (below SNV cutoff vs. an actual R-side error vs. an unmapped
OncoTree code) stays visible in data_clinical_sample.txt rather than the
sample silently having no SigMA columns at all. This mirrors the
project's "log and skip degrades gracefully, distinguish 'no data' from
'negative'" convention (see SIGMA_INTEGRATION_FEASIBILITY.md section 4,
point 3 on the ~25% of real panel samples that fall below snv_cutoff by
construction, not due to a bug).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
from loguru import logger

from config_loader import get_config
from filter_clinvar import check_bool
from sigma_cancer_type_map import get_sigma_call_params, is_mva_safe_for_panel_data
from sigma_filter import prepare_sigma_maf

config = get_config()

# Resolved relative to this module's own location, not the process's cwd -
# avoids the exact "cwd-relative external script path" bug class already
# found and fixed elsewhere in this codebase for cbio_validation's
# validateData.py invocation (see commit history on stabilize).
RUN_SIGMA_SCRIPT = Path(__file__).resolve().parent / "run_sigma.R"

# Varan's pipeline is hg19/GRCh37 end to end today (VEP_DATA cache, CNA
# file naming data_cna_hg19.*, etc. - see SIGMA_INTEGRATION_FEASIBILITY.md
# section 3). Not exposed as a conf.ini key yet since nothing else in
# Varan is build-agnostic either; trivial to promote to a [SigMA]
# REF_GENOME_BUILD key if/when Varan itself grows hg38 support.
REF_GENOME_BUILD = "hg19"

# SigMA's own run() requires a `catalog_name` argument (no usable
# default - confirmed from R/run.R source: it stop()s if NULL or unknown).
# Its "Parameter choices" wiki page recommends "cosmic_v2_inhouse" (COSMIC
# v2 plus WGS-discovered signatures not in the official catalog) for Sig3/
# MVA prediction specifically - the exact use case here. conf.ini's
# COSMIC_VERSION is deliberately kept in the short "v2"/"v3p2" form for a
# future GUI dropdown; this maps it to SigMA's actual catalog_name string.
_CATALOG_NAME_ALIASES = {
    "v2": "cosmic_v2_inhouse",
    "v3": "cosmic_v3p2_inhouse",
    "v3p2": "cosmic_v3p2_inhouse",
}


def resolve_catalog_name(cosmic_version: str) -> str:
    """Map conf.ini's [SigMA] COSMIC_VERSION to SigMA's `catalog_name`.

    Args:
        cosmic_version (str): conf.ini's [SigMA] COSMIC_VERSION value.

    Returns:
        str: A valid SigMA `catalog_name` value.

    """
    value = cosmic_version.strip()
    if value.lower().startswith("cosmic_"):
        # Power-user escape hatch: already a real SigMA catalog_name
        # (e.g. "cosmic_v2", "cosmic_v3p2" without "_inhouse") - pass
        # through unchanged instead of forcing the "_inhouse" variant.
        return value
    mapped = _CATALOG_NAME_ALIASES.get(value.lower())
    if mapped is None:
        logger.warning(
            f"[SigMA] COSMIC_VERSION='{cosmic_version}' is not one of the "
            f"recognized short forms ({sorted(_CATALOG_NAME_ALIASES)}) and "
            "doesn't start with 'cosmic_' - passing it to SigMA as-is; "
            "SigMA itself will raise an error if this isn't a real "
            "catalog_name.")
        return value
    return mapped


def blank_result(sample_id: str, status: str, **extra: object) -> dict:
    """Build a SigMA result row that still records SAMPLE_ID/status even
    when SigMA itself never actually ran for this sample.
    """
    row = {
        "SAMPLE_ID": sample_id,
        "SIGMA_STATUS": status,
        "SIGMA_TUMOR_TYPE": "",
        "SIGMA_DO_MVA": "",
        "SIGMA_TOTAL_SNVS": "",
        "SIGMA_SIGNATURE3_MVA": "",
        "SIGMA_SIGNATURE3_CALL": "",
        "SIGMA_CATEG": "",
    }
    row.update(extra)
    return row


def _classify_signature3(
    row: pd.Series,
    do_mva: bool,
    threshold_conf: str) -> str:
    """Resolve Varan's own SIGMA_SIGNATURE3_CALL from SigMA's raw output row.

    "default" (recommended - see conf.ini comment) uses SigMA's own
    pass_mva/pass_mva_strict cutoffs, already computed by SigMA::run()
    with do_assign=True. A numeric SIGNATURE3_POSITIVE_THRESHOLD instead
    thresholds Signature_3_mva directly, per SigMA's own "Parameter
    choices" wiki guidance that this is an explicitly-supported thing to
    tune, not a hack.
    """
    if not do_mva or "Signature_3_mva" not in row:
        return ""

    threshold_conf = threshold_conf.strip()
    if threshold_conf.lower() == "default":
        pass_strict = bool(row.get("pass_mva_strict", False))
        pass_loose = bool(row.get("pass_mva", False))
        if pass_strict:
            return "Positive_strict"
        if pass_loose:
            return "Positive"
        return "Negative"

    try:
        custom_threshold = float(threshold_conf)
    except ValueError:
        logger.warning(
            f"[SigMA] SIGNATURE3_POSITIVE_THRESHOLD='{threshold_conf}' is "
            "neither 'default' nor a number - falling back to SigMA's own "
            "pass_mva cutoff for this sample.")
        return "Positive" if row.get("pass_mva", False) else "Negative"

    mva_score = row.get("Signature_3_mva")
    if mva_score is None or pd.isna(mva_score):
        return ""
    return "Positive" if float(mva_score) >= custom_threshold else "Negative"


def run_sigma_for_sample(
    maf_df: pd.DataFrame,
    sample_id: str,
    oncotree_code: str,
    sigma_dir: Path) -> dict:
    """Run SigMA for one sample's freshly vcf2maf'd, unfiltered MAF.

    Args:
        maf_df (pd.DataFrame): The sample's unfiltered, annotated MAF
            (same DataFrame shape filter_main() reads from the `maf/`
            folder - i.e. read straight from vcf2maf's own output, before
            any clinical filtering).
        sample_id (str): The sample's SAMPLE_ID.
        oncotree_code (str): The sample's ONCOTREE_CODE (from sample.tsv).
        sigma_dir (Path): Intermediate working directory for this run's
            SigMA inputs/outputs (created by the caller), e.g.
            `<output_folder>/intermediate/sigma`.

    Returns:
        dict: A row of SIGMA_* columns plus SAMPLE_ID, always including
            SIGMA_STATUS ("OK", "BELOW_SNV_CUTOFF", "UNMAPPED_TUMOR_TYPE",
            or "ERROR"). Never raises - any failure is logged and
            reflected in SIGMA_STATUS instead, so one sample's SigMA
            problem never takes down the batch.

    """
    tumor_type_override = config.get("SigMA", "TUMOR_TYPE_OVERRIDE").strip()
    fallback_to_other = check_bool(config.get("SigMA", "FALLBACK_TO_OTHER"))

    if tumor_type_override:
        tumor_type = tumor_type_override
        do_mva = is_mva_safe_for_panel_data(tumor_type)
        if not do_mva:
            logger.warning(
                f"Sample {sample_id}: [SigMA] TUMOR_TYPE_OVERRIDE="
                f"'{tumor_type_override}' has no trained panel-data MVA "
                "classifier - forcing do_mva=False for this sample.")
    else:
        tumor_type, do_mva = get_sigma_call_params(
            oncotree_code, fallback_to_other)

    if tumor_type is None:
        logger.warning(
            f"Sample {sample_id}: ONCOTREE_CODE '{oncotree_code}' has no "
            "SigMA tumor_type mapping and [SigMA] FALLBACK_TO_OTHER is "
            "False - skipping SigMA for this sample.")
        return blank_result(sample_id, "UNMAPPED_TUMOR_TYPE")

    try:
        filtered_maf = prepare_sigma_maf(maf_df, sample_id)
    except Exception as err:  # noqa: BLE001 - never let one sample crash the batch
        logger.warning(
            f"Sample {sample_id}: prepare_sigma_maf() failed ({err!r}) - "
            "skipping SigMA for this sample.")
        return blank_result(sample_id, "ERROR", SIGMA_TUMOR_TYPE=tumor_type)

    snv_cutoff = int(config.get("SigMA", "SNV_CUTOFF"))
    n_variants = len(filtered_maf)
    if n_variants < snv_cutoff:
        logger.info(
            f"Sample {sample_id}: only {n_variants} SNV(s) survived SigMA's "
            f"own filtering (below SNV_CUTOFF={snv_cutoff}) - this is the "
            "expected, common case on panel-scale data (see "
            "SIGMA_INTEGRATION_FEASIBILITY.md), not an error. Skipping "
            "SigMA for this sample.")
        return blank_result(
            sample_id, "BELOW_SNV_CUTOFF",
            SIGMA_TUMOR_TYPE=tumor_type, SIGMA_DO_MVA=do_mva,
            SIGMA_TOTAL_SNVS=n_variants)

    maf_input_dir = sigma_dir / "maf_input"
    results_dir = sigma_dir / "results"
    maf_input_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    maf_input_path = maf_input_dir / f"{sample_id}.sigma_input.maf"
    filtered_maf.to_csv(maf_input_path, sep="\t", index=False)

    output_csv = results_dir / f"{sample_id}.sigma_output.csv"

    rscript_bin = config.get("Paths", "RSCRIPT") or "Rscript"
    catalog_name = resolve_catalog_name(config.get("SigMA", "COSMIC_VERSION"))
    data_platform = config.get("SigMA", "DATA_PLATFORM")
    check_msi = check_bool(config.get("SigMA", "CHECK_MSI"))
    lite_format = check_bool(config.get("SigMA", "LITE_FORMAT"))

    cmd = [
        rscript_bin, str(RUN_SIGMA_SCRIPT),
        "--maf", str(maf_input_path),
        "--sample-id", sample_id,
        "--tumor-type", tumor_type,
        "--data-platform", data_platform,
        "--do-mva", str(do_mva),
        "--check-msi", str(check_msi),
        "--catalog-name", catalog_name,
        "--lite-format", str(lite_format),
        "--snv-cutoff", str(snv_cutoff),
        "--ref-genome", REF_GENOME_BUILD,
        "--output", str(output_csv),
    ]

    logger.info(f"Sample {sample_id}: running SigMA (tumor_type={tumor_type}, "
                f"do_mva={do_mva}, data={data_platform})...")

    try:
        # Never shell=True - same fix already applied to
        # vcf2maf_constructor's vcf-query call elsewhere in this codebase.
        proc = subprocess.run(cmd, capture_output=True, check=False)
    except FileNotFoundError:
        logger.warning(
            f"Sample {sample_id}: could not run SigMA - '{rscript_bin}' was "
            "not found on PATH. Is R installed? (see [Paths] RSCRIPT in "
            "conf.ini). Skipping SigMA for this sample.")
        return blank_result(
            sample_id, "ERROR", SIGMA_TUMOR_TYPE=tumor_type,
            SIGMA_DO_MVA=do_mva, SIGMA_TOTAL_SNVS=n_variants)

    if proc.returncode != 0 or not output_csv.exists() or output_csv.stat().st_size == 0:
        stderr_text = proc.stderr.decode("utf-8", errors="replace").strip()
        logger.warning(
            f"Sample {sample_id}: SigMA run failed (exit code "
            f"{proc.returncode}). {stderr_text[-2000:] if stderr_text else ''}")
        return blank_result(
            sample_id, "ERROR", SIGMA_TUMOR_TYPE=tumor_type,
            SIGMA_DO_MVA=do_mva, SIGMA_TOTAL_SNVS=n_variants)

    try:
        result_df = pd.read_csv(output_csv)
    except Exception as err:  # noqa: BLE001
        logger.warning(
            f"Sample {sample_id}: could not parse SigMA output {output_csv} "
            f"({err!r}) - skipping SigMA for this sample.")
        return blank_result(
            sample_id, "ERROR", SIGMA_TUMOR_TYPE=tumor_type,
            SIGMA_DO_MVA=do_mva, SIGMA_TOTAL_SNVS=n_variants)

    if result_df.empty:
        logger.warning(
            f"Sample {sample_id}: SigMA output {output_csv} had no rows - "
            "skipping SigMA for this sample.")
        return blank_result(
            sample_id, "ERROR", SIGMA_TUMOR_TYPE=tumor_type,
            SIGMA_DO_MVA=do_mva, SIGMA_TOTAL_SNVS=n_variants)

    row = result_df.iloc[0]
    signature3_call = _classify_signature3(
        row, do_mva, config.get("SigMA", "SIGNATURE3_POSITIVE_THRESHOLD"))

    return {
        "SAMPLE_ID": sample_id,
        "SIGMA_STATUS": "OK",
        "SIGMA_TUMOR_TYPE": tumor_type,
        "SIGMA_DO_MVA": do_mva,
        "SIGMA_TOTAL_SNVS": row.get("total_snvs", n_variants),
        "SIGMA_SIGNATURE3_MVA": row.get("Signature_3_mva", ""),
        "SIGMA_SIGNATURE3_CALL": signature3_call,
        "SIGMA_CATEG": row.get("categ", ""),
    }
