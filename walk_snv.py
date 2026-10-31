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

"""Walk stage: SNV calls -> per-sample vcf2maf/VEP MAFs (with co-located
variant recovery) and, with -g/--sigma, SigMA.
"""

from __future__ import annotations

import os
import pickle
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from loguru import logger

import clonal_hematopoiesis
import sigma_runner
import specimen
import vcf_filter
from config_loader import config
from filter_clinvar import read_maf
from walk_common import (
    _path_setting,
    create_random_name_folder,
    output_filtered,
    WalkContext,
)


def get_snv_from_folder(inputfolder_snv: str) -> list[str]:
    """List all VCF files in a given folder.

    Args:
        inputfolder_snv (str): Path to the folder containing SNV VCF files.

    Returns:
        list[str]: List of filenames ending with 'vcf'.

    """
    return [str(file) for file in
            Path(inputfolder_snv).iterdir() if file.suffix == ".vcf"]


def get_sample_id_from_snv(snv_vcf: str) -> str:
    """Extract the sample ID from a VCF filename by converting it to a BAM filename.

    Args:
        snv_vcf (str): SNV VCF filename.

    Returns:
        str: Corresponding BAM filename.

    """
    if "MergedSmallVariants.genome.vcf" in snv_vcf:
        sample = snv_vcf.replace("_MergedSmallVariants.genome.vcf", ".bam")
    else:
        sample = "bam".join(snv_vcf.rsplit("vcf", 1))
    return sample


def snv_type_from_folder(input_pat: str,
                         snv_vcf_files: list,
                         output_folder: str,
                         id_map: dict[str, str] | None = None) -> dict:
    """Map sample IDs to their full SNV paths, handling duplicates.

    Args:
        input_pat (str): Path to the folder containing SNV files.
        snv_vcf_files (list): List of SNV VCF filenames.
        output_folder (str): The study's output folder, where sampleID_dup.log
            and noParsed_snv.log are written (read back into the report's
            Warnings section).
        id_map (dict, optional): VCF file name -> sample.tsv SAMPLE_ID (see
            walk.sample_ids_by_file()). A file not in it gets an id derived
            from its own name (get_sample_id_from_snv()).

    Returns:
        dict: Mapping from sample ID to full SNV file path.

    """
    c = 0
    sid_path = {}
    for case_folder in snv_vcf_files:
        try:
            snv_vcf = case_folder
            sample_id = (id_map or {}).get(Path(case_folder).name)                 or get_sample_id_from_snv(case_folder)
            if sample_id in sid_path:
                with (Path(output_folder) / "sampleID_dup.log").open("a", encoding="utf-8") as dup:
                    dup.write(sample_id + "\t" + "snv_vcf\n")
            else:
                sid_path[sample_id] = str(Path(input_pat) / snv_vcf)
        except Exception as err:
            logger.warning(f"Error while reading {case_folder}: {err}")
            with (Path(output_folder) / "noParsed_snv.log").open("a", encoding="utf-8") as log_noparsed:
                log_noparsed.write("[WARNING]" + case_folder + "\n")
        c = c + 1

    return sid_path


def vcf_filtering(sid_path: dict,
                  output_folder: str,
                  output_filtered: str) -> dict:
    """Apply filtering to VCF file.

    Args:
        sid_path (dict): Mapping from sample IDs to VCF file paths.
        output_folder (str): Path to the base output folder.
        output_filtered (str): Subfolder name for filtered output files.

    Returns:
        dict: Updated mapping from sample IDs
        to filtered VCF file paths.

    """
    sid_path_filtered = {}
    if output_filtered.strip() == "":
        output_filtered = "snv_filtered"
    output_path = Path(output_folder) / output_filtered
    output_path.mkdir(parents=True, exist_ok=True)

    for k, v in sid_path.items():
        _, vcf_file = os.path.split(v)
        out_filt = Path(output_folder) / output_filtered
        vcf_filtered = out_filt / (vcf_file.replace(".vcf","") + ".FILTERED.vcf")
        vcf_filter.main(v, str(vcf_filtered))
        sid_path_filtered[k] = str(vcf_filtered)
    return sid_path_filtered


def vcf2maf_constructor(v: str,
                        temporary: str,
                        output_folder: str,
                        sample_id: str | None = None) -> list:
    """Construct the vcf2maf command line based on input parameters and config file.

    Args:
        v (str): Path to the input VCF file.
        temporary (str): Temporary directory for intermediate files.
        output_folder (str): Base output folder.
        sample_id (str, optional): The sample.tsv SAMPLE_ID this VCF belongs
            to. When given, it becomes the MAF's Tumor_Sample_Barcode
            (--tumor-id) even if the VCF's own genotype column is named
            otherwise (e.g. a DRAGEN library id, then passed as
            --vcf-tumor-id) - so mutations line up with the clinical data.

    Returns:
        list: A list representing the command-line call to vcf2maf.

    """
    cache = config.get("Paths", "CACHE")
    cmd = ["vcf-query", "-l", v]
    try:
        tum_id = subprocess.check_output(cmd).decode("utf-8").strip()
    except Exception as err:
        logger.warning(
            f"Could not extract the tumor sample ID from {v} via vcf-query "
            f"({err}) - vcf2maf will be run with an empty --tumor-id.")
        tum_id = ""

    if _path_setting("VCF2MAF") == "" or _path_setting("REF_FASTA") == "" or _path_setting("VEP_PATH") == "" or _path_setting("VEP_DATA") == "":
        logger.critical("[Paths] section in conf.ini is not correctly compiled. "
        "Please check again!")
        msg = "Input error"
        raise Exception (msg)

    cl = ["perl"]
    cl.append(_path_setting("VCF2MAF"))
    cl.append("--input-vcf")
    cl.append(v)
    _, file_vcf = os.path.split(v)

    out_file = Path(output_folder) / "maf" / (file_vcf + ".maf")
    # No ClinVar annotation when [Paths] CLINV is blank - warned about once,
    # at setup (walk.validate_input), not for every VCF.
    if _path_setting("CLINV").strip() != "":
        cl.append("--vep-custom")
        cl.append(_path_setting("CLINV"))
    cl.append("--output-maf")
    cl.append(str(out_file))
    cl.append("--ref-fasta")
    cl.append(_path_setting("REF_FASTA"))
    cl.append("--tmp-dir")
    cl.append(str(temporary))
    cl.append("--retain-fmt")
    cl.append("GT,GQ,AD,DP,VF,AF")
    cl.append("--vep-path")
    cl.append(_path_setting("VEP_PATH"))
    cl.append("--vep-data")
    cl.append(_path_setting("VEP_DATA"))
    cl.append("--tumor-id")
    if sample_id and tum_id:
        cl.append(sample_id)
        if tum_id != sample_id:
            cl.append("--vcf-tumor-id")
            cl.append(tum_id)
    else:
        # Without the VCF's own column name vcf2maf couldn't find the
        # genotypes under SAMPLE_ID: keep its column name as the barcode.
        cl.append(tum_id)
    cl.append("--cache-version")
    cl.append(cache)
    # Same value the Snakefile declares as each SNV job's threads.
    cl.append("--vep-forks")
    cl.append(str(config.getint("Paths", "VEP_FORKS", fallback=4)))
    return cl


def run_vcf2maf(cl: list, sample: str) -> bool:
    """Run the vcf2maf conversion using the constructed command.

    Args:
        cl (list): The command-line list to execute.
        sample (str): The sample ID for logging purposes.

    Returns:
        bool: Whether vcf2maf produced a non-empty MAF.

    """
    logger.info(f"Starting vcf2maf conversion of {sample}... This may take several "
    "minutes (approx 5 minutes per sample).")
    logger.info(f"args={cl}")
    sout = subprocess.run(cl, capture_output=True, check=False)

    stderr = (sout.stderr or b"").decode("ascii", errors="replace").strip()
    if stderr:
        if "ERROR" not in stderr:
            logger.warning(stderr.replace("ERROR: ", ""))
        else:
            logger.error(stderr.replace("ERROR: ", ""))

    # stderr text alone is not a reliable failure signal: vcf2maf commonly
    # logs non-fatal INFO/warning text there even on success (hence only
    # escalating to ERROR above when "ERROR" literally appears in it) - so
    # a crash, OOM kill, or timeout that produces no matching stderr text
    # would otherwise go completely unnoticed. The exit code isn't checked
    # here on purpose (deliberately, not an oversight): vcf2maf can exit
    # non-zero on a run that still produced a perfectly good MAF, so it
    # isn't a trustworthy signal either - checking whether the MAF it was
    # told to produce actually exists and has content is the one check
    # that's actually reliable both ways.
    if "--output-maf" in cl:
        out_file = Path(cl[cl.index("--output-maf") + 1])
        if not out_file.exists() or out_file.stat().st_size == 0:
            logger.error(
                f"vcf2maf did not produce a non-empty MAF for sample {sample} "
                f"(expected at {out_file}) - this sample's mutations will be "
                "missing from the final study.")
            return False
    return True


def _read_vcf_header_and_data(vcf_path: str) -> tuple[list[str], list[str]]:
    """Split a VCF's lines into (header lines incl. #CHROM row, data rows)."""
    with open(vcf_path, encoding="utf-8") as f:
        lines = f.readlines()
    header_lines = [ln for ln in lines if ln.startswith("#")]
    data_lines = [ln for ln in lines if not ln.startswith("#")]
    return header_lines, data_lines


def split_colocated_variants(vcf_path: str, work_dir: str) -> list[str]:
    """vcf2maf keeps only the first VCF record for a given (CHROM, POS) when
    two or more records share it - a known limitation, seen in practice on
    a real Guardant case where the same homopolymer-adjacent indel position
    carried two independently valid (Guardant call==1) calls (an insertion
    and a deletion), and only the first survived into the MAF.

    This scans vcf_path for such groups and writes every record AFTER the
    first one in each group into its own single-variant VCF under
    work_dir, so it can be run through vcf2maf separately by the caller.
    The original vcf_path is left untouched
    (its first-of-group record still goes through the normal single
    vcf2maf run exactly as before - no change to today's behavior for the
    common, non-colliding case).

    Returns the list of extra VCF paths created (empty if no collision).
    """
    header_lines, data_lines = _read_vcf_header_and_data(vcf_path)
    by_pos: dict[tuple[str, str], list[str]] = {}
    for line in data_lines:
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        by_pos.setdefault((fields[0], fields[1]), []).append(line)

    extra_paths = []
    base = Path(vcf_path)
    for (chrom, pos), rows in by_pos.items():
        if len(rows) < 2:
            continue
        for i, row in enumerate(rows[1:], start=1):
            extra_path = Path(work_dir) / f"{base.stem}.colocated_{chrom}_{pos}_{i}{base.suffix}"
            with open(extra_path, "w", encoding="utf-8") as f:
                f.writelines(header_lines)
                f.write(row)
            extra_paths.append(str(extra_path))
    return extra_paths


def recover_colocated_variants(v: str, temporary: str, output_folder: str,
                               sample: str, tumor_id: str | None = None) -> None:
    """Call after the normal vcf2maf run for sample `sample` on VCF `v`:
    detects any variant(s) that shared a POS with another valid call in
    the same VCF (see split_colocated_variants) and, for each one, runs
    vcf2maf again on it alone, writing its own MAF into output_folder's
    `maf/` folder next to the sample's main one - recovering calls that
    vcf2maf's own same-POS handling would otherwise silently drop.

    No merging into the main MAF is done here on purpose: the later
    OncoKB-annotation + concatenation stage (filter_clinvar.py /
    concatenate.py) already globs every *.maf file under `maf/` and folds
    them all into data_mutations_extended.txt regardless of how many
    files exist per sample - writing this extra MAF there is enough for
    it to be picked up and annotated exactly like any other, with no risk
    of double-counting the sample's main variants.
    """
    extra_vcfs = split_colocated_variants(v, temporary)
    if not extra_vcfs:
        return
    logger.warning(
        f"{sample}: found {len(extra_vcfs)} variant(s) sharing a genomic "
        "position with another valid call in this sample's VCF - vcf2maf "
        "only keeps one call per position, so re-running it separately "
        "for each extra one; its MAF will be picked up by the normal "
        "OncoKB-annotation/concatenation stage like any other.")
    for extra_vcf in extra_vcfs:
        cl_extra = vcf2maf_constructor(extra_vcf, temporary, output_folder, tumor_id)
        run_vcf2maf(cl_extra, sample)


def extract_multiple_snv(multiple_vcf: str, input_dir: str) -> None:
    """Split a multi-sample SNV VCF file into single-sample files.

    Args:
        multiple_vcf (str): Path to the multi-sample VCF file.
        input_dir (str): Directory where single-sample VCFs are saved.

    """
    input_dir = Path(input_dir)
    vcf_dir = input_dir / "single_sample_vcf"

    vcf_dir.mkdir(parents=True, exist_ok=True)

    sample_id_path = input_dir / "sample_id.txt"
    try:
        with sample_id_path.open("w", encoding="utf-8") as f:
            subprocess.run(["vcf-query", "-l", str(multiple_vcf)], stdout=f, check=True)
    except subprocess.CalledProcessError as e:
        msg = f"Error during vcf-query: {e}"
        raise RuntimeError(msg) from e

    sample_id_file = input_dir / "sample_id.txt"

    with sample_id_file.open(encoding="utf-8") as file:
        lines = file.readlines()
        for sample_line in lines:
            sample = sample_line.strip()
            output_vcf = vcf_dir / f"{sample}.vcf"

            with output_vcf.open("w", encoding="utf-8") as out_f:
                try:
                    subprocess.run(
                        [
                            "vcf-subset",
                            "--exclude-ref",
                            "-c", sample,
                            str(multiple_vcf)],
                        stdout=out_f,
                        check=True)
                except subprocess.CalledProcessError as e:
                    msg = f"Error during vcf-subset for {sample}: {e}"
                    raise RuntimeError(msg)


def _walk_process_snv(ctx: WalkContext) -> None:
    """Stage: SNV calls -> per-sample vcf2maf (writes maf/*.maf), optionally
    followed by per-sample SigMA mutational-signature analysis.

    Depends only on WalkContext, independent of _walk_process_cnv/_fusion.
    Runs the three steps below in sequence, one sample at a time - the
    Snakemake path (see Snakefile) runs _walk_process_snv_sample() as one
    job per sample instead, in parallel, through the same three functions.
    """
    s_id_path_snv = _walk_snv_samples(ctx)
    job_keys = snv_job_keys(list(s_id_path_snv))
    for job_key, (k, v) in zip(job_keys, s_id_path_snv.items()):
        # varan.py's sequential path keeps going with the other samples
        # (the error is logged); under Snakemake the sample's own job fails.
        _walk_process_snv_sample(ctx, k, v, job_key, fail_on_missing_maf=False)
    _walk_snv_finish(ctx, job_keys)


def _walk_snv_samples(ctx: WalkContext) -> dict[str, str]:
    """SNV step 1 of 3: map each sample ID to the SNV VCF this run has to
    convert. Empty when there's nothing to convert: no SNV input, an
    analysis type (-t) without SNVs, or resume=True (MAFs reused as-is).
    """
    if not ctx.input_folder_snv.exists() or ctx.vcf_type in ["cnv", "fus", "tab"]:
        return {}

    logger.info("Managing SNV files...")
    s_id_path_snv = snv_type_from_folder(
        ctx.input_folder_snv, ctx.case_folder_arr, ctx.output_folder,
        ctx.snv_sample_ids)

    logger.info("Checking maf folder...")
    maf_path = Path(ctx.output_folder) / "maf"
    if maf_path.is_dir() and any(f.suffix == ".maf" for f in maf_path.iterdir()):
        logger.info("A non-empty maf folder already exists!")

    if ctx.resume:
        return {}
    return s_id_path_snv


def snv_job_keys(sample_ids: list[str]) -> list[str]:
    """One short, filename-safe key per SNV sample, in order: the sample
    IDs snv_type_from_folder() returns are SAMPLE_IDs, but a file it
    couldn't tie to one keeps a whole-path id (".../S1.snv.bam"), and
    either may hold characters unusable in file or Snakemake job names. Used by
    _walk_process_snv_sample()/_walk_snv_finish() and walk_stage.py."""
    keys = []
    for i, sample_id in enumerate(sample_ids, start=1):
        name = re.sub(r"[^A-Za-z0-9._-]", "_", Path(sample_id.replace("\\", "/")).name)
        keys.append(f"{i:04d}_{name}")
    return keys


def _walk_process_snv_sample(
    ctx: WalkContext, sample_id: str, vcf_path: str, job_key: str,
    fail_on_missing_maf: bool = True,
) -> None:
    """SNV step 2 of 3, for one sample: optional "d" VCF filter, vcf2maf
    (plus co-located variant recovery) and, with -g/--sigma, SigMA.

    Safe to run for several samples at once (one Snakemake job each): every
    call uses its own scratch subfolder and removes only that one -
    clear_scratch() would also remove the shared scratch/ parent, under
    other samples still running - and saves its SigMA row to its own file,
    combined later by _walk_snv_finish(). job_key is this sample's entry
    from snv_job_keys(), used to name that file.
    """
    if "d" in ctx.filters:
        logger.info(f"{sample_id}: filtering out VCF rows with dots in ALT column")
        vcf_path = vcf_filtering(
            {sample_id: vcf_path}, ctx.output_folder, output_filtered)[sample_id]

    temporary = create_random_name_folder(ctx.output_folder)
    # SAMPLE_ID as the MAF barcode only when sample.tsv really tied this
    # file to it (otherwise sample_id is just derived from the file name).
    tumor_id = sample_id if sample_id in ctx.snv_sample_ids.values() else None
    cl = vcf2maf_constructor(vcf_path, temporary, ctx.output_folder, tumor_id)
    if not run_vcf2maf(cl, sample_id) and fail_on_missing_maf:
        # One Snakemake job per sample: failing it makes a re-run convert
        # just this sample again, instead of the job counting as done with
        # the sample's mutations silently missing.
        msg = f"vcf2maf produced no MAF for sample {sample_id} - see the error above."
        raise RuntimeError(msg)
    recover_colocated_variants(
        vcf_path, temporary, ctx.output_folder, sample_id, tumor_id)

    if ctx.sigma:
        row = _run_sigma_for_snv_sample(ctx, cl, sample_id)
        rows_dir = _sigma_sample_rows_dir(ctx.output_folder)
        rows_dir.mkdir(parents=True, exist_ok=True)
        with (rows_dir / f"{job_key}.pkl").open("wb") as f:
            pickle.dump(row, f)

    # Only on success, like clear_scratch()'s callers: a failed sample's
    # scratch folder is left in place for inspection.
    shutil.rmtree(temporary, ignore_errors=True)


def _walk_snv_finish(ctx: WalkContext, job_keys: list[str]) -> None:
    """SNV step 3 of 3, once every sample is done: combine the per-sample
    SigMA rows into data_sigma.txt (in job_keys' order, i.e. sample order)
    and remove the now-empty scratch/ folder."""
    if ctx.sigma and job_keys:
        rows_dir = _sigma_sample_rows_dir(ctx.output_folder)
        sigma_results = []
        for job_key in job_keys:
            with (rows_dir / f"{job_key}.pkl").open("rb") as f:
                sigma_results.append(pickle.load(f))
        _write_sigma_intermediate(ctx.output_folder, sigma_results)
        shutil.rmtree(rows_dir, ignore_errors=True)

    logger.info("Clearing scratch folder...")
    scratch = Path(ctx.output_folder) / "scratch"
    if scratch.is_dir() and not any(scratch.iterdir()):
        scratch.rmdir()


def _sigma_sample_rows_dir(output_folder: str) -> Path:
    return Path(output_folder) / "intermediate" / "sigma" / "samples"


def _run_sigma_for_snv_sample(ctx: WalkContext, cl: list, sample_id: str) -> dict:
    """Run SigMA for one sample right after its vcf2maf conversion.

    Only called when ctx.sigma is True (the -g/--sigma CLI flag). Reads
    the MAF vcf2maf_constructor/run_vcf2maf were just told to produce (the
    same "--output-maf" path baked into `cl`, so this can never drift from
    what vcf2maf actually wrote to), resolves the sample's ONCOTREE_CODE
    from sample.tsv, and delegates the rest to
    sigma_runner.run_sigma_for_sample(). Never raises - any problem is
    logged and reflected in the returned row's SIGMA_STATUS instead, so
    one sample's SigMA failure never stops the batch (matches this
    codebase's existing log-and-skip philosophy, e.g. the VAF filter's
    missing-column handling in sigma_filter.prepare_sigma_maf).

    Args:
        ctx (WalkContext): Shared setup state.
        cl (list): The vcf2maf command-line list built by
            vcf2maf_constructor() for this sample - used only to recover
            the "--output-maf" path it specifies.
        sample_id (str): The sample's SAMPLE_ID.

    Returns:
        dict: A SigMA result row, see sigma_runner.run_sigma_for_sample().

    """
    try:
        maf_out_path = Path(cl[cl.index("--output-maf") + 1])
    except (ValueError, IndexError):
        logger.warning(
            f"Sample {sample_id}: could not determine the vcf2maf output MAF "
            "path - skipping SigMA for this sample.")
        return sigma_runner.blank_result(sample_id, "ERROR")

    if not maf_out_path.exists() or maf_out_path.stat().st_size == 0:
        logger.warning(
            f"Sample {sample_id}: no (or empty) MAF at {maf_out_path} after "
            "vcf2maf - vcf2maf likely failed for this sample. Skipping "
            "SigMA for this sample.")
        return sigma_runner.blank_result(sample_id, "ERROR")

    maf_df = read_maf(maf_out_path)

    oncotree_code = ""
    match = ctx.clin_file.loc[ctx.clin_file["SAMPLE_ID"].astype(str) == str(sample_id)]
    if not match.empty and "ONCOTREE_CODE" in match.columns:
        oncotree_code = str(match.iloc[0]["ONCOTREE_CODE"])

    if not oncotree_code or oncotree_code.lower() == "nan":
        logger.warning(
            f"Sample {sample_id}: no ONCOTREE_CODE found in sample.tsv - "
            "SigMA's tumor_type mapping will fall through to its "
            "fallback/'other' handling (see get_sigma_call_params).")

    kind = specimen.specimen_types(ctx.clin_file).get(str(sample_id).strip())
    if clonal_hematopoiesis.exclude_downstream():
        wbc = clonal_hematopoiesis.samples_from_table(ctx.clin_file, {}).wbc.get(str(sample_id))
        maf_df = clonal_hematopoiesis.drop_flagged(maf_df, sample_id, kind, wbc)

    sigma_dir = Path(ctx.output_folder) / "intermediate" / "sigma"
    return sigma_runner.run_sigma_for_sample(
        maf_df, sample_id, oncotree_code, sigma_dir, kind)


def _write_sigma_intermediate(output_folder: str, sigma_results: list[dict]) -> None:
    """Write this run's per-sample SigMA results to a shared intermediate
    file, so _walk_write_clinical_tables (a separate Snakemake stage/
    process, see WalkContext's docstring) can pick them up without any
    in-memory state - the same file-based handoff already used for CNA's
    intermediate/ artifacts.

    A no-op (no file written) when sigma_results is empty, e.g. because
    the SNV folder was empty - data_clinical_sample.txt then simply gets
    no SIGMA_* columns at all, exactly as if -g/--sigma had not been
    passed for this run.
    """
    if not sigma_results:
        return

    sigma_dir = Path(output_folder) / "intermediate" / "sigma"
    sigma_dir.mkdir(parents=True, exist_ok=True)
    out_path = sigma_dir / "data_sigma.txt"
    pd.DataFrame(sigma_results).to_csv(out_path, sep="\t", index=False)
    logger.info(
        f"Wrote SigMA results for {len(sigma_results)} sample(s) to {out_path}")
