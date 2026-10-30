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

"""Walk stage: RNA fusions and splice variants -> data_sv.txt.
"""

from __future__ import annotations

import contextlib
import os
from pathlib import Path

import pandas as pd
from loguru import logger

import tsv
from config_loader import config
from filter_clinvar import filter_oncokb
from populate_case_lists import read_profiled_samples
from thresholds import meets_threshold
from walk_common import WalkContext, get_combined_variant_output_from_folder


# data_sv.txt's columns as Varan writes them from CombinedVariantOutput.
# The Site*_Chromosome/Position columns (the fusion/splice breakpoints)
# stay blank when a source doesn't give them (fusion file, Guardant).
SV_COLUMNS = [
    "Sample_Id", "SV_Status", "Class", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol",
    "Normal_Paired_End_Read_Count", "Event_Info", "RNA_Support",
    "Site1_Chromosome", "Site1_Position", "Site2_Chromosome", "Site2_Position"]


def _chromosome(name: str) -> str:
    """Chromosome as cBioPortal's SV format writes it: no "chr" prefix."""
    name = str(name).strip()
    return name[3:] if name.lower().startswith("chr") else name


def _split_breakpoint(breakpoint: str) -> tuple[str, str]:
    """"chr7:55019278" -> ("7", "55019278"); ("", "") if not in that form."""
    chrom, _, pos = str(breakpoint).strip().partition(":")
    return (_chromosome(chrom), pos) if pos else ("", "")


def _warn_fusion_conflicts(
    fusion_folder: Path, fusion_files: list, combined_dict: dict[str, str]) -> None:
    """Log, per sample, the fusions the fusion file and the sample's
    CombinedVariantOutput disagree on - only for samples whose
    CombinedVariantOutput reports a fusion analysis (not "NA")."""
    from_file: dict[str, set] = {}
    for fusion_file in fusion_files:
        df = pd.read_csv(Path(fusion_folder) / fusion_file, sep="\t", dtype=str)
        if not {"Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol"} <= set(df.columns):
            continue
        for row in df.itertuples(index=False):
            from_file.setdefault(str(row.Sample_Id).strip(), set()).add(
                frozenset({row.Site1_Hugo_Symbol, row.Site2_Hugo_Symbol}))
    for sample_id, path in combined_dict.items():
        if not Path(path).is_file():
            continue
        if not _reports_fusions(Path(path)):
            if from_file.get(str(sample_id)):
                logger.warning(
                    f"{sample_id}: the fusion file lists fusions, but the "
                    "sample's CombinedVariantOutput reports no fusion analysis "
                    "(NA) - the fusion file's calls are used.")
            continue
        try:
            from_combined = {frozenset({f["Site1_Hugo_Symbol"], f["Site2_Hugo_Symbol"]})
                             for f in tsv.get_fusions(Path(path))}
        except Exception as err:
            logger.warning(f"{sample_id}: could not read CombinedVariantOutput fusions ({err}).")
            continue
        in_file = from_file.get(str(sample_id), set())
        if in_file != from_combined:
            def fmt(pairs: set) -> str:
                return ", ".join(sorted("-".join(sorted(p)) for p in pairs)) or "none"
            logger.warning(
                f"{sample_id}: fusion file and CombinedVariantOutput disagree - "
                f"file: {fmt(in_file)}; CombinedVariantOutput: {fmt(from_combined)}. "
                "The fusion file's calls are used.")


def _reports_fusions(combined_output: Path) -> bool:
    """Whether a CombinedVariantOutput's [Fusions] section holds an
    analysis (rows or an empty list) rather than "NA" (RNA not run)."""
    lines = combined_output.read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if "[Fusions]" in line or "[Data Fusions]" in line:
            return not (i + 2 < len(lines) and lines[i + 2].strip() == "NA")
    return False


def fill_fusion_from_temp(
    input_path: str,
    fusion_table_file: str,
    clin_file: pd.DataFrame,
    fusion_files: list,
    thr_fus: str) -> None:
    """Collect fusion events from temp directory into a unified table.

    Args:
        input_path (str): Temp input path.
        fusion_table_file (str): Output file path for combined fusions.
        clin_file (pd.DataFrame): Clinical file to match SAMPLE_IDs.
        fusion_files (list): List of fusion TSV filenames.
        thr_fus (str): Threshold expression for filtering read count (e.g.
            ">=15"), same conf.ini [FUSION] THRESHOLD_FUSION value and same
            eval pattern as fill_fusion_from_combined() - previously this
            path had its own separate hardcoded `min_read_count = 15`,
            silently dropping a Guardant-confirmed fusion (call=1) with
            fewer supporting molecules regardless of what THRESHOLD_FUSION
            was configured to, and regardless of any vendor's own
            confidence flag for that call.

    """
    logger.info(f"Found {len(fusion_files)} Fusion file(s)")

    fusion_input = Path(input_path) / "FUSIONS" / fusion_files[0]
    with fusion_input.open() as template:
        header = template.readline()

    with Path(fusion_table_file).open("w") as fusion_table:
        fusion_table.write(header)

        for fusion_file in fusion_files:
            ff = pd.read_csv(Path(input_path) / "FUSIONS" / fusion_file,
                              sep="\t", dtype=str)

            required_columns = {
                "Sample_Id",
                "SV_Status",
                "Site1_Hugo_Symbol",
                "Site2_Hugo_Symbol"}
            if not required_columns.issubset(ff.columns):
                logger.warning(f"{fusion_file} does not contain required columns")
                continue

            if ff.shape[0]==0:
                logger.info(f"No Fusions found in {fusion_file}")
                continue
            logger.info(f"Fusions found in {fusion_file}")
            # The read count column is optional in the fusion file template:
            # without a value THRESHOLD_FUSION can't be applied, so the
            # fusion is kept (and that's logged once).
            reads_col = "Normal_Paired_End_Read_Count"
            known_ids = clin_file["SAMPLE_ID"].astype(str).to_numpy()
            unfiltered = 0
            for fus in ff.itertuples(index=False):
                if str(fus.Sample_Id).strip() not in known_ids:
                    continue
                reads = str(getattr(fus, reads_col, "") or "").strip()
                if reads in ("", "nan"):
                    unfiltered += 1
                elif not meets_threshold(reads, thr_fus):
                    continue
                fusion_table.write("\t".join(map(str, fus)) + "\n")
            if unfiltered:
                logger.warning(
                    f"{fusion_file}: {unfiltered} fusion(s) without {reads_col} - "
                    f"kept without applying THRESHOLD_FUSION ({thr_fus}).")


def annotate_fusion(
    cancer: str,
    fusion_table_file: str,
    data_sv: pd.DataFrame,
    input_file: pd.DataFrame) -> str:
    """Annotate fusion events with OncoKB FusionAnnotator.

    Args:
        cancer (str): Cancer type code.
        fusion_table_file (str): Path to fusion TSV input file.
        data_sv (pd.DataFrame): DataFrame with fusion data.
        input_file (pd.DataFrame): Clinical data with SAMPLE_ID and optionally
        ONCOTREE_CODE.

    Returns:
        str: Path to annotated fusion output file.

    """
    if "ONCOTREE_CODE" in input_file.columns:
        fusion_table_df = data_sv.merge(
            input_file[["SAMPLE_ID", "ONCOTREE_CODE"]],
            how="inner",
            left_on="Sample_Id",
            right_on="SAMPLE_ID")
        fusion_table_df.to_csv(fusion_table_file, sep="\t", index=False)
        fusion_table_file_out = fusion_table_file.with_suffix(".ann.txt")
        os.system(
            f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file} "
            f"-o {fusion_table_file_out} -b {config.get('OncoKB', 'ONCOKB')}")

    else:
        fusion_table_file_out = fusion_table_file.with_suffix(".ann.txt")
        os.system(
            f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file} "
            f"-o {fusion_table_file_out} -t {cancer.upper()}  "
            f"-b {config.get('OncoKB', 'ONCOKB')}")

    return fusion_table_file_out


def fill_fusion_from_combined(
    fusion_table_file: str,
    combined_dict: dict[str, str],
    thr_fus: str) -> None:
    """Extract fusion events from combined variant output and write to TSV.

    Args:
        fusion_table_file (str): Output file path for fusion TSV.
        combined_dict (dict[str, str]): Map sampleID to CombinedVariantOutput file path.
        thr_fus (str): Threshold expression for filtering read count (e.g. ">5").

    Returns:
        None

    """
    logger.info("Writing data_sv.txt file...")

    fusion_table_path = Path(fusion_table_file)
    with fusion_table_path.open("w") as fusion_table:
        fusion_table.write("\t".join(SV_COLUMNS) + "\n")

        for k, v in combined_dict.items():
            fusions=[]
            try:
                fusions = tsv.get_fusions(Path(v))
            except Exception as err:
                logger.error("Something went wrong while reading Fusion section "
                f"for sample {k}: {err}")
            if len(fusions) == 0:
                continue

            for fus in fusions:
                if len(fusions) > 0:
                    site1_hugo_symbol = fus["Site1_Hugo_Symbol"]
                    site2_hugo_symbol = fus["Site2_Hugo_Symbol"]
                    if site2_hugo_symbol == "CASC1":
                        site2_hugo_symbol = "DNAI7"
                    if meets_threshold(fus["Normal_Paired_End_Read_Count"], thr_fus):
                        fusion_table.write("\t".join([
                            str(k).strip(), "SOMATIC", "FUSION",
                            str(site1_hugo_symbol), str(site2_hugo_symbol),
                            fus["Normal_Paired_End_Read_Count"],
                            fus["Event_Info"] + " Fusion", "Yes",
                            _chromosome(fus["Site1_Chromosome"]), fus["Site1_Position"],
                            _chromosome(fus["Site2_Chromosome"]), fus["Site2_Position"],
                        ]) + "\n")


def fill_splice_from_combined(
    fusion_table_file: str,
    combined_dict: dict[str, str],
    thr_splice: str) -> None:
    """Add splice variant events from combined variant output to data_sv.txt.

    Splice variants are written as intragenic structural variant rows (same
    gene on both Site1/Site2, Class="SPLICE") into the same data_sv.txt table
    fusions use, but deliberately called *after* the fusion pipeline's
    OncoKB annotation step (annotate_fusion/FusionAnnotator.py) rather than
    folded into it - that annotator expects a real two-gene fusion pair, and
    running a same-gene "fusion" through it would be meaningless at best.

    Rebuilds the file rather than blindly appending, so calling this twice
    on the same output folder (e.g. on resume) can't accumulate duplicate
    splice rows: any pre-existing Class="SPLICE" rows are dropped first,
    then this run's rows are added back.

    NOTE ON VERIFICATION: see tsv.get_splice_variants() - the column names
    are copied verbatim from real CombinedVariantOutput.tsv headers and
    cross-checked against Illumina's own release notes, but the row-split
    for a populated [Splice Variants] line has never been seen in a real
    file. Every parsed row is logged at INFO level below specifically so
    the first real one can be spot-checked against the source file.

    Args:
        fusion_table_file (str): data_sv.txt path (may not exist yet if the
            fusion pipeline found no fusions and removed it).
        combined_dict (dict[str, str]): Map sampleID to CombinedVariantOutput file path.
        thr_splice (str): Threshold expression for filtering read count (e.g. ">=15").

    Returns:
        None

    """
    logger.info("Adding splice variants to data_sv.txt file...")
    fusion_table_path = Path(fusion_table_file)
    # Appended by column name, not position: data_sv.txt may come from a
    # fusion file with its own columns, or carry OncoKB columns - every
    # splice row then simply gets blanks for the columns it doesn't fill.
    existing = (pd.read_csv(fusion_table_path, sep="\t", dtype=str, keep_default_na=False)
                if fusion_table_path.exists() else pd.DataFrame(columns=SV_COLUMNS))
    if "Class" in existing.columns:
        existing = existing[existing["Class"] != "SPLICE"]
    new_rows = []
    for k, v in combined_dict.items():
        splice_variants = []
        try:
            splice_variants = tsv.get_splice_variants(Path(v))
        except Exception as err:
            logger.error("Something went wrong while reading Splice Variants "
            f"section for sample {k}: {err}")
        for sv in splice_variants:
            gene = sv["Gene"]
            ssr = sv["Splice_Supporting_Reads"]
            try:
                if not meets_threshold(ssr, thr_splice):
                    continue
            except (ValueError, TypeError):
                logger.warning(
                    f"Non-numeric splice supporting read count for {gene} "
                    f"in sample {k}, skipping.")
                continue
            event_info = f"Exon {sv['Affected_Exon']} splice variant"
            logger.info(
                f"Splice variant found for sample {k}: gene={gene}, "
                f"exon={sv['Affected_Exon']}, breakpoints="
                f"{sv['Breakpoint_1']}/{sv['Breakpoint_2']}, "
                f"supporting_reads={ssr}, "
                f"reference_reads_transcript={sv['Reference_Reads_Transcript']}.")
            chrom1, pos1 = _split_breakpoint(sv["Breakpoint_1"])
            chrom2, pos2 = _split_breakpoint(sv["Breakpoint_2"])
            new_rows.append({
                "Sample_Id": str(k).strip(), "SV_Status": "SOMATIC", "Class": "SPLICE",
                "Site1_Hugo_Symbol": str(gene), "Site2_Hugo_Symbol": str(gene),
                "Normal_Paired_End_Read_Count": ssr, "Event_Info": event_info,
                "RNA_Support": "Yes", "Site1_Chromosome": chrom1, "Site1_Position": pos1,
                "Site2_Chromosome": chrom2, "Site2_Position": pos2})
    if not new_rows:
        logger.info(
            "No splice variants found for this batch - data_sv.txt left as "
            "is (not created if the fusion step above didn't create it "
            "either).")
        return
    merged = pd.concat([existing, pd.DataFrame(new_rows)], ignore_index=True).fillna("")
    merged.to_csv(fusion_table_path, sep="\t", index=False)
    logger.info(
        f"data_sv.txt now has {len(new_rows)} new splice variant row(s) "
        f"(plus {len(existing)} pre-existing non-splice row(s) carried "
        "over from the fusion step above, if any).")


def _walk_process_fusion(ctx: WalkContext) -> None:
    """Stage: RNA fusions + splice variants -> data_sv.txt.

    Depends only on WalkContext, independent of _walk_process_cnv/_snv.
    No-op if vcf_type excludes fusions, matching the original gate.
    """
    if ctx.vcf_type in ["cnv", "snv", "tab"]:
        return

    fusion_table_file = Path(ctx.output_folder) / "data_sv.txt"
    fusion_folder = Path(ctx.input_folder) / "FUSIONS"
    combined_dict = {}
    # Shared by both branches below - same conf.ini [FUSION] THRESHOLD_FUSION
    # value applies regardless of which fusion source produced the table.
    thr_fus = config.get("FUSION", "THRESHOLD_FUSION")

    # Fusion source: a fusion file (FUSIONS/*.tsv) has priority over the
    # CombinedVariantOutput [Fusions] sections - when both are there, the
    # samples on which they disagree are logged. Splice variants only ever
    # come from CombinedVariantOutput (below).
    combined_output_folder = Path(ctx.input_folder) / "CombinedOutput"
    if (
        combined_output_folder.exists()
        and any(f.is_file() for f in combined_output_folder.iterdir())):
        combined_dict = get_combined_variant_output_from_folder(
            ctx.input_folder, ctx.clin_file, ctx.isinputfile)
    fusion_files = (
        [f for f in fusion_folder.iterdir() if f.suffix == ".tsv"]
        if fusion_folder.exists() else [])

    if fusion_files:
        if combined_dict:
            logger.info(
                f"Fusions taken from {fusion_files[0].name}: a fusion file has "
                "priority over the CombinedVariantOutput [Fusions] sections "
                "(still used for splice variants).")
            _warn_fusion_conflicts(fusion_folder, fusion_files, combined_dict)
        else:
            logger.info(f"Fusions taken from {fusion_files[0].name} (the only fusion source).")
        fill_fusion_from_temp(
            ctx.input_folder, fusion_table_file, ctx.clin_file, fusion_files,
            thr_fus)
    elif combined_dict:
        logger.info(
            "Fusions taken from the CombinedVariantOutput [Fusions] sections "
            "(no fusion file given).")
        fill_fusion_from_combined(fusion_table_file, combined_dict, thr_fus)

    # Report the fusion outcome on its own - fusions and splice variants
    # (below) are two independent CombinedOutput data types that happen to
    # share the same data_sv.txt file. Neither one deletes the file itself
    # if it finds nothing: only the very end of this function does, once
    # both have had their chance to write to it.
    n_fusions = 0
    if fusion_table_file.exists():
        with fusion_table_file.open() as data_sv:
            n_fusions = max(len(data_sv.readlines()) - 1, 0)
    if n_fusions:
        logger.info(f"{n_fusions} fusion call(s) found for this batch.")
    else:
        logger.info("No fusion calls found for this batch.")

    if n_fusions:
        # Dedup always runs here, regardless of OncoKB annotation - a run
        # without --oncokb must not ship duplicate rows in data_sv.txt either,
        # since cBioPortal validation rejects those the same way either way.
        if ctx.oncokb:
            data_sv = pd.read_csv(fusion_table_file, sep="\t", dtype=str)
            input_file = pd.read_csv(ctx.clin_sample_path, sep="\t", dtype=str)
            fusion_table_file_out = annotate_fusion(
                ctx.cancer, fusion_table_file, data_sv, input_file)

            if "o" in ctx.filters:
                fus_file = pd.read_csv(fusion_table_file_out, sep="\t", dtype=str)
                fus_file = filter_oncokb(fus_file, "FUSION", "ONCOKB_FILTER_FUSION")
                fus_file.to_csv(fusion_table_file_out, index=False, sep="\t")

            data_sv_tmp = pd.read_csv(fusion_table_file_out, sep="\t", dtype=str)
            with contextlib.suppress(KeyError):
                data_sv_tmp = data_sv_tmp.drop(
                    ["SAMPLE_ID", "ONCOTREE_CODE"], axis=1)
        else:
            fusion_table_file_out = fusion_table_file
            data_sv_tmp = pd.read_csv(fusion_table_file_out, sep="\t", dtype=str)

        if "Normal_Paired_End_Read_Count" in data_sv_tmp.columns:
            data_sv_tmp["Normal_Paired_End_Read_Count"] = pd.to_numeric(data_sv_tmp["Normal_Paired_End_Read_Count"], errors='coerce')
            data_sv_tmp = data_sv_tmp.sort_values(by="Normal_Paired_End_Read_Count", ascending=False)
            col_subset = [col for col in data_sv_tmp.columns if col != "Normal_Paired_End_Read_Count"]

            data_sv_tmp = data_sv_tmp.drop_duplicates(subset=col_subset, keep='first')
            data_sv_tmp["Normal_Paired_End_Read_Count"] = data_sv_tmp["Normal_Paired_End_Read_Count"].astype(str).str.replace(r'\.0$', '', regex=True).replace('nan', '')

        else:
            data_sv_tmp = data_sv_tmp.drop_duplicates(keep='first')

        data_sv_tmp.to_csv(fusion_table_file_out, index=False, sep="\t")
        if fusion_table_file_out != fusion_table_file:
            Path(fusion_table_file_out).replace(fusion_table_file)

    # Runs after the fusion block above (annotation included) has fully
    # finished, so splice rows are never sent through the fusion-specific
    # OncoKB annotator. Logs its own found/not-found outcome independently
    # of the fusion one above (see fill_splice_from_combined).
    if combined_dict:
        thr_splice = config.get("SPLICE", "THRESHOLD_SPLICE")
        fill_splice_from_combined(fusion_table_file, combined_dict, thr_splice)

    # Only now, after both fusions and splice variants have had their
    # chance to write something, remove data_sv.txt if it ended up with no
    # rows at all - covers both "neither found anything" and "fusions were
    # found but all got filtered out by the OncoKB filter above". Kept
    # header-only instead when samples were tested for them (see
    # populate_case_lists.ensure_sv_file()): that's the study's "0 of N
    # tested" structural-variant profile.
    if fusion_table_file.exists():
        # Counted first, file closed before it may be removed below (an
        # open file can't be deleted on Windows).
        with fusion_table_file.open() as data_sv:
            header_only = len(data_sv.readlines()) == 1
        if header_only:
            tested = read_profiled_samples(ctx.output_folder).get("sv")
            if tested:
                logger.warning(
                    "No fusions or splice variants found for this batch - "
                    f"data_sv.txt kept header-only: {len(tested)} sample(s) "
                    "were tested for them.")
            else:
                fusion_table_file.unlink()
                logger.warning(
                    "No fusions or splice variants found for this batch - "
                    "data_sv.txt removed.")
