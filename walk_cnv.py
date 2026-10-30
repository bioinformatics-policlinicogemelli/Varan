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

"""Walk stage: CNV calls -> data_cna_hg19.seg(.fc.txt), data_cna.txt and the
BRCA exon-level CNV table.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

import tsv
import vcf2tab_cnv
from config_loader import VARAN_DIR, config
from filter_clinvar import check_bool, filter_oncokb
from walk_common import WalkContext, get_combined_variant_output_from_folder, ploidy


def get_cnv_from_folder(input_foldercnv: str) -> list:
    """Return a list of all VCF files in a directory.

    Args:
        input_foldercnv (str): Path to the CNV input folder.

    Returns:
        list: List of VCF filenames.

    """
    files = Path(input_foldercnv).iterdir()
    return[file.name for file in files if file.is_file() and file.name.endswith("vcf")]


def get_sample_id_from_cnv(cnv_vcf: str) -> str:
    """Extract sample ID from CNV.

    Args:
        cnv_vcf (str): VCF filename.

    Returns:
        str: Corresponding BAM filename.

    """
    if "_CopyNumberVariants.vcf" in cnv_vcf:
        sample=cnv_vcf.replace("_CopyNumberVariants.vcf", ".bam")
    else:
        sample="bam".join(cnv_vcf.rsplit("vcf", 1))
    return sample


def cnv_type_from_folder(input_path: str,
                         cnv_vcf_files: list,
                         output_folder: str,
                         oncokb: bool,
                         cancer: str,
                         multiple: bool,
                         filters: str = "",
                         id_map: dict[str, str] | None = None,
                         cnv_dir: str | None = None) -> dict:
    """Process CNV, converting them in CNA tables and performing annotation.

    Args:
    input_path : str
        Path to the input directory or the sample TSV file (if CNVKIT_algorithm=True).
    cnv_vcf_files : list
        List of CNV VCF file paths or file names (relative to CNV folder).
    output_folder : str
        Directory to store the output files.
    oncokb : bool
        If True, annotate results using OncoKB.
    cancer : str
        Default cancer type to use if ONCOTREE_CODE is not available.
    multiple : bool
        If True, indicates that VCFs are in "CNV/single_sample_vcf"; otherwise in "CNV".
    filters : str, optional
        Filter options as string, same as the SNV/fusion `filters` CLI
        option. Oncogenicity annotation (OncoKB) still runs whenever
        `oncokb` is True - that alone must not drop any row - but rows
        are only filtered down to Oncogenic/Likely Oncogenic when "o" is
        also present here, mirroring the SNV and fusion behavior.
        Defaults to "" (annotate but do not filter).
    cnv_dir : str, optional
        The folder the CNV VCFs (cnv_vcf_files) are in - e.g. a -m split's
        output inside the study's temp/. Default: <input_path>/CNV (or its
        single_sample_vcf/ subfolder with `multiple`).
    id_map : dict, optional
        CNV VCF file name -> sample.tsv SAMPLE_ID (see
        walk.sample_ids_by_file()). A file not in it gets an id derived
        from its own name (get_sample_id_from_cnv()).

    Returns:
    dict
        A mapping from sample ID to VCF file path used in processing.

    """
    counter = 0
    sid_path = {}

    # vcf_to_table()/vcf_to_table_fc() append to these two files whenever
    # they already exist on disk (that's how they accumulate rows across
    # the per-sample loop below). On a fresh output folder that's fine,
    # but -R/--resume (and the "temp/ already exists" resume path in
    # transform_input) can hand this function an output_folder left over
    # from an earlier CNV pass - without clearing these first, that stale
    # data just gets appended to again, silently tripling (or worse)
    # every sample's segments with stale seg.mean values mixed in.
    for stale in ("data_cna_hg19.seg", "data_cna_hg19.seg.fc.txt"):
        stale_path = Path(output_folder) / stale
        if stale_path.exists():
            stale_path.unlink()

    # Loaded ONCE for the whole batch instead of once per sample inside
    # vcf_to_table_fc(): TC availability (missing sample.tsv / missing TC
    # column) is a batch-wide fact, not a per-sample one, so it should log
    # ONCE - actual per-sample TC gaps are collected below and reported
    # either as one aggregate line (missing for every sample) or one line
    # per affected sample (missing for only some), matching how missing
    # SNV/CNV/CombinedOutput paths are already reported in transform_input.
    tc_lookup, tc_unavailable_reason = vcf2tab_cnv.load_tc_lookup(input_path)
    if tc_unavailable_reason:
        logger.warning(
            f"{tc_unavailable_reason}. Only unadjusted CN will be calculated.")

    for case_folder in cnv_vcf_files:
        try:
            cnv_vcf = case_folder
            # The id vcf2tab_cnv writes is everything before the first ".":
            # "<SAMPLE_ID>.bam" from sample.tsv's own SAMPLE_ID when id_map
            # ties this file to one, else derived from the file name.
            mapped = (id_map or {}).get(Path(case_folder).name)
            sample_id = f"{mapped}.bam" if mapped else get_sample_id_from_cnv(case_folder)

            if sample_id in sid_path:
                with (Path(output_folder) / "sampleID_dup.log").open("a") as dup_path:
                    dup_path.write(sample_id + "\t" + "cnv_vcf\n")
            else:
                if cnv_dir is not None:
                    sid_path[sample_id] = Path(cnv_dir) / cnv_vcf
                elif multiple:
                    sid_path[sample_id] = Path(
                        input_path) / "CNV" / "single_sample_vcf" / cnv_vcf
                else:
                    sid_path[sample_id] = Path(input_path) / "CNV" / cnv_vcf

                vcf2tab_cnv.vcf_to_table(
                    sid_path[sample_id], Path(
                        output_folder) / "data_cna_hg19.seg",
                    sample_id, "w")
                vcf2tab_cnv.vcf_to_table_fc(
                    tc_lookup,
                    sid_path[sample_id], Path(
                        output_folder) / "data_cna_hg19.seg.fc.txt",
                    sample_id, "w")

        except Exception as err:
            logger.warning(f"Error while reading {case_folder}: {err}")
            with (Path(output_folder) / "noParsed_cnv.log").open("a") as log_noparsed:
                log_noparsed.write("[WARNING] " + case_folder + "\n")

        counter += 1

    if tc_unavailable_reason is None and sid_path:
        missing_tc = [
            sid.split(".")[0] for sid in sid_path
            if sid.split(".")[0] not in tc_lookup]
        if len(missing_tc) == len(sid_path):
            logger.warning(
                "None of the samples have a TC value. Only unadjusted CN "
                "will be calculated.")
        else:
            for sid in missing_tc:
                logger.warning(
                    f"Sample '{sid}' does not have a TC value. Only "
                    "unadjusted CN will be calculated.")

    seg_path = Path(output_folder) / "data_cna_hg19.seg"
    segfc_path = Path(output_folder) / "data_cna_hg19.seg.fc.txt"

    check_data_cna(seg_path)
    check_data_cna(segfc_path)

    if seg_path.exists():
        logger.info("Writing data_cna_hg19.seg successfully completed!")

    if segfc_path.exists():
        logger.info("Writing data_cna_hg19.seg.fc.txt successfully completed!")

    ############################
    ### MANAGE DISCRETE TABLE ##
    ############################

        logger.info("Starting CNA evaluation (this step could take a while)...")

        df_table = pd.read_csv(
            Path(output_folder) / "data_cna_hg19.seg.fc.txt", 
            sep="\t", 
            header=0,
            dtype={"ID": str}
        )

        df_table=df_table.rename(columns={
            "discrete":"Copy_Number_Alteration",
            "ID":"Tumor_Sample_Barcode",
            "gene":"Hugo_Symbol"})
        df_table_filt = df_table[
            df_table["Copy_Number_Alteration"].isin([-2,2])]

        cnv_kit = config.get("Cna", "CNVKIT_algorithm")
        cnv_kit = check_bool(cnv_kit)

        if cnv_kit:
            intermediate_dir = Path(output_folder) / "intermediate"
            intermediate_dir.mkdir(parents=True, exist_ok=True)

            if not Path(input_path).is_file():
                input_file = pd.read_csv(
                    Path(input_path) / "sample.tsv", sep="\t",
                    dtype={"SAMPLE_ID": str})
            else:
                input_file = pd.read_csv(input_path, sep="\t",
                                          dtype={"SAMPLE_ID": str})

            if "TC" not in input_file.columns:
                input_file["TC"] = np.nan

            if len(input_file[input_file["TC"].isna()])>0:
                nan_sbj = input_file[input_file["TC"].isna()]
                nan_sbj = list(nan_sbj["SAMPLE_ID"])
                logger.warning(
                    f"Some samples have a NaN TC in the input sample.tsv: {nan_sbj}!")

            if "ONCOTREE_CODE" not in input_file.columns:
                input_file["ONCOTREE_CODE"] = cancer

            input_file["Tumor_Sample_Barcode"] = input_file["SAMPLE_ID"]

            annotate = df_table_filt[[
                "Tumor_Sample_Barcode", "Hugo_Symbol",
                "FC", "Copy_Number_Alteration"]].merge(
                    input_file[["Tumor_Sample_Barcode",
                                "ONCOTREE_CODE", "TC"]],
                    on="Tumor_Sample_Barcode")

            # Numbered so the build order of these intermediate artifacts -
            # what got read to produce what - is obvious just from listing
            # the intermediate/ folder, without having to read the code.
            temppath = intermediate_dir / "00_temp_cna_toannotate.txt"
            annotate.to_csv(temppath, index=False, sep="\t")

            if oncokb:

                out_can_ann=intermediate_dir / "01_CNV_ann"
                out_can_ann.mkdir(exist_ok=True)

                oncokb_key = config.get("OncoKB", "ONCOKB")

                temppath_df=pd.read_csv(temppath, sep="\t")
                
                for sample_df in temppath_df["Tumor_Sample_Barcode"].unique():
                    
                    out = Path(out_can_ann) / temppath.name.replace(
                    "toannotate.txt", f"annotated_{sample_df}.txt")

                    df_tmp=temppath_df[temppath_df["Tumor_Sample_Barcode"]==sample_df]
                    df_path=intermediate_dir / "tmp_ann.txt"
                    df_tmp.to_csv(df_path, sep="\t", index=False)


                    cmd = [
                        sys.executable, str(VARAN_DIR / "oncokb-annotator" / "CnaAnnotator.py"),
                        "-i", str(df_path),
                        "-o", str(out),
                        "-f", "individual",
                        "-b", oncokb_key,
                        "-t", df_tmp["ONCOTREE_CODE"].unique()[0],
                        "-z"
                        ]
                    
                    subprocess.run(cmd, check=True)

                    df_path.unlink(missing_ok=True)

                all_ann_files = list(out_can_ann.glob("*annotated_*.txt"))
                if not all_ann_files:
                    logger.warning(
                        "No sample in this run had a Copy Number Alteration "
                        "call of |2| (amplification or deep deletion) to "
                        "annotate with OncoKB - skipping CNA annotation and "
                        "leaving data_cna.txt unwritten for this run, "
                        "consistently with the no-CNVKIT_algorithm path.")
                    return sid_path
                merged_df = pd.concat(
                    [pd.read_csv(f, sep="\t") for f in all_ann_files],
                    ignore_index=True
                )
                out = intermediate_dir / "02_CNV_annotated_merged.txt"
                merged_df.to_csv(out, sep="\t", index=False)

                name = "03_annotated_oncokb_CNA_ndiscrete.txt"
                cna = pd.read_csv(out, sep="\t",
                                  dtype={"Copy_Number_Alteration":int})
                if "o" in filters:
                    cna = filter_oncokb(cna, "Cna", "ONCOKB_FILTER_CNV")
            else:
                out = temppath
                name = "01_CNA_ndiscrete.txt"
                cna = pd.read_csv(out, sep="\t",
                                  dtype={"Copy_Number_Alteration":int})

            logger.info("Analyzing cna sample(s)")

            cna["Copy_Number_Alteration"] = 0

            if cna["TC"].isna().all():
                logger.warning("TC column is empty or contains only NaNs! "
                               "This column is required when CNVKIT_algorithm = True!")
                return sid_path

            for sample_id, sample_df in cna.groupby("Tumor_Sample_Barcode"):
                tc_val = sample_df["TC"].iloc[0]

                if pd.isna(tc_val):
                    logger.warning(f"Skipping sample {sample_id} due to NaN TC value.")
                    continue

                try:
                    tc = int(tc_val)
                except (ValueError, TypeError):
                    logger.warning(f"Skipping sample {sample_id} due to invalid TC value: {tc_val}")
                    continue

                purity = tc / 100
                copy_nums = np.arange(6)
                thresholds = 2 ** (np.log2((1 - purity) + purity * (copy_nums + .5) / ploidy()))

                sample_mask = cna["Tumor_Sample_Barcode"] == sample_id

                # Intentional: thresholds[2] (CN2/CN3 boundary) and
                # thresholds[4] (CN4/CN5 boundary) are computed but not used
                # as bucket edges, so a single extra/missing copy (CN3, CN5)
                # gets folded into the neighboring Neutral/Gain bucket
                # instead of being called Gain/Amplification on its own.
                # This is a deliberate conservative choice, not an oversight:
                # a single-copy gain is often within noise at typical
                # tumor purity/depth, and clinical CNV calling guidelines
                # commonly recommend stricter thresholds (e.g. requiring
                # CN>3 before calling amplification) specifically to avoid
                # over-calling marginal single-copy gains as clinically
                # meaningful. Do not "fix" this by using all 6 thresholds
                # without checking with the wet-lab/clinical team first.
                cna.loc[sample_mask & (cna["FC"] < thresholds[0]), "Copy_Number_Alteration"] = -2
                cna.loc[sample_mask & (cna["FC"] >= thresholds[0]) & (cna["FC"] < thresholds[1]), "Copy_Number_Alteration"] = -1
                cna.loc[sample_mask & (cna["FC"] >= thresholds[1]) & (cna["FC"] < thresholds[3]), "Copy_Number_Alteration"] = 0
                cna.loc[sample_mask & (cna["FC"] >= thresholds[3]) & (cna["FC"] < thresholds[5]), "Copy_Number_Alteration"] = 1
                cna.loc[sample_mask & (cna["FC"] >= thresholds[5]), "Copy_Number_Alteration"] = 2

            cna.to_csv(intermediate_dir / name,
                       index=False, sep="\t")

            cna["Tumor_Sample_Barcode"] = cna[
                "Tumor_Sample_Barcode"].str.replace(
                    ".cnv.bam", "", regex=False)

            # astype(int): cBioPortal only accepts discrete CNA values
            # written as integers ("2", not "2.0") - same as the branch below.
            data_cna = _one_call_per_gene(cna).pivot_table(
                index="Hugo_Symbol",
                columns="Tumor_Sample_Barcode",
                values="Copy_Number_Alteration", fill_value=0).astype(int)
            data_cna.to_csv(Path(output_folder) / "data_cna.txt",
                            index=True, sep="\t")

        else:
            df_table_filt = df_table_filt.copy()
            df_table_filt["Tumor_Sample_Barcode"] = df_table_filt["Tumor_Sample_Barcode"].astype(str)
            df_table_filt.loc[:, "Tumor_Sample_Barcode"] = df_table_filt[
                "Tumor_Sample_Barcode"].str.replace(
                    ".cnv.bam", "", regex=True)

            data_cna = _one_call_per_gene(df_table_filt).pivot_table(
                index="Hugo_Symbol",
                columns="Tumor_Sample_Barcode",
                values="Copy_Number_Alteration",
                fill_value=0).astype(int)
            if not data_cna.empty:
                data_cna.to_csv(Path(output_folder) / "data_cna.txt",
                                index=True, sep="\t")

    return sid_path


def _one_call_per_gene(calls: pd.DataFrame) -> pd.DataFrame:
    """Keep one discrete call per (sample, gene), the most extreme one.

    data_cna.txt is built from data_cna_hg19.seg.fc.txt, one row per
    segment: a gene split into several segments (a partial gain/loss, or a
    caller writing both a gene-level and smaller records), or a SAMPLE_ID
    listed twice in sample.tsv, gives the same gene more than one call in
    one sample. Averaging them (pivot_table's default) and truncating to an
    integer could turn -1 and -2 into -1; instead the call with the
    largest absolute value is kept, and every duplicate is logged -
    loudly when the calls disagree in direction (gain vs. loss).
    """
    keys = ["Tumor_Sample_Barcode", "Hugo_Symbol"]
    value = "Copy_Number_Alteration"
    duplicated = calls[calls.duplicated(keys, keep=False)]
    for (sample_id, gene), group in duplicated.groupby(keys):
        values = sorted(int(v) for v in group[value])
        if values[0] < 0 < values[-1]:
            logger.warning(
                f"{sample_id} {gene}: conflicting CNA calls {values} (gain and "
                "loss) from different segments - keeping the most extreme one. "
                "Check this gene's segments.")
        else:
            logger.warning(
                f"{sample_id} {gene}: {len(values)} CNA calls {values} from "
                "different segments - keeping the most extreme one.")
    order = calls[value].astype(int).abs().sort_values(ascending=False, kind="stable").index
    return calls.loc[order].drop_duplicates(keys, keep="first")


def extract_multiple_cnv(multiple_vcf: str, input_dir: str) -> None:
    """Split a multi-sample CNV VCF file into single-sample VCFs using vcftools.

    Args:
        multiple_vcf (str): Path to multi-sample CNV VCF file.
        input_dir (str): Working directory where single VCFs will be saved.

    Side Effects:
        - Creates `single_sample_vcf/` directory and per-sample VCF files inside it.
        - Generates and removes a temporary `sample_id.txt` file.

    """
    input_dir = Path(input_dir)
    single_sample_vcf_dir = input_dir / "single_sample_vcf"
    single_sample_vcf_dir.mkdir(parents=True, exist_ok=True)

    sample_id_txt = input_dir / "sample_id.txt"

    cmd_extract = ["vcf-query", "-l", str(multiple_vcf)]
    try:
        with sample_id_txt.open("w") as f:
            subprocess.run(cmd_extract, stdout=f, check=True)
    except subprocess.CalledProcessError as e:
        msg = f"Error during vcf-query: {e}"
        raise RuntimeError(msg) from e

    with sample_id_txt.open() as f:
        for line in f:
            sample_id = line.strip()
            output_vcf = single_sample_vcf_dir / f"{sample_id}.vcf"
            cmd_vcf = [
                "vcftools",
                "--vcf", str(multiple_vcf),
                "--indv", sample_id,
                "--recode",
                "--recode-INFO-all",
                "--stdout"]
            try:
                with output_vcf.open("w") as out_f:
                    subprocess.run(cmd_vcf, stdout=out_f, check=True)
            except subprocess.CalledProcessError as e:
                logger.warning(f"Error during VCF creation for {sample_id}: {e}")


def check_data_cna(data_cna_path: str) -> None:
    """Check if CNA data file is empty; remove file if empty.

    Args:
        data_cna_path (str): Path to CNA input file.

    Returns:
        None

    """
    input_file = Path(data_cna_path).name
    logger.info(f"Checking {input_file} file...")
    try:
        with Path(data_cna_path).open() as data_cna:
            all_data_cna = data_cna.readlines()
            if len(all_data_cna) == 1:
                Path(data_cna_path).unlink()
                logger.warning(f"{input_file} is empty. File removed.")
    except FileNotFoundError:
        logger.warning(f"{input_file} does not exist!")
    except OSError as e:
        logger.error(f"Error reading {input_file}: {e}")


def write_exon_brca(
    output_file: str, combined_dict: dict[str, str]
) -> dict[str, dict[str, dict[str, str]]]:
    """Write a tab-separated file containing BRCA exon-level CNV data.

    This function reads exon-level CNV information from a set of input files,
    filters and formats the data, and writes it into a single output file
    named `exonic_BRCA.txt`. Each row in the output corresponds to a CNV call
    from a given sample related to BRCA genes.

    Args:
        output_file (str): Path to the output file where results will be written.
        combined_dict (dict[str, str]): Dictionary mapping sample IDs to input file
        paths. Each file is expected to contain a section with exon-level CNV data.

    Returns:
        dict[str, dict[str, dict[str, str]]]: The same events, keyed by sample
            then gene ({sample_id: {gene: {"chromosome", "start", "stop",
            "exon", "fold_change", "cnv_type"}}}), for callers that need the
            structured data (e.g. write_exon_brca_generic_assay) without
            re-parsing the source CombinedVariantOutput files.

    """
    logger.info("Checking for exonic BRCA data...")
    rows_to_write = []
    events_by_sample: dict[str, dict[str, dict[str, str]]] = {}

    header = (
        "Sample_Id\tGene\tChromosome\tStart\tStop\t"
        "Affected Exon(s)\tFold Change\tCNV Type\n")

    for k, v in combined_dict.items():
        exonic = []
        try:
            exonic = tsv.get_exons(Path(v))
        except Exception as err:
            logger.error(f"Something went wrong while reading Exon-Level CNVs for sample {k}: {err}")

        if not exonic:
            continue

        for ex in exonic:
            hugo_symbol = ex.get("Hugo_Symbol", "NA")
            chr_val = ex.get("Chromosome", "NA")
            start = ex.get("Start_Position", "NA")
            stop = ex.get("Stop_Position", "NA")
            exon = ex.get("Affected_Exon(s)", "NA")
            fc = ex.get("Fold_Change", "NA")
            cnv_type = ex.get("CNV_Type", "NA")

            line = f"{k}\t{hugo_symbol}\t{chr_val}\t{start}\t{stop}\t{exon}\t{fc}\t{cnv_type}\n"
            rows_to_write.append(line)

            events_by_sample.setdefault(k, {})[hugo_symbol] = {
                "chromosome": chr_val,
                "start": start,
                "stop": stop,
                "exon": exon,
                "fold_change": fc,
                "cnv_type": cnv_type,
            }

    if not rows_to_write:
        logger.warning(
            "No BRCA exon-level CNV data found across all samples. "
            "The exonic_BRCA.txt file will not be created."
        )
        return events_by_sample

    try:
        output_file_path = Path(output_file)
        with output_file_path.open("w") as exonic_table:
            exonic_table.write(header)
            exonic_table.writelines(rows_to_write)
        logger.info(f"Successfully created {output_file_path.name} with {len(rows_to_write)} records.")
    except Exception as e:
        logger.error(f"Failed to write output file: {e}")

    return events_by_sample


def write_exon_brca_generic_assay(
    output_folder: str,
    events_by_sample: dict[str, dict[str, dict[str, str]]],
) -> None:
    """Write a cBioPortal Generic Assay data file for BRCA exon-level CNVs.

    Derived from the same per-sample/per-gene events write_exon_brca already
    collected, so this file and exon_CNA_data.txt can't drift apart - one
    parse of the source CombinedVariantOutput data, two views of it:
    exon_CNA_data.txt keeps full per-event detail (coordinates, fold change,
    affected exon) for audit purposes; this file is the coarser
    LOSS/GAIN/NEUTRAL-per-gene matrix the Generic Assay format requires (one
    row per entity, one column per sample), so the calls show up as their
    own oncoprint/study-view track in cBioPortal instead of being invisible
    outside the two clinical-attribute summary columns. It's a label/value
    matrix with no genomic-coordinate field, so unlike a real CNA or SV
    profile it won't appear in position-based views (Genome View, Mutation
    Mapper) - only in the oncoprint track, study view and comparison plots.

    Kept as a separate profile/track (not merged into data_cna.txt) so it
    can never be confused with the whole-gene copy-number-alteration values
    already reported there - BRCA1/BRCA2 keep their normal gene-level CNA
    row, and this is an independent, additional track.

    Args:
        output_folder (str): Output directory (data_exon_brca_cna.txt is
            written there).
        events_by_sample (dict): {sample_id: {gene: {"cnv_type": "LOSS" |
            "GAIN", ...}}}, as returned by write_exon_brca.

    Returns:
        None

    """
    if not events_by_sample:
        return

    genes = ["BRCA1", "BRCA2"]
    samples = sorted(events_by_sample)

    lines = ["ENTITY_STABLE_ID\tNAME\t" + "\t".join(samples) + "\n"]
    for gene in genes:
        values = []
        for sample in samples:
            cnv_type = events_by_sample.get(sample, {}).get(gene, {}).get("cnv_type")
            values.append(cnv_type if cnv_type in {"LOSS", "GAIN"} else "NEUTRAL")
        lines.append(f"{gene}\t{gene} (exon-level)\t" + "\t".join(values) + "\n")

    out_path = Path(output_folder) / "data_exon_brca_cna.txt"
    try:
        with out_path.open("w") as f:
            f.writelines(lines)
        logger.info(f"Successfully created {out_path.name} for {len(samples)} sample(s).")
    except Exception as e:
        logger.error(f"Failed to write {out_path.name}: {e}")


def build_exon_df(output_folder: Path, exonic, sample_id: str) -> pd.DataFrame:
    """Build a df with BRCA1 and BRCA2 exon-level CNV info for a given sample.

    This function evaluates exon-level CNV data for BRCA1 and BRCA2 genes and classifies their status
    as "Positive" (if exons are reported) or "Negative" (if gene is present but no events). If the
    exon-level CNV section is missing entirely (i.e., `exonic` is None), all fields are set to "NA".

    Args:
        output_folder (Path): Path to the output directory (not used in this function).
        exonic (list[dict] or None): List of dictionaries representing exon-level CNVs,
            each dictionary must contain at least a "Hugo_Symbol" key. If None, the section
            was missing entirely from the input.
        sample_id (str): The sample identifier.

    Returns:
        pd.DataFrame: A DataFrame with one row and the following columns:
            - SAMPLE_ID
            - BRCA1: "Positive", "Negative", or "NaN"
            - BRCA2: "Positive", "Negative", or "NaN"
            - BRCA1_DETAILS: formatted string or "NaN"
            - BRCA2_DETAILS: formatted string or "NaN"

    """
    if exonic is None:
        row = {
            "SAMPLE_ID": sample_id,
            "EXONIC_BRCA1": "NaN",
            "EXONIC_BRCA2": "NaN",
            "BRCA1_DETAILS": "NaN",
            "BRCA2_DETAILS": "NaN"
        }
    else:
        brca1_exons = [entry for entry in exonic if entry.get("Hugo_Symbol") == "BRCA1"]
        brca2_exons = [entry for entry in exonic if entry.get("Hugo_Symbol") == "BRCA2"]

        if any(exon.get("CNV_Type") == "LOSS" for exon in brca1_exons):
            brca1_status = "LOSS"
        elif any(exon.get("CNV_Type") == "GAIN" for exon in brca1_exons):
            brca1_status = "GAIN"
        else:
            brca1_status = "NaN"

        if any(exon.get("CNV_Type") == "LOSS" for exon in brca2_exons):
            brca2_status = "LOSS"
        elif any(exon.get("CNV_Type") == "GAIN" for exon in brca2_exons):
            brca2_status = "GAIN"
        else:
            brca2_status = "NaN"

        brca1_details = "; ".join(str(e) for e in brca1_exons) if brca1_exons else "NaN"
        brca2_details = "; ".join(str(e) for e in brca2_exons) if brca2_exons else "NaN"

        row = {
            "SAMPLE_ID": sample_id,
            "EXONIC_BRCA1": brca1_status,
            "EXONIC_BRCA2": brca2_status,
            "BRCA1_DETAILS": brca1_details,
            "BRCA2_DETAILS": brca2_details
        }

    return pd.DataFrame([row])


def _walk_process_cnv(ctx: WalkContext) -> None:
    """Stage: CNV calls (data_cna*.txt) and BRCA exon-level CNV table.

    Depends only on WalkContext, independent of _walk_process_snv/_fusion.
    """
    if ctx.input_folder_cnv.exists() and ctx.vcf_type not in ["snv", "fus", "tab"]:
        logger.info("Managing CNV files...")
        cnv_type_from_folder(
            ctx.input_folder, ctx.case_folder_arr_cnv, ctx.output_folder,
            ctx.oncokb, ctx.cancer, ctx.multiple, ctx.filters, ctx.cnv_sample_ids,
            ctx.input_folder_cnv)

    combined_output_folder = Path(ctx.input_folder) / "CombinedOutput"
    if (combined_output_folder.exists()
        and any(f.is_file() for f in combined_output_folder.iterdir())):
        combined_dict = get_combined_variant_output_from_folder(
            ctx.input_folder, ctx.clin_file, ctx.isinputfile)

        exon_file_output = Path(ctx.output_folder) / "exon_CNA_data.txt"
        exon_events = write_exon_brca(exon_file_output, combined_dict)
        write_exon_brca_generic_assay(ctx.output_folder, exon_events)
