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

"""Module to convert VCF in MAF and to evaluate CNA.

This script supports:
- conversion and vep annotation
- cna extraction and filtering

"""
from __future__ import annotations

import ast
import contextlib
import os
import secrets
import shutil
import string
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

import sigma_runner
import tsv
import vcf2tab_cnv
import vcf_filter
from config_loader import get_config
from filter_clinvar import check_bool, filter_oncokb
from versioning import get_newest_version, get_version_list

config = get_config()

VCF2MAF = config.get("Paths", "VCF2MAF")
REF_FASTA = config.get("Paths", "REF_FASTA")
VEP_PATH = config.get("Paths", "VEP_PATH")
VEP_DATA = config.get("Paths", "VEP_DATA")
CLINV = config.get("Paths", "CLINV")
PLOIDY = int(config.get("Cna", "PLOIDY"))
SAMPLE_TYPE = (config.get("Sample_Type", "TYPE").strip().strip('"').strip("'").upper())
THRESHOLD_MSI_LIQUID = float(config.get("MSI", "THRESHOLD_MSI_LIQUID"))

output_filtered = "snv_filtered"


def create_random_name_folder(output_folder: str) -> str:
    """Create a temporary scratch folder with a random name inside output_folder.

    Scratch lives at `output_folder/scratch/<random>` instead of a shared,
    cwd-relative `scratch/` folder, so each run's temp files are namespaced
    under its own output folder - no more collision risk when several
    varan.py processes run in parallel from the same working directory.

    Args:
        output_folder (str): The study's output folder for this run.

    Returns:
        str: Path to the created temporary folder.

    """
    folder_name = "".join(
        secrets.choice(string.ascii_lowercase + string.digits) for _ in range(10))
    temporary = Path(output_folder) / "scratch" / folder_name

    try:
        temporary.mkdir(parents=True)
    except Exception as err:
        logger.critical("Something went wrong while creating the vep tmp folder")
        msg = "Error in create_random_name_folder: exiting from walk script!"
        raise(Exception(msg)) from err
    return(str(temporary))


def clear_scratch(folder: str | None = None) -> None:
    """Remove this run's scratch subfolder, then the (now-empty) scratch/
    parent directory itself, so nothing scratch-related is left in the
    study's output folder.

    Only ever removes `<output_folder>/scratch/<this run's random name>`
    and then `<output_folder>/scratch` itself - never scans/clears a shared,
    cwd-relative `scratch/` directory, since that would risk deleting other
    runs' still-in-use temporary files when multiple varan.py processes run
    in parallel from the same working directory. Removing the parent here is
    still safe: it's namespaced under this run's own output_folder, and this
    is the only place a scratch subfolder is ever created for that folder.

    Args:
        folder (str | None): Path to this run's own scratch subfolder, or None
            if this run never created one (e.g. resume=True).

    Returns:
        None

    """
    if folder is None:
        return
    to_rem = Path(folder)
    if to_rem.exists():
        shutil.rmtree(to_rem)
    scratch_parent = to_rem.parent
    if scratch_parent.name == "scratch" and scratch_parent.exists():
        shutil.rmtree(scratch_parent, ignore_errors=True)
        # create_random_name_folder()'s mkdir(parents=True) may have had to
        # create the study's output folder itself just to hold this
        # scratch/ subfolder - e.g. vendor-adapter preprocessing
        # (run_vendor_adapter() in varan.py) runs against the not-yet-
        # versioned output folder name, before create_output_folder() ever
        # creates the real "<name>_v1" folder. Removing scratch/ above can
        # then leave that now-empty output folder behind. Only ever remove
        # it when it's completely empty, so a real output folder that
        # already has pipeline output in it (e.g. walk.py's own mid-run
        # scratch cleanup) is never touched.
        output_folder = scratch_parent.parent
        if output_folder.exists() and not any(output_folder.iterdir()):
            output_folder.rmdir()


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
                         filters: str = "") -> dict:
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

    Returns:
    dict
        A mapping from sample ID to VCF file path used in processing.

    """
    counter = 0
    sid_path = {}

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
            sample_id = get_sample_id_from_cnv(case_folder)

            if sample_id in sid_path:
                with (Path(output_folder) / "sampleID_dup.log").open("a") as dup_path:
                    dup_path.write(sample_id + "\t" + "cnv_vcf\n")
            else:
                if multiple:
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

        except Exception:
            logger.warning(f"Error while reading {case_folder}")
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
        logger.info("Writing data_cna_hg19.seg succefully completed!")

    if segfc_path.exists():
        logger.info("Writing data_cna_hg19.seg.fc.txt succefully completed!")

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
                    f"Some subject have NaN TC in tsv input file: {nan_sbj}!")

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
                
                df_all=pd.DataFrame()
                for sample_df in temppath_df["Tumor_Sample_Barcode"].unique():
                    
                    out = Path(out_can_ann) / temppath.name.replace(
                    "toannotate.txt", f"annotated_{sample_df}.txt")

                    df_tmp=temppath_df[temppath_df["Tumor_Sample_Barcode"]==sample_df]
                    df_path=intermediate_dir / "tmp_ann.txt"
                    df_tmp.to_csv(df_path, sep="\t", index=False)


                    cmd = [
                        sys.executable, "./oncokb-annotator/CnaAnnotator.py",
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
                thresholds = 2 ** (np.log2((1 - purity) + purity * (copy_nums + .5) / PLOIDY))

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

            data_cna = cna.pivot_table(
                index="Hugo_Symbol",
                columns="Tumor_Sample_Barcode",
                values="Copy_Number_Alteration", fill_value=0)
            data_cna.to_csv(Path(output_folder) / "data_cna.txt",
                            index=True, sep="\t")

        else:
            df_table_filt = df_table_filt.copy()
            df_table_filt["Tumor_Sample_Barcode"] = df_table_filt["Tumor_Sample_Barcode"].astype(str)
            df_table_filt.loc[:, "Tumor_Sample_Barcode"] = df_table_filt[
                "Tumor_Sample_Barcode"].str.replace(
                    ".cnv.bam", "", regex=True)

            data_cna = df_table_filt.pivot_table(
                index="Hugo_Symbol",
                columns="Tumor_Sample_Barcode",
                values="Copy_Number_Alteration",
                fill_value=0).astype(int)
            if not data_cna.empty:
                data_cna.to_csv(Path(output_folder) / "data_cna.txt",
                                index=True, sep="\t")

    return sid_path


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
                         output_folder: str) -> dict:
    """Map sample IDs to their full SNV paths, handling duplicates.

    Args:
        input_pat (str): Path to the folder containing SNV files.
        snv_vcf_files (list): List of SNV VCF filenames.
        output_folder (str): The study's output folder, where sampleID_dup.log
            and noParsed_snv.log are written (read back into the report's
            Warnings section).

    Returns:
        dict: Mapping from sample ID (as BAM filename) to full SNV file path.

    """
    c = 0
    sid_path = {}
    for case_folder in snv_vcf_files:
        try:
            snv_vcf = case_folder
            sample_id = get_sample_id_from_snv(case_folder)
            if sample_id in sid_path:
                with (Path(output_folder) / "sampleID_dup.log").open("a") as dup:
                    dup.write(sample_id + "\t" + "snv_vcf\n")
            else:
                sid_path[sample_id] = str(Path(input_pat) / snv_vcf)
        except Exception:
            with (Path(output_folder) / "noParsed_snv.log").open("a") as log_noparsed:
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
                        output_folder: str) -> list:
    """Construct the vcf2maf command line based on input parameters and config file.

    Args:
        v (str): Path to the input VCF file.
        temporary (str): Temporary directory for intermediate files.
        output_folder (str): Base output folder.

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

    if VCF2MAF == "" or REF_FASTA == "" or VEP_PATH == "" or VEP_DATA == "":
        logger.critical("[Paths] section in conf.ini is not correctly compiled. "
        "Please check again!")
        msg = "Input error"
        raise Exception (msg)

    cl = ["perl"]
    cl.append(VCF2MAF)
    cl.append("--input-vcf")
    cl.append(v)
    _, file_vcf = os.path.split(v)

    out_file = Path(output_folder) / "maf" / (file_vcf + ".maf")
    if CLINV.strip() != "":
        cl.append("--vep-custom")
        cl.append(CLINV)
    else:
        logger.warning("CLINV section in [Paths] in conf.ini is not compiled. "
        "This step will be skipped")
    cl.append("--output-maf")
    cl.append(str(out_file))
    cl.append("--ref-fasta")
    cl.append(REF_FASTA)
    cl.append("--tmp-dir")
    cl.append(str(temporary))
    cl.append("--retain-fmt")
    cl.append("GT,GQ,AD,DP,VF,AF")
    cl.append("--vep-path")
    cl.append(VEP_PATH)
    cl.append("--vep-data")
    cl.append(VEP_DATA)
    cl.append("--tumor-id")
    cl.append(tum_id)
    cl.append("--cache-version")
    cl.append(cache)
    return cl


def run_vcf2maf(cl: list, sample: str) -> None:
    """Run the vcf2maf conversion using the constructed command.

    Args:
        cl (list): The command-line list to execute.
        sample (str): The sample ID for logging purposes.

    Returns:
        None

    """
    logger.info(f"Starting vcf2maf conversion of {sample}... This may take several "
    "minutes (approx 5 minutes per sample).")
    logger.info(f"args={cl}")
    sout = subprocess.run(cl, capture_output=True, check=False)

    if sout.stderr is not None:
        if "ERROR" not in sout.stderr.decode("ascii"):
            logger.warning(sout.stderr.decode("ascii").replace("ERROR: ",""))
        else:
            logger.error(sout.stderr.decode("ascii").replace("ERROR: ",""))

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


def create_folder(output_folder: str,
                  overwrite_output: bool, resume: bool) -> Path:
    """Create versioned output folder for storing results.

    Args:
        output_folder (str): Base output path.
        overwrite_output (bool): If True, delete existing version and create new.
        resume (bool): If True, resume from last existing version if available.

    Returns:
        Path: Path to the versioned output folder.

    """
    output_list = [
        Path(output_folder).parent / x
        for x in get_version_list(output_folder)
    ]

    if output_list and output_list[-1].exists():
        output = output_list[-1]
        logger.warning(
            f"It seems that a version of the folder '{output_folder}' already exists.")
        if overwrite_output:
            logger.info("Overwrite option set. Start removing folder")
            shutil.rmtree(output)
        elif resume:
            _,current_version = get_newest_version(output_folder)
            return Path(output_folder + current_version)

    if not output_list:
        version = "_v1"
        output_folder_version = Path(f"{output_folder}{version}")


    else:
        output_folder_version, _ = get_newest_version(output_folder)

    logger.info(f"Creating the output folder '{output_folder_version}'...")

    output_folder_version.mkdir(parents=True, exist_ok=True)

    maf_path = output_folder_version / "maf"
    maf_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"The folder '{output_folder_version}' was correctly created!")

    return output_folder_version


def get_table_from_folder(tsvpath: str) -> dict[str, list[str]]:
    """Extract sample-to-patient mapping from a TSV file.

    Args:
        tsvpath (str): Path to the input TSV file.

    Returns:
        dict[str, list[str]]: Mapping from SAMPLE_ID to list with
        corresponding PATIENT_ID.

    """
    table_dict = {}
    file = pd.read_csv(tsvpath, sep="\t", index_col=False, dtype=str)
    for _, row in file.iterrows():
        sample_id = str(row["SAMPLE_ID"])
        if ".bam" in sample_id:
           sample_id = sample_id.replace(".bam", "")
        if sample_id not in table_dict:
            table_dict[sample_id] = [str(row["PATIENT_ID"])]
    return table_dict


def flatten(nested_list: list[list]) -> list:
    """Flatten a nested list into a single flat list.

    Args:
        nested_list (list[list]): A list of lists to flatten.

    Returns:
        list: A single flattened list.

    """
    return [item for sublist in nested_list for item in sublist]


def check_multiple_file(input_file: str, multiple: bool) -> None:
    """Check sample.tsv and conf.ini for multiple sample input.

    Args:
        input_file (str): Path to the sample.tsv file.
        multiple (bool): Whether the input is in multiple-sample mode.

    Raises:
        Exception: If validation of configuration and sample.tsv fails.

    """
    conf_snv = config.get("Multiple", "SNV")
    conf_cnv = config.get("Multiple", "CNV")
    file_paths = pd.read_csv(input_file, sep="\t", header=0)
    snv_file_path = file_paths["snv_path"].isna().all()
    cnv_file_path = file_paths["cnv_path"].isna().all()

    #CASE 1: multiple = True; neither snv or cnv are filled in sample.tsv;
    #multiple snv or cnv are filled in conf.ini
    if (multiple
    and not (snv_file_path or cnv_file_path)
    and not (conf_snv == "" or conf_cnv == "")):
        logger.critical(
            "-m was selected and Muliple section in conf.ini was filled but the "
            "file doesn't looks like a multiVCF.")
        msg = "Input error"
        raise Exception(msg)
    #CASE 2: multiple = True; neither snv or cnv are filled in sample.tsv;
    #neither multiple snv or cnv are filled in conf.ini
    if (multiple
    and not (snv_file_path or cnv_file_path)
    and (conf_snv == "" or conf_cnv == "")):
        logger.critical("-m was selected but Muliple section in conf.ini wasn't "
        "filled and the file doesn't looks like a multiVCF.")
        msg = "Input error"
        raise Exception (msg)
    #CASE 3: multiple = False; snv or cnv are filled in sample.tsv;
    #multiple snv or cnv are filled in conf.ini
    if (not multiple
    and (snv_file_path or cnv_file_path)
    and not (conf_snv == "" or conf_cnv == "")):
        logger.critical("-m was not selected but both sample.tsv and Muliple "
        "section in conf.ini were filled.")
        msg = "Input error"
        raise Exception (msg)
    #CASE 4: multiple = False; snv or cnv are filled in sample.tsv;
    #neither multiple snv or cnv are filled in conf.ini
    if (not multiple
    and (snv_file_path or cnv_file_path)
    and (conf_snv == "" or conf_cnv == "")):
        logger.warning("SNV and/or CNV columns in sample.tsv were not filled.")


def check_multiple_folder(input_dir: str, multiple: bool) -> None:
    """Validate whether SNV/CNV VCFs are single or multi-sample based on -m option.

    Args:
        input_dir (str): Path to the base input directory.
        multiple (bool): Whether the input should be multi-sample.

    Raises:
        Exception: If actual file type does not match `multiple` expectation.

    """
    snv_mulitple, snv_single, cnv_mulitple, cnv_single = False, False, False, False
    snv_folder = Path(input_dir) / "SNV"
    try:
        vcf_files = (f.name for f in snv_folder.iterdir() if f.suffix == ".vcf")
        multiple_vcf_snv = next(vcf_files)
    except StopIteration:
        multiple_vcf_snv = ""
        snv_folder.mkdir(parents=True, exist_ok=True)

    snv_file = Path(input_dir) / "sample_id_snv.txt"
    cmd_snv = f"vcf-query -l {(snv_folder / multiple_vcf_snv)} > {snv_file}"
    os.system(cmd_snv)

    with snv_file.open() as file:
        lines = file.readlines()

    minimum_lines_for_multiple = 2

    if len(lines) >= minimum_lines_for_multiple and not multiple:
        snv_mulitple = True
        logger.error("-m option was not selected but the SNV file is multiple!")
    elif len(lines) < minimum_lines_for_multiple and multiple:
        snv_single = True
        logger.error("-m option was selected but the SNV file is not multiple!")

    snv_file.unlink()

    cnv_folder = Path(input_dir) / "CNV"
    try:
        vcf_files = (f.name for f in cnv_folder.iterdir() if f.suffix == ".vcf")
        multiple_vcf_cnv = next(vcf_files)
    except StopIteration:
        multiple_vcf_cnv = ""
        cnv_folder.mkdir(parents=True, exist_ok=True)

    cnv_file = Path(input_dir) / "sample_id_cnv.txt"
    cmd_cnv = f"vcf-query -l {(Path(cnv_folder) / multiple_vcf_cnv)} > {cnv_file}"
    os.system(cmd_cnv)

    with cnv_file.open() as file:
        lines = file.readlines()

    if len(lines) >= minimum_lines_for_multiple and not multiple:
        cnv_mulitple = True
        logger.error("-m option was not selected but the CNV file is multiple!")
    elif len(lines) < minimum_lines_for_multiple and multiple:
        cnv_single = True
        logger.error("-m option was selected but the CNV file is not multiple!")

    cnv_file.unlink()

    msg = "Input error"
    if cnv_mulitple or snv_mulitple:
        raise Exception(msg)
    if cnv_single or snv_single:
        raise Exception(msg)


def write_clinical_sample(
    clin_samp_path: str,
    output_folder: str,
    table_dict: dict,
    combined_output_used: bool = False,
) -> None:
    """Write the `data_clinical_sample.txt` file by merging sample and metrics data.

    The function reads the clinical sample file, optionally merges it with MSI/TMB
    metrics from a provided table, validates headers, and writes the output file
    with appropriate headers as required by cBioPortal format.

    Args:
        clin_samp_path (str): Path to the input clinical sample TSV file.
        output_folder (str): Directory where the final file will be saved.
        table_dict (dict): Dictionary containing metrics (MSI, TMB) by sample ID.
        combined_output_used (bool): True if `table_dict` was built from a real
            CombinedOutput folder (fill_from_combined); False if it was built as a
            fallback straight from sample.tsv (fill_from_file). Used to log which
            source MSI/TMB actually came from, and to avoid comparing sample.tsv
            against itself when no CombinedOutput was ever provided.

    Raises:
        NameError: If expected columns or header definitions are missing or incorrect.
        KeyError: If input template does not contain required fields.

    """
    logger.info("Writing data_clinical_sample.txt file...")
    conf_header_short = config.get("ClinicalSample", "HEADER_SAMPLE_SHORT")
    conf_header_long = config.get("ClinicalSample", "HEADER_SAMPLE_LONG")
    conf_header_type = config.get("ClinicalSample", "HEADER_SAMPLE_TYPE")

    data_clin_samp = pd.read_csv(clin_samp_path, sep="\t", header=0, dtype=str)
    data_clin_samp["ONCOTREE_CODE"] = data_clin_samp["ONCOTREE_CODE"].str.upper()

    try:
        data_clin_samp = data_clin_samp.drop(
            columns=["snv_path", "cnv_path", "comb_path"])
    except KeyError as err:
        logger.critical("snv_path, cnv_path or comb_path columns were removed or "
        "modified from template. Please use the correct template!")
        msg = "Exiting from script!"
        raise NameError(msg) from err

    data_clin_samp.columns = data_clin_samp.columns.str.upper()

    combout_df = pd.DataFrame.from_dict(table_dict).transpose().reset_index()
    combout_df = combout_df.rename(columns={"index": "SAMPLE_ID", 0: "PATIENT_ID"})

    final_data_sample = data_clin_samp
    limit = 2

    if len(combout_df.columns) > limit:
        combout_df = combout_df.rename(columns={
            1: "MSI",
            2: "TMB",
            3: "MSI_THR",
            4:"TMB_THR"})
        source = "CombinedOutput" if combined_output_used else "sample.tsv"
        logger.info(
            f"MSI/TMB values written to data_clinical_sample.txt were taken from {source}.")

        try:
            if combined_output_used:
                msi_sample = pd.to_numeric(data_clin_samp["MSI"], errors="coerce")
                msi_combined = pd.to_numeric(combout_df["MSI"], errors="coerce")
                tmb_sample = pd.to_numeric(data_clin_samp["TMB"], errors="coerce")
                tmb_combined = pd.to_numeric(combout_df["TMB"], errors="coerce")

                msi_mismatch = ((msi_sample != msi_combined)
                                 & msi_sample.notna() & msi_combined.notna()).any()
                tmb_mismatch = ((tmb_sample != tmb_combined)
                                 & tmb_sample.notna() & tmb_combined.notna()).any()

                if msi_mismatch or tmb_mismatch:
                    logger.warning(
                        "MSI and/or TMB values reported in sample.tsv differ from "
                        "the ones computed from CombinedOutput! CombinedOutput "
                        "values were used, sample.tsv values were discarded.")
            try:
                data_clin_samp = data_clin_samp.drop(
                    columns=["MSI", "TMB", "MSI_THR", "TMB_THR"])
            except KeyError as err:
                logger.critical(
                    "MSI_THR or TMB_THR columns were removed or modified from template."
                    " Please use the correct template!")
                msg = "Exiting from script!"
                raise NameError(msg) from err
        except KeyError as err:
            logger.warning("No MSI or TMB columns found in template.")
            msg = "Exiting from script!"
            raise(KeyError(msg)) from err

        final_data_sample = data_clin_samp.merge(
            combout_df, on=["PATIENT_ID", "SAMPLE_ID"])

    # SigMA columns only ever appear when -g/--sigma was actually passed for
    # this run: sigma_runner/_walk_process_snv only ever writes this file
    # when ctx.sigma is True and at least one sample went through the SigMA
    # pipeline (see _write_sigma_intermediate). No file -> no columns added,
    # matching this codebase's existing convention of only showing/recording
    # a filter's info when that filter was actually applied (see
    # ONCOKB_FILTER_CNV/ONCOKB_FILTER_FUSION in write_report.py).
    sigma_path = Path(output_folder) / "intermediate" / "sigma" / "data_sigma.txt"
    if sigma_path.exists():
        sigma_df = pd.read_csv(sigma_path, sep="\t", dtype=str)
        if "SAMPLE_ID" in sigma_df.columns:
            final_data_sample = final_data_sample.merge(
                sigma_df, on="SAMPLE_ID", how="left")
            logger.info(
                f"Merged SigMA results for {len(sigma_df)} sample(s) into "
                "data_clinical_sample.txt.")
        else:
            logger.warning(
                f"{sigma_path} exists but has no SAMPLE_ID column - "
                "skipping SigMA column merge into data_clinical_sample.txt.")

    basic_columns = ["SAMPLE_ID", "PATIENT_ID", "MSI", "TMB", "MSI_THR", "TMB_THR"]

    # A pipeline that has no native TMB (or MSI) of its own - e.g. Guardant,
    # which reports MSI as a text call but never reports TMB at all - must
    # not be forced to ship a clinical column full of "NA" placeholders.
    # Only applies to the sample.tsv-driven path (no CombinedOutput source
    # to fall back on): native TSO500 CombinedOutput ingestion always has
    # real MSI/TMB values.
    if not combined_output_used:
        for marker in ("MSI", "TMB"):
            thr_col = f"{marker}_THR"
            if marker in final_data_sample.columns and (
                final_data_sample[marker].apply(_is_blank_value).all()
                and final_data_sample[thr_col].apply(_is_blank_value).all()):
                logger.info(
                    f"No {marker} value or {marker}_THR was provided for any "
                    f"sample in this run, and there is no CombinedOutput to "
                    f"compute it from - dropping {marker}/{thr_col} from "
                    "data_clinical_sample.txt entirely instead of writing NA.")
                basic_columns.remove(marker)
                basic_columns.remove(thr_col)
                final_data_sample = final_data_sample.drop(columns=[marker, thr_col])

    other_columns = [c for c in final_data_sample.columns
                      if c not in ("SAMPLE_ID", "PATIENT_ID", "MSI", "TMB",
                                   "MSI_THR", "TMB_THR")]
    new_cols = [*basic_columns, *other_columns]

    final_data_sample = final_data_sample[new_cols]
    dataclin_columns = list(final_data_sample.columns)

    # Add header's fifth row
    default_row = pd.DataFrame([dataclin_columns], columns=dataclin_columns)
    final_data_sample = pd.concat([default_row, final_data_sample], ignore_index=True)

    # Add header's fourth row (SERIES OF 1s)
    header_numbers = pd.DataFrame(
        [[1] * len(final_data_sample.columns)], columns=dataclin_columns)
    final_data_sample = pd.concat(
        [header_numbers, final_data_sample], ignore_index=True)

    # Add header's third row (HEADER_SAMPLE_TYPE)
    if not conf_header_type:
        header_row = ["STRING"] * len(dataclin_columns)
        for marker in ("MSI", "TMB", "SIGMA_TOTAL_SNVS", "SIGMA_SIGNATURE3_MVA"):
            if marker in dataclin_columns:
                header_row[dataclin_columns.index(marker)] = "NUMBER"
        for marker in ("SIGMA_DO_MVA",):
            if marker in dataclin_columns:
                header_row[dataclin_columns.index(marker)] = "BOOLEAN"
        sample_header_type = pd.DataFrame([header_row], columns=dataclin_columns)
    else:
        types_list = conf_header_type.split(",")
        types_list = [x.strip() for x in types_list]
        for types in types_list:
            if types.upper() not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{types} is not a valid type. Please check the given "
                "input in conf.ini. Valid types: STRING, NUMBER, BOOLEAN")
                msg = "The type is not valid: exiting from walk script!"
                raise(NameError(msg))
        try:
            types_list = [x.upper() for x in types_list]
            sample_header_type = pd.DataFrame([types_list], columns=dataclin_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(types_list)}) in "
            "HEADER_SAMPLE_TYPE is different from the effective number of columns "
            f"({len(final_data_sample.columns)}).")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    final_data_sample = pd.concat(
        [sample_header_type, final_data_sample], ignore_index=True)

    # Add header's second row (HEADER_SAMPLE_LONG)
    if not conf_header_long:
        sample_header_long = default_row
    else:
        try:
            combined_headers = conf_header_long.split(",")
            sample_header_long = pd.DataFrame(
                [combined_headers], columns=dataclin_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(combined_headers)}) in "
            "HEADER_SAMPLE_LONG in conf.ini is different from the effective number of "
            f"columns ({len(final_data_sample.columns)}).")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err
    final_data_sample = pd.concat(
        [sample_header_long, final_data_sample], ignore_index=True)

    # Add header's first row (HEADER_SAMPLE_SHORT)
    if not conf_header_short:
        sample_header_short = default_row
    else:
        try:
            combined_headers = conf_header_short.split(",")
            sample_header_short = pd.DataFrame(
                [combined_headers], columns=dataclin_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(combined_headers)}) in "
            "HEADER_SAMPLE_SHORT in conf.ini is different from the effective number of "
            f"columns ({len(final_data_sample.columns)}).")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err
    final_data_sample = pd.concat(
        [sample_header_short, final_data_sample], ignore_index=True)
    final_data_sample.loc[0:3, "SAMPLE_ID"] = (
        final_data_sample.loc[0:3, "SAMPLE_ID"].apply(lambda x: f"#{x}"))

    output_folder = Path(output_folder)
    data_clin_txt = output_folder / "data_clinical_sample.txt"
    final_data_sample.to_csv(data_clin_txt, sep="\t", index=False, header=False)


def write_default_clinical_patient(output_folder: str, table_dict: dict) -> None:
    """Write a default `data_clinical_patient.txt` file.

    Args:
        output_folder (str): Path to the output directory.
        table_dict (dict): Dictionary containing sample info mapped to patient IDs.

    """
    logger.info("Writing data_clinical_patient.txt file...")

    output_folder = Path(output_folder)
    data_clin_samp = output_folder / "data_clinical_patient.txt"

    with data_clin_samp.open("w") as cil_sample:
        cil_sample.write("#Patient Identifier\tAge\tGender\n")
        cil_sample.write("#Patient identifier\tAge\tGender\n")
        cil_sample.write("#STRING\tNUMBER\tSTRING\n")
        cil_sample.write("#1\t1\t1\n")
        cil_sample.write("PATIENT_ID\tAGE\tGENDER\n")

        nested_list = list(table_dict.values())
        list_patients = set(flatten(nested_list))

        for v in list_patients:
            cil_sample.write(f"{v}\tNaN\tNaN\n")



def add_header_patient_type(
    patient_tsv: str,
    datapat_columns: list[str],
    conf_header_type: str,
    final_data_pat: pd.DataFrame) -> pd.DataFrame:
    """Add HEADER_PATIENT_TYPE row to the patient DataFrame.

    Args:
        patient_tsv (str): File name used in log messages for context.
        datapat_columns (list[str]): List of column names.
        conf_header_type (str): Comma-separated types from config.
        final_data_pat (pd.DataFrame): The data to which headers are appended.

    Returns:
        pd.DataFrame: DataFrame with the header row prepended.

    Raises:
        NameError: If types are invalid or column count mismatches.

    """
    if not conf_header_type:
        def_type = ["STRING", "NUMBER", "STRING"]
        header_type = def_type + ["STRING"] * (len(datapat_columns)-3)
        header_type_df = pd.DataFrame([header_type], columns=datapat_columns)
    else:
        types_list = conf_header_type.split(",")
        types_list = [x.strip().upper() for x in types_list]
        for types in types_list:
            if types not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{types} is not a valid type. Please check the given "
                "input in conf.ini. Valid types: STRING, NUMBER, BOOLEAN")
                msg = "The type is not valid: exiting from walk script!"
                raise(NameError(msg))
        try:
            header_type_df = pd.DataFrame([types_list], columns=datapat_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(types_list)}) in "
            "HEADER_PATIENT_TYPE is different from the effective number of "
            f"columns ({len(datapat_columns)}) in {patient_tsv}.")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    return pd.concat([header_type_df, final_data_pat], ignore_index=True)


def add_header_patient_short(
    patient_tsv: str,
    datapat_columns: list[str],
    conf_header_short: str,
    default_row: pd.DataFrame,
    final_data_pat: pd.DataFrame) -> pd.DataFrame:
    """Add HEADER_PATIENT_SHORT row to the patient DataFrame.

    Args:
        patient_tsv (str): File name for context in logging.
        datapat_columns (list[str]): List of column names.
        conf_header_short (str): Header from config, comma-separated.
        default_row (pd.DataFrame): Default row to use if config not provided.
        final_data_pat (pd.DataFrame): The target DataFrame.

    Returns:
        pd.DataFrame: Updated DataFrame.

    """
    if not conf_header_short:
        pat_header_short = default_row
    else:
        try:
            pat_header_short = pd.DataFrame(
                [conf_header_short.split(", ")], columns=datapat_columns)
        except ValueError as err:
            logger.critical(f"The number of column names "
            f"({len(conf_header_short.split(', '))}) in HEADER_PATIENT_SHORT "
            "in conf.ini is different from the effective number of columns "
            f"({len(datapat_columns)}) in {patient_tsv}.")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    return pd.concat([pat_header_short, final_data_pat], ignore_index=True)


def add_header_patient_long(
    patient_tsv: str,
    datapat_columns: list[str],
    conf_header_long: str,
    default_row: pd.DataFrame,
    final_data_pat: pd.DataFrame) -> pd.DataFrame:
    """Add HEADER_PATIENT_LONG row to the patient DataFrame.

    Args:
        patient_tsv (str): File name for context.
        datapat_columns (list[str]): List of column names.
        conf_header_long (str): Header from config, comma-separated.
        default_row (pd.DataFrame): Default row if config missing.
        final_data_pat (pd.DataFrame): DataFrame to prepend header to.

    Returns:
        pd.DataFrame: DataFrame with header row added.

    """
    if not conf_header_long:
        pat_header_long = default_row
    else:
        try:
            pat_header_long = pd.DataFrame(
                [conf_header_long.split(", ")], columns=datapat_columns)
        except ValueError as err:
            logger.critical(f"The number of column names "
            f"({len(conf_header_long.split(', '))}) in HEADER_PATIENT_LONG "
            "in conf.ini is different from the effective number of columns "
            f"({len(datapat_columns)}) in {patient_tsv}.")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    return pd.concat([pat_header_long, final_data_pat], ignore_index=True)


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
    single_sample_vcf_dir.mkdir(exist_ok=True)

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


def extract_multiple_snv(multiple_vcf: str, input_dir: str) -> None:
    """Split a multi-sample SNV VCF file into single-sample files.

    Args:
        multiple_vcf (str): Path to the multi-sample VCF file.
        input_dir (str): Directory where single-sample VCFs are saved.

    """
    input_dir = Path(input_dir)
    vcf_dir = input_dir / "single_sample_vcf"

    if not vcf_dir.exists():
        vcf_dir.mkdir()

    sample_id_path = input_dir / "sample_id.txt"
    cmd = f"vcf-query -l {multiple_vcf} > {sample_id_path}"
    os.system(cmd)

    sample_id_file = input_dir / "sample_id.txt"

    with sample_id_file.open() as file:
        lines = file.readlines()
        for sample_line in lines:
            sample = sample_line.strip()
            output_vcf = vcf_dir / f"{sample}.vcf"

            with output_vcf.open("w") as out_f:
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


def check_field_tsv(row: pd.Series, name: str) -> str:
    """Return a field from a TSV row. Exit if field is missing.

    Args:
        row (pd.Series): A row from a DataFrame.
        name (str): Name of the column to extract.

    Returns:
        str: The value in the specified column.

    """
    try:
        field=str(row[name])
    except KeyError as e:
        logger.critical(f"KeyError: {e} not found! Check if column name is correctly "
        "spelled or if there are tabs/spaces before or after the coloumn key: "
        f"\n{row.index}. \nThis error may also occur if the table columns have "
        "not been separated by tabs!")
        sys.exit(1)
    return field


def get_combined_variant_output_from_folder(
    input_folder: str, file: pd.DataFrame, isinputfile: bool) -> dict:
    """Create a dict mapping SAMPLE_ID to CombinedVariantOutput file paths.

    Args:
        input_folder (str): Base folder for default paths.
        file (pd.DataFrame): Input sample information.
        isinputfile (bool): Whether to use the path from the file or default.

    Returns:
        dict: Mapping from SAMPLE_ID to CombinedOutput file paths.

    """
    combined_dict = {}
    for _,row in file.iterrows():
        sample_id = check_field_tsv(row, "SAMPLE_ID")
        patient_id = check_field_tsv(row, "PATIENT_ID")
        if isinputfile:
            combined_path = Path(check_field_tsv(row, "comb_path"))
        else:
            combined_path = (
                Path(input_folder)
                / "CombinedOutput"
                / f"{patient_id}_CombinedVariantOutput.tsv")

        if combined_path.exists():
            pass
        else:
            logger.warning("comb_path in conf.ini does not exists")
        combined_dict[sample_id] = combined_path
    return combined_dict


def check_input_file(
    output_folder: str, file: str, copy_to: str, sample_id: str,
    relevant: bool = True) -> bool:
    """Copy an input file to a temporary folder, or report it's missing.

    Args:
        output_folder (str): Base output directory.
        file (str): File path to copy.
        copy_to (str): Subfolder name.
        sample_id (str): ID of the sample (for logging).
        relevant (bool): Whether this file type is actually part of the
            selected analysis type (-t). When False, a missing/unset path
            is expected and not worth a warning - only an unexpected state
            (path given but the file doesn't exist) still gets logged.

    Returns:
        bool: True when `file` is unset for this sample and `relevant` is
            True - the caller aggregates this flag across the whole sample
            batch (see check_folders/transform_input) so a config-wide gap
            (e.g. every Guardant sample's blank comb_path) logs ONE warning
            instead of repeating the identical line once per sample. A
            path that was given but doesn't exist on disk is still warned
            about immediately here, since that message differs per sample.

    """
    file_path = Path(file)
    destination = Path(output_folder) / "temp" / copy_to

    # `file` empty (e.g. Guardant's always-blank comb_path) must be checked
    # BEFORE file_path.exists(): Path("") == Path(".") and the current
    # directory always exists, so file_path.exists() would be True and this
    # would try to `cp .` the whole cwd - producing a confusing "cp: -r not
    # specified" shell error instead of the intended "no path set" warning,
    # without actually copying anything.
    if not file:
        return relevant

    if file_path.exists():
        os.system(f"cp {file_path} {destination}")
    else:
        logger.warning(f"{file_path} not found")
    return False


def check_folders(
    output_folder: str,
    snv_path: str,
    cnv_path: str,
    combout: str,
    sample_id: str,
    vcf_type: str | None = None) -> dict[str, bool]:
    """Check presence of SNV, CNV, and CombinedOutput files for a sample.

    Only flags a file type that's actually missing or unset for a sample
    where it's needed - a missing/unset path for a file type excluded by
    the selected analysis type (-t) is expected, not a problem (see
    check_input_file's `relevant` flag).

    Args:
        output_folder (str): Base output directory.
        snv_path (str): Path to SNV file.
        cnv_path (str): Path to CNV file.
        combout (str): Path to CombinedOutput file.
        sample_id (str): ID of the sample.
        vcf_type (str | None): Selected analysis type restriction (-t), or
            None if no restriction was given.

    Returns:
        dict[str, bool]: For each of "SNV"/"CNV"/"CombinedOutput", True if
            this sample has no path set for it (and it's relevant) - the
            caller aggregates these across the batch (see transform_input)
            to log one warning per file type instead of one per sample.

    """
    return {
        "SNV": check_input_file(
            output_folder, snv_path, "SNV", sample_id,
            relevant=vcf_type not in ["cnv", "fus", "tab"]),
        "CNV": check_input_file(
            output_folder, cnv_path, "CNV", sample_id,
            relevant=vcf_type not in ["snv", "fus", "tab"]),
        "CombinedOutput": check_input_file(
            output_folder, combout, "CombinedOutput", sample_id),
    }


def transform_input(
    tsv: str,
    clin_pzt: str,
    fusion_tsv: str,
    output_folder: str,
    multiple: bool,
    vcf_type: str | None = None) -> str:
    """Prepare temp folder structure and copy necessary files.

    Args:
        tsv (str): Sample TSV file path.
        clin_pzt (str): Clinical patient TSV file path.
        fusion_tsv (str): Fusion TSV file path.
        output_folder (str): Output directory.
        multiple (bool): Whether input is from multi-sample VCF.
        vcf_type (str | None): Selected analysis type restriction (-t), used
            to silence expected-missing-file warnings for file types the
            current run doesn't need (see check_folders).

    Returns:
        str: Path to the temporary folder created.

    """
    base_temp = Path(output_folder) / "temp"
    (base_temp).mkdir(parents=True, exist_ok=True)
    (base_temp / "SNV").mkdir(parents=True, exist_ok=True)
    (base_temp / "CNV").mkdir(parents=True, exist_ok=True)
    (base_temp / "CombinedOutput").mkdir(parents=True, exist_ok=True)
    (base_temp / "FUSIONS").mkdir(parents=True, exist_ok=True)

    sample_path = Path(output_folder) / "temp" / "sample.tsv"
    patient_path = Path(output_folder) / "temp" / "patient.tsv"
    fusion_path = Path(output_folder) / "temp" / "FUSIONS" / "fusion.tsv"

    shutil.copy(tsv, sample_path)

    if clin_pzt != "":
        shutil.copy(clin_pzt, patient_path)

    if fusion_tsv != "":
        shutil.copy(fusion_tsv, fusion_path)

    if multiple:
        snv_path = config.get("Multiple", "SNV")
        cnv_path = config.get("Multiple", "CNV")
        combout = config.get("Multiple", "COMBOUT")

        missing = check_folders(
            output_folder, snv_path, cnv_path, combout, "multiple", vcf_type)
        for copy_to, was_missing in missing.items():
            if was_missing:
                logger.warning(
                    f"No final_path set in conf.ini for {copy_to}!")

    else:
        tsv_file = pd.read_csv(tsv, sep="\t", dtype="string", keep_default_na=False)

        # Collected per file type instead of warning inline per sample: a
        # config-wide gap (e.g. every Guardant sample's blank comb_path,
        # since Guardant never has a CombinedOutput) would otherwise log
        # the identical line once per sample - see check_input_file.
        missing_by_type: dict[str, list[str]] = {
            "SNV": [], "CNV": [], "CombinedOutput": []}
        for _,row in tsv_file.iterrows():
            sample_id = row["SAMPLE_ID"]
            snv_path = row["snv_path"]
            cnv_path = row["cnv_path"]
            combout = row["comb_path"]

            missing = check_folders(
                output_folder, snv_path, cnv_path, combout, sample_id, vcf_type)
            for copy_to, was_missing in missing.items():
                if was_missing:
                    missing_by_type[copy_to].append(sample_id)

        total_samples = len(tsv_file)
        for copy_to, missing_samples in missing_by_type.items():
            if not missing_samples:
                continue
            if len(missing_samples) == total_samples:
                logger.warning(
                    f"No final_path set in conf.ini for {copy_to} on any "
                    "of the samples!")
            else:
                for sample_id in missing_samples:
                    logger.warning(
                        f"No final_path set in conf.ini for sample "
                        f"{sample_id}'s {copy_to}!")

    return str(Path(output_folder) / "temp")


def fill_fusion_from_temp(
    input_path: str,
    fusion_table_file: str,
    clin_file: pd.DataFrame,
    fusion_files: list) -> None:
    """Collect fusion events from temp directory into a unified table.

    Args:
        input_path (str): Temp input path.
        fusion_table_file (str): Output file path for combined fusions.
        clin_file (pd.DataFrame): Clinical file to match SAMPLE_IDs.
        fusion_files (list): List of fusion TSV filenames.

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
            min_read_count = 15
            for fus in ff.itertuples(index=False):
                if (str(fus.Sample_Id).strip() in clin_file["SAMPLE_ID"].astype(str).to_numpy() and
                int(fus.Normal_Paired_End_Read_Count) >= min_read_count):
                    fusion_table.write("\t".join(map(str, fus)) + "\n")


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
        header = (
            "Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\t"
            "Normal_Paired_End_Read_Count\tEvent_Info\tRNA_Support\n")
        fusion_table.write(header)

        for k, v in combined_dict.items():
            fusions=[]
            try:
                fusions = tsv.get_fusions(Path(v))
            except Exception:
                logger.error("Something went wrong while reading Fusion section "
                f"for sample {k}")
            if len(fusions) == 0:
                continue

            for fus in fusions:
                if len(fusions) > 0:
                    site1_hugo_symbol = fus["Site1_Hugo_Symbol"]
                    site2_hugo_symbol = fus["Site2_Hugo_Symbol"]
                    if site2_hugo_symbol == "CASC1":
                        site2_hugo_symbol = "DNAI7"
                    site1_Chromosome = fus["Site1_Chromosome"]
                    site2_Chromosome = fus["Site2_Chromosome"]
                    site1_Position = fus["Site1_Position"]
                    site2_Position = fus["Site2_Position"]

                    if eval("int(fus['Normal_Paired_End_Read_Count'])" + thr_fus):
                        fusion_table.write(
                            str(k).strip() + "\tSOMATIC\tFUSION\t" +
                            str(site1_hugo_symbol) + "\t" +
                            str(site2_hugo_symbol) + "\t" +
                            fus["Normal_Paired_End_Read_Count"] + "\t" +
                            fus["Event_Info"] + " Fusion\tYes\n")


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
    header = (
        "Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\t"
        "Normal_Paired_End_Read_Count\tEvent_Info\tRNA_Support\n")

    existing_lines = []
    if fusion_table_path.exists():
        with fusion_table_path.open() as f:
            existing_lines = f.readlines()
        if existing_lines:
            header = existing_lines[0]
            existing_lines = [
                line for line in existing_lines[1:]
                if len(line.split("\t")) <= 2 or line.split("\t")[2] != "SPLICE"]

    new_rows = []
    for k, v in combined_dict.items():
        splice_variants = []
        try:
            splice_variants = tsv.get_splice_variants(Path(v))
        except Exception:
            logger.error("Something went wrong while reading Splice Variants "
            f"section for sample {k}")

        for sv in splice_variants:
            gene = sv["Gene"]
            ssr = sv["Splice_Supporting_Reads"]
            try:
                if not eval("int(ssr)" + thr_splice):
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
                f"reference_reads_transcript={sv['Reference_Reads_Transcript']} "
                "- this is the first code path to see a real populated "
                "[Splice Variants] row, please spot-check this row against "
                "the source CombinedVariantOutput.tsv file.")
            new_rows.append(
                str(k).strip() + "\tSOMATIC\tSPLICE\t" +
                str(gene) + "\t" + str(gene) + "\t" +
                ssr + "\t" + event_info + "\tYes\n")

    if not new_rows and not existing_lines:
        logger.info(
            "No splice variants found for this batch - data_sv.txt left as "
            "is (not created if the fusion step above didn't create it "
            "either).")
        return

    with fusion_table_path.open("w") as fusion_table:
        fusion_table.write(header)
        fusion_table.writelines(existing_lines)
        fusion_table.writelines(new_rows)

    logger.info(
        f"data_sv.txt now has {len(new_rows)} new splice variant row(s) "
        f"(plus {len(existing_lines)} pre-existing non-splice row(s) carried "
        "over from the fusion step above, if any).")


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


def _is_blank_value(value: object) -> bool:
    """Return True if a MSI/TMB VALUE or THR cell should be treated as empty."""
    if value is None:
        return True
    if isinstance(value, float) and np.isnan(value):
        return True
    return str(value).strip().upper() in ("", "NAN", "NA")


def _resolve_biomarker_threshold(
    sample_id: str,
    marker_name: str,
    value: object,
    thr_given: object,
    compute_thr) -> str:
    """Resolve a MSI/TMB threshold label for one sample.

    Precedence (identical for MSI and TMB, per explicit user design so the
    input can carry a pipeline's own pre-computed threshold - e.g. Guardant
    - instead of forcing everything through conf.ini's thresholds):

    - THR pre-filled in the input -> keep it as-is, never overwritten by
      conf.ini. Logged as a warning (conf.ini threshold ignored) if VALUE is
      also populated, or as info (THR trusted alone) if VALUE is empty.
    - THR empty, VALUE populated -> compute from conf.ini's own threshold,
      exactly as before.
    - Neither populated -> "NA".

    Args:
        sample_id (str): Sample ID, for logging.
        marker_name (str): "MSI" or "TMB", for logging.
        value (object): The raw MSI/TMB numeric value cell (may be blank).
        thr_given (object): The raw MSI_THR/TMB_THR cell (may be blank).
        compute_thr (Callable[[str], str]): conf.ini-threshold-based
            classifier, called only when thr_given is blank.

    Returns:
        str: The resolved threshold label to write to data_clinical_sample.txt.

    """
    value_present = not _is_blank_value(value)
    thr_present = not _is_blank_value(thr_given)

    if thr_present:
        thr_given = str(thr_given).strip()
        if value_present:
            logger.warning(
                f"Sample {sample_id}: {marker_name}_THR is pre-filled "
                f"('{thr_given}') and {marker_name} value is also populated "
                f"- keeping the provided {marker_name}_THR and ignoring "
                f"conf.ini's {marker_name} threshold.")
        else:
            logger.info(
                f"Sample {sample_id}: {marker_name} value is empty but "
                f"{marker_name}_THR is pre-filled ('{thr_given}') - keeping "
                "it as-is.")
        return thr_given

    if value_present:
        return compute_thr(value)

    return "NA"


def fill_from_file(
    table_dict_patient: dict[str, list],
    file_input_clinical: pd.DataFrame,
    msi_thr: str,
    tmb_thr: dict[str, str]) -> dict[str, list]:
    """Populate clinical table dictionary with MSI and TMB values and statuses.

    Honors a pre-filled MSI_THR/TMB_THR column in the input (if present) over
    conf.ini's own thresholds - see _resolve_biomarker_threshold(). This
    replaces the previous behavior of always recomputing Stable/Unstable and
    the TMB category from conf.ini, which forced non-Illumina pipelines
    (e.g. Guardant) that already know their own threshold to encode it via
    a fake VALUE just to get the right label out.

    Args:
        table_dict_patient (dict): Patient clinical data dictionary.
        file_input_clinical (pd.DataFrame): Clinical input with SAMPLE_ID,
            MSI, TMB and, if present, MSI_THR/TMB_THR.
        msi_thr (str): Threshold condition string for MSI (used with eval).
        tmb_thr (dict): Map of TMB categories to threshold conditions.

    Returns:
        dict: Updated clinical data dictionary.

    """
    has_msi_thr_col = "MSI_THR" in file_input_clinical.columns
    has_tmb_thr_col = "TMB_THR" in file_input_clinical.columns

    def compute_msi(m: object) -> str:
        return "Stable" if eval("float(m)" + msi_thr) else "Unstable"

    def compute_tmb(t: object) -> str:
        for _k, _v in tmb_thr.items():
            if eval("float(t)" + _v):
                return _k
        logger.warning(f"TMB {t} out of range for {k}")
        return "Out of threshold ranges"

    for _, row in file_input_clinical.iterrows():
        k = row["SAMPLE_ID"]
        m = row["MSI"]
        t = row["TMB"]
        table_dict_patient[k].append(m)
        table_dict_patient[k].append(t)

        msi_thr_given = row["MSI_THR"] if has_msi_thr_col else None
        table_dict_patient[k].append(
            _resolve_biomarker_threshold(k, "MSI", m, msi_thr_given, compute_msi))

        tmb_thr_given = row["TMB_THR"] if has_tmb_thr_col else None
        table_dict_patient[k].append(
            _resolve_biomarker_threshold(k, "TMB", t, tmb_thr_given, compute_tmb))

    return table_dict_patient


def fill_from_combined(
    combined_dict: dict[str, str],
    table_dict_patient: dict[str, list],
    msi_sites_thr: str,
    msi_thr: str,
    tmb: dict[str, str]) -> dict[str, list]:
    """Extract MSI and TMB from combined output files and update clinical dict.

    Args:
        combined_dict (dict): Map sample ID to CombinedVariantOutput path.
        table_dict_patient (dict): Clinical data dictionary to update.
        msi_sites_thr (str): Threshold condition for MSI sites (eval string).
        msi_thr (str): Threshold condition for MSI stability (eval string).
        tmb (dict): Map TMB categories to threshold conditions.

    Returns:
        dict: Updated clinical data dictionary.

    """
    for k, v in combined_dict.items():
        try:
            tmv_msi = tsv.get_msi_tmb(Path(v), SAMPLE_TYPE)
        except Exception:
            logger.error(f"Something went wrong reading CombinedOutput for sample "
                         f"{k}! MSI/TMB set to NA for this sample.")
            # Keep table_dict_patient[k] the same length as every successfully
            # processed sample (4 appends below) - otherwise
            # pd.DataFrame.from_dict(table_dict_patient) in write_clinical_sample
            # raises "All arrays must be of the same length" for the whole batch.
            table_dict_patient[k].extend(["NA", "NA", "NA", "NA"])
            continue

        if (
            tmv_msi["MSI"][0][1] != "NA" and
            eval("float(tmv_msi['MSI'][0][1])" + msi_sites_thr)):
            table_dict_patient[k].append(tmv_msi["MSI"][1][1])
        else:
            table_dict_patient[k].append("NA")

        table_dict_patient[k].append(tmv_msi["TMB_Total"])

        if (tmv_msi["MSI"][1][1] != "NA" and table_dict_patient[k][-2] != "NA"):
                msi_value = float(tmv_msi["MSI"][1][1])

                if SAMPLE_TYPE == "LIQUID":
                    if msi_value >= THRESHOLD_MSI_LIQUID:
                        table_dict_patient[k].append("Unstable")
                    else:
                        table_dict_patient[k].append("Stable")

                elif SAMPLE_TYPE == "SOLID":
                    if eval(str(msi_value) + msi_thr):
                        table_dict_patient[k].append("Stable")
                    else:
                        table_dict_patient[k].append("Unstable")

        else:
            table_dict_patient[k].append("NA")

        if tmv_msi["TMB_Total"] != "NA":
            found = False
            for _k, _v in tmb.items():
                if eval(tmv_msi["TMB_Total"] + _v):
                    table_dict_patient[k].append(_k)
                    found = True
                    break
            if not found:
                logger.warning(
                    f"TMB {tmv_msi['TMB_Total']} out of range for {k}")
                table_dict_patient[k].append("Out of threshold ranges")
        else:
            table_dict_patient[k].append("NA")

    return table_dict_patient


def input_extraction_file(input_f: list) -> tuple:
    """Extract sample, patient, and fusion TSV file paths from input list.

    Args:
        input_f (list): List of file path strings (up to 3 elements).

    Returns:
        tuple: Paths to sample, patient, and fusion TSV files.

    """
    sample_tsv = input_f[0]
    patient_tsv, fusion_tsv = "",""

    limit = 2
    if len(input_f) > 1:
        patient_tsv = input_f[1].strip()
        if len(input_f) > limit:
            fusion_tsv=input_f[limit].strip()

    return sample_tsv, patient_tsv, fusion_tsv


def input_extraction_folder(input_path: str) -> tuple:
    """Extract sample, patient, and fusion TSV file paths from input folder.

    Args:
        input_path (str): Path to input folder.

    Returns:
        tuple: Paths to sample.tsv, patient.tsv, and fusion TSV files.

    """
    input_path = Path(input_path)
    patient_tsv, fusion_tsv = "",""
    sample_tsv = input_path / "sample.tsv"
    if (input_path / "patient.tsv").exists():
        patient_tsv = input_path / "patient.tsv"
    if (input_path / "FUSIONS" / "Fusions.tsv").exists():
        fusion_tsv = input_path / "FUSIONS" / "Fusions.tsv"
    return sample_tsv, patient_tsv, fusion_tsv


def validate_input(
    oncokb: bool,
    vcf_type: str | None,
    filters: str,
    cancer: str,
    input_path: str) -> None:
    """Validate the input configuration and required files for the pipeline.

    Args:
        oncokb (bool): Whether OncoKB annotation is enabled.
        vcf_type (str | None): VCF type or None.
        filters (str): Filter options as string.
        cancer (str): Cancer ID to validate against cBioPortal.
        input_path (str): Path to input file or folder.

    Raises:
        Exception: If required columns are missing or config fields are empty.
        SystemExit: If cancer ID is not recognized by cBioPortal.

    """
    if not isinputfile:
        sample_path = Path(input_path) / "sample.tsv"
        file_tsv = pd.read_csv(sample_path, sep="\t", dtype=str)
    elif isinputfile:
        file_tsv = pd.read_csv(input_path, sep="\t", dtype=str)

    column_list = [
        "SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path",
        "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR"]
    if not all(name in file_tsv.columns for name in column_list):
        logger.critical('Required columns: "SAMPLE_ID", "PATIENT_ID",'
        '"ONCOTREE_CODE", "snv_path", "cnv_path", "comb_path", "MSI", '
        '"TMB", "MSI_THR", "TMB_THR"')
        msg = "The input file is missing some important columns!"
        raise(ValueError(msg))

    if oncokb and config.get("OncoKB", "ONCOKB") == "":
        msg = "oncokb option was set but ONCOKB field in conf.ini is empty!"
        raise ValueError(msg)

    # Fail here, before the expensive VEP/vcf2maf step, rather than deep
    # inside filter_oncokb() after that work is already done: a conf.ini
    # that predates the SNV/CNV/FUSION split (a single ONCOKB_FILTER used
    # to cover all three) is otherwise only caught with a raw
    # configparser.NoOptionError once filtering actually runs.
    oncokb_filter_checks = []
    if oncokb and "o" in filters:
        oncokb_filter_checks.append(("Filters", "ONCOKB_FILTER_SNV"))
    if oncokb and "o" in filters and vcf_type not in ["snv", "fus", "tab"]:
        oncokb_filter_checks.append(("Cna", "ONCOKB_FILTER_CNV"))
    if oncokb and "o" in filters and vcf_type not in ["cnv", "snv", "tab"]:
        oncokb_filter_checks.append(("FUSION", "ONCOKB_FILTER_FUSION"))

    for section, key in oncokb_filter_checks:
        if not config.has_option(section, key):
            logger.critical(
                f"conf.ini is missing '{key}' under [{section}] - required "
                "for this run's OncoKB filtering. If this conf.ini predates "
                "the SNV/CNV/FUSION split (a single ONCOKB_FILTER used to "
                "cover all three), add ONCOKB_FILTER_SNV under [Filters], "
                "ONCOKB_FILTER_CNV under [Cna] and ONCOKB_FILTER_FUSION "
                "under [FUSION] - see the conf.ini template.")
            msg = f"Missing conf.ini option: [{section}] {key}"
            raise ValueError(msg)

    if SAMPLE_TYPE not in {"SOLID", "LIQUID"}:
        raise ValueError('Please select a sample type between "Solid" and "Liquid" in conf.ini.')

    if (vcf_type is None or "snv" in vcf_type) and (not VEP_PATH or not VEP_DATA):
        msg = "VEP_PATH and/or VEP_DATA field in conf.ini is empty!"
        raise ValueError(msg)

    if not REF_FASTA:
        msg = "REF_FASTA field in conf.ini is empty!"
        raise ValueError(msg)

    if "o" in filters and not oncokb:
        logger.warning("OncoKB filter was selected in filters options but -k option "
        "was not set. This filtering will be ignored.")

    cancer_cbio = pd.read_csv("cancer_list.txt", sep="\t")
    cancer_cbio = cancer_cbio["TYPE_OF_CANCER_ID"].to_numpy().tolist()

    if cancer not in cancer_cbio:
        logger.critical(f"The cancer_id '{cancer}' is not recognize by cBioPortal. "
        "Check the cancer_list.txt to find the correct cancer id")
        sys.exit(1)


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
        except Exception:
            logger.error(f"Something went wrong while reading Exon-Level CNVs for sample {k}")

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


def update_data_clinical_with_exon_info(
    combined_output: Path,
    combined_dict: dict,
    output_folder: Path):
    """Updates clinical_sample WITH exon-level CNV information for BRCA1 and BRCA2

    Args:
        combined_output (Path): Directory containing combined output files.
        combined_dict (dict): Mapping of {sample_id: path_to_tsv}.
        output_folder (Path): Output directory where data_clinical_sample.txt is located.

    """
    data_clin_path = Path(output_folder) / "data_clinical_sample.txt"

    if combined_output.is_dir() and any(combined_output.iterdir()):
        all_exon_dfs = []

        for sample_id, tsv_path in combined_dict.items():
            try:
                exonic = tsv.get_exons(Path(tsv_path))
            except Exception as e:
                logger.error(f"Error while reading exon-level CNVs for sample {sample_id}: {e}")
                exonic = []

            exon_df = build_exon_df(output_folder, exonic, sample_id=sample_id)
            all_exon_dfs.append(exon_df)

        if all_exon_dfs:
            all_exon_dfs = pd.concat(all_exon_dfs, ignore_index=True)

            with open(data_clin_path, "r") as f:
                header_lines = [next(f) for _ in range(4)]

            data_clin_df = pd.read_csv(
                data_clin_path, 
                sep="\t", 
                header=4, 
                dtype={"SAMPLE_ID": str}
            )
            merged_data_clin = pd.merge(data_clin_df, all_exon_dfs, on="SAMPLE_ID", how="left")

            new_columns = ["Exonic_BRCA1", "Exonic_BRCA2", "BRCA1_details", "BRCA2_details"]

            updated_headers = []
            updated_headers.append(header_lines[0].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n")
            updated_headers.append(header_lines[1].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n")
            updated_headers.append(header_lines[2].rstrip("\n") + "\t" + "\t".join(["STRING"] * len(new_columns)) + "\n")
            updated_headers.append(header_lines[3].rstrip("\n") + "\t" + "\t".join(["1"] * len(new_columns)) + "\n")

            with open(data_clin_path, "w") as f:
                f.writelines(updated_headers)
                merged_data_clin.to_csv(f, sep="\t", index=False)

        else:
            logger.warning("No exon CNV data available to update the clinical file.")


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


def update_data_clinical_with_hrd_info(
    combined_output: Path,
    combined_dict: dict,
    output_folder: Path) -> None:
    """Update data_clinical_sample.txt with HRD/GIS biomarker info.

    Genomic Instability Score, Tumor Fraction and Ploidy are only reported by
    DRAGEN for samples run with the TSO500 HRD feature enabled - for every
    other sample these three new columns are just "NA", the same way
    update_data_clinical_with_exon_info() layers BRCA exon-level info onto
    this same file without touching the rest of it.

    Args:
        combined_output (Path): Directory containing combined output files.
        combined_dict (dict): Mapping of {sample_id: path_to_tsv}.
        output_folder (Path): Output directory where data_clinical_sample.txt is located.

    """
    data_clin_path = Path(output_folder) / "data_clinical_sample.txt"

    if combined_output.is_dir() and any(combined_output.iterdir()):
        all_gis_rows = []

        for sample_id, tsv_path in combined_dict.items():
            try:
                gis = tsv.get_gis(Path(tsv_path))
            except Exception as e:
                logger.error(f"Error while reading GIS/HRD info for sample {sample_id}: {e}")
                gis = {"GIS": "NA", "Tumor_Fraction": "NA", "Ploidy": "NA"}

            all_gis_rows.append({
                "SAMPLE_ID": sample_id,
                "GENOMIC_INSTABILITY_SCORE": gis["GIS"],
                "TUMOR_FRACTION": gis["Tumor_Fraction"],
                "PLOIDY": gis["Ploidy"],
            })

        if all_gis_rows:
            all_gis_df = pd.DataFrame(all_gis_rows)

            with open(data_clin_path, "r") as f:
                header_lines = [next(f) for _ in range(4)]

            data_clin_df = pd.read_csv(
                data_clin_path,
                sep="\t",
                header=4,
                dtype={"SAMPLE_ID": str}
            )
            merged_data_clin = pd.merge(data_clin_df, all_gis_df, on="SAMPLE_ID", how="left")

            new_columns = ["GENOMIC_INSTABILITY_SCORE", "TUMOR_FRACTION", "PLOIDY"]

            updated_headers = [
                header_lines[0].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n",
                header_lines[1].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n",
                header_lines[2].rstrip("\n") + "\t" + "\t".join(["NUMBER"] * len(new_columns)) + "\n",
                header_lines[3].rstrip("\n") + "\t" + "\t".join(["1"] * len(new_columns)) + "\n",
            ]

            with open(data_clin_path, "w") as f:
                f.writelines(updated_headers)
                merged_data_clin.to_csv(f, sep="\t", index=False)
        else:
            logger.warning("No sample data available to update HRD/GIS clinical info.")


@dataclass
class WalkContext:
    """Shared state produced once by _walk_setup() and read by the four
    independent stages below (_walk_process_cnv/_snv/_fusion and
    _walk_write_clinical_tables).

    None of those four stages depends on another's *output* - CNV, SNV,
    fusion/splice and the clinical tables are each derived only from this
    setup context - they are simply called in sequence by walk_folder() for
    simplicity and readable logs. That independence is exactly what lets
    them be exposed as separate Snakemake rules (see Snakefile) with real
    per-stage resume/parallelism, instead of the previous single
    ~400-line walk_folder() that Snakemake could only shell out to as one
    opaque unit.

    isinputfile mirrors the module-level `isinputfile` global _walk_setup()
    sets: the four stage functions read that bare global directly (unchanged
    from the pre-split code), which only survives within a single process.
    A caller driving each stage from a separate process (e.g. walk_stage.py,
    for real per-stage Snakemake rules) must restore it from here
    (`walk.isinputfile = ctx.isinputfile`) before invoking a stage function.
    """

    output_folder: str
    input_folder: Path
    input_path: object
    patient_tsv: str
    fusion_tsv: str
    input_folder_snv: Path
    input_folder_cnv: Path
    case_folder_arr: dict
    case_folder_arr_cnv: dict | None
    clin_file: pd.DataFrame
    clin_sample_path: Path
    multiple: bool
    oncokb: bool
    cancer: str
    filters: str
    vcf_type: str | None
    resume: bool
    isinputfile: bool
    sigma: bool = False


def _walk_setup(
    input_path: list,
    multiple: bool,
    output_folder: str,
    oncokb: bool,
    cancer: str,
    overwrite_output: bool,
    resume: bool,
    vcf_type: str | None,
    filters: str,
    sigma: bool = False) -> WalkContext:
    """Stage 0: resolve input, create/resume the output folder, and check
    which of SNV/CNV/CombinedOutput are actually present.

    Every other stage (_walk_process_cnv/_snv/_fusion,
    _walk_write_clinical_tables) depends only on this stage's output, not
    on each other's - see WalkContext and walk_folder().
    """
    logger.info("Starting walk_folder script:")
    logger.info(
        f"walk_folder args [input:{input_path}, output_folder:{output_folder}, "
        f"overwrite:{overwrite_output}, resume:{resume}, vcf_type:{vcf_type}, "
        f"filters:{filters}, multiple:{multiple}]")

    if not Path(input_path[0]).exists():
        msg = f"No valid file/folder {input_path} found. Check your input path"
        raise FileNotFoundError(msg)

    global isinputfile

    path = Path(input_path[0])
    if path.is_dir():
        isinputfile = False
    elif path.is_file():
        isinputfile = True
    validate_input(oncokb, vcf_type, filters, cancer, input_path[0])

    if not resume or not (Path(output_folder) / "temp").exists():
        output_folder = create_folder(output_folder, overwrite_output, resume)
    else:
        get_version_list(output_folder)

    if not isinputfile:
        input_folder = input_path[0]
        input_path, patient_tsv, fusion_tsv = input_extraction_folder(input_folder)
        check_multiple_folder(input_folder, multiple)

    elif isinputfile:
        input_path, patient_tsv, fusion_tsv = input_extraction_file(input_path)
        check_multiple_file(input_path, multiple)
        input_folder = transform_input(
            input_path, patient_tsv, fusion_tsv, output_folder, multiple, vcf_type)

    else:
        logger.critical(f"The input {input_path} isn't a file nor a folder")
        msg = "Exiting from walk script!"
        raise(FileNotFoundError(msg))

    img_path = Path("docs") / "img" / "logo_VARAN.png"
    if img_path.exists():
        img_output_dir = Path(output_folder) / "img"
        img_output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(img_path, img_output_dir / "logo_VARAN.png")

    input_folder = Path(input_folder)
    input_folder_snv = (input_folder / "SNV").resolve()
    input_folder_cnv = (input_folder / "CNV").resolve()
    input_folder_comb_out = (input_folder / "CombinedOutput").resolve()

    if not (
        any(input_folder_snv.iterdir()) or
        any(input_folder_cnv.iterdir()) or
        any(input_folder_comb_out.iterdir())):
        msg = ("Empty input folder(s)! Either SNV, CNV or CombinedOutput folder should"
        "contain input files.")
        raise ValueError(msg)

    maf_path = Path(output_folder) / "maf"
    maf_zip_path = Path(output_folder) / "maf.zip"
    clin_sample_path = Path(input_folder) / "sample.tsv"

    try:
        clin_file = pd.read_csv(clin_sample_path, sep="\t", dtype=str)
    except Exception as err:
        logger.critical(f"Something went wrong while reading {clin_sample_path}!")
        msg = "Error in reading the input file! Please check again."
        raise OSError(msg) from err

    # MAF/clinical consistency only makes sense when this run actually
    # processes SNV data (vcf_type is None, or explicitly "snv") - the same
    # gate _walk_process_snv itself uses. A CNV-only/fusion-only/tab-only
    # resume (-t cnv/fus/tab -R) never produces a MAF in the first place, so
    # checking for one here was both pointless and a crash risk (maf_samples
    # was left undefined - see below - if neither maf/ nor maf.zip existed).
    if resume and vcf_type not in ["cnv", "fus", "tab"]:
        maf_samples = None
        if maf_path.exists():
            try:
                maf_samples = {
                    entry.name
                    for entry in maf_path.iterdir()
                    if entry.is_file()}
            except Exception as err:
                logger.critical(
                    "Error accessing the maf folder. Please check its integrity.")
                msg = "Error reading the maf folder."
                raise Exception(msg) from err

        elif maf_zip_path.exists():
            try:
                with zipfile.ZipFile(maf_zip_path, "r") as zipped_maf:
                    zipped_maf.extractall(maf_path)
                maf_samples = {
                    entry.name
                    for entry in maf_path.iterdir()
                    if entry.is_file()}
            except zipfile.BadZipFile as err:
                logger.critical("Corrupted maf.zip file.")
                msg = ("Unable to extract the ZIP file. Try decompressing it manually "
                "and run Varan again.")
                raise Exception(msg) from err

        if maf_samples is None:
            logger.warning(
                "Resuming an SNV-inclusive run, but no maf/ folder or maf.zip "
                "was found to resume from - proceeding as if this were the "
                "first run for the MAF step.")
        else:
            clin_samples = set(clin_file["SAMPLE_ID"])
            clin_in_maf = all(
                any(clin_sample in maf_sample for maf_sample in maf_samples)
                for clin_sample in clin_samples)
            maf_in_clin = all(
                any(clin_sample in maf_sample for clin_sample in clin_samples)
                for maf_sample in maf_samples)

            if not (clin_in_maf and maf_in_clin) and len(maf_samples) != 0:
                logger.critical("It seems you are resuming an existing study with a "
                "different set of input samples. Please verify the sample consistency!")
                msg = "Sample mismatch detected."
                raise FileNotFoundError(msg)

        zip_maf = config.get("Zip", "ZIP_MAF")
        zip_maf = check_bool(zip_maf)
        if maf_zip_path.exists() and not zip_maf:
            maf_zip_path.unlink()

    input_folder_snv = Path(input_folder_snv)
    input_folder_cnv = Path(input_folder_cnv)

    snv_folder_empty = len(list(input_folder_snv.iterdir())) == 0
    cnv_folder_empty = len(list(input_folder_cnv.iterdir())) == 0

    if vcf_type is None:
        if snv_folder_empty and not cnv_folder_empty:
            vcf_type = "cnv"
            logger.info("SNV path was empty, the analysis will exclude SNV")
        elif cnv_folder_empty and not snv_folder_empty:
            vcf_type = "snv"
            logger.info("CNV path was empty, the analysis will exclude CNV")
        elif snv_folder_empty and cnv_folder_empty:
            # Both empty is a valid CombinedOutput-only run (fusions, splice
            # variants, MSI/TMB all come from CombinedOutput, not SNV/CNV
            # VCFs). Leaving vcf_type as None here is required: setting it to
            # "cnv" or "snv" (as a naive single `if` chain would, since the
            # SNV check would win) makes _walk_process_fusion's
            # `vcf_type in ["cnv", "snv", "tab"]` gate skip fusion/splice
            # processing entirely, silently dropping data_sv.txt.
            logger.info(
                "SNV and CNV paths were both empty; proceeding with "
                "CombinedOutput-only processing.")

    case_folder_arr_cnv = None
    if input_folder_cnv.exists() and vcf_type not in ["snv", "fus", "tab"]:
        if multiple:
            multivcf = next(
                (i for i in input_folder_cnv.iterdir() if i.name.endswith(".vcf")),
                None)
            if multivcf:
                extract_multiple_cnv(multivcf, input_folder_cnv)
                input_folder_cnv = Path(input_folder_cnv) / "single_sample_vcf"
        logger.info("Checking CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(input_folder_cnv)
        logger.info("Everything ok!")

    input_folder_snv = Path(input_folder_snv)

    if input_folder_snv.exists() and vcf_type not in ["cnv", "fus", "tab"]:
        if multiple:
            multivcf = next(f for f in input_folder_snv.iterdir() if f.suffix == ".vcf")
            extract_multiple_snv(multivcf, input_folder_snv)
            input_folder_snv = input_folder_snv / "single_sample_vcf"

        logger.info("Checking SNV files...")

    case_folder_arr = get_snv_from_folder(input_folder_snv)
    logger.info("Everything ok!")

    return WalkContext(
        output_folder=output_folder,
        input_folder=Path(input_folder),
        input_path=input_path,
        patient_tsv=patient_tsv,
        fusion_tsv=fusion_tsv,
        input_folder_snv=Path(input_folder_snv),
        input_folder_cnv=Path(input_folder_cnv),
        case_folder_arr=case_folder_arr,
        case_folder_arr_cnv=case_folder_arr_cnv,
        clin_file=clin_file,
        clin_sample_path=clin_sample_path,
        multiple=multiple,
        oncokb=oncokb,
        cancer=cancer,
        filters=filters,
        vcf_type=vcf_type,
        resume=resume,
        isinputfile=isinputfile,
        sigma=sigma,
    )


def _walk_process_cnv(ctx: WalkContext) -> None:
    """Stage: CNV calls (data_cna*.txt) and BRCA exon-level CNV table.

    Depends only on WalkContext, independent of _walk_process_snv/_fusion.
    """
    if ctx.input_folder_cnv.exists() and ctx.vcf_type not in ["snv", "fus", "tab"]:
        logger.info("Managing CNV files...")
        cnv_type_from_folder(
            ctx.input_folder, ctx.case_folder_arr_cnv, ctx.output_folder,
            ctx.oncokb, ctx.cancer, ctx.multiple, ctx.filters)

    combined_output_folder = Path(ctx.input_folder) / "CombinedOutput"
    if (combined_output_folder.exists()
        and any(f.is_file() for f in combined_output_folder.iterdir())):
        combined_dict = get_combined_variant_output_from_folder(
            ctx.input_folder, ctx.clin_file, isinputfile)

        exon_file_output = Path(ctx.output_folder) / "exon_CNA_data.txt"
        exon_events = write_exon_brca(exon_file_output, combined_dict)
        write_exon_brca_generic_assay(ctx.output_folder, exon_events)


def _walk_process_snv(ctx: WalkContext) -> None:
    """Stage: SNV calls -> per-sample vcf2maf (writes maf/*.maf), optionally
    followed by per-sample SigMA mutational-signature analysis.

    Depends only on WalkContext, independent of _walk_process_cnv/_fusion.
    """
    temporary = None
    if ctx.input_folder_snv.exists() and ctx.vcf_type not in ["cnv", "fus", "tab"]:
        logger.info("Managing SNV files...")
        s_id_path_snv = snv_type_from_folder(
            ctx.input_folder_snv, ctx.case_folder_arr, ctx.output_folder)

        logger.info("Checking maf folder...")
        maf_path = Path(ctx.output_folder) / "maf"
        if maf_path.is_dir() and any(f.suffix == ".maf" for f in maf_path.iterdir()):
            logger.info("A non empty maf folder already exists!")

        if not ctx.resume:
            if "d" in ctx.filters:
                logger.info("Filtering out VCFs with dots in ALT column")
                s_id_path_snv = vcf_filtering(
                    s_id_path_snv, ctx.output_folder, output_filtered)

            temporary = create_random_name_folder(ctx.output_folder)
            sigma_results = []
            for k, v in s_id_path_snv.items():
                cl = vcf2maf_constructor(v, temporary, ctx.output_folder)
                run_vcf2maf(cl, k)

                if ctx.sigma:
                    sigma_results.append(_run_sigma_for_snv_sample(ctx, cl, k))

            if ctx.sigma:
                _write_sigma_intermediate(ctx.output_folder, sigma_results)

    logger.info("Clearing scratch folder...")
    clear_scratch(temporary)


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

    try:
        maf_df = pd.read_csv(maf_out_path, sep="\t", dtype=object)
    except Exception:
        maf_df = pd.read_csv(maf_out_path, sep="\t", dtype=object, skiprows=1)

    oncotree_code = ""
    match = ctx.clin_file.loc[ctx.clin_file["SAMPLE_ID"].astype(str) == str(sample_id)]
    if not match.empty and "ONCOTREE_CODE" in match.columns:
        oncotree_code = str(match.iloc[0]["ONCOTREE_CODE"])

    if not oncotree_code or oncotree_code.lower() == "nan":
        logger.warning(
            f"Sample {sample_id}: no ONCOTREE_CODE found in sample.tsv - "
            "SigMA's tumor_type mapping will fall through to its "
            "fallback/'other' handling (see get_sigma_call_params).")

    sigma_dir = Path(ctx.output_folder) / "intermediate" / "sigma"
    return sigma_runner.run_sigma_for_sample(
        maf_df, sample_id, oncotree_code, sigma_dir)


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

    combined_output_folder = Path(ctx.input_folder) / "CombinedOutput"
    if (
        combined_output_folder.exists()
        and any(f.is_file() for f in combined_output_folder.iterdir())):
        logger.info("Getting Fusions infos from CombinedOutput...")
        thr_fus = config.get("FUSION", "THRESHOLD_FUSION")
        combined_dict = get_combined_variant_output_from_folder(
            ctx.input_folder, ctx.clin_file, isinputfile)
        fill_fusion_from_combined(fusion_table_file, combined_dict, thr_fus)

    elif fusion_folder.exists() and any(fusion_folder.iterdir()):
        fusion_files = [f for f in fusion_folder.iterdir() if f.suffix == ".tsv"]
        if fusion_files:
            logger.info(f"Getting Fusions infos from {fusion_files[0].name} file.")
            fill_fusion_from_temp(
                ctx.input_folder, fusion_table_file, ctx.clin_file, fusion_files)

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
            os.system(f"mv {fusion_table_file_out} {fusion_table_file}")

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
    # found but all got filtered out by the OncoKB filter above".
    if fusion_table_file.exists():
        with fusion_table_file.open() as data_sv:
            if len(data_sv.readlines()) == 1:
                fusion_table_file.unlink()
                logger.warning(
                    "No fusions or splice variants found for this batch - "
                    "data_sv.txt removed.")


def _walk_write_clinical_tables(ctx: WalkContext) -> None:
    """Stage: data_clinical_patient.txt + data_clinical_sample.txt (MSI/TMB,
    exon, HRD info).

    Depends only on WalkContext, independent of _walk_process_cnv/_snv/_fusion
    - it reads sample.tsv/patient.tsv and CombinedOutput directly, not their
    output files.
    """
    table_dict_patient = get_table_from_folder(ctx.clin_sample_path)
    logger.info("Writing clinical files...")

    if Path(ctx.patient_tsv).name.strip() != "":
        logger.info("Writing data_clinical_patient.txt file...")

        input_file_path = Path(ctx.input_folder) / "patient.tsv"
        data_clin_pat = pd.read_csv(input_file_path, sep="\t", header=0, dtype=str)

        data_clin_pat.columns = data_clin_pat.columns.str.upper()
        datapat_columns = list(data_clin_pat.columns)

        # Get headers from conf.ini if they're present
        conf_header_short = config.get("ClinicalPatient", "HEADER_PATIENT_SHORT")
        conf_header_long = config.get("ClinicalPatient", "HEADER_PATIENT_LONG")
        conf_header_type = config.get("ClinicalPatient", "HEADER_PATIENT_TYPE")

        # Add header's fifth row
        default_row = pd.DataFrame([datapat_columns], columns=datapat_columns)
        final_data_pat = pd.concat([default_row, data_clin_pat], ignore_index=True)

        # Add header's fourth row (1s)
        header_numbers = pd.DataFrame(
            [[1] * len(datapat_columns)], columns=datapat_columns)
        final_data_pat = pd.concat([header_numbers, final_data_pat], ignore_index=True)

        # Add header's third row (HEADER_PATIENT_TYPE)
        final_data_pat = add_header_patient_type(
            ctx.patient_tsv, datapat_columns, conf_header_type,
            final_data_pat)

        # Add header's second row (HEADER_PATIENT_LONG)
        final_data_pat = add_header_patient_long(
            ctx.patient_tsv, datapat_columns, conf_header_long,
            default_row, final_data_pat)

        # Add header's first row (HEADER_PATIENT_SHORT)
        final_data_pat = add_header_patient_short(
            ctx.patient_tsv, datapat_columns, conf_header_short,
            default_row, final_data_pat)

        final_data_pat.loc[0:3, "PATIENT_ID"] = final_data_pat.loc[
            0:3, "PATIENT_ID"].apply(lambda x: f"#{x}")

        data_clin_txt = Path(ctx.output_folder) / "data_clinical_patient.txt"
        final_data_pat.to_csv(data_clin_txt, sep="\t", index=False, header=False)

    else:
        write_default_clinical_patient(ctx.output_folder, table_dict_patient)

    file_input_sample = pd.read_csv(
        ctx.clin_sample_path, sep="\t", index_col=False, dtype=str)

    msi_thr = config.get("MSI", "THRESHOLD_MSI")
    tmb_thr = ast.literal_eval(config.get("TMB", "THRESHOLD_TMB"))

    combined_output = Path(ctx.input_folder) / "CombinedOutput"
    if combined_output.exists() and len(list(combined_output.iterdir())) > 0:
        msi_sites_thr = config.get("MSI", "THRESHOLD_SITES")

        combined_dict = get_combined_variant_output_from_folder(
            ctx.input_folder, ctx.clin_file, isinputfile)
        new_table_dict_patient = fill_from_combined(
            combined_dict, table_dict_patient,
            msi_sites_thr, msi_thr, tmb_thr)
    else:
        new_table_dict_patient = fill_from_file(
            table_dict_patient, file_input_sample, msi_thr, tmb_thr)
        combined_dict = {}

    write_clinical_sample(ctx.clin_sample_path, ctx.output_folder, new_table_dict_patient,
                           combined_output_used=bool(combined_dict))

    update_data_clinical_with_exon_info(combined_output, combined_dict, ctx.output_folder)

    update_data_clinical_with_hrd_info(combined_output, combined_dict, ctx.output_folder)


def walk_folder(
    input_path: list,
    multiple: bool,
    output_folder: str,
    oncokb: bool,
    cancer: str,
    overwrite_output: bool = False,
    resume: bool = False,
    vcf_type: str | None = None,
    filters: str = "",
    sigma: bool = False) -> tuple:
    """Process input files/folders for SNV, CNV, fusions, and prepare output.

    Thin orchestrator over four independent stages - see WalkContext,
    _walk_setup, _walk_process_cnv, _walk_process_snv, _walk_process_fusion
    and _walk_write_clinical_tables. CNV, SNV, fusion/splice and the
    clinical tables each depend only on _walk_setup()'s output, never on
    each other, which is what makes them safe to expose as independent
    Snakemake rules (see Snakefile) instead of one opaque shell-out.

    Args:
        input_path (list): List with input path(s), either folder(s) or file(s).
        multiple (bool): Whether multiple samples per file are expected.
        output_folder (str): Path to the output directory.
        oncokb (bool): Whether to enable OncoKB annotation.
        cancer (str): Cancer ID used for validations and annotations.
        overwrite_output (bool, optional): Overwrite output folder if exists.
        Defaults to False.
        resume (bool, optional): Resume from previous run if True. Defaults to False.
        vcf_type (str | None, optional): Type of VCF to process (snv, cnv, fus, tab).
        Defaults to None.
        filters (str, optional): Filters to apply on VCF data. Defaults to "".
        sigma (bool, optional): Whether to run SigMA mutational-signature
        analysis per sample (see sigma_runner.run_sigma_for_sample). Defaults
        to False. All of SigMA's own parameters live in conf.ini's [SigMA]
        section - this is only the on/off switch.

    Returns:
        tuple: Returns output folder path, input path or file, and fusion TSV path.

    """
    ctx = _walk_setup(
        input_path, multiple, output_folder, oncokb, cancer,
        overwrite_output, resume, vcf_type, filters, sigma)

    _walk_process_cnv(ctx)
    _walk_process_snv(ctx)
    _walk_process_fusion(ctx)
    _walk_write_clinical_tables(ctx)

    logger.success("Walk script completed!\n")

    return ctx.output_folder, ctx.input_path, ctx.fusion_tsv
