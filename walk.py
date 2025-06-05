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
from configparser import ConfigParser
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

import tsv
import vcf2tab_cnv
import vcf_filter
from filter_clinvar import check_bool, filter_oncokb
from versioning import get_newest_version, get_version_list

config = ConfigParser()
config_file = config.read("conf.ini")

VCF2MAF = config.get("Paths", "VCF2MAF")
REF_FASTA = config.get("Paths", "REF_FASTA")
VEP_PATH = config.get("Paths", "VEP_PATH")
VEP_DATA = config.get("Paths", "VEP_DATA")
CLINV = config.get("Paths", "CLINV")
CNA = ast.literal_eval(config.get("Cna", "HEADER_CNV"))
PLOIDY = int(config.get("Cna", "PLOIDY"))
ONCOKB_FILTER = ast.literal_eval(config.get("Filters", "ONCOKB_FILTER"))

output_filtered = "snv_filtered"
tmp = "scratch"


def create_random_name_folder() -> str:
    """Create a temporary folder with a random name.

    Returns:
        str: Path to the created temporary folder.

    """
    folder_name = "".join(
        secrets.choice(string.ascii_lowercase + string.digits) for _ in range(10))
    tmp = "scratch"
    temporary = Path(tmp) / folder_name

    try:
        temporary.mkdir()
    except FileNotFoundError as err:
        logger.critical(f"Scratch folder '{tmp}' not found!")
        msg = "Error in create_random_name_folder: exiting from walk script!"
        raise(FileNotFoundError(msg)) from err
    except Exception as err:
        logger.critical("Something went wrong while creating the vep tmp folder")
        msg = "Error in create_random_name_folder: exiting from walk script!"
        raise(Exception(msg)) from err
    return(str(temporary))


def clear_scratch() -> None:
    """Remove all directories inside the `tmp` directory.

    Returns:
        None

    """
    tmp = "scratch"
    for root, dirs, _ in os.walk(tmp):
        for directory in dirs:
            to_rem=Path(root) / directory
            shutil.rmtree(to_rem)


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
        sample=cnv_vcf.replace("vcf", "bam")
    return sample


def reshape_cna(my_input: str,
                cna_df_path: str,
                cancer: str,
                output_dir: str) -> str:
    """Reshape and annotate CNA data and add ONCOTREE code.

    Args:
        my_input (str): Path to input file or directory.
        cna_df_path (str): Path to CNA dataframe.
        cancer (str): ONCOTREE code.
        output_dir (str): Output directory.

    Returns:
        str: Path to the intermediate annotated file.

    """
    my_input_path = Path(my_input)
    output_path = Path(output_dir) / "temp_cna.txt"
    if not os.path.is_file(my_input_path):
        input_file = pd.read_csv(my_input_path / "sample.tsv", sep="\t")
    else:
        input_file = pd.read_csv(my_input_path, sep="\t")

    cna_df = pd.read_csv(cna_df_path, sep="\t")

    cna_df=cna_df.rename({"ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"},
                         axis=1)
    input_file=input_file.rename({"SampleID":"Tumor_Sample_Barcode"}, axis=1)

    if "ONCOTREE_CODE" not in input_file.columns:
        input_file["ONCOTREE_CODE"] = cancer

    input_file["Tumor_Sample_Barcode"] = input_file["Tumor_Sample_Barcode"] + ".cnv.bam"

    annotate = cna_df[[
        "Tumor_Sample_Barcode", "Hugo_Symbol", "discrete",
        "Copy_Number_Alteration"]].merge(
    input_file[["Tumor_Sample_Barcode", "ONCOTREE_CODE"]],
    on="Tumor_Sample_Barcode",
)

    return output_path.name


def annotate_cna(path_cna: str, output_folder: str) -> None:
    """Annotate CNA file using OncoKB and filters oncogenic alterations.

    Args:
        path_cna (str): Path to CNA input file.
        output_folder (str): Output directory for results.

    Returns:
        None

    """
    out = path_cna.replace(".txt", "2.txt")
    oncokb_key = config.get("OncoKB", "ONCOKB")
    cmd = [
        "python3", "./oncokb-annotator/CnaAnnotator.py",
        "-i", path_cna,
        "-o", out,
        "-f", "individual",
        "-b", oncokb_key]
    subprocess.run(cmd, check=True)

    cna = pd.read_csv(out, sep="\t", dtype={"Copy_Number_Alteration":int})
    cna = cna[cna["ONCOGENIC"].isin(["Oncogenic", "Likely Oncogenic"])]

    data_cna = cna.pivot_table(
        index="Hugo_Symbol",
        columns="Tumor_Sample_Barcode",
        values="Copy_Number_Alteration",
        fill_value=0)
    outpath=Path(output_folder) / "data_cna.txt"
    data_cna.to_csv(outpath, index=True, sep="\t")


def cnv_type_from_folder(input_path: str,
                         cnv_vcf_files: list,
                         output_folder: str,
                         oncokb: bool,
                         cancer: str,
                         multiple: bool) -> dict:
    """Process CNV, converting them in CNA tables and performing annotation.

    Args:
    input_path : str
        Path to the input directory or the sample TSV file (if CNVkit=True).
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

    Returns:
    dict
        A mapping from sample ID to VCF file path used in processing.

    """
    c = 0
    sid_path = {}

    for case_folder in cnv_vcf_files:
        mode = "a" if Path("data_cna_hg19.seg").exists() else "w"
        try:
            cnv_vcf = case_folder
            sample_id = get_sample_id_from_cnv(case_folder)

            if sample_id in sid_path:
                with Path("sampleID_dup.log").open("w") as dup_path:
                    dup_path.write(sample_id + "\t" + "cnv_vcf")
            else:
                if multiple:
                    sid_path[sample_id] = Path(
                        input_path) / "CNV" / "single_sample_vcf" / cnv_vcf
                else:
                    sid_path[sample_id] = Path(input_path) / "CNV" / cnv_vcf

                vcf2tab_cnv.vcf_to_table(
                    sid_path[sample_id], Path(
                        output_folder) / "data_cna_hg19.seg",
                    sample_id, mode)
                vcf2tab_cnv.vcf_to_table_fc(
                    sid_path[sample_id], Path(
                        output_folder) / "data_cna_hg19.seg.fc.txt",
                    sample_id, mode)

        except Exception:
            with Path("noParsed_cnv.log").open("w") as log_noparsed:
                log_noparsed.write("[WARNING] " + case_folder + "\n")
                log_noparsed.close()

        c = c + 1

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
            Path(output_folder) / "data_cna_hg19.seg.fc.txt", sep="\t", header=0)

        df_table=df_table.rename(columns={
            "discrete":"Copy_Number_Alteration",
            "ID":"Tumor_Sample_Barcode",
            "gene":"Hugo_Symbol"})
        df_table_filt = df_table[
            df_table["Copy_Number_Alteration"].isin([-2,2])]

        cnv_kit = config.get("Cna", "CNVkit")
        cnv_kit = check_bool(cnv_kit)

        if cnv_kit:
            if not Path(input_path).is_file():
                input_file = pd.read_csv(
                    Path(input_path) / "sample.tsv", sep="\t")
            else:
                input_file = pd.read_csv(input_path, sep="\t")

            if "TC" not in input_file.columns:
                input_file["TC"] = np.nan

            if len(input_file[input_file["TC"].isna()])>0:
                nan_sbj = input_file[input_file["TC"].isna()]
                nan_sbj = list(nan_sbj["SAMPLE_ID"])
                logger.warning(
                    f"Some subject have NaN TC in tsv input file: {nan_sbj}!")

            if "ONCOTREE_CODE" not in input_file.columns:
                input_file["ONCOTREE_CODE"] = cancer

            input_file["Tumor_Sample_Barcode"] = input_file["SAMPLE_ID"] #+ ".cnv.bam"

            annotate = df_table_filt[[
                "Tumor_Sample_Barcode", "Hugo_Symbol",
                "seg.mean", "Copy_Number_Alteration"]].merge(
                    input_file[["Tumor_Sample_Barcode",
                                "ONCOTREE_CODE", "TC"]],
                    on="Tumor_Sample_Barcode")

            temppath = Path(output_folder) / "temp_cna_toannotate.txt"
            annotate.to_csv(temppath, index=False, sep="\t")

            if oncokb:
                out = temppath.name.replace(
                    "toannotate.txt", "annotated.txt")
                oncokb_key = config.get("OncoKB", "ONCOKB")
                cmd = [
                    "python3", "./oncokb-annotator/CnaAnnotator.py",
                    "-i", str(temppath),
                    "-o", out,
                    "-f", "individual",
                    "-b", oncokb_key]
                subprocess.run(cmd, check=True)

                name = "annotated_oncokb_CNA_ndiscrete.txt"
                cna = pd.read_csv(out, sep="\t",
                                  dtype={"Copy_Number_Alteration":int})
                cna = cna[cna["ONCOGENIC"].isin(ONCOKB_FILTER)]
            else:
                out = temppath
                name = "CNA_ndiscrete.txt"
                cna = pd.read_csv(out, sep="\t",
                                  dtype={"Copy_Number_Alteration":int})

            logger.info("Analyzing cna sample(s)")
            for _, row in cna.iterrows():
                try:
                    tc = int(
                        input_file[input_file[
                            "Tumor_Sample_Barcode"] == row[
                                "Tumor_Sample_Barcode"]]["TC"])
                except ValueError:
                    continue
                except Exception as err:
                    msg = "Something went wrong while reading TC!"
                    raise(Exception(msg)) from err

                purity = tc / 100
                copy_nums = np.arange(6)
                c = 2 ** (
                    np.log2((1 - purity)
                            + purity * (copy_nums + .5)
                            / PLOIDY ))

            # CNVkit filter
            if not cna["TC"].isna().all():
                cna["Copy_Number_Alteration"]=0
                cna.loc[(cna["seg.mean"]<c[0]
                         ), "Copy_Number_Alteration"]=-2
                cna.loc[(cna["seg.mean"]>=c[0]
                         )&(cna["seg.mean"]<c[1]
                            ), "Copy_Number_Alteration"]=-1
                cna.loc[(cna["seg.mean"]>=c[1]
                         )&(cna["seg.mean"]<c[3]
                            ), "Copy_Number_Alteration"]=0
                cna.loc[(cna["seg.mean"]>=c[3]
                         )&(cna["seg.mean"]<c[5]
                            ), "Copy_Number_Alteration"]=1
                cna.loc[cna["seg.mean"]>=c[5],
                        "Copy_Number_Alteration"]=2

            else:
                logger.warning("TC column is empty or does not exist in sample.tsv! "
                "This column is required when CNVkit = True! If TC values are not "
                "available the setting of CNVkit = False is recommended")
                return sid_path

            cna.to_csv(Path(output_folder) / name,
                       index=True, sep="\t")

            cna["Tumor_Sample_Barcode"] = cna[
                "Tumor_Sample_Barcode"].str.replace(
                    ".cnv.bam", "", regex=False)
            data_cna=cna.pivot_table(
                index="Hugo_Symbol",
                columns="Tumor_Sample_Barcode",
                values="Copy_Number_Alteration", fill_value=0)
            data_cna.to_csv(Path(output_folder) / "data_cna.txt",
                            index=True, sep="\t")

        else:
            df_table_filt = df_table_filt.copy()
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


def table_to_dict(df: pd.DataFrame) -> dict:
    """Convert a DataFrame into a dictionary grouped by the 'ID' column.

    Each key is a sample ID, and each value is a list
    of tuples containing:
    (chrom, _2, _3, _4, _5, gene, discrete).

    Args:
        df (pd.DataFrame): Input DataFrame with required columns.

    Returns:
        dict: Dictionary mapping sample IDs
        to lists of tuples with SNV data.

    """
    result = {}
    for row in df.itertuples(index=False):
        row_values = (row.chrom,
                      row._2, row._3, row._4,
                      row._5, row.gene, row.discrete)
        if row.ID not in result:
            result[row.ID] = []
        result[row.ID].append(row_values)
    return result


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
        sample = snv_vcf.replace("vcf", "bam")
    return sample


def snv_type_from_folder(input_pat: str,
                         snv_vcf_files: list) -> dict:
    """Map sample IDs to their full SNV paths, handling duplicates.

    Args:
        input_pat (str): Path to the folder containing SNV files.
        snv_vcf_files (list): List of SNV VCF filenames.

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
                with Path("sampleID_dup.log").open("w") as dup:
                    dup.write(sample_id + "\t" + "snv_vcf")
            else:
                sid_path[sample_id] = str(Path(input_pat) / snv_vcf)
        except Exception:
            with Path("noParsed_snv.log").open("a") as log_noparsed:
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
        vcf_filter.main(v, vcf_filtered.name)
        sid_path_filtered[k] = vcf_filtered.name
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
    cmd = "vcf-query -l " + v
    try:
        tum_id = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
    except Exception:
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
) -> None:
    """Write the `data_clinical_sample.txt` file by merging sample and metrics data.

    The function reads the clinical sample file, optionally merges it with MSI/TMB
    metrics from a provided table, validates headers, and writes the output file
    with appropriate headers as required by cBioPortal format.

    Args:
        clin_samp_path (str): Path to the input clinical sample TSV file.
        output_folder (str): Directory where the final file will be saved.
        table_dict (dict): Dictionary containing metrics (MSI, TMB) by sample ID.

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
        try:
            msi_notna = data_clin_samp["MSI"].notna().any()
            tmb_notna = data_clin_samp["TMB"].notna().any()

            if msi_notna or tmb_notna:
                msi_mismatch = (data_clin_samp["MSI"] != combout_df["MSI"]).any()
                tmb_mismatch = (data_clin_samp["TMB"] != combout_df["TMB"]).any()

                if msi_mismatch or tmb_mismatch:
                    logger.warning(
                        "MSI and/or TMB values are reported in sample.tsv and "
                        "CombinedOutput but they do not match! CombinedOutput "
                        "values were selected by default")
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

    basic_columns = ["SAMPLE_ID", "PATIENT_ID", "MSI", "TMB", "MSI_THR", "TMB_THR"]
    other_columns = [*final_data_sample.columns[2:-4]]
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
        header_row[dataclin_columns.index("MSI")] = "NUMBER"
        header_row[dataclin_columns.index("TMB")] = "NUMBER"
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
        raise RuntimeError(msg) from err

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
        sys.exit()
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
    output_folder: str, file: str, copy_to: str, sample_id: str) -> None:
    """Copy an input file to a temporary folder or log if missing.

    Args:
        output_folder (str): Base output directory.
        file (str): File path to copy.
        copy_to (str): Subfolder name.
        sample_id (str): ID of the sample (for logging).

    """
    file_path = Path(file)
    destination = Path(output_folder) / "temp" / copy_to

    if file_path.exists():
        os.system(f"cp {file_path} {destination}")
    elif not file:
        logger.warning(
            f"No final_path set in conf.ini for sample {sample_id}'s {copy_to}!")
    else:
        logger.warning(f"{file_path} not found")


def check_folders(
    output_folder: str,
    snv_path: str,
    cnv_path: str,
    combout: str,
    sample_id: str) -> None:
    """Check presence of SNV, CNV, and CombinedOutput files for a sample.

    Args:
        output_folder (str): Base output directory.
        snv_path (str): Path to SNV file.
        cnv_path (str): Path to CNV file.
        combout (str): Path to CombinedOutput file.
        sample_id (str): ID of the sample.

    """
    logger.info("Verifying the compilation of the conf.ini file for "
    f"sample {sample_id}")

    check_input_file(output_folder, snv_path, "SNV", sample_id)
    check_input_file(output_folder, cnv_path, "CNV", sample_id)
    check_input_file(output_folder, combout, "CombinedOutput", sample_id)


def transform_input(
    tsv: str,
    clin_pzt: str,
    fusion_tsv: str,
    output_folder: str,
    multiple: bool) -> str:
    """Prepare temp folder structure and copy necessary files.

    Args:
        tsv (str): Sample TSV file path.
        clin_pzt (str): Clinical patient TSV file path.
        fusion_tsv (str): Fusion TSV file path.
        output_folder (str): Output directory.
        multiple (bool): Whether input is from multi-sample VCF.

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

        check_folders(output_folder, snv_path, cnv_path, combout)

    else:
        tsv_file = pd.read_csv(tsv, sep="\t", dtype="string", keep_default_na=False)

        for _,row in tsv_file.iterrows():
            sample_id = row["SAMPLE_ID"]
            snv_path = row["snv_path"]
            cnv_path = row["cnv_path"]
            combout = row["comb_path"]

            check_folders(output_folder, snv_path, cnv_path, combout, sample_id)

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
            ff = pd.read_csv(fusion_input, sep="\t")

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
                if (fus.Sample_Id in clin_file["SAMPLE_ID"].to_numpy() and
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
                    # site1_Chromosome = fus["Site1_Chromosome"]
                    # site2_Chromosome = fus["Site2_Chromosome"]
                    # site1_Position = fus["Site1_Position"]
                    # site2_Position = fus["Site2_Position"]

                    if eval("int(fus['Normal_Paired_End_Read_Count'])" + thr_fus):
                        fusion_table.write(
                            k + "\tSOMATIC\tFUSION\t" +
                            str(site1_hugo_symbol) + "\t" +
                            str(site2_hugo_symbol) + "\t" +
                            fus["Normal_Paired_End_Read_Count"] + "\t" +
                            fus["Event_Info"] + " Fusion\tYes\n")


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


def fill_from_file(
    table_dict_patient: dict[str, list],
    file_input_clinical: pd.DataFrame,
    msi_thr: str,
    tmb_thr: dict[str, str]) -> dict[str, list]:
    """Populate clinical table dictionary with MSI and TMB values and statuses.

    Args:
        table_dict_patient (dict): Patient clinical data dictionary.
        file_input_clinical (pd.DataFrame): Clinical input with SAMPLE_ID, MSI, TMB.
        msi_thr (str): Threshold condition string for MSI (used with eval).
        tmb_thr (dict): Map of TMB categories to threshold conditions.

    Returns:
        dict: Updated clinical data dictionary.

    """
    for k, m, t in zip(
        file_input_clinical["SAMPLE_ID"],
        file_input_clinical["MSI"],
        file_input_clinical["TMB"]):
        table_dict_patient[k].append(m)
        table_dict_patient[k].append(t)

        if np.isnan(float(m)):
            table_dict_patient[k].append("NA")
        elif eval("float(m)" + msi_thr):
            table_dict_patient[k].append("Stable")
        else:
            table_dict_patient[k].append("Unstable")

        if not np.isnan(float(t)):
            for _k, _v in tmb_thr.items():
                if eval("float(t)" + _v):
                    table_dict_patient[k].append(_k)
                    break
        else:
            table_dict_patient[k].append("NA")
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
            tmv_msi = tsv.get_msi_tmb(Path(v))
        except Exception:
            logger.error("Something went wrong!")

        if (
            tmv_msi["MSI"][0][1] != "NA" and
            eval("float(tmv_msi['MSI'][0][1])" + msi_sites_thr)):
            table_dict_patient[k].append(tmv_msi["MSI"][1][1])
        else:
            table_dict_patient[k].append("NA")

        table_dict_patient[k].append(tmv_msi["TMB_Total"])

        if (
            tmv_msi["MSI"][0][1] != "NA" and
            tmv_msi["MSI"][1][1] != "NA" and
            table_dict_patient[k][1] != "NA"):
            if eval("float(tmv_msi['MSI'][1][1])" + msi_thr):
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
                    f"The TMB value {tmv_msi['TMB_Total']} is not within the conf.ini "
                    f"thresholds {list(tmb.values())}. For this value, the TMB_THR "
                    "will be set as 'Out of threshold ranges'.")
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
        sys.exit()


def walk_folder(
    input_path: list,
    multiple: bool,
    output_folder: str,
    oncokb: bool,
    cancer: str,
    overwrite_output: bool = False,
    resume: bool = False,
    vcf_type: str | None = None,
    filters: str = "") -> tuple:
    """Process input files/folders for SNV, CNV, fusions, and prepare output.

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

    Returns:
        tuple: Returns output folder path, input path or file, and fusion TSV path.

    """
    logger.info("Starting walk_folder script:")
    logger.info(
        f"walk_folder args [input:{input_path}, output_folder:{output_folder}, "
        f"overwrite:{overwrite_output}, resume:{resume}, vcf_type:{vcf_type}, "
        f"filters:{filters}, multiple:{multiple}]")

    config.read("conf.ini")

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


    ###############################
    ###      OUTPUT FOLDER      ###
    ###############################

    if not resume or not (Path(output_folder) / "temp").exists():
        output_folder = create_folder(output_folder, overwrite_output, resume)

    if not isinputfile:
        input_folder = input_path[0]
        input_path, patient_tsv, fusion_tsv = input_extraction_folder(input_folder)
        check_multiple_folder(input_folder, multiple)

    elif isinputfile:
        input_path, patient_tsv, fusion_tsv = input_extraction_file(input_path)
        check_multiple_file(input_path, multiple)
        input_folder = transform_input(
            input_path, patient_tsv, fusion_tsv, output_folder, multiple)

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

    maf_path = output_folder / "maf"
    maf_zip_path = output_folder / "maf.zip"
    clin_sample_path = input_folder / "sample.tsv"

    try:
        clin_file = pd.read_csv(clin_sample_path, sep="\t", dtype=str)
    except Exception as err:
        logger.critical(f"Something went wrong while reading {clin_sample_path}!")
        msg = "Error in reading the input file! Please check again."
        raise OSError(msg) from err

    if resume:
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

    if len(list(input_folder_snv.iterdir())) == 0 and vcf_type is None:
        vcf_type = "cnv"
        logger.info("SNV path was empty, the analysis will exclude SNV")

    if len(list(input_folder_cnv.iterdir())) == 0 and vcf_type is None:
        vcf_type = "snv"
        logger.info("CNV path was empty, the analysis will exclude CNV")

    if input_folder_cnv.exists() and vcf_type not in ["snv", "fus", "tab"]:
        if multiple:
            multivcf = next(
                (i for i in input_folder_cnv.iterdir() if i.name.endswith(".vcf")),
                None)
            if multivcf:
                extract_multiple_cnv(multivcf, input_folder_cnv)
                input_folder_cnv = input_folder_cnv / "single_sample_vcf"
        logger.info("Checking CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(input_folder_cnv)
        logger.info("Everything ok!")

    input_folder_snv = Path(input_folder_snv)

    if input_folder_snv.exists() and vcf_type not in ["cnv", "fus", "tab"]:
        tmp = Path("scratch")
        tmp.mkdir(parents=True, exist_ok=True)

        if multiple:
            multivcf = next(f for f in input_folder_snv.iterdir() if f.suffix == ".vcf")
            extract_multiple_snv(multivcf, input_folder_snv)
            input_folder_snv = input_folder_snv / "single_sample_vcf"

        logger.info("Checking SNV files...")

    case_folder_arr = get_snv_from_folder(input_folder_snv)
    logger.info("Everything ok!")

    ###############################
    ###       SNV AND CNV       ###
    ###############################

    input_folder_cnv = Path(input_folder_cnv)
    input_folder_snv = Path(input_folder_snv)

    if input_folder_cnv.exists() and vcf_type not in ["snv", "fus", "tab"]:
        logger.info("Managing CNV files...")
        sID_path_cnv = cnv_type_from_folder(input_folder, case_folder_arr_cnv, output_folder, oncokb, cancer, multiple)

    if input_folder_snv.exists() and vcf_type not in ["cnv", "fus", "tab"]:
        logger.info("Managing SNV files...")
        s_id_path_snv = snv_type_from_folder(input_folder_snv, case_folder_arr)

        logger.info("Checking maf folder...")
        maf_path = Path(output_folder) / "maf"
        if maf_path.is_dir() and any(f.suffix == ".maf" for f in maf_path.iterdir()):
            logger.info("A non empty maf folder already exists!")

        if not resume:
            if "d" in filters:
                logger.info("Filtering out VCFs with dots in ALT column")
                s_id_path_snv = vcf_filtering(
                s_id_path_snv, output_folder, output_filtered)

            temporary = create_random_name_folder()
            for k, v in s_id_path_snv.items():
                cl = vcf2maf_constructor(v, temporary, output_folder)
                run_vcf2maf(cl, k)

    logger.info("Clearing scratch folder...")
    clear_scratch()


    ###############################
    ###       GET FUSION        ###
    ###############################
    if vcf_type not in ["cnv","snv","tab"]:

        fusion_table_file = Path(output_folder) / "data_sv.txt"
        fusion_folder = Path(input_folder) / "FUSIONS"

        combined_output_folder = Path(input_folder) / "CombinedOutput"
        if (
            combined_output_folder.exists()
            and any(f.is_file() for f in combined_output_folder.iterdir())):
            logger.info("Getting Fusions infos from CombinedOutput...")
            thr_fus = config.get("FUSION", "THRESHOLD_FUSION")
            combined_dict = get_combined_variant_output_from_folder(
                input_folder, clin_file, isinputfile)
            fill_fusion_from_combined(fusion_table_file, combined_dict, thr_fus)

        elif fusion_folder.exists() and any(fusion_folder.iterdir()):
            fusion_files = [f for f in fusion_folder.iterdir() if f.suffix == ".tsv"]
            if fusion_files:
                logger.info(f"Getting Fusions infos from {fusion_files[0].name} file.")
                fill_fusion_from_temp(
                    input_folder, fusion_table_file, clin_file, fusion_files)

        if fusion_table_file.exists():
            with fusion_table_file.open() as data_sv:
                all_data_sv = data_sv.readlines()
                if len(all_data_sv) == 1:
                    fusion_table_file.unlink()
                    logger.warning("data.sv is empty. File removed.")

        if oncokb and fusion_table_file.exists():
            data_sv = pd.read_csv(fusion_table_file, sep="\t")
            input_file = pd.read_csv(clin_sample_path, sep="\t")
            fusion_table_file_out = annotate_fusion(
                cancer, fusion_table_file, data_sv, input_file)

            if "o" in filters:
                fus_file = pd.read_csv(fusion_table_file_out, sep="\t")
                fus_file = filter_oncokb(fus_file)
                fus_file.to_csv(fusion_table_file_out, index=False, sep="\t")

            data_sv_tmp = pd.read_csv(fusion_table_file_out, sep="\t")
            with contextlib.suppress(KeyError):
                data_sv_tmp = data_sv_tmp.drop(["SAMPLE_ID", "ONCOTREE_CODE"], axis=1)

            data_sv_tmp.to_csv(fusion_table_file_out, index=False, sep="\t")
            os.system(f"mv {fusion_table_file_out} {fusion_table_file}")


    ##############################
    ##       MAKES TABLE       ###
    ##############################

    table_dict_patient = get_table_from_folder(clin_sample_path)
    logger.info("Writing clinical files...")

    if Path(patient_tsv).name.strip() != "":
        logger.info("Writing data_clinical_patient.txt file...")

        input_file_path = Path(input_folder) / "patient.tsv"
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
            patient_tsv, datapat_columns, conf_header_type,
            final_data_pat)

        # Add header's second row (HEADER_PATIENT_LONG)
        final_data_pat = add_header_patient_long(
            patient_tsv, datapat_columns, conf_header_long,
            default_row, final_data_pat)

        # Add header's first row (HEADER_PATIENT_SHORT)
        final_data_pat = add_header_patient_short(
            patient_tsv, datapat_columns, conf_header_short,
            default_row, final_data_pat)

        final_data_pat.loc[0:3, "PATIENT_ID"] = final_data_pat.loc[
            0:3, "PATIENT_ID"].apply(lambda x: f"#{x}")

        data_clin_txt = Path(output_folder) / "data_clinical_patient.txt"
        final_data_pat.to_csv(data_clin_txt, sep="\t", index=False, header=False)

    else:
        write_default_clinical_patient(output_folder, table_dict_patient)

    fileinputclinical = pd.read_csv(
        clin_sample_path, sep="\t", index_col=False, dtype=str)

    msi_thr = config.get("MSI", "THRESHOLD_MSI")
    tmb_thr = ast.literal_eval(config.get("TMB", "THRESHOLD_TMB"))

    combined_output = Path(input_folder) / "CombinedOutput"
    if combined_output.exists() and len(list(combined_output.iterdir())) > 0:
        msi_sites_thr = config.get("MSI", "THRESHOLD_SITES")

        combined_dict = get_combined_variant_output_from_folder(
            input_folder, clin_file, isinputfile)
        new_table_dict_patient = fill_from_combined(
            combined_dict, table_dict_patient,
            msi_sites_thr, msi_thr, tmb_thr)
    else:
        new_table_dict_patient = fill_from_file(
            table_dict_patient, fileinputclinical, msi_thr, tmb_thr)

    write_clinical_sample(clin_sample_path, output_folder, new_table_dict_patient)
    logger.success("Walk script completed!\n")

    return output_folder, input_path, fusion_tsv
