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

import ast
import os
import random
import shutil
import string
import subprocess
import sys
import zipfile
from configparser import ConfigParser

import numpy as np
import pandas as pd
from loguru import logger

import tsv
import vcf2tab_cnv
import vcf_filter
from filter_clinvar import check_bool, filter_OncoKB
from versioning import get_newest_version, get_version_list

config = ConfigParser()
configFile = config.read("conf.ini")

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


def create_random_name_folder():
    nome_cartella = "".join(random.choices(string.ascii_lowercase + string.digits, k=10))
    temporary = os.path.join(tmp, nome_cartella)
    try:
        os.mkdir(temporary)
    except FileNotFoundError:
        logger.critical(f"Scratch folder '{tmp}' not found!")
        raise(FileNotFoundError("Error in create_random_name_folder: exiting from walk script!"))
    except Exception:
        logger.critical("Something went wrong while creating the vep tmp folder")
        raise(Exception("Error in create_random_name_folder: exiting from walk script!"))
    return(temporary)


def clear_scratch():
    for root, dirs, files in os.walk(tmp):
        for dir in dirs:
            shutil.rmtree(os.path.join(root,dir))


def get_cnv_from_folder(inputFolderCNV):
    files = os.listdir(inputFolderCNV)
    cnv_vcf_files = [file for file in files if file.endswith("vcf")]
    return cnv_vcf_files


def get_sampleID_from_cnv(cnv_vcf):
    if "_CopyNumberVariants.vcf" in cnv_vcf:
        sample=cnv_vcf.replace("_CopyNumberVariants.vcf", ".bam")
    else:
        sample=cnv_vcf.replace("vcf", "bam")
    return sample


def reshape_cna(input, cna_df_path, cancer, output_dir):
    if not os.path.isfile(input):
        input_file = pd.read_csv(os.path.join(input, "sample.tsv"), sep="\t")
    else:
        input_file = pd.read_csv(input, sep="\t")

    cna_df = pd.read_csv(cna_df_path, sep="\t")

    cna_df.rename({"ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"}, inplace=True, axis=1)
    input_file.rename({"SampleID":"Tumor_Sample_Barcode"}, inplace=True, axis=1)

    if "ONCOTREE_CODE" not in input_file.columns:
        input_file["ONCOTREE_CODE"] = cancer

    input_file["Tumor_Sample_Barcode"] = input_file["Tumor_Sample_Barcode"] + ".cnv.bam"

    annotate = pd.merge(cna_df[["Tumor_Sample_Barcode", "Hugo_Symbol", "discrete", "Copy_Number_Alteration"]], \
                        input_file[["Tumor_Sample_Barcode", "ONCOTREE_CODE"]], on="Tumor_Sample_Barcode")
    annotate.to_csv(os.path.join(output_dir, "temp_cna.txt"), index=False, sep="\t")

    return os.path.join(output_dir, "temp_cna.txt")


def annotate_cna(path_cna, output_folder):
    out = path_cna.replace(".txt", "2.txt")
    os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {path_cna}\
                        -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")

    cna = pd.read_csv(out, sep="\t", dtype={"Copy_Number_Alteration":int})
    cna = cna[cna["ONCOGENIC"].isin(["Oncogenic", "Likely Oncogenic"])]

    data_cna = cna.pivot_table(index="Hugo_Symbol", columns="Tumor_Sample_Barcode", values="Copy_Number_Alteration", fill_value=0)
    data_cna.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")


def cnv_type_from_folder(input, cnv_vcf_files, output_folder, oncokb, cancer, multiple):
    c = 0
    sID_path = dict()
    for case_folder in cnv_vcf_files:
        if os.path.exists("data_cna_hg19.seg"):
            MODE = "a"
        else:
            MODE = "w"
        try:
            cnv_vcf = case_folder
            sampleID = get_sampleID_from_cnv(case_folder)

            if sampleID in sID_path:
                dup = open("sampleID_dup".log, "w")
                dup.write(sampleID + "\t" + "cnv_vcf")
                dup.close()
            else:
                if multiple:
                    sID_path[sampleID] = os.path.join(os.path.join(input, "CNV", "single_sample_vcf"), cnv_vcf)
                else:
                    sID_path[sampleID] = os.path.join(os.path.join(input, "CNV"), cnv_vcf)

                vcf2tab_cnv.vcf_to_table(sID_path[sampleID], os.path.join(output_folder, "data_cna_hg19.seg"), sampleID, MODE)
                vcf2tab_cnv.vcf_to_table_fc(sID_path[sampleID], os.path.join(output_folder, "data_cna_hg19.seg.fc.txt"), sampleID, MODE)

        except Exception:
            log_noparsed = open("noParsed_cnv.log", "a")
            log_noparsed.write("[WARNING] " + case_folder + "\n")
            log_noparsed.close()

        c = c + 1

    seg_path = os.path.join(output_folder, "data_cna_hg19.seg")
    segfc_path = os.path.join(output_folder, "data_cna_hg19.seg.fc.txt")

    check_data_cna(seg_path)
    check_data_cna(segfc_path)

    if os.path.exists(seg_path):
        logger.info("Writing data_cna_hg19.seg succefully completed!")

    if os.path.exists(segfc_path):
        logger.info("Writing data_cna_hg19.seg.fc.txt succefully completed!")

    ############################
    ### MANAGE DISCRETE TABLE ##
    ############################

        logger.info("Starting CNA evaluation (this step could take a while)...")
        df_table = pd.read_csv(os.path.join(output_folder, "data_cna_hg19.seg.fc.txt"), sep="\t", header=0)
        result = table_to_dict(df_table)

        df_table.rename(columns={"discrete":"Copy_Number_Alteration", "ID":"Tumor_Sample_Barcode", "gene":"Hugo_Symbol"}, inplace=True)
        df_table_filt = df_table[df_table["Copy_Number_Alteration"].isin([-2,2])]

        CNVkit = config.get("Cna", "CNVkit")
        CNVkit = check_bool(CNVkit)
        if CNVkit:
            if not os.path.isfile(input):
                input_file = pd.read_csv(os.path.join(input, "sample.tsv"), sep="\t")
            else:
                input_file = pd.read_csv(input, sep="\t")

            try:
                input_file["TC"]
            except:
                input_file["TC"] = np.nan

            if len(input_file[input_file["TC"].isna()])>0:
                nan_sbj = input_file[input_file["TC"].isna()]
                nan_sbj = list(nan_sbj["SAMPLE_ID"])
                logger.warning(f"Some subject have NaN TC in tsv input file: {nan_sbj}!")

            if "ONCOTREE_CODE" not in input_file.columns:
                input_file["ONCOTREE_CODE"] = cancer

            input_file["Tumor_Sample_Barcode"] = input_file["SAMPLE_ID"] #+ ".cnv.bam"

            annotate = pd.merge(df_table_filt[["Tumor_Sample_Barcode","Hugo_Symbol", "seg.mean", "Copy_Number_Alteration"]], \
                            input_file[["Tumor_Sample_Barcode","ONCOTREE_CODE","TC"]], on="Tumor_Sample_Barcode")
            temppath = os.path.join(output_folder, "temp_cna_toannotate.txt")
            annotate.to_csv(temppath, index=False, sep="\t")

            if oncokb:
                out=temppath.replace("toannotate.txt", "annotated.txt")
                os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {temppath}\
                                -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")
                name = "annotated_oncokb_CNA_ndiscrete.txt"
                cna = pd.read_csv(out, sep="\t", dtype={"Copy_Number_Alteration":int})
                cna = cna[cna["ONCOGENIC"].isin(ONCOKB_FILTER)]
            else:
                out = temppath
                name = "CNA_ndiscrete.txt"
                cna = pd.read_csv(out, sep="\t", dtype={"Copy_Number_Alteration":int})

            cna["ESCAT"]="Unmatched"
            df_table["ESCAT"]="Unmatched"

            logger.info("Analyzing cna sample(s)")
            for _, row in cna.iterrows():
                try:
                    tc = int(input_file[input_file["Tumor_Sample_Barcode"] == row["Tumor_Sample_Barcode"]]["TC"])
                except ValueError:
                    continue
                except Exception:
                    raise(Exception("Something went wrong while reading TC!"))

                purity = tc / 100
                copy_nums = np.arange(6)
                c = 2 ** (np.log2((1 - purity) + purity * (copy_nums + .5) / PLOIDY ))

                # ESCAT classification
                escat_class(df_table, input_file, row)

            # CNVkit filter
            if not cna["TC"].isnull().all():
                cna["Copy_Number_Alteration"]=0
                cna.loc[(cna["seg.mean"]<c[0]), "Copy_Number_Alteration"]=-2
                cna.loc[(cna["seg.mean"]>=c[0])&(cna["seg.mean"]<c[1]), "Copy_Number_Alteration"]=-1
                cna.loc[(cna["seg.mean"]>=c[1])&(cna["seg.mean"]<c[3]), "Copy_Number_Alteration"]=0
                cna.loc[(cna["seg.mean"]>=c[3])&(cna["seg.mean"]<c[5]), "Copy_Number_Alteration"]=1
                cna.loc[cna["seg.mean"]>=c[5], "Copy_Number_Alteration"]=2

            else:
                logger.warning("TC column is empty or does not exist in sample.tsv! This column is required when CNVkit = True! If TC values are not available the setting of CNVkit = False is recommended")
                return sID_path

            df_table.to_csv(os.path.join(output_folder, "data_cna_hg19_escat.seg.fc.txt"), index=True, sep="\t")
            cna.to_csv(os.path.join(output_folder, name), index=True, sep="\t")

            cna["Tumor_Sample_Barcode"] = cna["Tumor_Sample_Barcode"].str.replace(".cnv.bam", "", regex=False)
            data_cna=cna.pivot_table(index="Hugo_Symbol", columns="Tumor_Sample_Barcode", values="Copy_Number_Alteration", fill_value=0)
            data_cna.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")

        else:
            df_table_filt = df_table_filt.copy()
            df_table_filt.loc[:, "Tumor_Sample_Barcode"] = df_table_filt["Tumor_Sample_Barcode"].str.replace(".cnv.bam", "", regex=True)

            data_cna = df_table_filt.pivot_table(index="Hugo_Symbol", columns="Tumor_Sample_Barcode", values="Copy_Number_Alteration", fill_value=0).astype(int)
            if not data_cna.empty:
                data_cna.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")

        return sID_path


def escat_class(df_table, input_file, row):
    oncocode = input_file[input_file["Tumor_Sample_Barcode"] == row["Tumor_Sample_Barcode"]]["ONCOTREE_CODE"].values
    gene = df_table[(df_table["Tumor_Sample_Barcode"] == row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"] == row["Hugo_Symbol"])]

    if oncocode =="NSCLC":
        if (row["Hugo_Symbol"]=="MET" and row["seg.mean"]>1) or (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"

    elif oncocode=="BREAST":
        if (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
            row["ESCAT"]="IA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"

    elif oncocode=="BOWEL":
        if (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"

    elif oncocode=="PROSTATE":
        if (row["Hugo_Symbol"]=="BRCA1" and row["seg.mean"]<1) or (row["Hugo_Symbol"]=="BRCA2" and row["seg.mean"]<1):
            row["ESCAT"]="IA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"
        elif (row["Hugo_Symbol"]=="PTEN" and row["seg.mean"]<1):
            row["ESCAT"]="IIA"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIA"
        elif (row["Hugo_Symbol"]=="ATM" and row["seg.mean"]<1) or (row["Hugo_Symbol"]=="PALB2" and row["seg.mean"]<1):
            row["ESCAT"]="IIB"
            df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"


def table_to_dict(df):
    result = {}
    for row in df.itertuples(index=False):
        row_values = (row.chrom, row._2, row._3, row._4, row._5, row.gene, row.discrete)
        if row.ID not in result:
            result[row.ID] = []
        result[row.ID].append(row_values)
    return result


def get_snv_from_folder(inputFolderSNV):
    files = os.listdir(inputFolderSNV)
    snv_vcf_files = [file for file in files if file.endswith("vcf")]
    return snv_vcf_files


def get_sampleID_from_snv(snv_vcf):
    if "MergedSmallVariants.genome.vcf" in snv_vcf:
        sample = snv_vcf.replace("_MergedSmallVariants.genome.vcf", ".bam")
    else:
        sample = snv_vcf.replace("vcf", "bam")
    return sample


def snv_type_from_folder(input, snv_vcf_files):
    c = 0
    sID_path = dict()
    for case_folder in snv_vcf_files:
        try:
            snv_vcf= case_folder
            sampleID = get_sampleID_from_snv(case_folder)
            if sampleID in sID_path:
                dup = open("sampleID_dup".log, "w")
                dup.write(sampleID + "\t" + "snv_vcf")
                dup.close()
            else:
                sID_path[sampleID] = os.path.join(input, snv_vcf)
        except Exception:
            log_noparsed = open("noParsed_snv.log", "a")
            log_noparsed.write("[WARNING]" + case_folder + "\n")
            log_noparsed.close()
        c = c + 1
    return sID_path


def vcf_filtering(sID_path, output_folder, output_filtered):
    sID_path_filtered = dict()
    if output_filtered.strip() == "":
        output_filtered = "snv_filtered"
    os.makedirs(os.path.join(output_folder, output_filtered), exist_ok=True)
    for k, v in sID_path.items():
        _, vcf_file = os.path.split(v)
        out_filt = os.path.join(output_folder, output_filtered)
        vcf_filtered = os.path.join(out_filt, vcf_file.replace(".vcf","") + ".FILTERED.vcf")
        vcf_filter.main(v, vcf_filtered)
        sID_path_filtered[k] = vcf_filtered
    return sID_path_filtered


def vcf2maf_constructor(k, v, temporary, output_folder):
    CACHE = config.get("Paths", "CACHE")
    cmd = "vcf-query -l "+v
    try:
        tum_id = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
    except Exception:
        tum_id = ""

    if VCF2MAF == "" or REF_FASTA == "" or VEP_PATH == "" or VEP_DATA == "":
        logger.critical("[Paths] section in conf.ini is not correctly compiled. Please check again!")
        raise Exception ("Input error")

    cl = ["perl"]
    cl.append(VCF2MAF)
    cl.append("--input-vcf")
    cl.append(v)
    _, file_vcf = os.path.split(v)

    out_file = os.path.join(output_folder, os.path.join("maf", file_vcf + ".maf"))
    if not CLINV.strip() == "":
        cl.append("--vep-custom")
        cl.append(CLINV)
    else:
        logger.warning("CLINV section in [Paths] in conf.ini is not compiled. This step will be skipped")
    cl.append("--output-maf")
    cl.append(out_file)
    cl.append("--ref-fasta")
    cl.append(REF_FASTA)
    cl.append("--tmp-dir")
    cl.append(temporary)
    cl.append("--retain-fmt")
    cl.append("GT,GQ,AD,DP,VF,AF")
    cl.append("--vep-path")
    cl.append(VEP_PATH)
    cl.append("--vep-data")
    cl.append(VEP_DATA)
    cl.append("--tumor-id")
    cl.append(tum_id)
    cl.append("--cache-version")
    cl.append(CACHE)
    return cl


def run_vcf2maf(cl, sample):
    logger.info(f"Starting vcf2maf conversion of {sample}... This may take several minutes (approx 5 minutes per sample).")
    logger.info(f"args={cl}")
    sout = subprocess.run(cl, capture_output=True, check=False)

    if sout.stderr != None:
        if "ERROR" not in sout.stderr.decode("ascii"):
            logger.warning(sout.stderr.decode("ascii").replace("ERROR: ",""))
        else:
            logger.error(sout.stderr.decode("ascii").replace("ERROR: ",""))


def create_folder(output_folder, overwrite_output, resume):
    output_list = get_version_list(output_folder)
    output_list=list(map(lambda x: os.path.join(os.path.dirname(output_folder),x),output_list))

    if output_list != [] and os.path.exists(output_list[-1]):
        output = output_list[-1]
        logger.warning(f"It seems that a version of the folder '{output_folder}' already exists.")
        if overwrite_output:
            logger.info("Overwrite option set. Start removing folder")
            shutil.rmtree(output)
        elif resume:
            _,current_version = get_newest_version(output_folder)
            return output_folder + current_version

    if output_list == []:
        version = "_v1"
        output_folder_version = output_folder + version
    else:
        output_folder_version,_ = get_newest_version(output_folder)
    logger.info(f"Creating the output folder '{output_folder_version}'...")

    os.makedirs(output_folder_version, exist_ok=True)
    maf_path = os.path.join(output_folder_version, "maf")
    os.makedirs(maf_path, exist_ok=True)
    logger.info(f"The folder '{output_folder_version}' was correctly created!")

    return output_folder_version


def get_table_from_folder(tsvpath):
    table_dict = dict()
    file = pd.read_csv(tsvpath, sep="\t", index_col=False, dtype=str)
    for _, row in file.iterrows():
        sampleID = str(row["SAMPLE_ID"])
        if ".bam" in sampleID:
           sampleID = sampleID.replace(".bam", "")
        if sampleID not in table_dict:
            table_dict[sampleID] = [str(row["PATIENT_ID"])]
    return table_dict


def flatten(nested_list):
    flat_list = []
    for sublist in nested_list:
        for item in sublist:
            flat_list.append(item)
    return flat_list


def check_multiple_file(input_file, multiple):
    conf_snv = config.get("Multiple", "SNV")
    conf_cnv = config.get("Multiple", "CNV")
    file_paths = pd.read_csv(input_file, sep="\t", header=0)
    snv_file_path = file_paths["snv_path"].isnull().all()
    cnv_file_path = file_paths["cnv_path"].isnull().all()

    #CASE 1: multiple = True; neither snv or cnv are filled in sample.tsv; multiple snv or cnv are filled in conf.ini
    if multiple and not (snv_file_path or cnv_file_path) and not (conf_snv == "" or conf_cnv == ""):
        logger.critical("-m was selected and Muliple section in conf.ini was filled but the file doesn't looks like a multiVCF.")
        raise Exception ("Input error")
    #CASE 2: multiple = True; neither snv or cnv are filled in sample.tsv; neither multiple snv or cnv are filled in conf.ini
    if multiple and not (snv_file_path or cnv_file_path) and (conf_snv == "" or conf_cnv == ""):
        logger.critical("-m was selected but Muliple section in conf.ini wasn't filled and the file doesn't looks like a multiVCF.")
        raise Exception ("Input error")
    #CASE 3: multiple = False; snv or cnv are filled in sample.tsv; multiple snv or cnv are filled in conf.ini
    if not multiple and (snv_file_path or cnv_file_path) and not (conf_snv == "" or conf_cnv == ""):
        logger.critical("-m was not selected but both sample.tsv and Muliple section in conf.ini were filled.")
        raise Exception ("Input error")
    #CASE 4: multiple = False; snv or cnv are filled in sample.tsv; neither multiple snv or cnv are filled in conf.ini
    if not multiple and (snv_file_path or cnv_file_path) and (conf_snv == "" or conf_cnv == ""):
        logger.warning("SNV and/or CNV columns in sample.tsv were not filled.")


def check_multiple_folder(input_dir, multiple):
    snv_mulitple, snv_single, cnv_mulitple, cnv_single = False, False, False, False
    snv_folder = os.path.join(input_dir, "SNV")
    try:
        multiple_vcf_snv = [file for file in os.listdir(snv_folder) if file.endswith(".vcf")][0]
    except Exception:
        multiple_vcf_snv = ""
        os.makedirs(snv_folder, exist_ok=True)

    snv_file = os.path.join(input_dir, "sample_id_snv.txt")
    cmd_snv = "vcf-query -l " + os.path.join(snv_folder, multiple_vcf_snv) + " > " + snv_file
    os.system(cmd_snv)
    with open (snv_file) as file:
        lines = file.readlines()
        if len(lines) >= 2 and not multiple:
            snv_mulitple = True
            logger.error("-m option was not selected but the SNV file is multiple!")
        elif len(lines) < 2 and multiple:
            snv_single = True
            logger.error("-m option was selected but the SNV file is not multiple!")
    os.remove(snv_file)

    cnv_folder = os.path.join(input_dir, "CNV")
    try:
        multiple_vcf_cnv = [file for file in os.listdir(cnv_folder) if file.endswith(".vcf")][0]
    except Exception:
        os.makedirs(cnv_folder, exist_ok=True)
        return

    cnv_file = os.path.join(input_dir, "sample_id_cnv.txt")
    cmd_cnv = "vcf-query -l " +  os.path.join(cnv_folder, multiple_vcf_cnv) + " > " + cnv_file
    os.system(cmd_cnv)
    with open (cnv_file) as file:
        lines = file.readlines()
        if len(lines) >= 2 and not multiple:
            cnv_mulitple = True
            logger.error("-m option was not selected but the CNV file is multiple!")
        elif len(lines) < 2 and multiple:
            cnv_single = True
            logger.error("-m option was selected but the CNV file is not multiple!")
    os.remove(cnv_file)

    if cnv_mulitple or snv_mulitple:
        raise Exception ("Input error")
    if cnv_single or snv_single:
        raise Exception ("Input error")


def write_clinical_sample(clin_samp_path, output_folder, table_dict):
    logger.info("Writing data_clinical_sample.txt file...")
    conf_header_short = config.get("ClinicalSample", "HEADER_SAMPLE_SHORT")
    conf_header_long = config.get("ClinicalSample", "HEADER_SAMPLE_LONG")
    conf_header_type = config.get("ClinicalSample", "HEADER_SAMPLE_TYPE")

    data_clin_samp = pd.read_csv(clin_samp_path, sep="\t", header=0, dtype=str)
    data_clin_samp["ONCOTREE_CODE"] = data_clin_samp["ONCOTREE_CODE"].str.upper()

    try:
        data_clin_samp.drop(["snv_path", "cnv_path", "comb_path"], axis=1, inplace=True)
    except KeyError:
        logger.critical("snv_path, cnv_path or comb_path columns were removed or modified from template. Please use the correct template!")
        raise NameError("Exiting from script!")

    data_clin_samp.columns = data_clin_samp.columns.str.upper()

    combout_df = pd.DataFrame.from_dict(table_dict).transpose().reset_index()
    combout_df = combout_df.rename(columns={"index": "SAMPLE_ID", 0: "PATIENT_ID"})

    final_data_sample = data_clin_samp

    if len(combout_df.columns) > 2:
        combout_df = combout_df.rename(columns={1: "MSI", 2: "TMB", 3: "MSI_THR", 4:"TMB_THR"})
        try:
            if data_clin_samp["MSI"].notna().any() or data_clin_samp["TMB"].notna().any():
                if (any(data_clin_samp["MSI"] != combout_df["MSI"]) or any(data_clin_samp["TMB"] != combout_df["TMB"])):
                    logger.warning("MSI and/or TMB values are reported in sample.tsv and CombinedOutput but they do not match! CombinedOutput values were selected by default")
            try:
                data_clin_samp.drop(columns=["MSI", "TMB", "MSI_THR", "TMB_THR"], inplace=True)
            except KeyError:
                logger.critical("MSI_THR or TMB_THR columns were removed or modified from template. Please use the correct template!")
                raise(NameError("Exiting from script!"))
        except KeyError:
            logger.critical("MSI or TMB columns were removed or modified from template. Please use the correct template!")
            raise(KeyError("Exiting from script!"))

        final_data_sample = pd.merge(data_clin_samp, combout_df, on=["PATIENT_ID", "SAMPLE_ID"])

    new_cols = ["SAMPLE_ID", "PATIENT_ID", "MSI", "TMB", "MSI_THR", "TMB_THR"] + final_data_sample.columns[2:-4].to_list()
    final_data_sample = final_data_sample[new_cols]

    dataclin_columns = list(final_data_sample.columns)

    # Add header's fifth row
    default_row = pd.DataFrame([dataclin_columns], columns=dataclin_columns)
    final_data_sample = pd.concat([default_row, final_data_sample], ignore_index=True)

    # Add header's fourth row (SERIES OF 1s)
    header_numbers = pd.DataFrame([[1] * len(final_data_sample.columns)], columns=dataclin_columns)
    final_data_sample = pd.concat([header_numbers, final_data_sample], ignore_index=True)

    # Add header's third row (HEADER_SAMPLE_TYPE)
    if not conf_header_type:
        header_row = ["STRING"] * len(dataclin_columns)
        header_row[dataclin_columns.index("MSI")] = "NUMBER"
        header_row[dataclin_columns.index("TMB")] = "NUMBER"
        sample_header_type = pd.DataFrame([header_row], columns=dataclin_columns)
    else:
        types_list = conf_header_type.split(",")
        types_list = list(map(lambda x: x.strip(), types_list))
        for types in types_list:
            if types.upper() not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{types} is not a valid type. Please check the given input in conf.ini. Valid types: STRING, NUMBER, BOOLEAN")
                raise(NameError("The type is not valid: exiting from walk script!"))
        try:
            types_list = list(map(lambda x: x.upper(), types_list))

            sample_header_type = pd.DataFrame([types_list], columns=dataclin_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(types_list)}) in HEADER_SAMPLE_TYPE is different from the effective number of columns ({len(final_data_sample.columns)}).")
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_sample = pd.concat([sample_header_type, final_data_sample], ignore_index=True)

    # Add header's second row (HEADER_SAMPLE_LONG)
    if not conf_header_long:
        sample_header_long = default_row
    else:
        try:
            combined_headers = conf_header_long.split(",")
            sample_header_long = pd.DataFrame([combined_headers], columns=dataclin_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(combined_headers)}) in HEADER_SAMPLE_LONG in conf.ini is different from the effective number of columns ({len(final_data_sample.columns)}).")
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_sample = pd.concat([sample_header_long, final_data_sample], ignore_index=True)

    # Add header's first row (HEADER_SAMPLE_SHORT)
    if not conf_header_short:
        sample_header_short = default_row
    else:
        try:
            combined_headers = conf_header_short.split(",")
            sample_header_short = pd.DataFrame([combined_headers], columns=dataclin_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(combined_headers)}) in HEADER_SAMPLE_SHORT in conf.ini is different from the effective number of columns ({len(final_data_sample.columns)}).")
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_sample = pd.concat([sample_header_short, final_data_sample], ignore_index=True)
    final_data_sample.loc[0:3, "SAMPLE_ID"] = final_data_sample.loc[0:3, "SAMPLE_ID"].apply(lambda x: f"#{x}")

    data_clin_txt = os.path.join(output_folder, "data_clinical_sample.txt")
    final_data_sample.to_csv(data_clin_txt, sep="\t", index=False, header=False)


def write_default_clinical_patient(output_folder, table_dict):
    logger.info("Writing data_clinical_patient.txt file...")
    data_clin_samp = os.path.join(output_folder, "data_clinical_patient.txt")
    cil_sample = open(data_clin_samp, "w")
    cil_sample.write("#Patient Identifier\tAge\tGender\n")
    cil_sample.write("#Patient identifier\tAge\tGender\n")
    cil_sample.write("#STRING\tNUMBER\tSTRING\n")
    cil_sample.write("#1\t1\t1\n")
    cil_sample.write("PATIENT_ID\tAGE\tGENDER\n")

    nested_list = list(table_dict.values())
    list_patients = flatten(nested_list)
    list_patients = set(list_patients)
    for v in list_patients:
        cil_sample.write(v + "\tNaN\tNaN\n")
    cil_sample.close()


def add_header_patient_type(patient_tsv, datapat_columns, conf_header_type, final_data_pat):
    if not conf_header_type:
        def_type = ["STRING", "NUMBER", "STRING"]
        header_type = def_type + ["STRING"] * (len(datapat_columns)-3)
        header_type_df = pd.DataFrame([header_type], columns=datapat_columns)
    else:
        types_list = conf_header_type.split(",")
        types_list = list(map(lambda x: x.strip().upper(), types_list))
        for types in types_list:
            if types not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{types} is not a valid type. Please check the given input in conf.ini. " +
                                    "Valid types: STRING, NUMBER, BOOLEAN")
                raise(NameError("The type is not valid: exiting from walk script!"))
        try:
            header_type_df = pd.DataFrame([types_list], columns=datapat_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(types_list)}) in HEADER_PATIENT_TYPE " +
                                f"is different from the effective number of columns ({len(datapat_columns)}) in {patient_tsv}.")
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_pat = pd.concat([header_type_df, final_data_pat], ignore_index=True)
    return final_data_pat


def add_header_patient_short(patient_tsv, datapat_columns, conf_header_short, default_row, final_data_pat):
    if not conf_header_short:
        pat_header_short = default_row
    else:
        try:
            pat_header_short = pd.DataFrame([conf_header_short.split(", ")], columns=datapat_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(conf_header_short.split(', '))}) in HEADER_PATIENT_SHORT in conf.ini " +
                                f"is different from the effective number of columns ({len(datapat_columns)}) in {patient_tsv}.")
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_pat = pd.concat([pat_header_short, final_data_pat], ignore_index=True)
    return final_data_pat


def add_header_patient_long(patient_tsv, datapat_columns, conf_header_long, default_row, final_data_pat):
    if not conf_header_long:
        pat_header_long = default_row
    else:
        try:
            pat_header_long = pd.DataFrame([conf_header_long.split(", ")], columns=datapat_columns)
        except ValueError:
            logger.critical(f"The number of column names ({len(conf_header_long.split(', '))}) in HEADER_PATIENT_LONG " +
                                f"in conf.ini is different from the effective number of columns ({len(datapat_columns)}) in {patient_tsv}.")
            raise(NameError("Different number of columns: exiting from walk script!"))
    final_data_pat = pd.concat([pat_header_long, final_data_pat], ignore_index=True)
    return final_data_pat


def extract_multiple_cnv(multiple_vcf, input_dir):
    single_sample_vcf_dir = os.path.join(input_dir, "single_sample_vcf")
    if not os.path.exists(os.path.join(input_dir, "single_sample_vcf")):
        os.mkdir(os.path.join(input_dir, "single_sample_vcf"))
    cmd = "vcf-query -l " + multiple_vcf + " > " + os.path.join(input_dir, "sample_id.txt")
    os.system(cmd)

    with open (os.path.join(input_dir, "sample_id.txt")) as file:
        lines = file.readlines()
        for sample in lines:
            sample = sample.strip()
            single_sample_vcf = os.path.join(single_sample_vcf_dir, f"{sample}.vcf")
            cmd = f"vcftools --vcf {multiple_vcf} --indv {sample} --recode --recode-INFO-all --stdout > {single_sample_vcf}"
            os.system(cmd)

def extract_multiple_snv(multiple_vcf,input_dir):
    if not os.path.exists(os.path.join(input_dir, "single_sample_vcf")):
        os.mkdir(os.path.join(input_dir, "single_sample_vcf"))

    cmd = "vcf-query -l " + multiple_vcf + " > " + os.path.join(input_dir, "sample_id.txt")
    os.system(cmd)

    with open (os.path.join(input_dir, "sample_id.txt")) as file:
        lines = file.readlines()
        for sample in lines:
            cmd="vcf-subset --exclude-ref -c " + sample.strip() + " " + multiple_vcf + " > " + os.path.join(os.path.join(input_dir, "single_sample_vcf"), sample.strip() + ".vcf")
            os.system(cmd)


def check_field_tsv(row, name):
    try:
        field=str(row[name])
    except KeyError as e:
        logger.critical(f"KeyError: {e} not found! Check if column name is correctly spelled or if there are tabs/spaces before or after the coloumn key: \n{row.index}. \nThis error may also occur if the table columns have not been separated by tabs!")
        sys.exit()
    return field


def get_combinedVariantOutput_from_folder(inputFolder, file, isinputfile):
    combined_dict = dict()
    for _,row in file.iterrows():
        sampleID = check_field_tsv(row, "SAMPLE_ID")
        patientID = check_field_tsv(row, "PATIENT_ID")
        if isinputfile:
            combined_path = check_field_tsv(row, "comb_path")
        else:
            combined_path = os.path.join(inputFolder, "CombinedOutput", patientID + "_CombinedVariantOutput.tsv")

        if os.path.exists(combined_path):
            pass
        else:
            logger.warning("comb_path in conf.ini does not exists")
        combined_dict[sampleID] = combined_path
    return combined_dict


def check_input_file(output_folder, file, copy_to, sample_id):
    if os.path.exists(file):
        os.system("cp " + file + " " + os.path.join(output_folder, "temp", copy_to))
    elif not file:
        logger.warning(f"No final_path set in conf.ini for sample {sample_id}'s {copy_to}!")
    else:
        logger.warning(f"{file} not found")


def check_folders(output_folder, snv_path, cnv_path, combout, sample_id):
    logger.info(f"Verifying the compilation of the conf.ini file for sample {sample_id}")
    check_input_file(output_folder, snv_path, "SNV", sample_id)
    check_input_file(output_folder, cnv_path, "CNV", sample_id)
    check_input_file(output_folder, combout, "CombinedOutput", sample_id)


def transform_input(tsv, clin_pzt, fusion_tsv, output_folder, multiple):
    os.makedirs(os.path.join(output_folder, "temp"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "SNV"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "CNV"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "CombinedOutput"), exist_ok=True)
    os.makedirs(os.path.join(output_folder, "temp", "FUSIONS"), exist_ok=True)

    os.system("cp " + tsv + " " + os.path.join(output_folder, "temp", "sample.tsv"))
    if clin_pzt!="":
        os.system("cp " + clin_pzt + " " + os.path.join(output_folder, "temp", "patient.tsv"))
    if fusion_tsv!="":
        os.system("cp " + fusion_tsv + " " + os.path.join(output_folder, "temp", "FUSIONS", "fusion.tsv"))

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

    return os.path.join(output_folder, "temp")


def fill_fusion_from_temp(input, fusion_table_file, clin_file, fusion_files):
    nfusion=len(fusion_files)
    logger.info(f"Found {nfusion} Fusion file(s)")

    fusion_input = os.path.join(input, "FUSIONS", fusion_files[0])
    with open(fusion_input) as template:
        header = template.readline()

    with open(fusion_table_file, "w") as fusion_table:

        fusion_table.write(header)

        for fusion_file in fusion_files:
            ff = pd.read_csv(fusion_input, sep="\t")

            if not set(["Sample_Id","SV_Status","Site1_Hugo_Symbol","Site2_Hugo_Symbol"]).issubset(ff.columns):
                logger.warning(f"{fusion_file} does not contain required columns")
                continue

            if ff.shape[0]==0:
                logger.info(f"No Fusions found in {fusion_file}")
                continue
            logger.info(f"Fusions found in {fusion_file}")
            for fus in ff.itertuples(index=False):
                if fus.Sample_Id in clin_file["SAMPLE_ID"].values:
                    if int(fus.Normal_Paired_End_Read_Count) >= 15:
                        fusion_table.write("\t".join(map(str, fus)) + "\n")


def annotate_fusion(cancer, fusion_table_file, data_sv, input_file):
    if "ONCOTREE_CODE" in input_file.columns:
        fusion_table_df = data_sv.merge(input_file[["SAMPLE_ID", "ONCOTREE_CODE"]], how="inner", left_on="Sample_Id", right_on="SAMPLE_ID")
        fusion_table_df.to_csv(fusion_table_file, sep="\t", index=False)

        fusion_table_file_out = fusion_table_file.replace(".txt", "ann.txt")
        os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                        -o {fusion_table_file_out} -b {config.get('OncoKB', 'ONCOKB')}")

    else:
        fusion_table_file_out = fusion_table_file.replace(".txt", "ann.txt")
        os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                            -o {fusion_table_file_out} -t {cancer.upper()}  -b {config.get('OncoKB', 'ONCOKB')}")

    return fusion_table_file_out


def fill_fusion_from_combined(fusion_table_file, combined_dict, THR_FUS):
    logger.info("Writing data_sv.txt file...")

    with open(fusion_table_file, "w") as fusion_table:

        header = "Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\tNormal_Paired_End_Read_Count\tEvent_Info\tRNA_Support\n"
        fusion_table.write(header)

        for k, v in combined_dict.items():
            fusions=[]
            try:
                fusions = tsv.get_fusions(v)
            except Exception:
                logger.error(f"Something went wrong while reading Fusion section for sample {k}")
            if len(fusions) == 0:
                continue

            for fus in fusions:
                if len(fusions) > 0:
                    Site1_Hugo_Symbol = fus["Site1_Hugo_Symbol"]
                    Site2_Hugo_Symbol = fus["Site2_Hugo_Symbol"]
                    if Site2_Hugo_Symbol == "CASC1":
                        Site2_Hugo_Symbol = "DNAI7"
                    Site1_Chromosome = fus["Site1_Chromosome"]
                    Site2_Chromosome = fus["Site2_Chromosome"]
                    Site1_Position = fus["Site1_Position"]
                    Site2_Position = fus["Site2_Position"]

                    if eval("int(fus['Normal_Paired_End_Read_Count'])" + THR_FUS):
                        fusion_table.write(k + "\tSOMATIC\tFUSION\t" + str(Site1_Hugo_Symbol) + "\t" + \
                        str(Site2_Hugo_Symbol) + "\t" + fus["Normal_Paired_End_Read_Count"] + "\t" + \
                        fus["Event_Info"] + " Fusion\t" + "Yes\n")


def check_data_cna(data_cna_path):
    input_file = os.path.basename(data_cna_path)
    logger.info(f"Checking {input_file} file...")
    try:
        with open(data_cna_path) as data_cna:
            all_data_cna = data_cna.readlines()
            if (len(all_data_cna) == 1):
                os.remove(data_cna_path)
                logger.warning(f"{input_file} is empty. File removed.")
    except Exception:
        logger.warning(f"{input_file} does not exist!")


def fill_from_file(table_dict_patient, file_input_clinical, MSI_THR, TMB_THR):
    for k, m, t in zip(file_input_clinical["SAMPLE_ID"], file_input_clinical["MSI"], file_input_clinical["TMB"]):
        table_dict_patient[k].append(m)
        table_dict_patient[k].append(t)

        if np.isnan(float(m)):
            table_dict_patient[k].append("NA")
        elif eval("float(m)" + MSI_THR):
            table_dict_patient[k].append("Stable")
        else:
            table_dict_patient[k].append("Unstable")

        if not np.isnan(float(t)):
            for _k, _v in TMB_THR.items():
                if eval("float(t)" + _v):
                    table_dict_patient[k].append(_k)
                    break
        else:
            table_dict_patient[k].append("NA")
    return table_dict_patient


def fill_from_combined(combined_dict, table_dict_patient, MSI_SITES_THR, MSI_THR, TMB):
    for k, v in combined_dict.items():
        try:
            tmv_msi = tsv.get_msi_tmb(v)
        except Exception:
            logger.error("Something went wrong!")

        if not tmv_msi["MSI"][0][1]=="NA" and eval("float(tmv_msi['MSI'][0][1])" + MSI_SITES_THR):
            table_dict_patient[k].append(tmv_msi["MSI"][1][1])
        else:
            table_dict_patient[k].append("NA")

        table_dict_patient[k].append(tmv_msi["TMB_Total"])

        if not tmv_msi["MSI"][0][1]=="NA" and not tmv_msi["MSI"][1][1] =="NA" and not table_dict_patient[k][1]=="NA":
            if eval("float(tmv_msi['MSI'][1][1])" + MSI_THR):
                table_dict_patient[k].append("Stable")
            else:
                table_dict_patient[k].append("Unstable")
        else:
            table_dict_patient[k].append("NA")

        if not tmv_msi["TMB_Total"]=="NA":
            found = False
            for _k, _v in TMB.items():
                if eval(tmv_msi["TMB_Total"] + _v):
                    table_dict_patient[k].append(_k)
                    found = True
                    break
            if found == False:
                logger.warning(f"The TMB value {tmv_msi['TMB_Total']} is not within the conf.ini thresholds {list(TMB.values())}. For this value, the TMB_THR will be set as 'Out of threshold ranges'.")
                table_dict_patient[k].append("Out of threshold ranges")
        else:
            table_dict_patient[k].append("NA")

    return table_dict_patient


def input_extraction_file(input):
    sample_tsv=input[0]
    patient_tsv, fusion_tsv = "",""

    if len(input)>1:
        patient_tsv=input[1].strip()
        if len(input)>2:
            fusion_tsv=input[2].strip()

    return sample_tsv, patient_tsv, fusion_tsv


def input_extraction_folder(input):
    patient_tsv, fusion_tsv = "",""
    sample_tsv = os.path.join(input, "sample.tsv")
    if os.path.exists(os.path.join(input, "patient.tsv")):
        patient_tsv = os.path.join(input, "patient.tsv")
    if os.path.exists(os.path.join(input, "FUSIONS", "Fusions.tsv")):
        fusion_tsv = os.path.join(input, "FUSIONS", "Fusions.tsv")
    return sample_tsv, patient_tsv, fusion_tsv


def validate_input(oncokb, vcf_type, filters, cancer, input):
    #check that the input file is correctly made
    if not isinputfile:
        file_tsv = pd.read_csv(os.path.join(input, "sample.tsv"), sep="\t", dtype=str)
    elif isinputfile:
        file_tsv = pd.read_csv(input, sep="\t", dtype=str)
    column_list = ["SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path", "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR"]
    if not all(name in file_tsv.columns for name in column_list):
        logger.critical('Necessary columns: "SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path", "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR"')
        raise(Exception("The input file is missing some important columns!"))

    #check that oncokb key is filled in conf.ini when oncokb annotation is selected
    if oncokb:
        assert config.get("OncoKB", "ONCOKB")!="", \
               "oncokb option was set but ONCOKB field in conf.ini is empty!"

    #check that vep info in conf.ini are set when snv analysis is request
    if vcf_type==None or "snv" in vcf_type:
        assert VEP_PATH !="" and VEP_DATA !="", \
               "VEP_PATH and/or VEP_DATA field in conf.ini is empty!"
        assert REF_FASTA !="", \
               "REF_FASTA field in conf.ini is empty!"

    #verify that oncokb filter function only when oncokb annotation is set
    if "o" in filters and oncokb==False:
        logger.warning("OncoKB filter was selected in filters options but -k option was not set. This filtering will be ignored.")

    #check if cancer id is compatible with cbioportal
    cancer_cbio=pd.read_csv("cancer_list.txt", sep="\t")
    cancer_cbio=cancer_cbio["TYPE_OF_CANCER_ID"].values.tolist()

    if cancer not in cancer_cbio:
        logger.critical(f"cancer_id '{cancer}' is not recognize by cBioPortal. Check the cancer_list.txt to find the correct cancer id")
        sys.exit()


def walk_folder(input, multiple, output_folder, oncokb, cancer, overwrite_output=False, resume=False, vcf_type=None, filters=""):

    logger.info("Starting walk_folder script:")
    logger.info(f"walk_folder args [input:{input}, output_folder:{output_folder}, overwrite:{overwrite_output}, resume:{resume}, vcf_type:{vcf_type}, filters:{filters}, multiple:{multiple}]")

    config.read("conf.ini")

    assert os.path.exists(input[0]), \
            f"No valid file/folder {input} found. Check your input path"

    global isinputfile

    if os.path.isdir(input[0]):
        isinputfile = False
    elif os.path.isfile(input[0]):
        isinputfile = True
    validate_input(oncokb, vcf_type, filters, cancer, input[0])


    ###############################
    ###      OUTPUT FOLDER      ###
    ###############################

    if not resume or not os.path.exists(os.path.join(output_folder, "temp")):
        output_folder = create_folder(output_folder, overwrite_output, resume)

    if not isinputfile:
        input_folder = input[0]
        input, patient_tsv, fusion_tsv = input_extraction_folder(input_folder)
        check_multiple_folder(input_folder, multiple)

    elif isinputfile:
        input, patient_tsv, fusion_tsv = input_extraction_file(input)
        check_multiple_file(input, multiple)
        input_folder = transform_input(input, patient_tsv, fusion_tsv, output_folder, multiple)

    else:
        logger.critical(f"The input {input} isn't a file nor a folder")
        raise(FileNotFoundError("Exiting from walk script!"))

    img_path = os.path.join("docs", "img", "logo_VARAN.png")
    if os.path.exists(img_path):
        img_output_dir = os.path.join(output_folder, "img")
        os.makedirs(img_output_dir, exist_ok=True)
        shutil.copy(img_path, os.path.join(img_output_dir, "logo_VARAN.png"))

    inputFolderSNV = os.path.abspath(os.path.join(input_folder, "SNV"))
    inputFolderCNV = os.path.abspath(os.path.join(input_folder, "CNV"))
    inputFolderCombOut = os.path.abspath(os.path.join(input_folder, "CombinedOutput"))

    assert (len(os.listdir(inputFolderSNV))>0 or len(os.listdir(inputFolderCNV))>0 or len(os.listdir(inputFolderCombOut))>0), \
        "No valid input file was found for neither SNV, CNV or CombinedOutput! Check your input file/folder and input options."

    maf_path = os.path.join(output_folder, "maf")
    maf_zip_path = os.path.join(output_folder, "maf.zip")
    clin_sample_path = os.path.join(input_folder, "sample.tsv")
    try:
        clin_file = pd.read_csv(clin_sample_path, sep="\t", dtype=str)
    except Exception:
        logger.critical(f"Something went wrong while reading {clin_sample_path}!")
        raise(Exception("Error in reading the input file! Please check again."))

    if resume:
        if os.path.exists(maf_path):
            try:
                maf_samples = set(os.listdir(maf_path))
            except:
                logger.critical("Error accessing the maf folder. Please check its integrity.")
                raise(Exception("Error reading the maf folder."))
        elif os.path.exists(maf_zip_path):
            try:
                with zipfile.ZipFile(maf_zip_path, "r") as zipped_maf:
                    zipped_maf.extractall(maf_path)
                maf_samples = set(os.listdir(maf_path))
            except:
                logger.critical("Error extracting the maf ZIP file. Please check its integrity.")
                raise(Exception("Unable to extract the ZIP file. Try decompressing it manually and run Varan again."))


        clin_samples = set(clin_file["SAMPLE_ID"])
        clin_in_maf = all(any(clin_sample in maf_sample for maf_sample in maf_samples) for clin_sample in clin_samples)
        maf_in_clin = all(any(clin_sample in maf_sample for clin_sample in clin_samples) for maf_sample in maf_samples)

        if not (clin_in_maf and maf_in_clin) and len(maf_samples) != 0:
            logger.critical("It seems you are resuming an existing study with a different set of input samples. Please verify the sample consistency!")
            raise FileNotFoundError("Sample mismatch detected.")

        ZIP_MAF = config.get("Zip", "ZIP_MAF")
        ZIP_MAF = check_bool(ZIP_MAF)
        if os.path.exists(maf_zip_path) and not ZIP_MAF:
            os.remove(maf_zip_path)

    if len(os.listdir(inputFolderSNV))==0 and vcf_type == None:
        vcf_type = "cnv"
        logger.info("SNV path was empty, the analysis will exclude SNV")
    if len(os.listdir(inputFolderCNV))==0 and vcf_type == None:
        vcf_type = "snv"
        logger.info("CNV path was empty, the analysis will exclude CNV")

    if os.path.exists(inputFolderCNV) and vcf_type not in ["snv", "fus", "tab"]:
        if multiple:
            multivcf = [i for i in os.listdir(inputFolderCNV) if i.endswith(".vcf")][0]
            extract_multiple_cnv(os.path.join(inputFolderCNV, multivcf), inputFolderCNV)
            inputFolderCNV= os.path.join(inputFolderCNV, "single_sample_vcf")
        logger.info("Checking CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(inputFolderCNV)
        logger.info("Everything ok!")

    if os.path.exists(inputFolderSNV) and vcf_type not in ["cnv", "fus", "tab"]:
        os.makedirs(tmp, exist_ok=True)
        if multiple:
            multivcf = [i for i in os.listdir(inputFolderSNV) if i.endswith(".vcf")][0]
            extract_multiple_snv(os.path.join(inputFolderSNV, multivcf), inputFolderSNV)
            inputFolderSNV = os.path.join(inputFolderSNV, "single_sample_vcf")
        logger.info("Checking SNV files...")
    case_folder_arr = get_snv_from_folder(inputFolderSNV)
    logger.info("Everything ok!")

    ###############################
    ###       SNV AND CNV       ###
    ###############################

    if os.path.exists(inputFolderCNV) and vcf_type not in ["snv", "fus", "tab"]:
        logger.info("Managing CNV files...")
        sID_path_cnv = cnv_type_from_folder(input_folder, case_folder_arr_cnv, output_folder, oncokb, cancer, multiple)

    if os.path.exists(inputFolderSNV) and vcf_type not in ["cnv", "fus", "tab"]:
        logger.info("Managing SNV files...")
        sID_path_snv = snv_type_from_folder(inputFolderSNV, case_folder_arr)

        logger.info("Checking maf folder...")
        maf_path = os.path.join(output_folder, "maf")
        if os.path.isdir(maf_path) and len([i for i in os.listdir(maf_path) if i.endswith(".maf")])>0:
            logger.info("A non empty maf folder already exists!")

        if not resume:
            if "d" in filters:
                logger.info("Filtering out VCFs with dots in ALT column")
                sID_path_snv = vcf_filtering(sID_path_snv, output_folder, output_filtered)
            temporary = create_random_name_folder()
            for k, v in sID_path_snv.items():
                cl = vcf2maf_constructor(k, v, temporary, output_folder)
                run_vcf2maf(cl, k)

    logger.info("Clearing scratch folder...")
    clear_scratch()


    ###############################
    ###       GET FUSION        ###
    ###############################
    if vcf_type not in ["cnv","snv","tab"]:

        fusion_table_file = os.path.join(output_folder, "data_sv.txt")
        fusion_folder = os.path.join(input_folder, "FUSIONS")

        if os.path.exists(os.path.join(input_folder, "CombinedOutput")) and len(os.listdir(os.path.join(input_folder, "CombinedOutput")))>0:
            logger.info("Getting Fusions infos from CombinedOutput...")
            THR_FUS = config.get("FUSION", "THRESHOLD_FUSION")
            combined_dict = get_combinedVariantOutput_from_folder(input_folder, clin_file, isinputfile)
            fill_fusion_from_combined(fusion_table_file, combined_dict, THR_FUS)

        elif os.path.exists(os.path.abspath(fusion_folder)) and os.listdir(fusion_folder):
            fusion_files=[file for file in os.listdir(fusion_folder) if "tsv" in file]
            logger.info(f"Getting Fusions infos from {fusion_files[0]} file...")
            if fusion_files != []:
                fill_fusion_from_temp(input_folder, fusion_table_file, clin_file, fusion_files)

        if os.path.exists(fusion_table_file):
            with open(fusion_table_file) as data_sv:
                all_data_sv = data_sv.readlines()
                if (len(all_data_sv) == 1):
                    os.remove(fusion_table_file)
                    logger.warning("data.sv is empty. File removed.")

        if oncokb and os.path.exists(fusion_table_file):
            data_sv = pd.read_csv(fusion_table_file, sep="\t")
            input_file = pd.read_csv(clin_sample_path, sep="\t")
            fusion_table_file_out = annotate_fusion(cancer, fusion_table_file, data_sv, input_file)

            if "o" in filters:
                fus_file = pd.read_csv(fusion_table_file_out, sep="\t")
                fus_file = filter_OncoKB(fus_file)
                fus_file.to_csv(fusion_table_file_out, index=False, sep="\t")

            data_sv_tmp = pd.read_csv(fusion_table_file_out, sep="\t")
            try:
                data_sv_tmp.drop(["SAMPLE_ID", "ONCOTREE_CODE"], inplace=True, axis=1)
            except KeyError:
                pass

            data_sv_tmp.to_csv(fusion_table_file_out, index=False, sep="\t")
            os.system(f"mv {fusion_table_file_out} {fusion_table_file}")


    ##############################
    ##       MAKES TABLE       ###
    ##############################

    table_dict_patient = get_table_from_folder(clin_sample_path)

    logger.info("Writing clinical files...")
    if patient_tsv.strip()!="":
        logger.info("Writing data_clinical_patient.txt file...")

        input_file_path = os.path.join(input_folder, "patient.tsv")
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
        header_numbers = pd.DataFrame([[1] * len(datapat_columns)], columns=datapat_columns)
        final_data_pat = pd.concat([header_numbers, final_data_pat], ignore_index=True)

        # Add header's third row (HEADER_PATIENT_TYPE)
        final_data_pat = add_header_patient_type(patient_tsv, datapat_columns, conf_header_type, final_data_pat)

        # Add header's second row (HEADER_PATIENT_LONG)
        final_data_pat = add_header_patient_long(patient_tsv, datapat_columns, conf_header_long, default_row, final_data_pat)

        # Add header's first row (HEADER_PATIENT_SHORT)
        final_data_pat = add_header_patient_short(patient_tsv, datapat_columns, conf_header_short, default_row, final_data_pat)

        final_data_pat.loc[0:3, "PATIENT_ID"] = final_data_pat.loc[0:3, "PATIENT_ID"].apply(lambda x: f"#{x}")

        data_clin_txt = os.path.join(output_folder, "data_clinical_patient.txt")
        final_data_pat.to_csv(data_clin_txt, sep="\t", index=False, header=False)

    else:
        write_default_clinical_patient(output_folder, table_dict_patient)

    fileinputclinical = pd.read_csv(clin_sample_path, sep="\t", index_col=False, dtype=str)

    MSI_THR = config.get("MSI", "THRESHOLD_MSI")
    TMB_THR = ast.literal_eval(config.get("TMB", "THRESHOLD_TMB"))

    if os.path.exists(os.path.join(input_folder, "CombinedOutput")) and \
        len(os.listdir(os.path.join(input_folder, "CombinedOutput"))) > 0:
        MSI_SITES_THR = config.get("MSI", "THRESHOLD_SITES")

        combined_dict = get_combinedVariantOutput_from_folder(input_folder, clin_file, isinputfile)
        new_table_dict_patient = fill_from_combined(combined_dict, table_dict_patient, MSI_SITES_THR, MSI_THR, TMB_THR)
    else:
        new_table_dict_patient = fill_from_file(table_dict_patient, fileinputclinical, MSI_THR, TMB_THR)

    write_clinical_sample(clin_sample_path, output_folder, new_table_dict_patient)
    logger.success("Walk script completed!\n")

    return output_folder, input, fusion_tsv
