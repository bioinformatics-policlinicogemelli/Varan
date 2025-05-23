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
import shutil
import sys
from configparser import ConfigParser
from pathlib import Path

import pandas as pd
from loguru import logger

import concatenate

config = ConfigParser()
configFile = config.read("conf.ini")

def check_bool(key_value):
    """
    Converts a string representation of a boolean value from the configuration file into a Python boolean.

    Accepted true values: "True", "true", "T"  
    Accepted false values: "False", "false", "F", "" (empty string)

    Parameters:
        key_value (str): A string value from the config file that should represent a boolean.

    Returns:
        bool: The corresponding boolean value (True or False).

    Raises:
        ValueError: If the input string does not match any of the accepted boolean representations.
    """

    bool_key_value = key_value
    if bool_key_value.strip() in ["True", "true", "T"]:
        bool_key_value = True
    elif bool_key_value.strip() in ["False", "false", "F", ""]:
        bool_key_value = False
    else:
        logger.critical(f'Please insert a boolean value in {bool_key_value} section of conf.ini: accepted values are ["True", "true", "T", "False", "false", "F", ""]')
        raise ValueError("Check again the compilation conf.ini")
    return bool_key_value


def print_unique_clin_sig(df):
    """
    Prints the unique values found in the 'CLIN_SIG' column of the input DataFrame.

    This function is useful for inspecting the different clinical significance annotations
    present in a mutation annotation file (MAF).

    Parameters:
        df (pd.DataFrame): A pandas DataFrame containing at least the 'CLIN_SIG' column.

    Returns:
        None
    """
    unique_clin_sig = df["CLIN_SIG"].unique()
    print(unique_clin_sig)


def filter_OncoKB(df):
    """
    Filters the input DataFrame based on the 'ONCOGENIC' column using allowed values
    specified in the 'ONCOKB_FILTER' field of the config file.

    This function is intended to retain only mutations with oncogenicity annotations
    matching those defined as relevant in the configuration.

    Parameters:
        df (pd.DataFrame): A pandas DataFrame that includes an 'ONCOGENIC' column,
                           typically a MAF (Mutation Annotation Format) table.

    Returns:
        pd.DataFrame: A filtered DataFrame containing only rows where 'ONCOGENIC'
                      matches the allowed values from config.
    """

    oncokb_filter=ast.literal_eval(config.get("Filters", "ONCOKB_FILTER"))

    return df[df["ONCOGENIC"].isin(oncokb_filter)]


def check_CLIN_SIG(row):
    """
    Checks whether any of the comma-separated clinical significance annotations in the 
    'CLIN_SIG' field of a row are NOT included in the allowed values specified in the 
    configuration file.

    This function is used as a row-wise filter to identify variants with at least one 
    non-accepted clinical significance annotation.

    Parameters:
        row (pd.Series): A row of a pandas DataFrame containing a 'CLIN_SIG' column,
                         with comma-separated clinical significance annotations (e.g., "benign,likely_benign").

    Returns:
        bool: True if at least one annotation in 'CLIN_SIG' is NOT in the allowed list
              from the config ('CLIN_SIG' field under [Filters]); False otherwise.
    """

    clin_sig=ast.literal_eval(config.get("Filters", "CLIN_SIG"))
    output=[]
    for _e in str(row["CLIN_SIG"]).split(","):
        if _e not in clin_sig:
            output.append(True)
        else:
            output.append(False)
    return any(output)


def check_consequences(row):
    """
    Checks whether any of the comma-separated consequence annotations in the 'Consequence'
    field of a row are included in the list of accepted consequences specified in the
    configuration file.

    This function is typically used to filter variants based on their predicted functional
    consequences (e.g., missense_variant, stop_gained).

    Parameters:
        row (pd.Series): A row of a pandas DataFrame containing a 'Consequence' column
                         with comma-separated annotations.

    Returns:
        bool: True if at least one annotation is in the list of allowed consequences
              from the config ('CONSEQUENCES' field under [Filters]); False otherwise.
    """

    consequences=ast.literal_eval(config.get("Filters", "CONSEQUENCES"))
    output=[]
    for _e in str(row["Consequence"]).split(","):
        if _e in consequences:
            output.append(True)
        else:
            output.append(False)
    return any(output)


def check_polyphen(row):
    """
    Verifica se la predizione PolyPhen riportata in una riga del DataFrame è tra quelle accettate,
    come specificato nella sezione 'POLYPHEN' del file di configurazione.

    Il valore di 'PolyPhen' viene splittato sul carattere '(', e viene considerata solo la parte
    principale della predizione (es. 'probably_damaging') per il confronto.

    Parametri:
        row (pd.Series): Una riga del DataFrame che contiene la colonna 'PolyPhen'.

    Restituisce:
        bool: True se il valore di 'PolyPhen' è tra quelli accettati, False altrimenti.
    """
    consequences=ast.literal_eval(config.get("Filters", "POLYPHEN"))
    output=[]
    if str(row["PolyPhen"]).split("(")[0] in consequences:
        output.append(True)
    else:
        output.append(False)
    return any(output)


def check_sift(row):
    consequences=ast.literal_eval(config.get("Filters", "SIFT"))
    output=[]
    if str(row["SIFT"]).split("(")[0] in consequences:
        output.append(True)
    else:
        output.append(False)
    return any(output)


def write_csv_with_info(df, file_path):
    """
    Writes a pandas DataFrame to a CSV file with a custom version header.

    This function first writes a version line (`#version 2.4`) at the top of the file,
    and then appends the DataFrame content in tab-separated format, including headers
    and without the index.

    Parameters:
        df (pd.DataFrame): The DataFrame to be written to file.
        file_path (str): The path to the output file where the content will be saved.

    Returns:
        None
    """
    file_path = Path(file_path)
    with file_path.open("w") as f:
        f.write("#version 2.4\n")
        df.to_csv(file_path, sep="\t", index=False, header=True)


def filter_main(input,folder, output_folder, oncokb, filters, cancer, resume, overwrite=False):
    """
    Main function to apply a series of filters to MAF (Mutation Annotation Format) files, 
    optionally annotate them using OncoKB, and save the filtered results.

    The function performs the following steps:
    1. Checks whether the output folder for OncoKB-annotated files exists and handles it
       according to the `overwrite` and `resume` flags.
    2. Optionally runs OncoKB annotation on input MAF files using the provided cancer type(s).
    3. Applies various filters defined in a configuration file and saves the resulting
       filtered MAF files in a dedicated output directory.

    Parameters:
        input (str): Path to the sample information file (used to retrieve ONCOTREE_CODEs).
        folder (str): Path to the main input folder containing a "maf" subfolder with MAF files.
        output_folder (str): Path to the folder where output files will be written.
        oncokb (bool): Whether to annotate MAF files using OncoKB.
        filters (str): String of single-letter codes indicating which filters to apply:
                       - 'i': remove entries with low/unknown impact
                       - 'p': keep only variants with FILTER == "PASS"
                       - 'o': keep only oncogenic variants (from OncoKB)
                       - 'v': apply VAF thresholds (with optional distinction for novel variants using 'n')
                       - 'a': apply population AF filter
                       - 'c': filter using ClinVar significance
                       - 'q': filter using variant consequence annotations
                       - 'y': filter using PolyPhen predictions
                       - 's': filter using SIFT predictions
        cancer (str): Default cancer code used for OncoKB annotation if not found in input file.
        resume (bool): Whether to resume processing if partial results exist.
        overwrite (bool, optional): If True, overwrite existing OncoKB annotation results.
                                    Default is False.

    Returns:
        None

    Raises:
        Exception: If no input MAF files are found.
        SystemExit: If the output folder exists and overwrite/resume conditions are not satisfied.
    """

    logger.info("Starting filter_main script:")
    logger.info(f"filter_main args [maf_folder:{folder}, output_folder:{output_folder}, filters:{filters}, cancer:{cancer}, overwrite:{overwrite}]")

    maf_oncokb_path = Path(output_folder) / "MAF_OncoKB"
    if maf_oncokb_path.exists() and any(maf_oncokb_path.iterdir()):
        if overwrite:
            logger.warning("It seems that the folder 'MAF_OncoKB' already exists. Start removing process...")
            shutil.rmtree(os.path.join(output_folder,"MAF_OncoKB"))
        elif not resume:
            logger.critical("The folder 'MAF_OncoKB' already exists. To overwrite an existing folder add the -w option!")
            logger.critical("Exit without completing the task!")
            sys.exit()

    maf_folder = Path(folder) / "maf"
    file_list = concatenate.get_files_by_ext(maf_folder, "maf")

    if len(file_list) == 0:
        logger.critical(f"The maf folder {maf_folder} seems to be empty: check if SNV folder exists and contains the VCF files. If you don't have SNV vcf, use -t option to select specific analysis.")
        raise(Exception("Exiting from filter_clinvar script!"))

    if oncokb:
        maf_oncokb_path.mkdir(parents=True, exist_ok=True)
        extension = "_OncoAnnotated.maf"
        input_file = pd.read_csv(input, sep="\t")

        for f in file_list:
                    if extension in f:
                        continue
                    file = Path(f).name
                    file_no = file.replace(".maf", "") + extension
                    file_path = maf_oncokb_path / file_no
                    if "ONCOTREE_CODE" in input_file.columns:
                        for _, row in input_file.iterrows():
                            if row["SAMPLE_ID"] in file_no:
                                cancer_onco = row["ONCOTREE_CODE"] or cancer
                                os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f} "
                                        f"-o {file_path} -t {cancer_onco.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
                    else:
                        os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f} "
                                f"-o {file_path} -t {cancer.upper()} -b {config.get('OncoKB', 'ONCOKB')}")

    file_list = concatenate.get_files_by_ext(maf_folder, "maf")
    out_filter = output_folder/"MAF_filtered"

    if oncokb and "o" in filters:
        file_list = concatenate.get_files_by_ext(maf_oncokb_path, "maf")
        out_filter=output_folder/"MAF_Onco_filtered"

    if filters!="" and filters!="d":

        logger.info("Start filtering vcf...")
        out_filter.mkdir(parents=True, exist_ok=True)

        for file in file_list:
            file_to_filter=pd.read_csv(file, sep="\t", comment="#", dtype=object)

            if "i" in filters:
                file_to_filter = file_to_filter[~file_to_filter["IMPACT"].isin(ast.literal_eval(config.get("Filters", "IMPACT")))]

            if "p" in filters:
                file_to_filter=file_to_filter[file_to_filter["FILTER"]=="PASS"]

            if oncokb and "o" in filters:
                oncokb_filter=ast.literal_eval(config.get("Filters", "ONCOKB_FILTER"))
                file_to_filter= file_to_filter[file_to_filter["ONCOGENIC"].isin(oncokb_filter)]

            if "v" in filters:
                t_VAF_min=float(config.get("Filters", "t_VAF_min"))
                t_VAF_max=float(config.get("Filters", "t_VAF_max"))

                temp = file_to_filter.dropna(subset=["t_AF"])

                if len(temp) == 0:
                    vaf_colname = "t_VF"
                    file_to_filter = file_to_filter.dropna(subset=[vaf_colname])

                else:
                    vaf_colname = "t_AF"
                    file_to_filter = file_to_filter.dropna(subset=[vaf_colname])

                file_to_filter[vaf_colname] = pd.to_numeric(file_to_filter[vaf_colname])

                if "n" in filters:
                    t_VAF_min_novel = float(config.get("Filters", "t_VAF_min_novel"))
                    file_to_filter = file_to_filter.dropna(subset=["dbSNP_RS"])

                    file_to_filter[vaf_colname] = pd.to_numeric(file_to_filter[vaf_colname])
                    file_to_filter_nov = file_to_filter[(file_to_filter[vaf_colname] > t_VAF_min_novel) & (file_to_filter[vaf_colname] <= t_VAF_max) & (file_to_filter["dbSNP_RS"]=="novel")]
                    file_to_filter_notnov = file_to_filter[(file_to_filter[vaf_colname] > t_VAF_min) & (file_to_filter[vaf_colname] <= t_VAF_max) & (file_to_filter["dbSNP_RS"]!="novel")]
                    file_to_filter = pd.concat([file_to_filter_nov, file_to_filter_notnov], axis=0)

                else:
                    file_to_filter = file_to_filter[(file_to_filter[vaf_colname] > t_VAF_min) & (file_to_filter[vaf_colname] <= t_VAF_max)]

            if "a" in filters:
                af=config.get("Filters", "AF")
                drop_NA = config.get("Filters", "drop_NA_AF")
                drop_NA = check_bool(drop_NA)
                file_to_filter["AF"] = pd.to_numeric(file_to_filter["AF"], errors="coerce")
                na_file_to_filter = file_to_filter[~file_to_filter.index.isin(file_to_filter["AF"].dropna().index)]

                file_to_filter = file_to_filter.dropna(subset=["AF"])
                file_to_filter = file_to_filter[eval(f"file_to_filter['AF'] {af}")]

                if not drop_NA:
                        file_to_filter = pd.concat([file_to_filter, na_file_to_filter], ignore_index=True)

            if "c" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_CLIN_SIG,axis=1)]

            if "q" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_consequences,axis=1)]

            if "y" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_polyphen,axis=1)]

            if "s" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_sift,axis=1)]

            file_out = out_filter / Path(file).name
            file_to_filter.to_csv(file_out, sep="\t", index=False)

    logger.success("Filter script completed!\n")
