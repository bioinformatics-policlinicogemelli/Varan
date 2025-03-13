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

import os
import ast
import concatenate
import pandas as pd
from configparser import ConfigParser
from loguru import logger
import sys
import shutil
import numpy as np

config = ConfigParser()
configFile = config.read("conf.ini")

def check_bool(key_value):
    bool_key_value = key_value
    if bool_key_value.strip() in ["True", "true", "T"]:
        bool_key_value = True
    elif bool_key_value.strip() in ["False", "false", "F", ""]:
        bool_key_value = False
    else:
        logger.critical(f"Please insert a boolean value in {bool_key_value} section of conf.ini: accepted values are [\"True\", \"true\", \"T\", \"False\", \"false\", \"F\", \"\"]")
        raise ValueError("Check again the compilation conf.ini")
    return bool_key_value


def print_unique_clin_sig(df):
    unique_clin_sig = df['CLIN_SIG'].unique()
    print(unique_clin_sig)
    
def filter_OncoKB(df):
    oncokb_filter=ast.literal_eval(config.get('Filters', 'ONCOKB_FILTER'))
    
    df_filtered=df[df["ONCOGENIC"].isin(oncokb_filter)]
    return df_filtered


def check_CLIN_SIG(row):
    clin_sig=ast.literal_eval(config.get('Filters', 'CLIN_SIG'))
    output=[]
    for _e in str(row["CLIN_SIG"]).split(","):
        if _e not in clin_sig:
            output.append(True)
        else:
            output.append(False)
    return any(output)
        
def check_consequences(row):
    consequences=ast.literal_eval(config.get('Filters', 'CONSEQUENCES'))
    output=[]
    for _e in str(row['Consequence']).split(","):
        if _e in consequences:
            output.append(True)
        else:
            output.append(False)
    return any(output)

def check_polyphen(row):
    consequences=ast.literal_eval(config.get('Filters', 'POLYPHEN'))
    output=[]
    if str(row['PolyPhen']).split("(")[0] in consequences:
        output.append(True)
    else:
        output.append(False)
    return any(output)

def check_sift(row):
    consequences=ast.literal_eval(config.get('Filters', 'SIFT'))
    output=[]
    if str(row['SIFT']).split("(")[0] in consequences:
        output.append(True)
    else:
        output.append(False)
    return any(output)


def write_csv_with_info(df, file_path):
    f = open(file_path, 'w')
    f.write('#version 2.4\n')
    f.close()
    df.to_csv(file_path, sep='\t', index=False, header=True, mode='w')


def filter_main(input,folder, output_folder, oncokb, filters, cancer, resume, overwrite=False):
    
    logger.info("Starting filter_main script:")
    logger.info(f"filter_main args [maf_folder:{folder}, output_folder:{output_folder}, filters:{filters}, cancer:{cancer}, overwrite:{overwrite}]")
    
    if os.path.exists(os.path.join(output_folder,'MAF_OncoKB')) and len(os.listdir(os.path.join(output_folder,'MAF_OncoKB')))>0:
        if overwrite:
            logger.warning(f"It seems that the folder 'MAF_OncoKB' already exists. Start removing process...")        
            shutil.rmtree(os.path.join(output_folder,'MAF_OncoKB'))
        elif not resume:
            logger.critical(f"The folder 'MAF_OncoKB' already exists. To overwrite an existing folder add the -w option!")
            logger.critical(f"Exit without completing the task!")
            sys.exit()

    file_list = concatenate.get_files_by_ext(os.path.join(folder, "maf"), 'maf')

    if len(file_list)==0:
        logger.critical(f"The maf folder {os.path.join(folder, 'maf')} seems to be empty: check if SNV folder exists and contains the VCF files. If you don't have SNV vcf, use -t option to select specific analysis.")
        raise(Exception("Exiting from filter_clinvar script!"))
    
    if oncokb:
        output_onco=os.path.join(output_folder, 'MAF_OncoKB')
        os.makedirs(output_onco, exist_ok=True)
        extension="_OncoAnnotated.maf"
        input_file=pd.read_csv(input, sep="\t")
    
        for f in file_list:
            if extension in f:
                continue
            _, file = os.path.split(f)
            file_No = file.replace('.maf','') + extension
            file_path = os.path.join(output_onco, file_No)
            if "ONCOTREE_CODE" in input_file.columns:
                
                for _ ,row in input_file.iterrows():
                    if row["SAMPLE_ID"] in file_No:
                        cancer_onco=row["ONCOTREE_CODE"]
                        if cancer_onco == "":
                            cancer_onco = cancer
                        os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f}\
                                -o {file_path} -t {cancer_onco.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
            else:               
                os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f}\
                            -o {file_path} -t {cancer.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
        
    file_list = concatenate.get_files_by_ext(os.path.join(folder,"maf"), 'maf')
    out_filter=os.path.join(output_folder, 'MAF_filtered')
    
    if oncokb and "o" in filters:   
        file_list = concatenate.get_files_by_ext(output_onco, 'maf')
        out_filter=os.path.join(output_folder, 'MAF_Onco_filtered')
    
    if filters!="" and filters!="d":
        
        logger.info("Start filtering MAF...")
        os.makedirs(out_filter, exist_ok=True)
        
        for file in file_list:
            file_to_filter=pd.read_csv(file, sep="\t", dtype=object)
            if "C" in filters:
                file_to_filter = file_to_filter.dropna(subset=["t_alt_count"])
                file_to_filter["t_alt_count"] = file_to_filter["t_alt_count"].astype(int)
                file_to_filter = file_to_filter[(file_to_filter["Variant_Classification"].isin(["Missense_Mutation", "Nonsense_Mutation", "Splice_Region", "Splice_Site"])) & (file_to_filter["t_alt_count"] > 3) |
                 ~(file_to_filter["Variant_Classification"].isin(["Missense_Mutation", "Nonsense_Mutation", "Splice_Region", "Splice_Site"])) & (file_to_filter["t_alt_count"] > 10)]
            
            if "i" in filters:
                file_to_filter = file_to_filter[~file_to_filter["IMPACT"].isin(ast.literal_eval(config.get('Filters',"IMPACT")))]    
            
            if "p" in filters:
                file_to_filter=file_to_filter[file_to_filter["FILTER"]=="PASS"]
            
            if oncokb and "o" in filters:
                oncokb_filter=ast.literal_eval(config.get('Filters', 'ONCOKB_FILTER'))
                file_to_filter= file_to_filter[file_to_filter["ONCOGENIC"].isin(oncokb_filter)]
            
            if "v" in filters:
                t_VAF_min=float(config.get('Filters', 't_VAF_min'))
                t_VAF_max=float(config.get('Filters', 't_VAF_max'))
                
                temp = file_to_filter.dropna(subset=["t_AF"])
        
                if len(temp) == 0:
                    vaf_colname = "t_VF"
                    file_to_filter.dropna(subset=[vaf_colname], inplace=True)

                else:
                    vaf_colname = "t_AF"
                    file_to_filter.dropna(subset=[vaf_colname], inplace=True)

                file_to_filter[vaf_colname] = pd.to_numeric(file_to_filter[vaf_colname])
                
                if "n" in filters:
                    t_VAF_min_novel=float(config.get('Filters', 't_VAF_min_novel'))
                    file_to_filter.dropna(subset=["dbSNP_RS"], inplace=True)

                    file_to_filter[vaf_colname] = pd.to_numeric(file_to_filter[vaf_colname])
                    file_to_filter_nov = file_to_filter[(file_to_filter[vaf_colname] > t_VAF_min_novel) & (file_to_filter[vaf_colname] <= t_VAF_max) & (file_to_filter["dbSNP_RS"]=="novel")]
                    file_to_filter_notnov = file_to_filter[(file_to_filter[vaf_colname] > t_VAF_min) & (file_to_filter[vaf_colname] <= t_VAF_max) & (file_to_filter["dbSNP_RS"]!="novel")]
                    file_to_filter = pd.concat([file_to_filter_nov, file_to_filter_notnov], axis=0)
                    
                else:
                    file_to_filter = file_to_filter[(file_to_filter[vaf_colname] > t_VAF_min) & (file_to_filter[vaf_colname] <= t_VAF_max)]
            
            if "a" in filters:
                af=config.get('Filters', 'AF')
                drop_NA = config.get('Filters', 'drop_NA_AF')
                drop_NA = check_bool(drop_NA)
                file_to_filter['AF'] = pd.to_numeric(file_to_filter['AF'], errors='coerce')
                na_file_to_filter = file_to_filter[~file_to_filter.index.isin(file_to_filter["AF"].dropna().index)]
                
                file_to_filter.dropna(subset=["AF"], inplace=True)
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
                
            file_to_filter.to_csv(os.path.join(out_filter, os.path.basename(file)),sep="\t",index=False)  

    logger.success("Filter script completed!\n")