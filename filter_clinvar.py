#####################################
# NAME: filter_clinvar.py
# AUTHOR: Luciano Giaco'
# Date: 23/01/2023
version = "1.0"
# ===================================

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

def print_unique_clin_sig(df):
    unique_clin_sig = df['CLIN_SIG'].unique()
    print(unique_clin_sig)
    
def filter_OncoKB(df):
    oncokb_filter=ast.literal_eval(config.get('Filters', 'ONCOKB_FILTER'))
    
    df_filtered=df[df["ONCOGENIC"].isin(oncokb_filter)]
    return df_filtered

def filter_benign(df):
    benign_filter = ~df['CLIN_SIG'].str.contains(config.get('Filters', 'BENIGN')
            , case=False
            , na=False
            , regex=True)
    return df[benign_filter]

def check_CLIN_SIG(row):
    clin_sig=ast.literal_eval(config.get('Filters', 'CLIN_SIG'))
    output=[]
    for _e in str(row["CLIN_SIG"]).split(","):
        if _e in clin_sig:
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
    for _e in str(row['PolyPhen']).split(","):
        if any(_e in s for s in consequences):
            output.append(True)
        else:
            output.append(False)
    return any(output)

def check_sift(row):
    consequences=ast.literal_eval(config.get('Filters', 'SIFT'))
    output=[]
    for _e in str(row['SIFT']).split(","):
        if any(_e in s for s in consequences):
            output.append(True)
        else:
            output.append(False)
    return any(output)

def keep_risk_factors(df):
    benign_filter = ~df['CLIN_SIG'].str.contains(config.get('Filters', 'BENIGN')
        , case=False
        , na=False
        , regex=True)
    df=df[benign_filter]
    df = df[df.apply(check_CLIN_SIG,axis=1)|df.apply(check_consequences,axis=1)]
    return df


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
        logger.warning(f"The maf folder {os.path.join(folder, 'maf')} seems to be empty! Filtering cannot be done.")
        logger.critical("Empty maf folder: Filter script exited before completing!")
        raise(Exception("Exiting from filter_clinvar script!"))
    
    if oncokb:
        output_onco=os.path.join(output_folder, 'MAF_OncoKB')
        os.makedirs(output_onco, exist_ok=True)
        extension="_OncoAnnotated.maf"
        if not os.path.isfile(input[0]):
            input_file=pd.read_csv(os.path.join(input[0], "sample.tsv"), sep="\t")
            
        else:
            input_file=pd.read_csv(input,sep="\t")
    
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
        
        logger.info("Start filtering vcf...")
        
        os.makedirs(out_filter, exist_ok=True)
        
        for file in file_list:
            file_to_filter=pd.read_csv(file,sep="\t")
            
            if "i" in filters:
                file_to_filter = file_to_filter[~file_to_filter["IMPACT"].isin(ast.literal_eval(config.get('Filters',"IMPACT")))]    
            
            if "p" in filters:
                file_to_filter=file_to_filter[file_to_filter["FILTER"]=="PASS"]
            
            if "b" in filters:
                benign_filter = ~file_to_filter['CLIN_SIG'].str.contains(config.get('Filters', 'BENIGN')
                    , case=False
                    , na=False
                    , regex=True)   
                file_to_filter=file_to_filter[benign_filter]
            
            if oncokb and "o" in filters:
                oncokb_filter=ast.literal_eval(config.get('Filters', 'ONCOKB_FILTER'))
                file_to_filter= file_to_filter[file_to_filter["ONCOGENIC"].isin(oncokb_filter)]
            
            if "v" in filters:
                t_VAF_min=float(config.get('Filters', 't_VAF_min'))
                t_VAF_max=float(config.get('Filters', 't_VAF_max'))
                
                file_to_filter.dropna(subset=["t_AF"], inplace=True)

                if "n" in filters:
                    t_VAF_min_novel=float(config.get('Filters', 't_VAF_min_novel'))
                    file_to_filter.dropna(subset=["dbSNP_RS"], inplace=True)

                    if file_to_filter[(file_to_filter["dbSNP_RS"]=="novel")]:
                        file_to_filter = file_to_filter[(file_to_filter['t_AF'] > t_VAF_min_novel) & (file_to_filter['t_AF'] <= t_VAF_max)]
                    else:
                        file_to_filter = file_to_filter[(file_to_filter['t_AF'] > t_VAF_min) & (file_to_filter['t_AF'] <= t_VAF_max)]
                    
                else:
                    file_to_filter = file_to_filter[(file_to_filter['t_AF'] > t_VAF_min) & (file_to_filter['t_AF'] <= t_VAF_max)]
            # usato per le varianti introniche
            if "a" in filters:    
                    af=config.get('Filters', 'AF')
                    
                    file_to_filter.dropna(subset=["AF"], inplace=True)
                    file_to_filter = file_to_filter[(eval("file_to_filter['AF']" + af))]
                
            if "c" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_CLIN_SIG,axis=1)]
                    
            if "q" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_consequences,axis=1)]
                
            if "y" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_polyphen,axis=1)]
            
            if "s" in filters:
                file_to_filter = file_to_filter[file_to_filter.apply(check_sift,axis=1)]
                 
                
            #logger.info(f"Filtered file: {file}")             
            
            file_to_filter.to_csv(os.path.join(out_filter, os.path.basename(file)),sep="\t",index=False)  

    logger.success("Filter script completed!\n")