#####################################
# NAME: filter_clinvar.py
# AUTHOR: Luciano Giaco'
# Date: 23/01/2023
version = "1.0"
# ===================================

import os
import ast
import argparse
import concatenate
import pandas as pd
from configparser import ConfigParser
from loguru import logger
import sys
import shutil

config = ConfigParser()
configFile = config.read("conf.ini")


def print_unique_clin_sig(df):
    unique_clin_sig = df['CLIN_SIG'].unique()
    print(unique_clin_sig)

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

def filter_vf(df):
    t_vaf=float(config.get('Filters', 't_VAF'))
    gnomAD=float(config.get('Filters', 'gnomAD'))
    df = df[(df['t_VF'] > t_vaf) | (df['t_VF'].isnull())]
    df = df[(df['gnomAD_AF'] <gnomAD) | (df['gnomAD_AF'].isnull())]
    return df

def filter_main(input,folder, output_folder, vus,oncokb,cancer, overwrite=False, log=False):
    if not log:
        logger.remove()
        logfile="filter_main_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")#,mode="w")
    
    logger.info("Starting filter_main script:")
    logger.info(f"filter_main args [maf_folder:{folder}, output_folder:{output_folder}, vus:{vus}, overwrite:{overwrite}]")

    if os.path.exists(os.path.join(output_folder,'MAF_OncoKB')) and len(os.listdir(os.path.join(output_folder,'MAF_OncoKB')))>0:
        if overwrite:
            logger.warning(f"It seems that the folder 'MAF_OncoKB' already exists. Start removing process...")        
            shutil.rmtree(os.path.join(output_folder,'MAF_OncoKB'))
        else:
            logger.critical(f"The folder 'MAF_OncoKB' already exists. To overwrite an existing folder add the -w option!")
            logger.critical(f"Exit without completing the task!")
            exit()

    # if os.path.exists(os.path.join(output_folder,'NoBenign')) and len(os.listdir(os.path.join(output_folder,'NoBenign')))>0:
    #     if overwrite:
    #         logger.warning(f"It seems that the folder 'NoBenign' already exists. Start removing process...")        
    #         shutil.rmtree(os.path.join(output_folder,'NoBenign'))
    #     else:
    #         logger.critical(f"The folder 'NoBenign' already exists. To overwrite an existing folder add the -w option!")
    #         logger.critical(f"Exit without completing the task!")
    #         raise(Exception('Exiting from filter_clinvar script!'))
    
    # if os.path.exists(os.path.join(output_folder,'NoVus')) and len(os.listdir(os.path.join(output_folder,'NoVus')))>0:
    #     if overwrite:
    #         logger.warning(f"It seems that the folder 'NoVus' already exists. Start removing process...")       
    #         shutil.rmtree(os.path.join(output_folder,'NoVus'))
    #     else:
    #         logger.critical(f"The folder 'NoVus' already exists. To overwrite an existing folder add the -w option!")
    #         logger.critical(f"Exit without completing the task!")
    #         raise(Exception('Exiting from filter_clinvar script!'))

    file_list = concatenate.get_files_by_ext(folder, 'maf')
    if len(file_list)==0:
        logger.warning(f"The maf folder {os.path.join(folder, 'maf')} seems to be empty! Filtering cannot be done.")
        logger.critical("Empty maf folder: Filter script exited before completing!")
        raise(Exception("Exiting from filter_clinvar script!"))
    
    out_folders=[]
    extensions=[]
    
    if oncokb:
        output_onco=os.path.join(output_folder, 'MAF_OncoKB')
        os.mkdir(output_onco)
        extension="_OncoAnnotated.maf"
       
        if not os.path.isfile(input):
            tsv_file=[file for file in os.listdir(input) if "tsv" in file][0]
            input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
            
        else:
            input_file=pd.read_csv(input,sep="\t")
    
        
        for f in file_list:
            root, file = os.path.split(f)
            file_No = file.replace('.maf','') + extension
            file_path = os.path.join(output_onco, file_No)
            if "ONCOTREE_CODE" in input_file.columns:
                
                for _ ,row in input_file.iterrows():
                    if row["SampleID"] in file_No:
                        cancer_onco=row["ONCOTREE_CODE"]
                        
                        os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f}\
                                -o {file_path} -t {cancer_onco.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
            else:               
                os.system(f"python3 oncokb-annotator/MafAnnotator.py -i {f}\
                            -o {file_path} -t {cancer.upper()} -b {config.get('OncoKB', 'ONCOKB')}")
        


    if oncokb:   
        file_list = concatenate.get_files_by_ext(output_onco, 'maf')
        
        os.mkdir(os.path.join(output_folder, 'MAF_Onco_filtered'))
        for file in file_list:
            file_to_filter=pd.read_csv(file,sep="\t")
            file_to_filter= file_to_filter[file_to_filter["ONCOGENIC"].isin(["Oncogenic","Likely Oncogenic"])]
           
            file_to_filter=file_to_filter[file_to_filter["FILTER"]=="PASS"]
            
            file_to_filter=file_to_filter[~file_to_filter["IMPACT"].isin(["LOW","MODIFIER"])]
            file_to_filter=file_to_filter[file_to_filter["t_AF"]>=0.02]
            
            file_to_filter.to_csv(os.path.join(output_folder,"MAF_Onco_filtered",file.split("/")[-1]),sep="\t",index=False)   
    else:
        file_list = concatenate.get_files_by_ext(folder, 'maf')

    # out_folders.append(os.path.join(output_folder, 'NoBenign'))    
    # extensions.append("_NoBenign.maf")

    # out_folders.append(os.path.join(output_folder, 'NoBenign'))
    # extensions.append("_NoBenign.maf")
    # if vus:
    #     out_folders.append(os.path.join(output_folder, 'NoVus'))
    #     extensions.append('_NoVus.maf')

    # for out_folder,extension in zip(out_folders,extensions):

    #     if os.path.exists(out_folder):
    #         pass
    #     else:
    #         logger.info(f"Creating folder {out_folder}...")
    #         os.mkdir(out_folder)

    #     for f in file_list:
    #         root, file = os.path.split(f)
    #         file_No = file.replace('.maf','') + extension
    #         file_path = os.path.join(out_folder, file_No)

    #         if os.path.isfile(file_path):
    #             logger.warning(f"Skipping {file_path}: already filtered!")
    #             continue
    #         else:
    #             logger.info(f"Filtering file {f}")
                
    #             data = pd.read_csv(f, sep='\t', comment="#")
                
    #             if out_folder == os.path.join(output_folder, 'NoVus'):
    #                 try:
    #                     df = keep_risk_factors(data)
    #                 except KeyError as e:
    #                     logger.error(f"{e} key value not found: check your maf file!")
    #                     continue
    #                 except Exception as e:
    #                     logger.error(f"Something went wrong!")
    #                     continue
    #             else:
    #                 try:
    #                     df = filter_benign(data)
    #                 except KeyError as e:
    #                     logger.error(f"{e} key value not found: check your maf file!")
    #                     continue
    #                 except Exception as e:
    #                     logger.error(f"Something went wrong! Cannot create {file_No}")
    #                     continue
      
    #             filtered_data = filter_vf(df)
    #             logger.info(f"Filtered file: {file_path}")
    #             write_csv_with_info(filtered_data, file_path)
    
    logger.success("Filter script completed!\n")

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':
    
    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='filter out the variants with benign and likely benign from maf file \
                                                The filter is on clinVar annotation')

    parser.add_argument('-i', '--input', required=True,
                                            help='Path folder containing the maf files')
    parser.add_argument('-f', '--folder', required=True,
                                            help='Path folder containing the maf files')
    parser.add_argument('-o', '--output_folder', required=True,
                                            help='Output folder')
    parser.add_argument('-v', '--vus', required=False,
                                            action='store_true',
                                            help='Filter out VUS variants')
    parser.add_argument('-w', '--overWrite', required=False,action='store_true',
                                                help='Overwrite output folder if it exists')
    parser.add_argument('-k', '--oncoKB', required=False,action='store_true',help='OncoKB annotation')
    parser.add_argument('-c', '--Cancer', required=False,
                        help='Cancer Name')

    try:
        args = parser.parse_args()
    except Exception as err:
        logger.remove()
        logfile="filter_clinvar_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True,catch=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
        logger.critical(f"error: {err}", file=sys.stderr)
    
    folder = args.folder
    output_folder = args.output_folder
    vus = args.vus
    overwrite = args.overWrite
    input=args.input
    oncokb=args.oncoKB
    cancer=args.Cancer
    filter_main(input,folder, output_folder, vus,oncokb,cancer, overwrite=False, log=False)