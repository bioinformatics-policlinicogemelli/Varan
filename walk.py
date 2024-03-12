#####################################
# NAME: walk.py
# Date: 10/01/2023
version = "1.0"
# ===================================

from operator import index
import os
import ast
import argparse
import subprocess
import vcf2tab_cnv
import vcf_filter
import tsv
import pandas as pd
from configparser import ConfigParser
import shutil
from loguru import logger 
import random
import string
import sys
import traceback

config = ConfigParser()

configFile = config.read("conf.ini")
OUTPUT_FILTERED = config.get('Paths', 'OUTPUT_FILTERED')
OUTPUT_MAF = config.get('Paths', 'OUTPUT_MAF')
VCF2MAF = config.get('Paths', 'VCF2MAF')
REF_FASTA = config.get('Paths', 'REF_FASTA')
TMP = config.get('Paths', 'TMP')
VEP_PATH = config.get('Paths', 'VEP_PATH')
VEP_DATA = config.get('Paths', 'VEP_DATA')
CLINV = config.get('Paths', 'CLINV')
CNA=ast.literal_eval(config.get('Cna', 'HEADER_CNV'))

def create_random_name_folder():
    nome_cartella = ''.join(random.choices(string.ascii_lowercase + string.digits, k=10))
    temporary = os.path.join(TMP, nome_cartella)
    try:
        os.mkdir(temporary)
    except FileNotFoundError:
        logger.critical(f"Scratch folder '{TMP}' not found! Check TMP field in conf.ini")
        raise(FileNotFoundError("Error in create_random_name_folder: exiting from walk script!"))
    except Exception:
        logger.critical("Something went wrong while creating the vep tmp folder")
        raise(Exception("Error in create_random_name_folder: exiting from walk script!"))
    return(temporary)

def clear_scratch():
    for root, dirs, files in os.walk(TMP):
        for dir in dirs:
            shutil.rmtree(os.path.join(root,dir))
        

def clear_temp(folder):
    shutil.rmtree(os.path.join(folder,"temp"))

def get_cnv_from_folder(inputFolderCNV):
    files= os.listdir(inputFolderCNV)
    cnv_vcf_files=[file for file in files if file.endswith("vcf")]
    # check=list(map(lambda x: check_cna_vcf(x,inputFolderCNV),cnv_vcf_files))
    # incorrect_files=[]
    # for i, check_res in enumerate(check):
    #     if not check_res:
    #         incorrect_files.append(cnv_vcf_files[i])
    # if len(incorrect_files)!=0:
    #     logger.critical(f"It seems that the files \n{incorrect_files} \nare not CNV! Please check your CNV input data and try again.")
    #     raise Exception("Error in get_cnv_from_folder script: exiting from walk script!")       
    # logger.info(f"#{len(cnv_vcf_files)} vcf files found in CNV folder")
    return cnv_vcf_files

def get_sampleID_from_cnv(cnv_vcf):
    if "_CopyNumberVariants.vcf" in cnv_vcf:
        sample=cnv_vcf.replace("_CopyNumberVariants.vcf",".bam")
    else:
        sample=cnv_vcf.replace("vcf","bam")
    return sample




def reshape_cna(input,cna_df_path,cancer,output_dir):
   
    if not os.path.isfile(input):
        tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
        input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")    
    else:
        input_file=pd.read_csv(input,sep="\t")

    cna_df=pd.read_csv(cna_df_path,sep="\t")
    
    cna_df.rename({"ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"},inplace=True,axis=1)
    input_file.rename({"SampleID":"Tumor_Sample_Barcode"}, inplace=True,axis=1)

    if not "ONCOTREE_CODE" in input_file.columns:
        input_file["ONCOTREE_CODE"]=cancer
  
    input_file["Tumor_Sample_Barcode"]=input_file["Tumor_Sample_Barcode"]+".cnv.bam"
    #cna_df["Tumor_Sample_Barcode"]=cna_df["Tumor_Sample_Barcode"].str.replace(".bam","")
    annotate= pd.merge(cna_df[["Tumor_Sample_Barcode","Hugo_Symbol","discrete","Copy_Number_Alteration"]] ,input_file[["Tumor_Sample_Barcode","ONCOTREE_CODE"]],on="Tumor_Sample_Barcode")
   
    
  
    annotate.to_csv(os.path.join(output_dir,"temp_cna.txt"),index=False,sep="\t")
    
    return os.path.join(output_dir,"temp_cna.txt")


def annotate_cna(path_cna,output_folder):

    out=path_cna.replace(".txt","2.txt")
    os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {path_cna}\
                        -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")
                
    
    cna=pd.read_csv(out,sep="\t",dtype={"Copy_Number_Alteration":int})
    cna=cna[cna["ONCOGENIC"].isin(["Oncogenic","Likely Oncogenic"])]
    
    
    
    data_cna=cna.pivot_table(index="Hugo_Symbol",columns="Tumor_Sample_Barcode",values="Copy_Number_Alteration",fill_value=0)

    data_cna.to_csv(os.path.join(output_folder,"data_cna.txt"),index=True,sep="\t")
    




def cnv_type_from_folder(input,cnv_vcf_files,output_folder,oncokb,cancer):
    c = 0
    sID_path = dict()

    for case_folder in cnv_vcf_files:
        if os.path.exists('data_cna_hg19.seg'):
            MODE = 'a'
        else:
            MODE = 'w'
        try:
            cnv_vcf=case_folder 
            sampleID = get_sampleID_from_cnv(case_folder)
            if sampleID in sID_path:
                dup = open('sampleID_dup'.log, 'w')
                dup.write(sampleID+'\t'+'cnv_vcf')
                dup.close()
            else:
                sID_path[sampleID] = os.path.join(os.path.join(input,"CNV"),cnv_vcf)
                vcf2tab_cnv.vcf_to_table(sID_path[sampleID], os.path.join(output_folder,'data_cna_hg19.seg'), sampleID, MODE)
                vcf2tab_cnv.vcf_to_table_fc(sID_path[sampleID], os.path.join(output_folder,'data_cna_hg19.seg.fc.txt'), sampleID, MODE)
               
        except Exception:
            log_noparsed = open('noParsed_cnv.log', 'a')
            log_noparsed.write('[WARNING] '+case_folder+'\n')
            log_noparsed.close()
        
        c = c +1
    logger.info("Writing data_cna_hg19.seg succefully completed!")
    logger.info("Writing data_cna_hg19.seg.fc.txt succefully completed!")
    
    ############################
    ### MANAGE DISCRETE TABLE ##
    ############################
    
    df_table = pd.read_csv(os.path.join(output_folder,'data_cna_hg19.seg.fc.txt'),sep="\t",header=0)
    result = tabella_to_dict(df_table)
    
    
    if oncokb:
        # if not os.path.isfile(input):
        #     tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
        #     input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
        
        # else:
        #     input_file=pd.read_csv(input,sep="\t")
        
        # for index, row in df_table.iterrows():
        #     tc= int(input_file[input_file["SampleID"]+".cnv.bam"==row["ID"]]["TC"])
        #     df_table.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc
            
        
        df_table_filtered=df_table#[(df_table["discrete"]>=3)|(df_table["discrete"]<=0.8)] 
        
        #df_table_filtered["Copy_Number_Alteration"]=0
        # df_table_filtered.loc[(df_table_filtered["discrete"]>3)&(df_table_filtered["discrete"]<5), "Copy_Number_Alteration"]=1
        # df_table_filtered.loc[df_table_filtered["discrete"]>5, "Copy_Number_Alteration"]=2
        # df_table_filtered.loc[(df_table_filtered["discrete"]>0)&(df_table_filtered["discrete"]<0.8), "Copy_Number_Alteration"]=-1
        # df_table_filtered.loc[df_table_filtered["discrete"]<=0, "Copy_Number_Alteration"]=-2
        
       
        df_table_filtered["Copy_Number_Alteration"]=df_table_filtered["discrete"]
        df_table_filtered.to_csv(os.path.join(output_folder,"data_cna_hg19.seg.filtered.txt"),sep="\t",index=False)
        temp_cna=reshape_cna(input,os.path.join(output_folder,"data_cna_hg19.seg.filtered.txt"),cancer,output_folder)
        annotate_cna(temp_cna,output_folder)
   
    return sID_path

def tabella_to_dict(df):
    result = {}
    for index, row in df.iterrows():
        row_values = (row['chrom'], row['loc.start'], row['loc.end'], row['num.mark'], row['seg.mean'], row['gene'], row['discrete'])
        if row['ID'] not in result:
            result[row['ID']] = []
        result[row['ID']].append(row_values)
    return result

def get_snv_from_folder(inputFolderSNV):
    files= os.listdir(inputFolderSNV)
    
    snv_vcf_files=[file for file in files if file.endswith("vcf")]
    # check=list(map(lambda x: check_snv_vcf(x,inputFolderSNV),snv_vcf_files))
    # incorrect_files=[]
    # for i, check_res in enumerate(check):
    #     if not check_res:
    #         incorrect_files.append(snv_vcf_files[i])
    # if len(incorrect_files)!=0:
    #     logger.critical(f"It seems that the files \n{incorrect_files} \nare not SNV! Please check your SNV input data and try again.")
    #     raise Exception("Error in get_snv_from_folder: exiting from walk script")   
    # logger.info(f"#{len(snv_vcf_files)} vcf files found in SNV folder")
    return snv_vcf_files

def get_sampleID_from_snv(snv_vcf):
    if "MergedSmallVariants.genome.vcf" in snv_vcf:
        sample=snv_vcf.replace("_MergedSmallVariants.genome.vcf",".bam")
    else:
        sample=snv_vcf.replace("vcf","bam")
    return sample

def snv_type_from_folder(input,snv_vcf_files):
    c = 0
    sID_path = dict()
    for case_folder in snv_vcf_files:
        try:
            snv_vcf= case_folder
            sampleID =get_sampleID_from_snv(case_folder)
            if sampleID in sID_path:
                dup = open('sampleID_dup'.log, 'w')
                dup.write(sampleID+'\t'+'snv_vcf')
                dup.close()
            else:
                sID_path[sampleID] = os.path.join(input,snv_vcf)
        except Exception:
            log_noparsed = open('noParsed_snv.log', 'a')
            log_noparsed.write('[WARNING]'+case_folder+'\n')
            log_noparsed.close()
        c = c +1
    return sID_path

def vcf_filtering(sID_path,output_folder):
    sID_path_filtered = dict()
    for k, v in sID_path.items():
        root, vcf_file = os.path.split(v)
        out_filt=os.path.join(output_folder,OUTPUT_FILTERED) #TEST
        vcf_filtered = os.path.join(out_filt, vcf_file.replace('.vcf','')+'.FILTERED.vcf')
        logger.info(f'[FILTERING] {v}')
        vcf_filter.main(v, vcf_filtered)
        logger.info(f'filtered file {vcf_filtered} created!')
        sID_path_filtered[k] = vcf_filtered
    return sID_path_filtered

def vcf2maf_constructor(k, v, temporary,output_folder):
    tum_id=k.replace(".hard-filtered.bam","")
    
    cl = ['perl']
    cl.append(VCF2MAF)
    cl.append('--input-vcf')
    cl.append(v)
    root, file_vcf = os.path.split(v)
    out_file = os.path.join(output_folder,os.path.join(OUTPUT_MAF, file_vcf+'.maf'))
    cl.append('--vep-custom')
    cl.append(CLINV)
    cl.append('--output-maf')
    cl.append(out_file)
    cl.append('--ref-fasta')
    cl.append(REF_FASTA)
    cl.append('--tmp-dir')
    cl.append(temporary)
    cl.append('--retain-fmt')
    cl.append('GT,GQ,AD,DP,VF,AF')
    cl.append('--vep-path')
    cl.append(VEP_PATH)
    cl.append('--vep-data')
    cl.append(VEP_DATA)
    cl.append('--tumor-id') 
    cl.append(tum_id)
    cl.append(" --cache-version 111")
    
    return cl

def run_vcf2maf(cl):
    logger.info('Starting vcf2maf conversion...')
    logger.info(f'args={cl}')
    sout = subprocess.run(cl, capture_output=True)
 
    if sout.stderr!=None:
        if 'ERROR' not in sout.stderr.decode('ascii'):
            logger.warning(sout.stderr.decode('ascii').replace('ERROR: ',''))
        else:
            logger.error(sout.stderr.decode('ascii').replace('ERROR: ',''))
    
def create_folder(output_folder,overwrite_output):
    if overwrite_output:
        if os.path.exists(output_folder):
            logger.warning(f"It seems that the folder '{output_folder}' already exists. Start removing process...")
            shutil.rmtree(output_folder)
    
    if not os.path.exists(output_folder):
        logger.info(f"Creating the output folder '{output_folder}' in {os.getcwd()}...")
        os.mkdir(output_folder)
        maf_path = os.path.join(output_folder, 'maf')
        os.mkdir(maf_path)
        filtered_path = os.path.join(output_folder, 'snv_filtered')
        os.mkdir(filtered_path)
        logger.info(f"The folder '{output_folder}' was correctly created!")

    else:
        logger.critical(f"The folder '{output_folder}' already exists. To overwrite an existing folder add the -w option!")
        raise(Exception('Error in create_folder script: exiting from walk script!'))

def get_table_from_folder(tsvpath):
    table_dict = dict()
    file=pd.read_csv(tsvpath,sep="\t",index_col=False, dtype=str)
    for _, row in file.iterrows():
        sampleID=str(row["SampleID"])
        if ".bam" in sampleID:
           sampleID=sampleID.replace(".bam","")
        if sampleID not in table_dict.keys():
            table_dict[sampleID]=[str(row["PatientID"])]  
    return table_dict
    
def flatten(nested_list):
    flat_list = []
    for sublist in nested_list:
        for item in sublist:
            flat_list.append(item)
    return flat_list

def write_clinical_patient(output_folder, table_dict):
    logger.info("Writing data_clinical_patient.txt file...")
    data_clin_samp = os.path.join(output_folder,'data_clinical_patient.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tAge\tGender\n')
    cil_sample.write('#Patient identifier\tAge\tGender\n')
    cil_sample.write('#STRING\tNUMBER\tSTRING\n')
    cil_sample.write('#1\t1\t1\n')
    cil_sample.write('PATIENT_ID\tAGE\tGENDER\n')

    nested_list = list(table_dict.values())
    list_patients = flatten(nested_list)
    list_patients = set(list_patients)
    for v in list_patients:
        cil_sample.write(v+"\tNaN\tNa\n")
    cil_sample.close()

def write_clinical_sample(output_folder, table_dict):
    logger.info("Writing data_clinical_sample.txt file...")
    data_clin_samp = os.path.join(output_folder, 'data_clinical_sample.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tSample Identifier\tMSI\tTMB\tMSI_THR\tTMB_THR\n')
    cil_sample.write('#Patient identifier\tSample Identifier\tMicro Satellite Instability\tMutational Tumor Burden\tMSI_THR\tTMB_THR\n')
    cil_sample.write('#STRING\tSTRING\tNUMBER\tNUMBER\tSTRING\tSTRING\n')
    cil_sample.write('#1\t1\t1\t1\t1\t1\n')
    cil_sample.write('PATIENT_ID\tSAMPLE_ID\tMSI\tTMB\tMSI_THR\tTMB_THR\n')
    for k, v in table_dict.items():
        cil_sample.write(v[0]+'\t'+k+'.bam\t'+v[1]+'\t'+v[2]+'\t'+v[3]+'\t'+v[4]+'\n')
    cil_sample.close()
    
def write_clinical_sample_empty(output_folder, table_dict):
    logger.info("Writing data_clinical_sample.txt file...")
    data_clin_samp = os.path.join(output_folder, 'data_clinical_sample.txt')
    cil_sample = open(data_clin_samp, 'w')
    cil_sample.write('#Patient Identifier\tSample Identifier\n')
    cil_sample.write('#Patient identifier\tSample Identifier\n')
    cil_sample.write('#STRING\tSTRING\n')
    cil_sample.write('#1\t1\n')
    cil_sample.write('PATIENT_ID\tSAMPLE_ID\n')
    for k, v in table_dict.items():
        cil_sample.write(v[0]+'\t'+k+'.bam'+'\n')
    cil_sample.close()

def check_cna_vcf(file,inputFolderCNV):
    vcf=pd.read_csv(os.path.join(inputFolderCNV,file),comment="#",sep="\t",names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"])
    if vcf.loc[0]["FORMAT"]=="FC":
        return True
    else:
        return False
    
def check_snv_vcf(file,inputFolderSNV):
    vcf=pd.read_csv(os.path.join(inputFolderSNV,file),comment="#",sep="\t",names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"])
    if vcf.loc[0]["FORMAT"].startswith("GT"):
        return True
    else:
        return False

def get_combinedVariantOutput_from_folder(inputFolder, tsvpath):
    combined_dict = dict()
    try:
        file = pd.read_csv(tsvpath,sep="\t",dtype=str)
    except Exception as e:
        logger.critical(f"Something went wrong while reading {tsvpath}!")
        raise(Exception("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))
    for _,row in file.iterrows():
        try:
            patientID=str(row["PatientID"])
        except KeyError as e: 
            logger.critical(f"KeyError: {e} not found! Check if column name is correctly spelled or if there are tabs/spaces before or after the coloumn key: \n{row.index}. \nThis error may also occur if the table columns have not been separated by tabs!")
            raise(KeyError("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))
        try:
            sampleID=str(row["SampleID"])
        except KeyError as e: 
            logger.critical(f"KeyError: {e} not found! Check if column name is correctly spelled or if there are tabs/spaces before or after the coloumn key: \n{row.index}. \nThis error may also occur if the table columns have not been separated by tabs!")
            raise(KeyError("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))
        combined_file = patientID+'_CombinedVariantOutput.tsv'
        combined_path = os.path.join(inputFolder,"CombinedOutput",combined_file)
        if os.path.exists(combined_path):
            pass 
        else:
            logger.warning(f"{combined_path} not exists")
        combined_dict[sampleID] = combined_path
    return combined_dict

def transform_input(tsv,output_folder):
    os.mkdir(os.path.join(output_folder,"temp"))
    os.mkdir(os.path.join(output_folder,"temp","SNV"))
    os.mkdir(os.path.join(output_folder,"temp","CNV"))
    os.mkdir(os.path.join(output_folder,"temp","CombinedOutput"))

    os.system("cp "+tsv +" "+os.path.join(output_folder,"temp"))

    tsv_file=pd.read_csv(tsv,sep="\t",dtype="string")

    for _,row in tsv_file.iterrows():
        res_folder="/data/data_storage/novaseq_results/Dragen_2.5"
        #
        snv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+".hard-filtered.vcf")     #"_MergedSmallVariants.genome.vcf")
        cnv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+".cnv.vcf") # "_CopyNumberVariants.vcf")
        combout=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["PatientID"]+"_CombinedVariantOutput.tsv")
    
        if os.path.exists(combout):
            os.system("cp "+combout + " "+os.path.join(output_folder,"temp","CombinedOutput"))
        else:
            logger.warning(f"{combout} not found")
        if os.path.exists(snv_path):
            os.system("cp "+snv_path + " "+os.path.join(output_folder,"temp","SNV"))
        else:
            logger.warning(f"{snv_path} not found")
        if os.path.exists(cnv_path):
            os.system("cp "+cnv_path + " "+os.path.join(output_folder,"temp","CNV"))
        else:
            logger.warning(f"{cnv_path} not found")
    
    return os.path.join(output_folder,"temp")

def walk_folder(input, output_folder,oncokb,cancer, overwrite_output=False, vcf_type=None ,filter_snv=False, log=False):
    if not log:
        logger.remove()
        logfile="Walk_folder_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")#,mode="w")
    
    logger.info("Starting walk_folder script:")
    logger.info(f"walk_folder args [input:{input}, output_folder:{output_folder}, Overwrite:{overwrite_output}, vcf_type:{vcf_type}, filter_snv:{filter_snv}]")
 
    config.read('conf.ini')
    
    ###############################
    ###       OUTPUT FOLDER     ###
    ###############################
 
    create_folder(output_folder,overwrite_output)
    if os.path.isfile(input):
        
        input=transform_input(input,output_folder)
            
    inputFolderSNV=os.path.abspath(os.path.join(input,"SNV"))
    inputFolderCNV=os.path.abspath(os.path.join(input,"CNV"))
 
    if os.path.exists(inputFolderCNV) and vcf_type!="snv":
        logger.info("Check CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(inputFolderCNV)
        logger.info("Everything ok!")
    if os.path.exists(inputFolderSNV) and vcf_type!="cnv":
        logger.info("Check SNV files...")
        
        case_folder_arr = get_snv_from_folder(inputFolderSNV)
        logger.info("Everything ok!")

    ###############################
    ###       SNV AND CNV       ###
    ###############################
    
    if os.path.exists(inputFolderCNV) and vcf_type!="snv":
        sID_path_cnv = cnv_type_from_folder(input,case_folder_arr_cnv,output_folder,oncokb,cancer)
    if os.path.exists(inputFolderSNV) and vcf_type!="cnv":
        sID_path_snv = snv_type_from_folder(inputFolderSNV,case_folder_arr)
        
        logger.info("Starting vcf2maf conversion...")
        if filter_snv == True:
            logger.info("filter option on")
            temporary = create_random_name_folder()
            sID_path_filtered = vcf_filtering(sID_path_snv,output_folder)
            for k, v in sID_path_filtered.items():
                cl = vcf2maf_constructor(k, v, temporary,output_folder)
                run_vcf2maf(cl)
        else:
            temporary = create_random_name_folder()
            for k, v in sID_path_snv.items():
                cl = vcf2maf_constructor(k, v, temporary,output_folder)
                run_vcf2maf(cl)
    
    logger.info("Clearing scratch folder...")
    clear_scratch()
    

    ###############################
    ###       GET FUSION        ###
    ###############################

    try:
        tsvfiles=[file for file in os.listdir(input) if file.endswith("tsv")][0]
    except IndexError:
        logger.critical(f"It seems that no tsv file is in your folder!")
        raise(IndexError("Exiting from walk script!"))
    except FileNotFoundError:
        logger.critical(f"No input directory '{input}' was found: try check your path")
        raise(FileNotFoundError("Exiting from walk script!"))
    except Exception as e:
        logger.critical(f"Something went wrong! {print(traceback.format_exc())}")
        raise(Exception("Error while reading clinical_info.tsv: exiting from walk script!"))

    try:
        tsvpath=os.path.join(input,tsvfiles)    
        combined_dict = get_combinedVariantOutput_from_folder(input,tsvpath)
        
        fusion_table_file = os.path.join(output_folder,'data_sv.txt')
            
        for k, v in combined_dict.items():
            fusions=[]
            logger.info(f"Reading Fusion info in CombinedOutput file {v}...")
            try:
                fusions = tsv.get_fusions(v)
            except Exception as e:
                logger.error(f"Something went wrong while reading Fusion section of file {v}")
            if len(fusions)==0:
                logger.info(f"No Fusions found in {v}")
                continue
            else:
                logger.info(f"Fusions found in {v}")
            if not os.path.exists(fusion_table_file):
                logger.info(f"Creating data_sv.txt file...")
                fusion_table = open(fusion_table_file, 'w')
                header = 'Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\tNormal_Paired_End_Read_Count\tEvent_Info\tRNA_Support\n'
                fusion_table.write(header)
            else:
                fusion_table = open(fusion_table_file, 'a')
            for fus in fusions:
                if len(fusions) > 0:
                    Site1_Hugo_Symbol = fus['Site1_Hugo_Symbol']
                    Site2_Hugo_Symbol = fus['Site2_Hugo_Symbol']
                    if Site2_Hugo_Symbol == 'CASC1':
                        Site2_Hugo_Symbol = 'DNAI7'
                    Site1_Chromosome = fus['Site1_Chromosome']
                    Site2_Chromosome = fus['Site2_Chromosome']
                    Site1_Position = fus['Site1_Position']
                    Site2_Position = fus['Site2_Position']

                    if int(fus['Normal_Paired_End_Read_Count'])>=15 :
                        fusion_table.write(k+'.bam\tSOMATIC\tFUSION\t'+\
        str(Site1_Hugo_Symbol)+'\t'+str(Site2_Hugo_Symbol)+'\t'+fus['Normal_Paired_End_Read_Count']+\
        '\t'+fus['Event_Info']+' Fusion\t'+'Yes\n')
            
        
        if oncokb and  os.path.exists(fusion_table_file):
            
            data_sv=pd.read_csv(fusion_table_file,sep="\t")
            if not os.path.isfile(input):
                tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
                input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
            
            else:
                input_file=pd.read_csv(input,sep="\t")

            
            if "ONCOTREE_CODE" in input_file.columns:
            
                input_file["SampleID"]=input_file["SampleID"]+".bam"
                fusion_table_df=data_sv.merge(input_file,how="inner",left_on="Sample_Id",right_on="SampleID")
                fusion_table_df.to_csv(fusion_table_file,sep="\t",index=False)
                fusion_table_file_out=fusion_table_file.replace(".txt","ann.txt")
                os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                        -o {fusion_table_file_out} -b {config.get('OncoKB', 'ONCOKB')}")
                

            else:   
                    
                fusion_table_file_out=fusion_table_file.replace(".txt","ann.txt")       
                os.system(f"python3 oncokb-annotator/FusionAnnotator.py -i {fusion_table_file}\
                            -o {fusion_table_file_out} -t {cancer.upper()}  -b {config.get('OncoKB', 'ONCOKB')}")

            fileonco=pd.read_csv(fusion_table_file_out,sep="\t")
            
        # fileonco=fileonco[fileonco["ONCOGENIC"].isin(["Oncogenic","Likely Oncogenic"])]
            fileonco.to_csv(fusion_table_file,sep="\t",index=False)           
            
            
    except:
        logger.warning("No fusion found")


    ##############################
    ##       MAKES TABLE       ###
    ##############################
    
    tsv_file=[file for file in os.listdir(input) if file.endswith("tsv")][0]
    tsvpath=os.path.join(input,tsv_file)    

    table_dict_patient = get_table_from_folder(tsvpath)

    logger.info("Writing clinical files...")
    write_clinical_patient(output_folder, table_dict_patient)
    
    try:
        combined_dict = get_combinedVariantOutput_from_folder(input,tsvpath)
        MSI_THR=config.get('MSI', 'THRESHOLD')
        TMB=ast.literal_eval(config.get('TMB', 'THRESHOLD'))
        
        for k, v in combined_dict.items():
            logger.info(f"Reading Tumor clinical parameters info in CombinedOutput file {v}...")
            try:
                tmv_msi = tsv.get_msi_tmb(v)
            except Exception as e:
                logger.error(f"Something went wrong!")
            logger.info(f"Tumor clinical parameters Values found: {tmv_msi}")
            
            if not tmv_msi['MSI'][0][1]=="NA" and float(tmv_msi['MSI'][0][1]) >= 40:
                table_dict_patient[k].append(tmv_msi['MSI'][1][1])   
            else:
                table_dict_patient[k].append('NA')
            table_dict_patient[k].append(tmv_msi['TMB_Total'])
            if not tmv_msi['MSI'][0][1]=="NA":
                if not tmv_msi['MSI'][1][1] =="NA":
                    if float(tmv_msi['MSI'][1][1]) < float(MSI_THR):
                        table_dict_patient[k].append("Stable")   
                    else:
                        table_dict_patient[k].append('Unstable')

                    found = False
                else:
                    table_dict_patient[k].append('NA')
            else:
                table_dict_patient[k].append('NA')
            found=False
            for _k, _v in TMB.items():
                if not tmv_msi["TMB_Total"]=="NA":
                    if float(tmv_msi["TMB_Total"])<float(_v):
                        table_dict_patient[k].append(_k)
                        found=True
                        break
                else:
                    table_dict_patient[k].append("NA")
            if found==False:
                table_dict_patient[k].append(list(TMB.keys())[-1])
        write_clinical_sample(output_folder, table_dict_patient)
    except:
        print("No combined output found")
        write_clinical_sample_empty(output_folder, table_dict_patient)
        
    logger.info("Deleting temp folder")
    #clear_temp(output_folder)
    
    logger.success("Walk script completed!\n")

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':
         
    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='Argument of walk script')

    parser.add_argument('-i', '--input', required=True,
                                            help='input folder with data')
    parser.add_argument('-t', '--vcf_type', required=False,
                                            choices=['snv', 'cnv'],
                                            help='Select the vcf file to parse')
    parser.add_argument('-f', '--filter_snv', required=False,
                                            action='store_true',
                                            help='Filter out from the vcf the variants wit dot (.) in Alt column')
    parser.add_argument('-o', '--output_folder', required=True,
                                            help='Output folder')
    parser.add_argument('-k', '--oncoKB', required=False,action='store_true',help='OncoKB annotation')
    parser.add_argument('-w', '--overWrite', required=False,action='store_true',
                                                help='Overwrite output folder if it exists')
    parser.add_argument('-c', '--Cancer', required=False,
                        help='Cancer Name')
    try:
        args = parser.parse_args()
    except Exception as err:
        logger.remove()
        logfile="walk_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True,catch=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
        logger.critical(f"error: {err}", file=sys.stderr)
    
    
    input = args.input
    vcf_type = args.vcf_type
    filter_snv = args.filter_snv
    output_folder = args.output_folder
    overwrite_output=args.overWrite
    onco=args.oncoKB
    cancer=args.Cancer
    
    walk_folder(input, output_folder,onco,cancer, overwrite_output, vcf_type, filter_snv, log=False)
