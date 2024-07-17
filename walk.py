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
import numpy as np
from filter_clinvar import filter_OncoKB
import pandas as pd

config = ConfigParser()

configFile = config.read("conf.ini")

INPUT_PATH=config.get('Paths', 'INPUT_PATH')
SNV=config.get('Paths', 'SNV')
CNV=config.get('Paths', 'CNV')
COMBOUT=config.get('Paths', 'COMBOUT')
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
    




def cnv_type_from_folder(input,cnv_vcf_files,output_folder,oncokb,cancer, multiple):
    
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
                if multiple:
                    sID_path[sampleID] = os.path.join(os.path.join(input,"CNV", "single_sample_vcf"),cnv_vcf)
                else:
                    sID_path[sampleID] = os.path.join(os.path.join(input,"CNV"),cnv_vcf)
                
                vcf2tab_cnv.vcf_to_table(sID_path[sampleID], os.path.join(output_folder,'data_cna_hg19.seg'), sampleID, MODE)
                vcf2tab_cnv.vcf_to_table_fc(sID_path[sampleID], os.path.join(output_folder,'data_cna_hg19.seg.fc.txt'), sampleID, MODE)
               
        except Exception:
            log_noparsed = open('noParsed_cnv.log', 'a')
            log_noparsed.write('[WARNING] '+case_folder+'\n')
            log_noparsed.close()
        
        c = c + 1
    logger.info("Writing data_cna_hg19.seg succefully completed!")
    logger.info("Writing data_cna_hg19.seg.fc.txt succefully completed!")
    
    ############################
    ### MANAGE DISCRETE TABLE ##
    ############################
    
    df_table = pd.read_csv(os.path.join(output_folder,'data_cna_hg19.seg.fc.txt'),sep="\t",header=0)
    result = tabella_to_dict(df_table)
        
    df_table.rename(columns={"discrete":"Copy_Number_Alteration","ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"},inplace=True)
    df_table_filt=df_table[df_table["Copy_Number_Alteration"].isin([-2,2])]
        
    if oncokb:
        
        ################################# OPZIONE 1 #########################################################
        # if not os.path.isfile(input):
        #     tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
        #     input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
        
        # else:
        #     input_file=pd.read_csv(input,sep="\t")
        
        # for index, row in df_table.iterrows():
        #     tc= int(input_file[input_file["SampleID"]+".cnv.bam"==row["ID"]]["TC"])
        #     df_table.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc
            
        
        # df_table_filtered=df_table#[(df_table["discrete"]>=3)|(df_table["discrete"]<=0.8)] 
        
        #df_table_filtered["Copy_Number_Alteration"]=0
        # df_table_filtered.loc[(df_table_filtered["discrete"]>3)&(df_table_filtered["discrete"]<5), "Copy_Number_Alteration"]=1
        # df_table_filtered.loc[df_table_filtered["discrete"]>5, "Copy_Number_Alteration"]=2
        # df_table_filtered.loc[(df_table_filtered["discrete"]>0)&(df_table_filtered["discrete"]<0.8), "Copy_Number_Alteration"]=-1
        # df_table_filtered.loc[df_table_filtered["discrete"]<=0, "Copy_Number_Alteration"]=-2
        
       
        #df_table_filtered["Copy_Number_Alteration"]=df_table_filtered["discrete"]
        #df_table_filtered.to_csv(os.path.join(output_folder,"data_cna_hg19.seg.filtered.txt"),sep="\t",index=False)
        #temp_cna=reshape_cna(input,os.path.join(output_folder,"data_cna_hg19.seg.filtered.txt"),cancer,output_folder)
        #annotate_cna(temp_cna,output_folder)
   
        ################################# OPZIONE 2 #########################################################
   
        # rename discrete 
        
        #df_table.rename(columns={"discrete":"Copy_Number_Alteration","ID":"Tumor_Sample_Barcode","gene":"Hugo_Symbol"},inplace=True)
        #df_table_filt=df_table[df_table["Copy_Number_Alteration"].isin([-2,2])]
        if not os.path.isfile(input):
            tsv_file=[file for file in os.listdir(input) if file.endswith(".tsv")][0]
            input_file=pd.read_csv(os.path.join(input,tsv_file),sep="\t")
        else:
            input_file=pd.read_csv(input,sep="\t")
        
        #TODO inserire filtro per TC (mettere specifiche su valori di TC da considerare)
        
        #check for nan TC values
        
        try:
            input_file["TC"]
        except:
            input_file["TC"]=np.nan
        
        if len(input_file[input_file["TC"].isna()])>0:
            nan_sbj=input_file[input_file["TC"].isna()]
            nan_sbj=list(nan_sbj["SampleID"])
            logger.warning(f"Some subject have Nan TC in tsv input file: {nan_sbj}!")
        
        if not "ONCOTREE_CODE" in input_file.columns:
            input_file["ONCOTREE_CODE"]=cancer

        input_file["Tumor_Sample_Barcode"]=input_file["SampleID"]+".cnv.bam"
        
        annotate= pd.merge( df_table_filt[["Tumor_Sample_Barcode","Hugo_Symbol","seg.mean","Copy_Number_Alteration"]] ,input_file[["Tumor_Sample_Barcode","ONCOTREE_CODE","TC"]],on="Tumor_Sample_Barcode")
        temppath=os.path.join(output_folder,"temp_cna_toannotate.txt")
        annotate.to_csv(temppath,index=False,sep="\t")

        
        out=temppath.replace("toannotate.txt","annotated.txt")
        os.system(f"python3 ./oncokb-annotator/CnaAnnotator.py -i {temppath}\
                        -o {out} -f individual -b {config.get('OncoKB', 'ONCOKB')}")
              
        cna=pd.read_csv(out,sep="\t",dtype={"Copy_Number_Alteration":int})
        cna=cna[cna["ONCOGENIC"].isin(["Oncogenic","Likely Oncogenic"])]
        cna["ESCAT"]="Unmatched"
        df_table["ESCAT"]="Unmatched"
        
        for index, row in cna.iterrows():
            
            logger.info("Analyzing cna sample " + row["Tumor_Sample_Barcode"])            
            try:
                tc= int(input_file[input_file["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]]["TC"])
            except ValueError:
                continue
            except Exception:
                raise(Exception(f"Something went wrong while reading TC!"))

            #decomm se cutoff su copy number discreto
            #cna.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc 
            
            purity = tc/100
            ploidy = 2
            copy_nums = np.arange(6)
            c=2**(np.log2((1 - purity) + purity * (copy_nums + .5) / ploidy ))
            
            # ESCAT classification
            oncocode=input_file[input_file["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]]["ONCOTREE_CODE"].values
            gene=df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]

            if oncocode=="NSCLC":  
                if (row["Hugo_Symbol"]=="MET" and row["seg.mean"]>1):
                    row["ESCAT"]="IIB"
                    df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
                elif (row["Hugo_Symbol"]=="ERBB2" and row["seg.mean"]>1):
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
                if (row["Hugo_Symbol"]=="BRCA1" and row["seg.mean"]<1):
                    row["ESCAT"]="IA"
                    df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"
                elif (row["Hugo_Symbol"]=="BRCA2" and row["seg.mean"]<1):
                    row["ESCAT"]="IA"
                    df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IA"
                elif (row["Hugo_Symbol"]=="PTEN" and row["seg.mean"]<1):
                    row["ESCAT"]="IIA"
                    df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIA"
                elif (row["Hugo_Symbol"]=="ATM" and row["seg.mean"]<1):
                    row["ESCAT"]="IIB"
                    df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
                elif (row["Hugo_Symbol"]=="PALB2" and row["seg.mean"]<1):
                    row["ESCAT"]="IIB"
                    df_table[(df_table["Tumor_Sample_Barcode"]==row["Tumor_Sample_Barcode"]) & (df_table["Hugo_Symbol"]==row["Hugo_Symbol"])]["ESCAT"]="IIB"
                    
        #nuovi filtri (filtro cnvkit - conservativo)
        
        if not cna.empty:
            cna["Copy_Number_Alteration"]=0
            cna.loc[(cna["seg.mean"]<c[0]), "Copy_Number_Alteration"]=-2
            cna.loc[(cna["seg.mean"]>=c[0])&(cna["seg.mean"]<c[1]), "Copy_Number_Alteration"]=-1
            cna.loc[(cna["seg.mean"]>=c[1])&(cna["seg.mean"]<c[3]), "Copy_Number_Alteration"]=0
            cna.loc[(cna["seg.mean"]>=c[3])&(cna["seg.mean"]<c[5]), "Copy_Number_Alteration"]=1
            cna.loc[cna["seg.mean"]>=c[5], "Copy_Number_Alteration"]=2
        
        else:
            logger.warning("After filtering cna dataframe is empty!")
               
        #vecchi filtri (decomm se filtro copy number discreto)
        
        #decomm se cutoff su copy number discreto
        #cna.at[index,"discrete"]= ((200*float(row["seg.mean"]))-2*(100-tc))/tc 
            
        # cna["Copy_Number_Alteration"]=0
        # cna.loc[(cna["discrete"]>3)&(cna["discrete"]<5), "Copy_Number_Alteration"]=1
        # cna.loc[cna["discrete"]>5, "Copy_Number_Alteration"]=2
        # cna.loc[(cna["discrete"]>0.25)&(cna["discrete"]<1), "Copy_Number_Alteration"]=-1
        # cna.loc[cna["discrete"]<=0.25, "Copy_Number_Alteration"]=-2

        df_table.to_csv(os.path.join(output_folder,"data_cna_hg19_escat.seg.fc.txt"),index=True,sep="\t")
        
        cna.to_csv(os.path.join(output_folder,"CNA_ANNOTATI_NDISCRETO.txt"),index=True,sep="\t")
        
        cna["Tumor_Sample_Barcode"]=cna["Tumor_Sample_Barcode"].str.replace(".cnv.bam","")
        data_cna=cna.pivot_table(index="Hugo_Symbol",columns="Tumor_Sample_Barcode",values="Copy_Number_Alteration",fill_value=0)
        data_cna.to_csv(os.path.join(output_folder,"data_cna.txt"),index=True,sep="\t")
        
    else:
        
        df_table_filt["Tumor_Sample_Barcode"]=df_table_filt["Tumor_Sample_Barcode"].str.replace(".cnv.bam","")
        
        data_cna=df_table_filt.pivot_table(index="Hugo_Symbol",columns="Tumor_Sample_Barcode",values="Copy_Number_Alteration",fill_value=0)
        data_cna.to_csv(os.path.join(output_folder,"data_cna.txt"),index=True,sep="\t")
    
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
        _, vcf_file = os.path.split(v)
        
        out_filt=os.path.join(output_folder,OUTPUT_FILTERED) #TEST
        vcf_filtered = os.path.join(out_filt, vcf_file.replace('.vcf','')+'.FILTERED.vcf')
        logger.info(f'[FILTERING] {v}')
        vcf_filter.main(v, vcf_filtered)
        logger.info(f'filtered file {vcf_filtered} created!')
        sID_path_filtered[k] = vcf_filtered
    return sID_path_filtered

def vcf2maf_constructor(k, v, temporary,output_folder):
    
    cmd="grep $'\t'"+k.split(".")[0]+" "+v
    try:
        tum_id = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
        tum_id=[i for i in tum_id.split("\t") if k.split(".")[0] in i][0]
    except Exception:
        tum_id=""
    
    cl = ['perl']
    cl.append(VCF2MAF)
    cl.append('--input-vcf')
    cl.append(v)
    _, file_vcf = os.path.split(v)
    out_file = os.path.join(output_folder,os.path.join(OUTPUT_MAF, file_vcf+'.maf'))
    if not CLINV =="":
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
    
def create_folder(output_folder,overwrite_output, resume):
    if overwrite_output:
        if os.path.exists(output_folder):
            logger.warning(f"It seems that the folder '{output_folder}' already exists. Start removing process...")
            shutil.rmtree(output_folder)
    
    if not os.path.exists(output_folder):
        logger.info(f"Creating the output folder '{output_folder}' in {os.getcwd()}...")
        os.mkdir(output_folder)
        maf_path = os.path.join(output_folder, 'maf')
        os.mkdir(maf_path)
        #TODO capire se serve tenere cartella snv_filtered
        filtered_path = os.path.join(output_folder, 'snv_filtered')
        os.mkdir(filtered_path)
        logger.info(f"The folder '{output_folder}' was correctly created!")
    
    elif resume:
        maf_path = os.path.join(output_folder, 'maf')
        filtered_path = os.path.join(output_folder, 'snv_filtered')

    else:
        logger.critical(f"The folder '{output_folder}' already exists. To overwrite an existing folder add the -w option!")
        raise(Exception('Error in create_folder script: exiting from walk script!'))

def get_table_from_folder(tsvpath):
    table_dict = dict()
    file = pd.read_csv(tsvpath, sep="\t", index_col=False, dtype=str)
    for _, row in file.iterrows():
        sampleID = str(row["SampleID"])
        if ".bam" in sampleID:
           sampleID = sampleID.replace(".bam","")
        if sampleID not in table_dict.keys():
            table_dict[sampleID] = [str(row["PatientID"])]  
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

def write_clinical_sample(input, output_folder, table_dict):

    logger.info("Writing data_clinical_sample.txt file...")

    tsv_file = [file for file in os.listdir(input) if file.endswith(".tsv")][0]
    
    conf_header_short = config.get('ClinicalSample', 'HEADER_SAMPLE_SHORT')
    conf_header_long = config.get('ClinicalSample', 'HEADER_SAMPLE_LONG')
    conf_header_type = config.get('ClinicalSample', 'HEADER_SAMPLE_TYPE')

    input_file_path = os.path.join(input, tsv_file)
    with open(input_file_path, 'r') as input_file:
        file_header_string = input_file.readline()

    file_header = file_header_string.split("\t")
    file_header = list(map(lambda x: x.strip(), file_header))
    snv_index = file_header.index('snv_path')
    cnv_index = file_header.index('cnv_path')
    comb_index = file_header.index('comb_path')

    file_header = list(map(lambda x: x.upper(), file_header))
    header_string = "\t".join(file_header) + "\n"

    data_clin_samp = os.path.join(output_folder, 'data_clinical_sample.txt')
    
    if not conf_header_short:
        sample_header_short = file_header
    else:
        sample_header_short = conf_header_short
    header_list_short = "\t".join(sample_header_short)
    header_list_short = header_list_short + "\n" if not header_list_short.endswith("\n") else header_list_short

    if not conf_header_long:
        sample_header_long = sample_header_short
    else:
        sample_header_long = [field.upper() for field in conf_header_long]
    header_list_long = "\t".join(sample_header_long)
    header_list_long = header_list_long + "\n" if not header_list_long.endswith("\n") else header_list_long

    if not conf_header_type:
        sample_header_type = ("STRING\t" * (len(file_header) - 1) + "STRING\n")
    else: 
        sample_header_type = "\t".join(sample_header_type)
        header_type_list = header_type_list + "\n" if not header_type_list.endswith("\n") else header_type_list
  
    with open(data_clin_samp, "w") as cil_sample:
        cil_sample.write("#" + header_list_short)
        cil_sample.write("#" + header_list_long)
        cil_sample.write("#" + sample_header_type)
        cil_sample.write("#" + "1\t" * (len(file_header) - 1) + "1\n") 
        cil_sample.write(header_string)

    with open(input_file_path, 'r') as input_file:
        with open(data_clin_samp, 'a') as cil_sample:
            cil_sample.writelines(input_file.readlines()[1:])

    import pdb; pdb.set_trace()

    # Leggi i file come dataframe
    cil_sample_df = pd.read_csv(data_clin_samp, sep="\t", header=None, dtype=str)
    table_df = pd.DataFrame.from_dict(table_dict).transpose().reset_index()
    table_df.columns = ["SAMPLEID", "PATIENTID", "MSI", "TMB", "MSI_THR", "TMB_THR"]

    # Rimuovi le colonne indesiderate
    columns_to_drop = [snv_index, cnv_index, comb_index]
    cil_sample_df.drop(columns=columns_to_drop, axis=1, inplace=True)

    # Dividi il dataframe in intestazione e corpo
    cil_sample_header = cil_sample_df.iloc[:4, :]
    cil_sample_body = cil_sample_df.iloc[4:, :]

    # Imposta i nomi delle colonne per il corpo del dataframe
    cil_sample_body.columns = cil_sample_body.iloc[0]
    cil_sample_body = cil_sample_body[1:]  # Rimuovi la riga usata come intestazione

    # Unisci il corpo del dataframe con table_df
    cil_sample_body = pd.merge(cil_sample_body, table_df, on=["PATIENTID", "SAMPLEID"])

    # Pulisci il dataframe di intestazione
    cil_sample_header.columns = ["SAMPLEID", "PATIENTID"]
    extra_values = {
        "MSI": ["MSI", "MSI", "NUMBER", 1],
        "TMB": ["TMB", "TMB", "NUMBER", 1],
        "MSI_THR": ["MSI_THR", "MSI_THR", "STRING", 1],
        "TMB_THR": ["TMB_THR", "TMB_THR", "STRING", 1],
    }
    for col, values in extra_values.items():
        cil_sample_header[col] = values

    # Aggiungi la riga dei nomi delle colonne al dataframe di intestazione
    column_names = cil_sample_body.columns.tolist()
    cil_sample_header = cil_sample_header.append(pd.Series(column_names, index=cil_sample_header.columns), ignore_index=True)

    # Unisci il dataframe di intestazione e corpo
    result = pd.concat([cil_sample_header, cil_sample_body], ignore_index=True)

    # Scrivi il dataframe risultante in un file di testo
    result.to_csv(data_clin_samp, sep="\t", index=False, header=False)


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
        cil_sample.write(v[0]+'\t'+k+'\n')
    cil_sample.close()

def check_cna_vcf(file, inputFolderCNV, multivcf):
    vcf=pd.read_csv(os.path.join(inputFolderCNV,file),comment="#",sep="\t",names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"])
        
    nsample=os.popen(f'bcftools query -l {os.path.join(inputFolderCNV,file)}').read().split("\n")
    nsample=[sample for sample in nsample if not sample ==""]
    if len(nsample)>1 and not multivcf:
        logger.critical("VCF contains multiple samples")
        exit(1)

    if vcf.loc[0]["FORMAT"]=="FC":
        return True
    else:
        return False
    
def check_snv_vcf(file,inputFolderSNV, multivcf):
    vcf=pd.read_csv(os.path.join(inputFolderSNV,file),comment="#",sep="\t",names=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample"])
    
    nsample=os.popen(f'bcftools query -l {os.path.join(inputFolderSNV,file)}').read().split("\n")
    nsample=[sample for sample in nsample if not sample ==""]
    if len(nsample)>1 and not multivcf:
        logger.critical("VCF contains multiple samples")
        exit(1)
    
    if vcf.loc[0]["FORMAT"].startswith("GT"):
        return True
    else:
        return False
    
import os

def extract_multiple_cnv(multiple_vcf,input_dir):
    single_sample_vcf_dir = os.path.join(input_dir, "single_sample_vcf")  
    if not os.path.exists(os.path.join(input_dir,"single_sample_vcf")):
        os.mkdir(os.path.join(input_dir,"single_sample_vcf"))

    cmd= "vcf-query -l " + multiple_vcf + " > " + os.path.join(input_dir, "sample_id.txt")
    os.system(cmd)

    with open (os.path.join(input_dir, "sample_id.txt"), "r") as file:
        lines = file.readlines()
        for sample in lines:
            sample = sample.strip()
            single_sample_vcf = os.path.join(single_sample_vcf_dir, f"{sample}.vcf")
            cmd = f"vcftools --vcf {multiple_vcf} --indv {sample} --recode --recode-INFO-all --stdout > {single_sample_vcf}"
            os.system(cmd)

    
def extract_multiple_snv(multiple_vcf,input_dir):    
    if not os.path.exists(os.path.join(input_dir,"single_sample_vcf")):
        os.mkdir(os.path.join(input_dir,"single_sample_vcf"))

    cmd= "vcf-query -l " + multiple_vcf + " > " + os.path.join(input_dir, "sample_id.txt")
    os.system(cmd)

    with open (os.path.join(input_dir, "sample_id.txt"), "r") as file:
        lines = file.readlines()
        for sample in lines:
            print(sample.strip())
            cmd="vcf-subset --exclude-ref -c " + sample.strip() + " " + multiple_vcf + " > " + os.path.join(os.path.join(input_dir, "single_sample_vcf"), sample.strip() + ".vcf")
            os.system(cmd)

def check_field_tsv(row, name):
    try:
        field=str(row[name])
    except KeyError as e: 
        logger.critical(f"KeyError: {e} not found! Check if column name is correctly spelled or if there are tabs/spaces before or after the coloumn key: \n{row.index}. \nThis error may also occur if the table columns have not been separated by tabs!")
        raise(KeyError("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))
    return field

def get_combinedVariantOutput_from_folder(inputFolder, tsvpath):
    
    combined_dict = dict()
    
    try:
        file = pd.read_csv(tsvpath,sep="\t",dtype=str)
    except Exception as e:
        logger.critical(f"Something went wrong while reading {tsvpath}!")
        raise(Exception("Error in get_combinedVariantOutput_from_folder script: exiting from walk script!"))
    
    for _,row in file.iterrows():
        
        #patientID = check_field_tsv(row, "PatientID")
        sampleID = check_field_tsv(row, "SampleID")
        combined_path = check_field_tsv(row, "comb_path")
        # combined_file = patientID+COMBOUT #da verificare
        # combined_path = os.path.join(inputFolder,"CombinedOutput",combined_file)
        
        if os.path.exists(combined_path):
            pass 
        else:
            logger.warning(f"{combined_path} not exists")
        combined_dict[sampleID] = combined_path
    return combined_dict

def check_input_file(output_folder, file, copy_to):
    if os.path.exists(file):
        os.system("cp " + file + " " + os.path.join(output_folder, "temp", copy_to))
    elif not file:
        logger.warning(f"No final_path set in conf.ini!")
    else:
        logger.warning(f"{file} not found")
        
def check_folders(output_folder, snv_path, cnv_path, combout):
    check_input_file(output_folder, snv_path, "SNV")
    check_input_file(output_folder, cnv_path, "CNV")
    check_input_file(output_folder, combout, "CombinedOutput")
        
def transform_input(tsv,output_folder, multiple):
    os.makedirs(os.path.join(output_folder,"temp"), exist_ok=True)
    os.makedirs(os.path.join(output_folder,"temp","SNV"), exist_ok=True)
    os.makedirs(os.path.join(output_folder,"temp","CNV"), exist_ok=True)
    os.makedirs(os.path.join(output_folder,"temp","CombinedOutput"), exist_ok=True)

    os.system("cp "+ tsv +" "+os.path.join(output_folder,"temp"))

    if multiple:
        snv_path=config.get('Multiple', 'SNV')
        cnv_path=config.get('Multiple', 'CNV')
        combout=config.get('Multiple', 'COMBOUT')
        
        check_folders(output_folder, snv_path, cnv_path, combout)
    
    else:   
        tsv_file=pd.read_csv(tsv,sep="\t",dtype="string", keep_default_na=False)
            
        for _,row in tsv_file.iterrows():

            res_folder=INPUT_PATH
            #res_folder="/data/data_storage/novaseq_results"
            snv_path=row["snv_path"]
            cnv_path=row["cnv_path"]
            combout=row["comb_path"]
            
            # snv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+SNV)     #"_MergedSmallVariants.genome.vcf")
            # #snv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+"_MergedSmallVariants.genome.vcf")
            # cnv_path=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["SampleID"],row["SampleID"]+CNV) # "_CopyNumberVariants.vcf")
            # combout=os.path.join(res_folder,row["RunID"],"Results",row["PatientID"],row["PatientID"]+COMBOUT)

            check_folders(output_folder, snv_path, cnv_path, combout)
                
    return os.path.join(output_folder,"temp")

def walk_folder(input, multiple, output_folder,oncokb,cancer, overwrite_output=False, resume=False, vcf_type=None ,filters="", log=False):
    
    if not log:
        logger.remove()
        logfile="Walk_folder_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")
    
    logger.info("Starting walk_folder script:")
    logger.info(f"walk_folder args [input:{input}, output_folder:{output_folder}, Overwrite:{overwrite_output}, resume:{resume}, vcf_type:{vcf_type}, filters:{filters}, multiple:{multiple}]")
 
    config.read('conf.ini')

    assert os.path.exists(input), \
        f"No valid file/folder {input} found. Check your input path"
    
    ###############################
    ###       OUTPUT FOLDER     ###
    ###############################

    if not resume or not os.path.exists(os.path.join(output_folder,"temp")):
        create_folder(output_folder, overwrite_output, resume)
        if os.path.isfile(input):
            input=transform_input(input, output_folder, multiple)
       
    input=os.path.join(output_folder,"temp") 
    inputFolderSNV=os.path.abspath(os.path.join(input,"SNV"))
    inputFolderCNV=os.path.abspath(os.path.join(input,"CNV"))
    inputFolderCombOut=os.path.abspath(os.path.join(input,"CombinedOutput"))
    
    assert (len(os.listdir(inputFolderCNV))>0 or len(os.listdir(inputFolderCNV))>0 or len(os.listdir(inputFolderCombOut))>0), \
        "No valid input file was found for neither SNV, CNV or CombinedOutput! Check your input file/folder and input options."
       
    if os.path.exists(inputFolderCNV) and not type in ["snv", "fus", "tab"]:
        if multiple:
            multivcf = os.listdir(inputFolderCNV)[0]
            extract_multiple_cnv(os.path.join(inputFolderCNV,multivcf),inputFolderCNV)
            inputFolderCNV= os.path.join(inputFolderCNV,"single_sample_vcf")
        logger.info("Check CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(inputFolderCNV)
        logger.info("Everything ok!")

    if os.path.exists(inputFolderSNV) and not type in ["cnv", "fus", "tab"]:
        os.makedirs(TMP, exist_ok=True)
        if multiple:
            multivcf = os.listdir(inputFolderSNV)[0]
            extract_multiple_snv(os.path.join(inputFolderSNV,multivcf),inputFolderSNV)
            inputFolderSNV= os.path.join(inputFolderSNV,"single_sample_vcf")
        logger.info("Check SNV files...")
    case_folder_arr = get_snv_from_folder(inputFolderSNV)
    logger.info("Everything ok!")

    ###############################
    ###       SNV AND CNV       ###
    ###############################
    if os.path.exists(inputFolderCNV) and not type in ["snv","fus","tab"]:
        sID_path_cnv = cnv_type_from_folder(input,case_folder_arr_cnv,output_folder,oncokb,cancer, multiple)
    
    if os.path.exists(inputFolderSNV) and not type in ["cnv","fus","tab"]:
        sID_path_snv = snv_type_from_folder(inputFolderSNV,case_folder_arr)
        
        logger.info("Check maf folder...")
        maf_path=os.path.join(output_folder,"maf")
        if os.path.isdir(maf_path) and len([i for i in os.listdir(maf_path) if i.endswith('.maf')])>0:
            logger.info("a non empty maf folder already exists!")
        
        if not resume:
            logger.info("Starting vcf2maf conversion...")
            if "d" in filters:
                logger.info("filtering out vcfs with dots in ALT column")
            #     temporary = create_random_name_folder()
                sID_path_snv = vcf_filtering(sID_path_snv,output_folder)
            #     for k, v in sID_path_filtered.items():
            #         cl = vcf2maf_constructor(k, v, temporary,output_folder)
            #         run_vcf2maf(cl)
            # else:
            temporary = create_random_name_folder()
            for k, v in sID_path_snv.items():
                cl = vcf2maf_constructor(k, v, temporary,output_folder)
                run_vcf2maf(cl)
    
    logger.info("Clearing scratch folder...")
    clear_scratch()
    
   
    ###############################
    ###       GET FUSION        ###
    ###############################
    
    if not COMBOUT=="":
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

        # try:
        tsvpath=os.path.join(input, tsvfiles)    
        combined_dict = get_combinedVariantOutput_from_folder(input, tsvpath)
        
        fusion_table_file = os.path.join(output_folder, 'data_sv.txt')
            
    #     for k, v in combined_dict.items():
    #         fusions=[]
    #         logger.info(f"Reading Fusion info in CombinedOutput file {v}...")
    #         try:
    #             fusions = tsv.get_fusions(v)
    #         except Exception as e:
    #             logger.error(f"Something went wrong while reading Fusion section of file {v}")
    # except FileNotFoundError:
    #     logger.critical(f"No input file '{input}' was found: try check your path")
    #     raise(FileNotFoundError("Exiting from walk script!"))

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
                #TODO far leggere l'header da template fusioni
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
                        fusion_table.write(k+'\tSOMATIC\tFUSION\t'+\
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
            
            if "o" in filters:
                fus_file=pd.read_csv(fusion_table_file_out,sep="\t")
                fus_file=filter_OncoKB(fus_file)
                fus_file.to_csv(fusion_table_file_out,index=False,sep="\t")
            os.system(f"mv {fusion_table_file_out}  {fusion_table_file}") 
            

            
            
        #     fileonco=pd.read_csv(fusion_table_file_out,sep="\t")
            
        # # fileonco=fileonco[fileonco["ONCOGENIC"].isin(["Oncogenic","Likely Oncogenic"])]
        #     fileonco.to_csv(fusion_table_file,sep="\t",index=False)           
            
    
    # except FileNotFoundError:
    #     logger.critical(f"No input file '{input}' was found: try check your path")
    #     raise(FileNotFoundError("Exiting from walk script!"))
    
    ##############################
    ##       MAKES TABLE       ###
    ##############################
    
    tsv_file=[file for file in os.listdir(input) if file.endswith("tsv")][0]
    tsvpath=os.path.join(input,tsv_file)    

    table_dict_patient = get_table_from_folder(tsvpath)

    logger.info("Writing clinical files...")
    write_clinical_patient(output_folder, table_dict_patient)
    #fileinputclinical=pd.read_csv(tsvpath,sep="\t",index_col=False, dtype=str)
    
    if len(os.listdir(os.path.join(output_folder, "temp", "CombinedOutput"))) > 0:
        combined_dict = get_combinedVariantOutput_from_folder(input, tsvpath)
        MSI_THR = config.get('MSI', 'THRESHOLD')
        TMB = ast.literal_eval(config.get('TMB', 'THRESHOLD'))
        
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
        
        # run=str(fileinputclinical[fileinputclinical["SampleID"]==k]["RunID"].values[0])
        # tc=str(fileinputclinical[fileinputclinical["SampleID"]==k]["TC"].values[0])
        # oncotree=str(fileinputclinical[fileinputclinical["SampleID"]==k]["ONCOTREE_CODE"].values[0])
        
        # table_dict_patient[k].append(run) 
        # table_dict_patient[k].append(oncotree)
        # table_dict_patient[k].append(tc)

    write_clinical_sample(input, output_folder, table_dict_patient)

    #DECOMMENTARE UNA VOLTA FINITI TEST
    #shutil.rmtree(os.path.join(output_folder,"temp"))
    # except:
    #     print("No combined output found")
    #     write_clinical_sample_empty(output_folder, table_dict_patient)
 
    logger.info("Deleting temp folder")
    #clear_temp(output_folder)
    
    logger.success("Walk script completed!\n")

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':
         
    parser = MyArgumentParser(add_help=True, exit_on_error=False, usage=None, description='Argument of walk script')

    parser.add_argument('-i', '--input', required=True,
                                            help='input folder with data')
    parser.add_argument('-t', '--vcf_type', required=False,
                                            choices=['snv', 'cnv'],
                                            help='Select the vcf file to parse')
    parser.add_argument('-f', '--filters', required=False,
                                            action='store_true',
                                            help='Select filter for SNV [p -> filter==PASS , b-> Benign , v-> vaf, o-> Oncokb , g -> gnomAD, q > Consequence, y-> polyphens -> clin_sig, n -> novel]',default="")
    parser.add_argument('-o', '--output_folder', required=True,
                                            help='Output folder')
    parser.add_argument('-k', '--oncoKB', required=False,action='store_true',help='OncoKB annotation')
    parser.add_argument('-w', '--overWrite', required=False,action='store_true',
                                                help='Overwrite output folder if it exists')
    parser.add_argument('-c', '--Cancer', required=False,
                        help='Cancer Name')
    parser.add_argument('-m', '--multiple', required=False, action='store_true', help='Multiple sample VCF?')

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
    filters = args.filters
    output_folder = args.output_folder
    overwrite_output=args.overWrite
    onco=args.oncoKB
    cancer=args.Cancer
    multiple=args.multiple
    
    walk_folder(input, multiple, output_folder, onco, cancer, overwrite_output, vcf_type, filters, log=False)

           