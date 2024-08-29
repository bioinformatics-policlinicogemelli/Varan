version = "1.0"
# ===============================================

import os 
import sys
import argparse
from loguru import logger 
from configparser import ConfigParser
from walk import walk_folder 
from filter_clinvar import filter_main
from concatenate import concatenate_main
from ValidateFolder import validateFolderlog
from Make_meta_and_cases import meta_case_main
from Update_script import update_main 
from Delete_script import delete_main 
from ExtractSamples_script import extract_main
import shutil

config = ConfigParser()

# configFile = config.read("conf.ini")
# vaf_default = config.get('Filters', 't_VAF')
# vaf_hotspot = config.get('Filters', 't_VAF')
# vaf_novel = config.get('Filters', 't_VAF_NOVEL')

def varan(input, cancer, output_folder, oncoKB, filters, vcf_type=None, overwrite_output=False, resume=False, multiple=False, update=False, extract=False, remove=False, log=False):
    
    # if not log:
    #     logger.remove()
    #     logfile="Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
    #     logger.level("INFO", color="<green>")
    #     logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
    #     logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")
    #     logger.info("Welcome to VARAN") 

    logger.info(f"Varan args [input:{input}, output_folder:{output_folder}, filters:{filters}, cancer:{cancer}, vcf_type:{vcf_type},",
                f"overwrite_output:{overwrite_output}, resume:{resume}, multiple:{multiple}],",
                            f"update:{update}, extract:{extract}, remove:{remove}")

    if not any([update, extract, remove]) :       
            
            ###########################
            #        1.  WALK         #
            ###########################
            
            logger.info("Starting preparation study folder")
            
            if len(input) == 1 or input[1].strip == "": 
                patient_tsv = ""
            else: patient_tsv = input[1]
                
            input = input[0]
            output_folder = walk_folder(input, patient_tsv, multiple, output_folder, oncoKB, cancer, overwrite_output, resume, vcf_type, filters)


            ###########################
            #       2. FILTER         #
            ###########################
        
            logger.info("Starting filter")    

            #filter_main(output_folder, output_folder, vus,overwrite_output,log)
            if args.vcf_type =="snv" or (vcf_type==None and os.path.exists(os.path.join(input,"SNV"))):
                filter_main(input,output_folder, output_folder, oncoKB, filters, cancer, resume)
            elif os.path.exists(os.path.exists(os.path.join(input,"maf"))) and not vcf_type=="cnv":
                filter_main(input,output_folder, output_folder, oncoKB, filters, cancer, resume)

                
            ############################
            #      3. CONCATENATE      #
            ############################

            if  os.path.exists(os.path.join(output_folder,"maf")) and not vcf_type=="cnv":
                logger.info("Concatenate mutation file")
                #folders=[]
                # if vus:
                #     folders.append("NoVus")
                if oncoKB and "o" in filters:
                    folder="MAF_Onco_filtered"
                    
                elif "o" not in filters and not filters=="d":
                    folder="MAF_filtered"
                
                else:
                    folder="maf"
                
                #for folder in folders:
                input_folder=os.path.join(output_folder,folder)
                output_file=os.path.join(input_folder,"data_mutations_extended.txt")
                concatenate_main(input_folder,"maf",output_file,log)
            
                if oncoKB and "o" in filters:
                    logger.info("Extracting data_mutations_extended from OncoKB folder") 
                    os.system("cp "+os.path.join(output_folder,os.path.join("MAF_Onco_filtered","data_mutations_extended.txt"))+" "+ output_folder )
                
                elif "o" not in filters and not filters=="d":
                    os.system("cp "+os.path.join(output_folder,os.path.join("MAF_filtered","data_mutations_extended.txt"))+" "+ output_folder )
                
                else:
                    os.system("cp "+os.path.join(output_folder,os.path.join("maf","data_mutations_extended.txt"))+" "+ output_folder )
                
                
                # elif vus:
                #     logger.info("Extracting data_mutations_extended from NoVUS folder") 
                #     os.system("cp "+os.path.join(output_folder,os.path.join("NoVus","data_mutations_extended.txt"))+" "+ output_folder )
                # else:
                #     logger.info("Extracting data_mutations_extended from NoBenign folder") 
                #     os.system("cp "+os.path.join(output_folder,os.path.join("NoBenign","data_mutations_extended.txt"))+" "+ output_folder )
                
            
            ###########################################
            #      4. MAKE AND POPULATE TABLES        #
            ###########################################

            logger.info("It's time to create tables!")
            meta_case_main(cancer,output_folder)

            
            ############################
            #      5. VALIDATION       #
            ############################

            logger.info("Starting Validation Folder")
            validateFolderlog(output_folder,log)
            logger.success("The end! The study is ready to be uploaded on cBioportal")

            ##########################################################################################
            ############ DA RIVEDERE #################################################################
            # shutil.make_archive(os.path.join(output_folder,"snv_filtered"),"zip",os.path.join(output_folder,"snv_filtered"))
            # shutil.rmtree(os.path.join(output_folder,"snv_filtered"))
            # shutil.make_archive(os.path.join(output_folder,"maf"),"zip",os.path.join(output_folder,"maf"))
            # shutil.rmtree(os.path.join(output_folder,"maf"))
            ##########################################################################################
            ##########################################################################################

    ############################
    #         UPDATE           #
    ############################

    if update: 
        logger.info("Starting Update study")
        oldpath=args.Path
        new=args.NewPath
        new_name = args.Name
        output_folder=args.output_folder
        update_main(oldpath, new, output_folder, new_name, overwrite_output)

    ############################
    #         DELETE           #
    ############################

    if remove:
        logger.info("Starting Delete sample(s) from study")
        oldpath=args.Path
        removepath=args.SampleList
        new_name = args.Name
        output_folder=args.output_folder
        delete_main(oldpath, removepath, output_folder, new_name, overwrite_output)

    ############################
    #         EXTRACT          #
    ############################

    if extract:
        logger.info("Starting Extract sample(s) from study")
        oldpath = args.Path
        removepath = args.SampleList
        new_name = args.Name
        output_folder = args.output_folder
        extract_main(oldpath, removepath, output_folder, new_name, overwrite_output)

#################################################################################################################

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':
	
    logger.remove()
    logfile="Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
    logger.level("INFO", color="<green>")
    logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}", colorize=True, catch=True)
    logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
    logger.info("Welcome to VARAN")
    log=True


    parser = MyArgumentParser(add_help=True, exit_on_error=False, usage=None, description='Argument of Varan script')
    
    # WALK BLOCK
    parser.add_argument('-c', '--Cancer', required=False, 
                        help='Cancer Name')
    # parser.add_argument('-i', '--input', required=False,
    #                     help='input folder tsv with data or tsv with path of data')
    parser.add_argument('-i', '--input', nargs='+', required=False, type=str,
                        help='list with input folder/sample file tsv (required) and patient tsv')
    parser.add_argument('-t', '--vcf_type', required=False, choices=['snv', 'cnv'],
                        help='Select the vcf file to parse')
    parser.add_argument('-w', '--overWrite', required=False, action='store_true',
                        help='Overwrite output folder if it exists')
    parser.add_argument('-R', '--resume', required=False, action='store_true',
                        help='Resume an already started analysis')
    
    # ANNOTATION BLOCK
    parser.add_argument('-k', '--oncoKB', required=False, action='store_true', 
                        help='OncoKB annotation')
    parser.add_argument('-m', '--multiple', required=False, action='store_true', 
                        help='Multiple sample VCF?')
    
    # FILTER BLOCK
    parser.add_argument('-f', '--Filter', required=False, 
                        help='Select filter for SNV [d -> filter, p -> filter==PASS , b-> Benign , v-> vaf, o-> Oncokb , g -> gnomAD, q > Consequence, y-> polyphens, c -> clin_sig, n -> novel]',default="")
    
    # UPDATE BLOCK
    parser.add_argument('-u', '--Update', required=False,action='store_true',
                        help='Add this argument if you want to concatenate two studies')
    parser.add_argument('-n', '--NewPath', required=False,
                        help='Path of new study folder to add')
        
    # DELETE BLOCK
    parser.add_argument('-r', '--Remove', required=False,action='store_true',
                        help='Add this argument if you want to remove samples from a study')
    
    # EXTRACT BLOCK
    parser.add_argument('-e', '--Extract', required=False,action='store_true',
                        help='Add this argument if you want to extract samples from a study')
    
    # COMMON BLOCK 
    parser.add_argument('-o', '--output_folder', required=False, default="",
                        help='Output folder')
    parser.add_argument('-s', '--SampleList', required=False,
                        help='Path of file with list of SampleIDs')
    parser.add_argument('-p', '--Path', required=False,
                        help='Path of original study folder')
    parser.add_argument('-N', '--Name', required=False,default="",
                        help='Add this argument if you want to give a custom name to the extract study')
    
    
    try:
        args = parser.parse_args()
        
        cancer = args.Cancer
        input = args.input
        filters=args.Filter
        output_folder = args.output_folder
        vcf_type=args.vcf_type
        overwrite_output=args.overWrite
        resume=args.resume
        oncoKB=args.oncoKB
        multiple=args.multiple
        
        update=args.Update
        extract=args.Extract
        remove=args.Remove
                    
        if sum([args.Update, args.Extract, args.Remove]) > 1:
            logger.critical("Please select only one option between Update, Extract and Remove")
            sys.exit()

        if not any([args.Update, args.Extract, args.Remove]) and (args.input==None or args.input[0].strip()==""):
            logger.critical("Error Argument: Valid Input is required")
            sys.exit()
        
        if not any([args.Update, args.Extract, args.Remove]) and args.output_folder=="":
            logger.critical("Error Argument: Output is required")
            sys.exit()
        
        if not any([args.Update, args.Extract, args.Remove]) and args.Cancer==None:
            logger.critical("Error Argument: Cancer name is required")
            sys.exit()

        if args.Update and (args.Path==None or args.NewPath==None):
            logger.critical("To update a study, you need to specify both original and new folder paths")
            sys.exit()

        if (any([args.Remove,args.Extract]) and args.Path==None) or (any([args.Remove,args.Extract]) and args.SampleList==None):
            logger.critical("To remove/extract samples from a study, you need to specify both original folder path and samples' list")
            sys.exit()
            
        if (args.output_folder=="" and args.Name!=""):
            logger.critical("To use -N option it's required to set also -o")
            sys.exit()
        
        if resume:
            if overwrite_output:
                logger.critical("Both resume and overwrite options are selected. Please select only one!")
                sys.exit()
         
        varan(input, cancer, output_folder, oncoKB, filters, vcf_type, overwrite_output, resume, multiple, update, extract, remove, log)
    
    except Exception as err:
        logger.critical(f"error: {err}", file=sys.stderr)