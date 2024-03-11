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

def varan(input, cancer, output_folder,oncoKB, filter_snv=False, vcf_type=None, overwrite_output=False, vus=False, update=False, extract=False, remove=False, log=False):
    
    if not log:
        logger.remove()
        logfile="Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")
        logger.info("Welcome to VARAN") 

    logger.info(f"Varan args [input:{input}, output_folder:{output_folder}, filter_snv:{filter_snv}, cancer:{cancer}, \
                            vcf_type:{vcf_type}, overwrite_output:{overwrite_output}, vus:{vus}], \
                            update:{update}, extract:{extract}, remove:{remove}")

    if not any([update ,extract , remove]) :       
            
            ###########################
            #        1.  WALK         #
            ###########################
            
            logger.info("Starting preparation study folder")
            walk_folder(input, output_folder,oncoKB,args.Cancer,overwrite_output, vcf_type,filter_snv,log)


            ###########################
            #       2. FILTER         #
            ###########################
        
            logger.info("Starting filter")    
            #filter_main(output_folder, output_folder, vus,overwrite_output,log)
            if args.vcf_type =="snv" or (args.vcf_type==None and os.path.exists(os.path.join(args.input,"SNV"))):
                filter_main(input,output_folder, output_folder, args.vus,oncoKB,args.Cancer)
            elif os.path.exists(os.path.exists(os.path.join(args.input,"maf"))) and not args.vcf_type=="cnv":
                filter_main(input,output_folder, output_folder, args.vus,oncoKB,args.Cancer,False)
                                
                
            ############################
            #      3. CONCATENATE      #
            ############################
            if  os.path.exists(os.path.join(output_folder,"maf"))and not args.vcf_type=="cnv":
                logger.info("Concatenate mutation file")
                folders=[]
                # if vus:
                #     folders.append("NoVus")
                if oncoKB:
                    folders.append("MAF_Onco_filtered")
                
                for folder in folders:
                    input_folder=os.path.join(output_folder,folder)
                    output_file=os.path.join(input_folder,"data_mutations_extended.txt")
                    concatenate_main(input_folder,"maf",output_file,log)
            
                if oncoKB:
                    logger.info("Extracting data_mutations_extended from OncoKB folder") 
                    os.system("cp "+os.path.join(output_folder,os.path.join("MAF_Onco_filtered","data_mutations_extended.txt"))+" "+ output_folder )
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
            meta_case_main(cancer,vus,output_folder,log)

            
            ############################
            #      5. VALIDATION       #
            ############################

            logger.info("Starting Validation Folder")
            validateFolderlog(output_folder,log)
            logger.success("The end! The study is ready to be uploaded on cBioportal")
            shutil.make_archive(os.path.join(output_folder,"snv_filtered"),"zip",os.path.join(output_folder,"snv_filtered"))
            shutil.rmtree(os.path.join(output_folder,"snv_filtered"))
            #
            # if os.path.exists(os.path.join(output_folder,"NoVus")):
            #     shutil.make_archive(os.path.join(output_folder,"NoVus"),"zip",os.path.join(output_folder,"NoVus"))
            #     shutil.rmtree(os.path.join(output_folder,"NoVus"))
            # #
            # #
            # if os.path.exists(os.path.join(output_folder,"NoBenign")):
            #     shutil.make_archive(os.path.join(output_folder,"NoBenign"),"zip",os.path.join(output_folder,"NoBenign"))
            #     shutil.rmtree(os.path.join(output_folder,"NoBenign"))
            #
            shutil.make_archive(os.path.join(output_folder,"maf"),"zip",os.path.join(output_folder,"maf"))
            shutil.rmtree(os.path.join(output_folder,"maf"))


    ############################
    #         UPDATE           #
    ############################

    if update: 
        logger.info("Starting Update study")
        oldpath=args.Path
        new=args.NewPath
        output_folder=args.output_folder
        update_main(oldpath,new,output_folder,log)

    ############################
    #         DELETE           #
    ############################

    if remove:
        logger.info("Starting Delete samples from study")
        oldpath=args.Path
        removepath=args.SampleList
        output_folder=args.output_folder
        delete_main(oldpath,removepath,output_folder,log)

    ############################
    #         EXTRACT          #
    ############################

    if extract:
        logger.info("Starting Extract samples from study")
        oldpath=args.Path
        removepath=args.SampleList
        output_folder=args.output_folder
        extract_main(oldpath,removepath,output_folder,log)

#################################################################################################################

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':
	
    logger.remove()
    logfile="Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
    logger.level("INFO", color="<green>")
    logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True, catch=True)
    logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
    logger.info("Welcome to VARAN")
    log=True


    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='Argument of Varan script')
    
    # WALK BLOCK
    parser.add_argument('-c', '--Cancer', required=False,
                        help='Cancer Name')
    parser.add_argument('-i', '--input', required=False,
                                            help='input folder tsv with data or tsv with path of data')
    parser.add_argument('-f', '--filter_snv', required=False,
                                            action='store_true',
                                            help='Filter out from the vcf the variants wit dot (.) in Alt column')
    parser.add_argument('-t', '--vcf_type', required=False, 
                                            choices=['snv', 'cnv'],
                                            help='Select the vcf file to parse')
    parser.add_argument('-w', '--overWrite', required=False,action='store_true',
                                                help='Overwrite output folder if it exists')

    parser.add_argument('-k', '--oncoKB', required=False,action='store_true',help='OncoKB annotation')
    # FILTER_CLINVAR BLOCK

    parser.add_argument('-v', '--vus', required=False,
                                            action='store_true',
                                            help='Filter out VUS variants')
    
    # UPDATE BLOCK

    parser.add_argument('-u', '--Update', required=False,action='store_true',
                                                help='Add this argument if you want to concatenate two studies')

    parser.add_argument('-n', '--NewPath', required=False,help='Path of new study folder to add')
        
    # DELETE BLOCK

    parser.add_argument('-r', '--Remove', required=False,action='store_true',
                                                help='Add this argument if you want to remove samples from a study')
    
    # EXTRACT BLOCK

    parser.add_argument('-e', '--Extract', required=False,action='store_true',
                                                help='Add this argument if you want to extract samples from a study')
    
    # COMMON BLOCK 

    parser.add_argument('-o', '--output_folder', required=True,
                                            help='Output folder')
    parser.add_argument('-s', '--SampleList', required=False,
                            help='Path of file with list of SampleIDs')
    parser.add_argument('-p', '--Path', required=False,
                                                help='Path of original study folder')
    
    try:
        args = parser.parse_args()

        cancer = args.Cancer
        input = args.input
        filter_snv=args.filter_snv
        output_folder = args.output_folder
        vcf_type=args.vcf_type
        overwrite_output=args.overWrite
        vus=args.vus
        oncoKB=args.oncoKB
        
        update=args.Update
        extract=args.Extract
        remove=args.Remove
        
        if not any([args.Update ,args.Extract , args.Remove]) and args.input==None:
            logger.critical("Error Argument: Input is required")
            raise argparse.ArgumentError("Input is required")
        
        if not any([args.Update ,args.Extract , args.Remove]) and args.Cancer==None:
            logger.critical("Error Argument: Cancer name is required")
            raise argparse.ArgumentError("Cancer is required")

        if args.Update and (args.Path==None or args.NewPath==None):
            logger.critical("To update a study, you need to specify both original and new folder paths")
            raise argparse.ArgumentError("To update a study, you need to specify both old and new paths")

        if (any([args.Remove,args.Extract]) and args.Path==None) or (any([args.Remove,args.Extract]) and args.SampleList==None):
            logger.critical("To remove/extract samples from a study, you need to specify both original folder path and list samples")
            raise argparse.ArgumentError("To remove/extract samples from a study, you need to specify both original folder path and list samples")
        
        varan(input, cancer, output_folder,oncoKB, filter_snv, vcf_type, overwrite_output, vus, update, extract, remove, log)
    
    except Exception as err:
        logger.critical(f"error: {err}", file=sys.stderr)