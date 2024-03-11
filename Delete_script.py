import os
import sys
import argparse
from Delete_functions import *
from ValidateFolder import validateFolderlog
import loguru
from loguru import logger

def delete_main(oldpath,removepath,destinationfolder,log=False):
    
    if not log:
        logger.remove()
        logfile="delete_main_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")
    
    logger.info("Starting delete_main script:")
    logger.info(f"delete_main args [oldpath:{oldpath}, removepath:{removepath}, destinationfolder:{destinationfolder}]")	
    
    if os.path.exists(oldpath):
        logger.info("Original folder found")
    
    if os.path.exists(removepath):
        logger.info("Sample list to remove found")
    
    output=os.path.join(destinationfolder,"filtered_data")
    if os.path.exists(output):
        logger.critical("Filtered_data folder already exists. Please change destination folder (--Destination arg)" )
        logger.critical("Exit")
        sys.exit()
    else:
        logger.info("Creating a new folder: filtered_data")    
        os.mkdir(output)
        output_caseslists=os.path.join(output,"case_lists")
        os.mkdir(output_caseslists)   

    logger.info("Great! Everything is ready to start")

    os.system("cp "+oldpath+"/*meta* "+output)
    sampleIds=open(removepath,"r").readlines()
    sampleIds=[sample.strip() for sample in sampleIds]

    
    o_clinical_patient=os.path.join(oldpath,"data_clinical_patient.txt")
    if os.path.exists(o_clinical_patient):
        delete_clinical_patient(oldpath,sampleIds,output)
    else:
        logger.warning("data_clinical_patient.txt not found in current folder. Skipping")
    #
    o_clinical_sample=os.path.join(oldpath,"data_clinical_sample.txt")
    if os.path.exists(o_clinical_sample) :
        delete_clinical_samples(o_clinical_sample,sampleIds,output)
    else:
        logger.warning("data_clinical_sample.txt not found in current folder. Skipping")
    
    o_cna_hg19=os.path.join(oldpath,"data_cna_hg19.seg")
    if os.path.exists(o_cna_hg19):
        delete_cna_hg19(o_cna_hg19,sampleIds,output)
    else:
        logger.warning("data_cna_hg19.seg not found in current folder. Skipping")
    
    #
    o_cna=os.path.join(oldpath,"data_cna.txt")
    if os.path.exists(o_cna):
        delete_cna(o_cna,sampleIds,output)
    else:
        logger.warning("data_cna.txt not found in current folder. Skipping")
    
    #
    o_mutations=os.path.join(oldpath,"data_mutations_extended.txt")
    if os.path.exists(o_mutations):
        delete_mutations(o_mutations,sampleIds,output)
    else:
        logger.warning("data_mutations_extended.txt not found in current folder. Skipping")
    #
    o_sv=os.path.join(oldpath,"data_sv.txt")
    if os.path.exists(o_sv):
        delete_sv(o_sv,sampleIds,output)
    else:
        logger.warning("data_sv.txt not found in current folder. Skipping")
    #
    
    o_cases_cna=os.path.join(oldpath,"case_lists/cases_cna.txt")
    if os.path.exists(o_cases_cna):
        delete_caselist_cna(o_cases_cna,sampleIds,output_caseslists)
    else:
        logger.warning("cases_cna.txt not found in 'case_lists' folder. Skipping")
    
    o_cases_sequenced=os.path.join(oldpath,"case_lists/cases_sequenced.txt")
    if os.path.exists(o_cases_sequenced):
        delete_caselist_sequenced(o_cases_sequenced,sampleIds,output_caseslists)
    else:
        logger.warning("cases_sequenced.txt not found in 'case_lists' folder. Skipping")
    #  
    o_cases_sv=os.path.join(oldpath,"case_lists/cases_sv.txt")
    if os.path.exists(o_cases_sv):
        delete_caselist_sv(o_cases_sv,sampleIds,output_caseslists)
    else:
        logger.warning("cases_sv.txt not found in 'case_lists' folder. Skipping")



    logger.info("Starting Validation Folder...")

    validateFolderlog(output)
    
    logger.success("The process ended without errors")
    logger.success("Please, check DeleteScript.log to verify that everything went as expected.")
    logger.success("Successfully deleted sample(s)!")
    
class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)     
    
if __name__ == '__main__':
    
    logfile="Delete_Script_{time:HH:mm:ss.SS}.log"
    
    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='Parser of Delete script for cBioportal')
     
    parser.add_argument('-o', '--OldDataPath', required=True,
                        help='Folder containing old existing data files')
    parser.add_argument('-s', '--SampleToRemove', required=True,
                        help='Path of file with SampleIDs to remove')
    parser.add_argument('-d', '--Destination', required=True,
                        help='Path of new folder to store updated data',default="./")
    
    try:
        args = parser.parse_args()
    except Exception as err:
        logger.remove()
        logfile="delete_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True,catch=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
        logger.critical(f"error: {err}", file=sys.stderr)

    oldpath=args.OldDataPath
    removepath=args.SampleToRemove
    destinationfolder=args.Destination
    
    delete_main(oldpath,removepath,destinationfolder,log=False)