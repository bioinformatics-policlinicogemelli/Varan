import os
import sys
import argparse
from Update_functions import *
import loguru
from loguru import logger
from ValidateFolder import validateFolderlog


def update_main(path,newpath,output,log=False):
    
    if not log:
        logger.remove()
        logfile="update_main_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")
    
    logger.info("Starting update_main script:")
    logger.info(f"update_main args [oldpath:{path}, newpath:{newpath}, output_folder:{output}]")	

    if os.path.exists(path):    
        logger.info("Old folder found")

    if os.path.exists(newpath):
        logger.info("New folder found")
    

    output=os.path.join(output,"updated_data")
    
    if os.path.exists(output):
        logger.critical("Updated_data folder for this study already exists. Please change destination folder (--Destination arg)" )
        logger.critical("Exit")
        sys.exit()
    else:    
        logger.info("Creating a new folder: updated_data")
        os.mkdir(output)
        output_caseslists=os.path.join(output,"case_lists")
        os.mkdir(output_caseslists)   

    logger.info("Great! Everything is ready to start")      

    os.system("cp "+ path+"/*meta* "+output)
    
    o_clinical_sample=os.path.join(path,"data_clinical_sample.txt")
    n_clinical_sample=os.path.join(newpath,"data_clinical_sample.txt")
    if os.path.exists(o_clinical_sample) and os.path.exists(n_clinical_sample):
        update_clinical_samples(o_clinical_sample,n_clinical_sample,output)
    else:
        logger.warning("data_clinical_sample.txt not found in current folders. Skipping")

    o_clinical_patient=os.path.join(path,"data_clinical_patient.txt")
    n_clinical_patient=os.path.join(newpath,"data_clinical_patient.txt")
    if os.path.exists(o_clinical_patient) and os.path.exists(n_clinical_patient):
        update_clinical_patient(o_clinical_patient,n_clinical_patient,output)
    else:
        logger.warning("data_clinical_patient.txt not found in current folders. Skipping")
    #
    o_cna_hg19=os.path.join(path,"data_cna_hg19.seg")
    n_cna_hg19=os.path.join(newpath,"data_cna_hg19.seg")
    if os.path.exists(o_cna_hg19) and os.path.exists(n_cna_hg19):
        update_cna_hg19(o_cna_hg19,n_cna_hg19,output)
    else:
        logger.warning("data_cna_hg19_seg.txt not found in current folders. Skipping")
    #
    o_cna=os.path.join(path,"data_cna.txt")
    n_cna=os.path.join(newpath,"data_cna.txt")
    if os.path.exists(o_cna) and os.path.exists(n_cna):
        update_cna(o_cna,n_cna,output)
    else:
        logger.warning("data_cna.txt not found in current folders. Skipping")
    #
    o_mutations=os.path.join(path,"data_mutations_extended.txt")
    n_mutations=os.path.join(newpath,"data_mutations_extended.txt")
    if os.path.exists(o_mutations) and os.path.exists(n_mutations):
        update_mutations(o_mutations,n_mutations,output)
    else:
        logger.warning("data_mutations_extended.txt not found in current folders. Skipping")
    
    o_sv=os.path.join(path,"data_sv.txt")
    n_sv=os.path.join(newpath,"data_sv.txt")
    if os.path.exists(o_sv) and os.path.exists(n_sv):
        update_sv(o_sv,n_sv,output)
    else:
        logger.warning("data_sv.txt not found in current folders. Skipping")
  
    o_cases_cna=os.path.join(path,"case_lists/cases_cna.txt")
    n_cases_cna=os.path.join(newpath,"case_lists/cases_cna.txt")
    if os.path.exists(o_cases_cna) and os.path.exists(n_cases_cna):
        update_caselist_cna(o_cases_cna,n_cases_cna,output_caseslists)
    else:
        logger.warning("cases_cna.txt not found in 'case_lists' folders. Skipping")

    #
    o_cases_sequenced=os.path.join(path,"case_lists/cases_sequenced.txt")
    n_cases_sequenced=os.path.join(newpath,"case_lists/cases_sequenced.txt")
    if os.path.exists(o_cases_sequenced) and os.path.exists(n_cases_sequenced):
        update_caselist_sequenced(o_cases_sequenced,n_cases_sequenced,output_caseslists)
    else:
        logger.warning("cases_sequenced.txt not found in 'case_lists' folders. Skipping")

    o_cases_sv=os.path.join(path,"case_lists/cases_sv.txt")
    n_cases_sv=os.path.join(newpath,"case_lists/cases_sv.txt")
    if os.path.exists(o_cases_sv) and os.path.exists(n_cases_sv):
        update_caselist_sv(o_cases_sv,n_cases_sv,output_caseslists)
    else:
        logger.warning("cases_sv.txt not found in 'case_lists' folders. Skipping")

    logger.info("Starting Validation Folder...")

    validateFolderlog(output)

    logger.success("The process ended without errors")
    logger.success("Successfully updated study!")
        
class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':
        
    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='Parser of Update script for cBioportal')

    parser.add_argument('-o', '--OldDataPath', required=True,
                        help='Folder containing old existing data files')
    parser.add_argument('-n', '--NewDataPath', required=True,
                        help='Path containing new data files to add to old ones')
    parser.add_argument('-d', '--Destination', required=True,
                        help='Path of new folder to store updated data',default="./")
    
    try:    
        args = parser.parse_args()
    except Exception as err:
        logger.remove()
        logfile="update_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True,catch=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
        logger.critical(f"error: {err}", file=sys.stderr)

    path=args.OldDataPath
    newpath=args.NewDataPath
    output=args.Destination

    update_main(path,newpath,output,log=False)
