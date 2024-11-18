import os
import argparse
from loguru import logger
import sys
from configparser import ConfigParser
import subprocess
import shutil
from walk import write_filters_in_report

config = ConfigParser()
configFile = config.read("conf.ini")

ZIP_MAF = config.get('Zip', 'ZIP_MAF')
ZIP_SNV_FILTERED = config.get('Zip', 'ZIP_SNV_FILTERED')
        
        
def cBio_validation(output_folder):
    config = ConfigParser()
    config.read('conf.ini')
    PORT = config.get('Validation', 'PORT')
    logger.info(f"Starting online validation. Connecting to {PORT}")
    
    try:
        process1 = subprocess.Popen(["python3", "importer/validateData.py", "-s", output_folder, "-u", PORT, "-e", os.path.join(output_folder, "report.txt"), "--html_table", os.path.join(output_folder, "report_validate.html"), "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout1, stderr1 = process1.communicate()
        warn = stderr1
        if process1.returncode == 1:
            raise subprocess.CalledProcessError(process1.returncode, process1.args, output=stdout1, stderr=stderr1)
        logger.info(stdout1)
        
        return process1.returncode

    except subprocess.CalledProcessError as e:
        logger.error("Connection to localhost failed. This may be due to an incorrect port selection" +\
                     " or invalid Docker settings.")
        logger.info("Starting offline validation...")

        process2 = subprocess.Popen(['python3', 'importer/validateData.py', '-s', output_folder, '-n', "-e", os.path.join(output_folder, "report.txt"), "--html_table", os.path.join(output_folder, "report_validate.html"), "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        _, stderr2 = process2.communicate()
        warn = stderr2
        if process2.returncode == 1:
            logger.error(f"Error: {stderr2.strip()} Check the report files in the study folder for more info!")
        
        check_process_status(process1, warn)

        if process1.returncode == 1 and 'process2' in locals():
            check_process_status(process2, warn)

        return process2.returncode


def check_process_status(process, warn_msg):
    if process.returncode in [2, 3]:
        logger.warning(f"{warn_msg.strip()} Check report files in the study folder for details!")
    elif process.returncode == 0:
        logger.success("The validation proceeded without errors and warnings! The study is ready to be uploaded!")


def clean_multi(input_folder, folder, file):
    mypath = os.path.join(input_folder, folder, file)
    if os.path.exists(mypath) and os.path.isdir(mypath):
        shutil.rmtree(mypath)
    if os.path.exists(mypath) and os.path.isfile(mypath):
        os.remove(mypath)


def validateFolderlog(folder):
    """
    Validates the contents of the folder against required files for cBioPortal data upload.
    
    This function checks the contents of the folder against a set of required files for different categories
    (e.g., Patient, Study, CNA, Fusion, SNV) that are necessary for uploading data to cBioPortal. It logs any
    missing files and provides a success message if all required files are present.

    Args:
        folder (str): Path to the folder to be validated.
        logfile (str): Path to the log file where validation messages will be logged.

    Notes:
        - The function checks the presence of required files within the specified 'folder' and its subdirectories.
        - Required file paths are defined for each category in the 'required_files' dictionary.
        - If any required file is missing, a warning message is logged along with the missing file names.
        - If all required files are present for all categories, a success message is logged.

    Example:
        >>> validateFolderlog('data_folder/', 'validation_log.txt')
    """
    list_files=[]
    for file in os.listdir(folder):
        if os.path.isdir(os.path.join(folder,file)):
            subdir=file
            sudbirfiles=os.listdir(os.path.join(folder,subdir))
            for subdirfile in sudbirfiles:
                list_files.append(os.path.join(subdir,subdirfile))
        else:
            list_files.append(file)

 
    # Define required files for each category
    required_files = {
        "Patient": [
            "data_clinical_patient.txt",
            "meta_clinical_patient.txt",
        ],
        "Study": [
            "data_clinical_sample.txt",
            "meta_study.txt",
            "meta_clinical_sample.txt",
        ],
        "CNA": [
            "case_lists/cases_cna.txt",
            "data_cna.txt",
            "data_cna_hg19.seg",
            "meta_cna.txt",
            "meta_cna_hg19_seg.txt",
        ],
        "Fusion": [
            "case_lists/cases_sv.txt",
            "data_sv.txt",
            "meta_sv.txt",
        ],
        "SNV": [
            "case_lists/cases_sequenced.txt",
            "data_mutations_extended.txt",
            "meta_mutations_extended.txt",
        ],
    }
    
    result_all = {}
    for category, required_files_list in required_files.items():
        missing_files = [elem for elem in required_files_list if elem not in list_files]
        result_all[category] = len(missing_files) == 0
    
        if not result_all[category]:
            logger.warning(f"Missing required file(s) for {category} from {folder}")
            logger.warning(f"Missing file(s):")
            for missing in missing_files:
                logger.warning("- " + missing)
        
            
    if all(result_all.values()):
        logger.success("Folder contains all required files for cBioportal")


def validateFolder(folder,log=False):
    """
    Validates the contents of a folder against required files for cBioPortal data upload.

    This function checks the contents of a folder against a set of required files for different categories
    (e.g., Patient, Study, CNA, Fusion, SNV) that are necessary for uploading data to cBioPortal. It prints
    any missing files and associated warnings.

    Args:
        folder (str): Path to the folder to be validated.

    Returns:
        None

    Notes:
        - The function checks the presence of required files within the specified 'folder' and its subdirectories.
        - Required file paths are defined for each category in the 'required_files' dictionary.
        - If any required file is missing, a warning message is printed along with the missing file names.

    Example:
        >>> validateFolder('data_folder/')
        
    """
    
    if not log:
        logger.remove()
        logfile="validateFolder_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")#,mode="w")
	
    logger.info("Starting validateFolder script:")
    logger.info(f"validateFolder args [folder:{folder}]")
    
    list_files=[]
    for file in os.listdir(folder):
        if os.path.isdir(os.path.join(folder,file)):
            subdir=file
            sudbirfiles=os.listdir(os.path.join(folder,subdir))
            for subdirfile in sudbirfiles:
                list_files.append(os.path.join(subdir,subdirfile))
        else:
            list_files.append(file)

 
    # Define required files for each category
    required_files = {
        "Patient": [
            "data_clinical_patient.txt",
            "meta_clinical_patient.txt",
        ],
        "Study": [
            "data_clinical_sample.txt",
            "meta_study.txt",
            "meta_clinical_sample.txt",
        ],
        "CNA": [
            "case_lists/cases_cna.txt",
            "data_cna.txt",
            "data_cna_hg19.seg",
            "meta_cna.txt",
            "meta_cna_hg19_seg.txt",
        ],
        "Fusion": [
            "case_lists/cases_sv.txt",
            "data_sv.txt",
            "meta_sv.txt",
        ],
        "SNV": [
            "case_lists/cases_sequenced.txt",
            "data_mutations_extended.txt",
            "meta_mutations_extended.txt",
        ],
    }
    
    result_all = {}
    for category, required_files_list in required_files.items():
        missing_files = [elem for elem in required_files_list if elem not in list_files]
        result_all[category] = len(missing_files) == 0
    
        
        if not result_all[category]:
            logger.warning(f"Missing required files for {category} from {folder}")
            print(f"Missing file(s):")
            for missing in missing_files:
                print("- ", missing)
                
def validateOutput(folder, input, multi, block2=False):
    
    validateFolderlog(folder)
    val = cBio_validation(folder)

    if val != 1:
        if not block2:
            write_filters_in_report(folder)    
            
            if os.path.exists(os.path.join(folder, "maf")) and ZIP_MAF:
                logger.info("Zipping maf folder...")    
                shutil.make_archive(os.path.join(folder, "maf"), "zip", os.path.join(folder, "maf"))
                logger.info("Deleting unzipped maf folder...")
                shutil.rmtree(os.path.join(folder, "maf"))

            if os.path.exists(os.path.join(folder, "snv_filtered")) and ZIP_SNV_FILTERED:
                logger.info("Zipping snv_filtered folder...") 
                shutil.make_archive(os.path.join(folder, "snv_filtered"), "zip", os.path.join(folder, "snv_filtered"))
                logger.info("Deleting unzipped snv_filtered folder...")
                shutil.rmtree(os.path.join(folder,"snv_filtered"))

            if os.path.exists(os.path.join(folder, "temp")):
                logger.info("Deleting temp folder")
                shutil.rmtree(os.path.join(folder, "temp"))

            if multi:
                clean_multi(input[0], "CNV", "single_sample_vcf")
                clean_multi(input[0], "CNV", "sample_id.txt")
                clean_multi(input[0], "SNV", "single_sample_vcf")
                clean_multi(input[0], "SNV", "sample_id.txt")

        logger.success("The study is ready to be uploaded on cBioportal")
    
    else:
        raise Exception("Validation Failed!")