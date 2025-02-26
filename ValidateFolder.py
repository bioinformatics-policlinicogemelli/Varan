import os
import argparse
from loguru import logger
import sys
from configparser import ConfigParser
import subprocess
import shutil
from write_report import *
from Create_graphs import *
from filter_clinvar import check_bool
import zipfile

config = ConfigParser()
configFile = config.read("conf.ini")

def cBio_validation(output_folder):
    config = ConfigParser()
    config.read('conf.ini')
    PORT = config.get('Validation', 'PORT')
    logger.info(f"Starting online validation. Connecting to {PORT}...")
    
    try:
        process1 = subprocess.Popen(["python3", "importer/validateData.py", "-s", output_folder, "-u", PORT, "-html", os.path.join(output_folder, "report_validate.html"), "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout1, stderr1 = process1.communicate()
        warn = stderr1
        if process1.returncode == 1:
            raise subprocess.CalledProcessError(process1.returncode, process1.args, output=stdout1, stderr=stderr1)
        check_process_status(process1, warn)
        return process1.returncode

    except subprocess.CalledProcessError as e:
        logger.warning("Connection to localhost failed. This may be due to an incorrect port selection." +\
                     " This warning is expected if Varan Docker version is in use.")
        logger.info("Starting offline validation... Be warned that files succeeding this validation may still fail to load (correctly).")

        process2 = subprocess.Popen(["python3", "importer/validateData.py", "-s", output_folder, "-n", "-html", os.path.join(output_folder, "report_validate.html"), "-v"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout2, stderr2 = process2.communicate()
        warn = stderr2
        if process2.returncode == 1:
            logger.error(f"Error: {stderr2.strip()} Check the report file in the study folder for more info!")
        
        if process1.returncode == 1 and 'process2' in locals():
            check_process_status(process2, warn)
        return process2.returncode


def check_process_status(process, warn_msg):
    if process.returncode in [2, 3]:
        logger.warning(f"{warn_msg.strip()} Check the report file in the study folder for details!")
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
                list_files.append(os.path.join(subdir, subdirfile))
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
            logger.warning(f"Missing file(s) for {category} from {folder}:")
            for missing in missing_files:
                logger.warning("- " + missing)
        
            
    if all(result_all.values()):
        logger.success("Folder contains all required files for cBioportal")

def validateOutput(folder, input, multi, block2=False, cancer=None, oncoKB=None, filters=None):
    validateFolderlog(folder)
    val = cBio_validation(folder)

    cases_path = os.path.join(folder, "case_lists")
    if os.path.exists(cases_path) and not os.listdir(cases_path):
        shutil.rmtree(cases_path)

    if val != 1:
        number_for_graph = int(create_barplots(folder))
        if not block2:  
            write_report_main(folder, cancer, oncoKB, filters, number_for_graph)

            maf_path = os.path.join(folder, "maf")
            snv_path = os.path.join(folder, "snv_filtered")
            temp_path = os.path.join(folder, "temp")
            
            if os.path.exists(maf_path):
                ZIP_MAF = config.get('Zip', 'ZIP_MAF')
                ZIP_MAF = check_bool(ZIP_MAF)
                if not os.listdir(maf_path):
                    shutil.rmtree(maf_path)

                elif ZIP_MAF:
                    logger.info("Zipping maf folder...")    
                    shutil.make_archive(maf_path, "zip", maf_path)
                    logger.info("Deleting unzipped maf folder...")
                    shutil.rmtree(maf_path)

            ZIP_SNV_FILTERED = config.get('Zip', 'ZIP_SNV_FILTERED')
            ZIP_SNV_FILTERED = check_bool(ZIP_SNV_FILTERED)
            if os.path.exists(snv_path) and ZIP_SNV_FILTERED:
                logger.info("Zipping snv_filtered folder...") 
                shutil.make_archive(snv_path, "zip", snv_path)
                logger.info("Deleting unzipped snv_filtered folder...")
                shutil.rmtree(snv_path)

            if os.path.exists(temp_path):
                shutil.rmtree(temp_path)

            if multi and input != None:
                clean_multi(input[0], "CNV", "single_sample_vcf")
                clean_multi(input[0], "CNV", "sample_id.txt")
                clean_multi(input[0], "SNV", "single_sample_vcf")
                clean_multi(input[0], "SNV", "sample_id.txt")

        logger.success("The study is ready to be uploaded on cBioportal")
        return number_for_graph
    
    else:
        raise Exception("Validation Failed!")


def copy_maf(oldpath, output, COPY_MAF, ZIP_MAF):

    COPY_MAF = config.get('Zip', 'COPY_MAF')
    COPY_MAF = check_bool(COPY_MAF)
    if COPY_MAF:
        clin_sample = pd.read_csv(os.path.join(output, "data_clinical_sample.txt"), sep="\t", header=4)
        sample_IDs = clin_sample["SAMPLE_ID"]

        maf_dir = os.path.join(oldpath, 'maf')
        maf_zip_path = os.path.join(oldpath, "maf.zip")
        output_maf_dir = os.path.join(output, 'maf')

        if os.path.exists(maf_zip_path):
            with zipfile.ZipFile(maf_zip_path, 'r') as zip_maf:
                zip_maf.extractall(maf_dir)

        if os.path.exists(maf_dir):
            os.makedirs(output_maf_dir, exist_ok=True)

            for sample in sample_IDs:
                for file_type in ["_MergedSmallVariants.genome.FILTERED.vcf.maf", "_MergedSmallVariants.genome.vcf.maf"]:
                    old_file = os.path.join(maf_dir, f"{sample}{file_type}")
                    new_file = os.path.join(output_maf_dir, f"{sample}{file_type}")
                    
                    if os.path.exists(old_file):
                        shutil.copy2(old_file, new_file)

            ZIP_MAF = config.get('Zip', 'ZIP_MAF')
            ZIP_MAF = check_bool(ZIP_MAF)
            if ZIP_MAF:
                maf_zip_path = os.path.join(output, "maf.zip")
                with zipfile.ZipFile(maf_zip_path, 'w', zipfile.ZIP_DEFLATED) as zip_ref:
                    for root, _, files in os.walk(output_maf_dir):
                        for file in files:
                            file_path = os.path.join(root, file)
                            zip_ref.write(file_path, os.path.relpath(file_path, output_maf_dir))
                shutil.rmtree(output_maf_dir)

