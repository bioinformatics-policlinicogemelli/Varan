import os
from Update_functions import *
from loguru import logger
from ValidateFolder import validateFolderlog
from versioning import *
from Make_meta_and_cases import meta_case_main
import shutil
import sys


def update_main(oldpath, newpath, output, study_id, overwrite):
    
    logger.info("Starting update_main script:")
    logger.info(f"update_main args [oldpath:{oldpath}, newpath:{newpath}, output_folder:{output}]")	

    logger.info("Checking inputs...")
    if not os.path.isdir(oldpath):
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit()
    if not os.path.isdir(newpath):
        logger.critical(f"{newpath} is not a valid folder!")
        sys.exit()
    
    if output!="":
        no_out=False
        if os.path.exists(oldpath):    
            logger.info("Old folder found")
        if os.path.exists(newpath):
            logger.info("New folder found")
    else:
        no_out=True
        output=re.split(r'_v[0-9]$',oldpath)[0]
    
    old_versions = get_version_list(output)

    if len(old_versions) > 0 and os.path.exists(old_versions[-1]):
        if overwrite:
            logger.info(f"Overwrite option set. Start removing folder")
            shutil.rmtree(old_versions[-1])

    output = create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")

    output_caseslists = os.path.join(output,"case_lists")
    os.mkdir(output_caseslists)   

    logger.info("Great! Everything is ready to start")      

    os.system("cp "+ oldpath+"/*meta* "+output)
    
    check_files(oldpath, newpath, output, "data_clinical_sample.txt")
    check_files(oldpath, newpath, output, "data_clinical_patient.txt")
    check_files(oldpath, newpath, output, "data_cna_hg19.seg")
    check_files(oldpath, newpath, output, "data_cna.txt")
    check_files(oldpath, newpath, output, "data_mutations_extended.txt")
    check_files(oldpath, newpath, output, "data_sv.txt")

    check_files_cases(oldpath, newpath, output_caseslists,"cases_cna.txt")
    check_files_cases(oldpath, newpath, output_caseslists,"cases_sequenced.txt")
    check_files_cases(oldpath, newpath, output_caseslists,"cases_sv.txt")

    logger.info("Starting Validation Folder...")
    
    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)

    meta_case_main(cancer, output, study_info, study_id)
 
    validateFolderlog(output)

    # if len(old_versions)>=1:
    #     compare_version(newpath, oldpath, "update", output)

    compare_version(output, newpath, "update")

    logger.success("The process ended without errors")
    logger.success("Successfully updated study!")