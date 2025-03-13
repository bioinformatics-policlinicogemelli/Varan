#Copyright 2025 bioinformatics-policlinicogemelli

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

import os
from Update_functions import *
from loguru import logger
from ValidateFolder import validateOutput, copy_maf
from versioning import *
from Make_meta_and_cases import meta_case_main
import shutil
import sys
from write_report import *
from filter_clinvar import check_bool

config = ConfigParser()
configFile = config.read("conf.ini")


def update_main(oldpath, newpath, output, study_id, overwrite):
    
    logger.info("Starting update_main script:")
    logger.info(f"update_main args [oldpath:{oldpath}, newpath:{newpath}, output_folder:{output}]")	
    oldpath = oldpath.rstrip("/")
    
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
        output=re.split(r'_v[0-9]+$', oldpath)[0]

    old_versions = get_version_list(output)

    if len(old_versions) > 0 and os.path.exists(old_versions[-1]):
        if overwrite:
            logger.info(f"Overwrite option set. Start removing folder")
            shutil.rmtree(old_versions[-1])

    output = create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")

    img_path = os.path.join(oldpath, "img", "logo_VARAN.png")
    if os.path.exists(img_path):
        img_output_dir = os.path.join(output, "img")
        os.makedirs(img_output_dir, exist_ok=True)
        shutil.copy(img_path, os.path.join(img_output_dir, "logo_VARAN.png"))  

    output_caseslists = os.path.join(output,"case_lists")
    os.mkdir(output_caseslists)   

    logger.info("Great! Everything is ready to start")
    os.system("cp " + oldpath + "/*meta* " + output)
    
    file_names = ["data_clinical_sample.txt", "data_clinical_patient.txt", "data_cna_hg19.seg", "data_cna_hg19.seg.fc.txt", "data_cna.txt", "data_mutations_extended.txt", "data_sv.txt"]
    for file in file_names:
        try:
            check_files(oldpath, newpath, output, file)
        except pd.errors.ParserError as e:
            line_number = int(re.search(r'line (\d+)', str(e)).group(1))
            logger.critical(f"error: Wrong column number in line {line_number} of {file} file")
            raise(IndexError("Exiting from Update script!"))


    check_files_cases(oldpath, newpath, output_caseslists,"cases_cna.txt")
    check_files_cases(oldpath, newpath, output_caseslists,"cases_sequenced.txt")
    check_files_cases(oldpath, newpath, output_caseslists,"cases_sv.txt")
    
    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)

    meta_case_main(cancer, output, study_info, study_id)

    COPY_MAF = config.get('Zip', 'COPY_MAF')
    COPY_MAF = check_bool(COPY_MAF)
    ZIP_MAF = config.get('Zip', 'ZIP_MAF')
    ZIP_MAF = check_bool(ZIP_MAF)
    copy_maf(oldpath, output, COPY_MAF, ZIP_MAF)
    copy_maf(newpath, output, COPY_MAF, ZIP_MAF)
    
    logger.info("Starting Validation Folder...")
    number_for_graph = validateOutput(output, None, False, True, None, None, None)

    logger.info("Starting writing report_VARAN.html...")
    write_report_update(oldpath, newpath, output, number_for_graph)

    logger.success("The process ended without errors")
    logger.success("Successfully updated study!")



