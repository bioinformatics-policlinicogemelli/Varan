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

"""Module: update_main.py.

Main script for updating a cBioPortal-style study folder.

Includes:
- Version control and cleanup
- File checking and metadata extraction
- Study report generation and folder validation

"""

import sys
from configparser import ConfigParser
from pathlib import Path

from loguru import logger

from filter_clinvar import check_bool
from Make_meta_and_cases import meta_case_main
from Update_functions import (
    check_files_cases,
    copy_logo,
    copy_metadata_files,
    prepare_output_folder,
    safe_check_file,
)
from ValidateFolder import copy_maf, validateOutput
from versioning import extract_info_from_meta
from write_report import write_report_update


config = ConfigParser()
configFile = config.read("conf.ini")


def update_main(oldpath: str, newpath: str,
                output: str, study_id: str, overwrite: bool) -> None:
    """Update a versioned study folder using new data.

    Args:
        oldpath (str): Path to the previous version of the study folder.
        newpath (str): Path to the new data folder.
        output (str): Base output path (can be empty to auto-detect).
        study_id (str): Identifier for the study.
        overwrite (bool): Whether to overwrite the last version if it exists.

    Returns:
        None

    """
    logger.info("Starting update_main script:")
    logger.info(
    f"update_main args [oldpath:{oldpath}, newpath:{newpath}, output_folder:{output}]"
)

    oldpath = oldpath.rstrip("/")

    logger.info("Checking inputs...")
    if not Path(oldpath).is_dir():
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit()
    if not Path(newpath).is_dir():
        logger.critical(f"{newpath} is not a valid folder!")
        sys.exit()

    output, no_out, output_caseslists = prepare_output_folder(
        oldpath, output, overwrite)

    copy_logo(oldpath, output)

    logger.info("Great! Everything is ready to start")

    # Ensure paths are valid
    oldpath = Path(oldpath)
    output = Path(output)

    copy_metadata_files(Path(oldpath), Path(output))

    file_names = ["data_clinical_sample.txt", "data_clinical_patient.txt",
                  "data_cna_hg19.seg", "data_cna_hg19.seg.fc.txt", "data_cna.txt",
                  "data_mutations_extended.txt", "data_sv.txt"]

    for file in file_names:
        safe_check_file(oldpath, newpath, output, file)

    check_files_cases(oldpath, newpath, output_caseslists,"cases_cna.txt")
    check_files_cases(oldpath, newpath, output_caseslists,"cases_sequenced.txt")
    check_files_cases(oldpath, newpath, output_caseslists,"cases_sv.txt")

    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)

    meta_case_main(cancer, output, study_info, study_id)

    copy_maf_flag = check_bool(config.get("Zip", "COPY_MAF"))
    zip_maf_flag = check_bool(config.get("Zip", "ZIP_MAF"))
    copy_maf(oldpath, output, copy_maf_flag, zip_maf_flag)
    copy_maf(newpath, output, copy_maf_flag, zip_maf_flag)

    logger.info("Starting Validation Folder...")
    number_for_graph = validateOutput(output, None, False, True, None, None, None)

    logger.info("Starting writing report_VARAN.html...")
    write_report_update(oldpath, newpath, output, number_for_graph)

    logger.success("The process ended without errors")
    logger.success("Successfully updated study!")
