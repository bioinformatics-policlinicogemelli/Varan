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
from datetime import datetime
from pathlib import Path

from loguru import logger

from filter_clinvar import check_bool
from log_utils import attach_file_log_sink
from Make_meta_and_cases import meta_case_main
from populate_case_lists import profiled_from_case_lists
from study_files import STUDY_FILES
from Update_functions import (
    copy_metadata_files,
    prepare_output_folder,
    safe_check_file,
)
from ValidateFolder import copy_maf, remove_meta, validate_output
from versioning import extract_info_from_meta
from write_report import write_report_failure, write_report_update
from config_loader import config



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
    start_time = datetime.now().astimezone().strftime("%d/%m/%Y, %H:%M:%S")
    logger.info("Starting update_main script:")
    logger.info(
    f"update_main args [oldpath:{oldpath}, newpath:{newpath}, "
    f"output_folder:{output}]")

    oldpath = oldpath.rstrip("/")

    logger.info("Checking inputs...")
    if not Path(oldpath).is_dir():
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit(1)
    if not Path(newpath).is_dir():
        logger.critical(f"{newpath} is not a valid folder!")
        sys.exit(1)

    output, no_out, output_caseslists = prepare_output_folder(
        oldpath, output, overwrite)

    # The run's log file lives inside its own output (study) folder - see
    # attach_file_log_sink()'s docstring for why any failure to write it is
    # swallowed rather than crashing an otherwise-fine run.
    attach_file_log_sink(
        Path(output), "Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log")

    logger.info("Great! Everything is ready to start")

    oldpath = Path(oldpath)
    output = Path(output)

    # `stage`/`cancer` let the failure report (if the except below fires)
    # say more than just "it failed" - `cancer` starts as None since it
    # isn't known until extract_info_from_meta() runs partway through.
    stage = "Copying metadata files"
    cancer = None
    try:
        copy_metadata_files(Path(oldpath), Path(newpath), Path(output))

        stage = "Merging data files"
        for study_file in STUDY_FILES:
            safe_check_file(oldpath, newpath, output, study_file.name)

        remove_meta(output)

        stage = "Creating tables"
        cancer, study_info = extract_info_from_meta(oldpath)
        study_info.append(oldpath)
        study_info.append(no_out)

        # Samples tested without events only survive through the source
        # studies' case lists - the merged data files don't show them.
        meta_case_main(
            cancer, output, study_info, study_id,
            profiled=profiled_from_case_lists([oldpath, newpath]))

        stage = "Copying MAF files"
        copy_maf_flag = check_bool(config.get("Zip", "COPY_MAF"))
        zip_maf_flag = check_bool(config.get("Zip", "ZIP_MAF"))
        copy_maf(oldpath, output, copy_maf_flag, zip_maf_flag)
        copy_maf(newpath, output, copy_maf_flag, zip_maf_flag)

        stage = "Validation"
        logger.info("Starting Validation Folder...")
        number_for_graph = validate_output(output, None, False, True, None, None, None)

        stage = "Writing report_VARAN.html"
        logger.info("Starting writing report_VARAN.html...")
        write_report_update(oldpath, newpath, output, number_for_graph, start_time)
    except Exception as err:
        logger.critical(f"Update failed during '{stage}': {err}")
        write_report_failure(
            output, err, cancer=cancer, start_time=start_time, stage=stage)
        raise

    logger.success("The process ended without errors")
    logger.success("Successfully updated study!")
