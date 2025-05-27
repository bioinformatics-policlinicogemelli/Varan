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

"""Module to delete samples from a study and produce a new version.

This script:
- Validates paths and inputs.
- Removes all information associated with the listed samples from clinical,
  mutation, CNA, and structural variant files.
- Creates a new versioned output folder.
- Copies relevant metadata.
- Validates the result and generates a summary HTML report.

Requires configuration in a `conf.ini` file and various helper modules:
- Delete_functions.py
- filter_clinvar.py
- Make_meta_and_cases.py
- ValidateFolder.py
- versioning.py
- write_report.py

"""
import re
import shutil
import sys
from configparser import ConfigParser
from pathlib import Path

from loguru import logger

from Delete_functions import check_sample_list, delete_all_data
from filter_clinvar import check_bool
from Make_meta_and_cases import meta_case_main
from ValidateFolder import copy_maf, validate_output
from versioning import (create_newest_version_folder, extract_info_from_meta,
                        get_version_list)
from write_report import write_report_remove

config = ConfigParser()
config_file = config.read("conf.ini")


def delete_main(oldpath: str, removepath: str, output: str,
                study_id: str, overwrite: bool) -> None:
    """Remove samples from a stdy and generate a cleaned new version.

    Args:
        oldpath (str): Path to the original dataset folder.
        removepath (str): Path to a text file listing samples to remove.
        output (str): Destination path for the new cleaned version folder.
        study_id (str): Identifier of the study (used in metadata).
        overwrite (bool): If True, overwrites existing version folder.

    Returns:
        None

    """
    logger.info(f"delete_main args [old_path:{oldpath}, remove_path:{removepath}, "
    f"destination_folder:{output}]")
    logger.info("Checking input...")

    oldpath = Path(oldpath.rstrip("/"))

    if not oldpath.is_dir():
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit()

    if output != "":
        no_out = False
        if oldpath.exists():
            logger.info("Original folder found")
        if Path(removepath).exists():
            logger.info("Sample list to remove found")
    else:
        no_out = True
        output = re.split(r"_v[0-9]+$", str(oldpath))[0]

    check_sample_list(removepath, oldpath)
    old_versions = get_version_list(output)

    last_version = Path(old_versions[-1]) if old_versions else None

    if last_version and last_version.exists() and overwrite:
        logger.info("Overwrite option set. Start removing folder")
        shutil.rmtree(last_version)

    output = create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")

    output_caseslists = Path(output) / "case_lists"
    output_caseslists.mkdir(parents=True, exist_ok=True)

    logger.info("Great! Everything is ready to start")

    for file in Path(oldpath).glob("*meta*"):
        if file.is_file():
            shutil.copy(file, Path(output))

    with Path(removepath).open() as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.warning(f"The file {removepath} contains more than a column. "
            "It may not be in the correct format!")

    with Path(removepath).open("r") as f:
        sample_ids = [line.strip() for line in f]

    delete_all_data(oldpath, sample_ids, output)

    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)

    meta_case_main(cancer, output, study_info, study_id)

    zip_maf = check_bool(config.get("Zip", "ZIP_MAF"))
    maf_copy = check_bool(config.get("Zip", "COPY_MAF"))
    copy_maf(oldpath, output, maf_copy, zip_maf)

    logger.info("Starting Validation Folder...")
    number_for_graph = validate_output(output, None, False, True, None, None, None)

    logger.info("Starting writing report_VARAN.html...")
    write_report_remove(oldpath, output, number_for_graph)

    logger.success("The process ended without errors")
    logger.success("Successfully removed sample(s)!")
