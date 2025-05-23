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

"""Module: extract_main.

Main script for extracting specific samples from a versioned study folder
in a cBioPortal-style workflow.

Provides:
  - Validation of input paths and overwrite handling
  - Sample list parsing and extraction of clinical, CNA, mutation, and SV data
  - Versioned output folder creation and metadata copy
  - Report generation and folder validation
"""

import re
import shutil
import sys
from configparser import ConfigParser
from pathlib import Path

from loguru import logger

from ExtractSamples_functions import (
    check_sample_list,
    extract_all_data,
)
from filter_clinvar import check_bool
from Make_meta_and_cases import meta_case_main
from ValidateFolder import copy_maf, validateOutput
from versioning import (
    create_newest_version_folder,
    extract_info_from_meta,
    get_version_list,
)
from write_report import write_report_extract

config = ConfigParser()
configFile = config.read("conf.ini")

def extract_main(oldpath: str,
                 extract_path: str,
                 output: str,
                 study_id: str,
                 overwrite: str) -> None:
    """Extract specified samples from an existing study folder.

    Validates inputs, creates a versioned output folder,
    copies metadata, and extracts clinical, CNA, mutation,
    and structural variation data for samples listed in
    `extract_path`.

    Args:
        oldpath (str): Path to the original versioned study folder.
        extract_path (str): Path to file containing one sample ID per line.
        output (str): Base output path. If empty, inferred from oldpath.
        study_id (str): Identifier for the study.
        overwrite (bool): Whether to overwrite the latest version folder if it exists.

    Returns:
        None

    """
    logger.info(
        f"extract_main args [old_path:{oldpath}, "
        f"extract_path:{extract_path}, output_folder:{output}]")
    logger.info("Checking input...")
    oldpath = oldpath.rstrip("/")

    if not Path(oldpath).is_dir():
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit()

    if output!="":
        no_out=False
        if Path(oldpath).exists():
            logger.info("Original folder found")
        if Path(extract_path).exists():
            logger.info("Sample list to extract found")
    else:
        no_out=True
        output=re.split(r"_v[0-9]+$",oldpath)[0]

    check_sample_list(extract_path, oldpath)
    old_versions=get_version_list(output)

    if overwrite and len(old_versions)>0 and Path(old_versions[-1]).exists():
        logger.info("Overwrite option set. Start removing folder")
        shutil.rmtree(old_versions[-1])

    output=create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")

    output_caseslists=Path(output) / "case_lists"
    output_caseslists.mkdir(parents=True, exist_ok=True)

    logger.info("Great! Everything is ready to start")
    meta_files = Path(oldpath).glob("*meta*")
    for file in meta_files:
        shutil.copy(file, output)

    with Path(extract_path).open() as f:
        sample_ids = [line.strip() for line in f]

    extract_all_data(oldpath, sample_ids, output)

    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)
    meta_case_main(cancer, output, study_info, study_id)

    zip_maf_set = config.get("Zip", "ZIP_MAF")
    zip_maf_set = check_bool(zip_maf_set)
    copy_maf_set = config.get("Zip", "COPY_MAF")
    copy_maf_set = check_bool(copy_maf_set)
    copy_maf(oldpath, output, copy_maf_set, zip_maf_set)

    logger.info("Starting Validation Folder...")
    number_for_graph = validateOutput(output, None, False, True, None, None, None)

    logger.info("Starting writing report_VARAN.html...")
    write_report_extract(oldpath, output, number_for_graph)

    logger.success("The process ended without errors")
    logger.success("Successfully extracted sample(s)!")
