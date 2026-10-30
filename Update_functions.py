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

"""Update (-u): output folder preparation, meta file carry-over and the
merge of two studies' data files. The per-file merge rules live in
study_files.py, shared with extract and remove.
"""

import re
import shutil
from pathlib import Path

import pandas as pd
from loguru import logger

from study_files import merge_study_file

from versioning import (
    create_newest_version_folder,
    get_version_list,
)



def safe_check_file(oldpath: Path, newpath: Path, output: Path, file: str) -> None:
    """Safely checks a file for consistency between old and new study folders.

    This function wraps study_files.merge_study_file() in a try-except block.

    Args:
        oldpath (Path): Path to the old version of the study folder.
        newpath (Path): Path to the new version of the study folder.
        output (Path): Path to the output (updated) study folder.
        file (str): Name of the file to check (e.g., "data_clinical_sample.txt").

    Raises:
        IndexError: Raised when a parsing error is detected in the file structure.

    """
    try:
        merge_study_file(oldpath, newpath, output, file)
    except pd.errors.ParserError as e:
        line_match = re.search(r"line (\d+)", str(e))
        line_number = line_match.group(1) if line_match else "unknown"
        logger.critical(
            f"Wrong column number in line {line_number} of {file} file")
        msg = "Exiting from Update script!"
        raise IndexError(msg) from e

def copy_metadata_files(oldpath: Path, newpath: Path, output: Path) -> None:
    """Copy meta_*.txt files forward from both source studies.

    Copies from oldpath, then adds any meta file that only exists in
    newpath - e.g. meta_exon_brca_cna.txt when oldpath predates that data
    type but newpath (a freshly walked batch being merged in) already has
    it. Only scanning oldpath silently dropped any such new-in-this-batch
    meta file. Which source wins on a same-named file doesn't matter: every
    file copied here is just a placeholder for remove_meta's cleanup pass
    and meta_case_main's proper regeneration afterward, both of which run
    later in update_main based on the actual merged data.

    Args:
        oldpath (Path): Path to previous version of the study folder.
        newpath (Path): Path to the incoming data folder being merged in.
        output (Path): Path to output (new version) study folder.

    Returns:
        None

    """
    meta_files = {f.name: f for f in newpath.glob("*meta*")}
    meta_files.update({f.name: f for f in oldpath.glob("*meta*")})
    if meta_files:
        for file in meta_files.values():
            shutil.copy(file, output)
    else:
        logger.warning("No meta files found in either study folder!")

def prepare_output_folder(oldpath: str,
                          output: str,
                          overwrite: bool) -> tuple[str, bool, str]:
    """Prepare output folder for the updated study version.

    This function handles:
    - Deriving the base output path if not provided.
    - Removing the last versioned folder if overwrite is enabled.
    - Creating a new versioned output folder.
    - Creating a 'case_lists' subfolder inside the output directory.

    Args:
        oldpath (str): Path to the previous version of the study folder.
        output (str): Base output path for the new version.
        overwrite (bool): Whether to overwrite the latest version folder.

    Returns:
        tuple[str, bool, str]: A tuple containing:
            - The path to the newly created output folder.
            - Boolean flag indicating whether output was auto-inferred.
            - Path to the 'case_lists' subfolder inside the output folder.

    """
    if output != "":
        no_out = False
    else:
        no_out = True
        output = re.split(r"_v[0-9]+$", oldpath)[0]

    old_versions = get_version_list(output)
    if old_versions and Path(old_versions[-1]).exists() and overwrite:
        logger.info("Overwrite option set. Start removing folder")
        shutil.rmtree(old_versions[-1])

    output = create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")

    output_caseslists = Path(output) / "case_lists"
    Path(output_caseslists).mkdir()

    return output, no_out, output_caseslists
