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


"""Utility functions for handling versioned output folders.

Includes:
- Version parsing and comparison logic.
- Folder creation for versioned study outputs.
- Metadata extraction from cBioPortal-style meta_study.txt files.

"""

import re
from pathlib import Path

from loguru import logger


def extract_version_str(foldername: str) -> str:
    """Extract the version string (e.g., "_v1") from a folder name.

    Args:
        foldername (str): The folder name to extract version from.

    Returns:
        str: The version suffix (e.g., "_v1").

    """
    version = extract_version_int(foldername)
    return "_v" + str(version)


def extract_version_int(foldername: str) -> list:
    """Extract the version number as an integer from a folder name.

    Args:
        foldername (str): Folder name ending with version.

    Returns:
        list: Extracted original foldername version.

    """
    if isinstance(foldername, Path):
        foldername = foldername.name
    elif not isinstance(foldername, str):
        return None

    match = re.search(r"_v(\d+)$", foldername)
    if match:
        return int(match.group(1))
    return None


def get_version_list(output_folder: str) -> list:
    """List all versioned folders related to the given output folder.

    Args:
        output_folder (str): Path to a versioned output folder.

    Returns:
        list: Sorted list of versioned folder names.

    """
    foldername = re.split(r"_v[0-9]+$",Path(output_folder).name)[0]
    outputfolderpath = Path(output_folder).parent
    if outputfolderpath == "":
        outputfolderpath = Path.cwd()

    outputfolderpath = outputfolderpath.resolve()

    old_versions = [
    f.name for f in Path(outputfolderpath).resolve().iterdir()
    if f.is_dir() and re.match(rf"^{re.escape(foldername)}_v\d+$", f.name)
]

    global old_version_exists
    old_version_exists = bool(old_versions)

    version_n = [extract_version_int(version) for version in old_versions]
    version_n = [elem for elem in version_n if isinstance(elem, int)]
    sorted_version = sorted(set(version_n), key=int)
    return [foldername + "_v" + str(x) for x in sorted_version]


def get_newest_version(output_folder: str) -> tuple:
    """Compute the next available versioned folder name.

    Args:
        output_folder (str): Base output folder name.

    Returns:
        tuple: Tuple of (new versioned folder path, previous version suffix).

    """
    foldername = re.split(r"_v[0-9]+$", Path(output_folder).name)[0]
    outputfolderpath = Path(output_folder).parent
    if outputfolderpath == "":
        outputfolderpath = Path.cwd()

    old_versions = [file.name for file in outputfolderpath.iterdir()
                    if re.match(rf"^{re.escape(foldername)}_v[0-9]+$", file.name)]

    logger.info(
    f"{len(old_versions)} version(s) of the selected output folder found: "
    f"{old_versions}")

    old_versions_number = list(map(extract_version_int, old_versions))
    if old_versions_number == []:
        v = "_v0"
        version = "_v1"
    else:
        v = max(old_versions_number)
        version = "_v" + str(v + 1)

    output_folder_version = foldername + version

    return Path(output_folder).parent / output_folder_version, f"_v{v}"


def create_newest_version_folder(outputfolder: str) -> str:
    """Create a new folder with the next version suffix.

    Args:
        outputfolder (str): Base folder path.

    Returns:
        str: Path to the newly created versioned folder.

    """
    if len(get_version_list(outputfolder)) == 0:
        version = "_v1"
        outputfolder_newest_version = Path(outputfolder + version)
    else:
        outputfolder_newest_version, _= get_newest_version(outputfolder)
    Path(outputfolder_newest_version).mkdir()

    return outputfolder_newest_version


def extract_info_from_meta(folder: str) -> tuple:
    """Extract cancer type, study ID, and study name from meta_study.txt.

    Args:
        folder (str): Path to folder containing "meta_study.txt".

    Returns:
        tuple: Cancer type, and a list containing.

    """
    file_meta = Path(folder) / "meta_study.txt"
    with file_meta.open() as meta:
        for line in meta:
            if line.startswith("type_of_cancer"):
                cancer = line.split(" ")[1].strip()
            if line.startswith("cancer_study_identifier"):
                study_id = re.split(r"_v[0-9]+$",line.split(" ")[1])[0].strip()
            if line.startswith("name"):
                study_name = re.split(r"V[0-9]+$",
                                      line.split(":")[1].split("(")[0]
                                      .strip())[0].strip()
    return cancer, [study_id, study_name]
