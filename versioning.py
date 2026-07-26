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
import subprocess
from pathlib import Path

from loguru import logger


def get_git_version() -> str:
    """Return Varan's version from the nearest git tag, or "unknown".

    Lives here (not in write_report.py, where it used to be) so varan.py can
    print the version banner - which happens before argument parsing, and
    therefore before conf.ini's path is known - without importing a module
    that reads conf.ini at import time.

    Returns:
        str: The tag name with a leading "v" stripped, or "unknown" if this
            isn't a git checkout or no tag exists.

    """
    try:
        repo_dir = Path(__file__).resolve().parent
        version = subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0"],
            cwd=repo_dir,
            stderr=subprocess.DEVNULL,
        ).decode().strip()

        return version.lstrip("v")
    except Exception:
        return "unknown"


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


def create_newest_version_folder(outputfolder: str, max_retries: int = 5) -> str:
    """Create a new folder with the next version suffix.

    Scans for the next free version number and creates it right after, with
    no lock in between - so two Varan runs started against the same output
    folder at nearly the same time (a real risk on a shared cluster) can
    both compute the same next version number. Rather than let the second
    one crash with a raw FileExistsError, this re-scans and retries a
    handful of times: by the time it retries, the first run's mkdir() has
    already landed, so the re-scan picks the next number after it.

    Args:
        outputfolder (str): Base folder path.
        max_retries (int): How many times to re-scan and retry after losing
            a race to another process, before giving up.

    Returns:
        str: Path to the newly created versioned folder.

    """
    for attempt in range(max_retries):
        if len(get_version_list(outputfolder)) == 0:
            version = "_v1"
            outputfolder_newest_version = Path(outputfolder + version)
        else:
            outputfolder_newest_version, _= get_newest_version(outputfolder)
        try:
            Path(outputfolder_newest_version).mkdir()
            return outputfolder_newest_version
        except FileExistsError:
            logger.warning(
                f"'{outputfolder_newest_version}' already exists - another "
                "process (likely a concurrent Varan run against the same "
                f"output) claimed it first. Re-checking and retrying with "
                f"the next version number (attempt {attempt + 1}/{max_retries}).")

    msg = (
        f"Could not claim a new version folder for '{outputfolder}' after "
        f"{max_retries} attempts - too many concurrent runs against the same "
        "output folder. Please retry, or use a different output folder.")
    logger.critical(msg)
    raise FileExistsError(msg)


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
