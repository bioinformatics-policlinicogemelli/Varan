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

"""Main validation and post-processing script for cBioPortal data preparation.

This script provides a set of utility functions used to validate study folders.
It verifies file structure, performs validation through cBioPortal's
`validateData.py`, generates plots and reports, handles optional compression
of SNV/MAF data, and performs clean-up of intermediate files.

Modules imported include:
- `subprocess`, `os`, `shutil`, and `zipfile` for file and process management.
- `loguru` for logging.
- Custom imports for report generation, filtering, and plotting.

Functions include:
- `cbio_validation`: Wraps validation script and parses output.
- `clean_multi`: Deletes intermediate files or directories.
- `validate_folder_log`: Checks for required files in the study folder.
- `validate_output`: Orchestrates full validation and reporting workflow.
- `copy_maf`: Copies and optionally zips MAF files based on config.

Configuration is read from `conf.ini` using `ConfigParser`.

"""

from __future__ import annotations

import collections
import os
import shutil
import subprocess
import sys
import zipfile
from configparser import ConfigParser
from pathlib import Path

import pandas as pd
from loguru import logger

from Create_graphs import create_barplots
from filter_clinvar import check_bool
from write_report import write_report_main

config = ConfigParser()
config_file = config.read("conf.ini")

def cbio_validation(output_folder: str) -> None:
    """Execute cBioPortal's validateData.py script on the specified output folder.

    Args:
        output_folder (str): The folder containing the study files to be validated.

    Behavior:
        - Runs a subprocess to invoke the cBioPortal validator.
        - Logs different outcomes based on the validator's return code.
        - Generates an HTML report in the output folder.

    Returns:
        None

    """
    config = ConfigParser()
    config.read("conf.ini")

    logger.info("Starting validation... ")
    logger.warning("Warning: succeeding files may still fail to load (correctly).")

    output_path = Path(output_folder)
    report_path = output_path / "report_validate.html"

    try:
        result = subprocess.run(
            [
                sys.executable,
                "importer/validateData.py",
                "-s", str(output_path),
                "-n",
                "-html", str(report_path),
                "-v",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    except subprocess.SubprocessError as e:
        logger.error(f"Validation subprocess failed to start: {e}")
        return

    if result.returncode == 1:
        logger.error(
            f"Error: {result.stderr.strip()} "
            f"Check the report file in the study folder for more info!",
        )
    elif result.returncode in [2, 3]:
        logger.warning(
            f"{result.stderr.strip()} "
            f"Check the report file in the study folder for details!",
        )
    elif result.returncode == 0:
        logger.success(
            "The validation proceeded without errors and warnings! "
            "The study is ready to be uploaded!",
        )

def clean_multi(input_folder: str, folder: str, file: str) -> None:
    """Remove a specified file or directory from a nested path.

    Args:
        input_folder (str): Base directory path.
        folder (str): Subfolder name inside input_folder.
        file (str): File or subdirectory name to be removed.

    Behavior:
        - Deletes file or folder if it exists at the constructed path.

    """
    mypath = Path(input_folder) / folder / file

    if mypath.exists() and mypath.is_dir():
        shutil.rmtree(mypath)
    elif mypath.exists() and mypath.is_file():
        mypath.unlink()


def validate_folder_log(folder: str) -> None:
    """Validate the folder against required files for cBioPortal data upload.

    This function checks the contents of the folder against a set of required
    files for different categories (e.g., Patient, Study, CNA, Fusion, SNV)
    that are necessary for uploading data to cBioPortal. It logs any missing
    files and provides a success message if all required files are present.

    Args:
        folder (str): Path to the folder to be validated.
        logfile (str): Path to the log file where validation messages will be logged.

    Notes:
        - The function checks the presence of required files within the specified
        'folder' and its subdirectories.
        - Required file paths are defined for each category in the 'required_files'
        dictionary.
        - If any required file is missing, a warning message is logged along with the
        missing file names.
        - If all required files are present for all categories, a success message is
        logged.

    Example:
        >>> validate_folder_log('data_folder/', 'validation_log.txt')

    """
    folder_path = Path(folder)
    list_files = []

    for item in folder_path.iterdir():
        if item.is_dir():
            subfiles = [f"{item.name}/{subitem.name}" for subitem in item.iterdir()]
            list_files.extend(subfiles)
        else:
            list_files.append(item.name)

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
    for category, required in required_files.items():
        missing = [f for f in required if f not in list_files]
        result_all[category] = len(missing) == 0

        if missing:
            logger.warning(f"Missing file(s) for {category} from {folder}:")
            for f in missing:
                logger.warning(f"- {f}")

    if all(result_all.values()):
        logger.success("Folder contains all required files for cBioportal")


def validate_output(
    folder: str,
    varan_input: tuple | list | None,
    multi: bool,
    block2: bool = False,
    cancer: str | None = None,
    oncokb: str | None = None,
    filters: dict | None = None,
) -> int:
    """Perform complete validation and post-processing for a cBioPortal study folder.

    Args:
        folder (str): Path to the folder to be validated.
        varan_input (tuple or None): Original input folder structure, used for cleanup.
        multi (bool): Whether multiple samples are present.
        block2 (bool): Flag to skip optional report generation steps.
        cancer (str or None): Cancer type to include in the report.
        oncokb (str or None): OncoKB configuration, if any.
        filters (dict or None): Filter configuration used in report generation.

    Behavior:
        - Validates required files.
        - Runs cBioPortal validator.
        - Generates plots and report.
        - Optionally zips and removes SNV/MAF/temp folders.
        - Cleans up temporary input files.

    Returns:
        int: Number of variants used in report (from plot function).

    Raises:
        Exception: If validation fails.

    """
    validate_folder_log(folder)
    val = cbio_validation(folder)

    cases_path = Path(folder) / "case_lists"
    if cases_path.exists() and not any(cases_path.iterdir()):
        shutil.rmtree(cases_path)
    if val != 1:
        number_for_graph = int(create_barplots(folder))
        if not block2:
            write_report_main(folder, cancer, filters, number_for_graph, oncokb)

            maf_path = Path(folder) / "maf"
            snv_path = Path(folder) / "snv_filtered"
            temp_path = Path(folder) / "temp"

            if maf_path.exists():
                zip_maf = config.get("Zip", "ZIP_MAF")
                zip_maf = check_bool(zip_maf)
                if not any(maf_path.iterdir()):
                    shutil.rmtree(maf_path)

                elif zip_maf:
                    logger.info("Zipping maf folder...")
                    shutil.make_archive(maf_path, "zip", maf_path)
                    logger.info("Deleting unzipped maf folder...")
                    shutil.rmtree(maf_path)

            zip_snv_filtered = config.get("Zip", "ZIP_SNV_FILTERED")
            zip_snv_filtered = check_bool(zip_snv_filtered)

            if snv_path.exists() and zip_snv_filtered:
                logger.info("Zipping snv_filtered folder...")
                shutil.make_archive(str(snv_path), "zip", str(snv_path))
                logger.info("Deleting unzipped snv_filtered folder...")
                shutil.rmtree(snv_path)

            if temp_path.exists():
                shutil.rmtree(temp_path)

            if multi and varan_input is not None:
                clean_multi(varan_input[0], "CNV", "single_sample_vcf")
                clean_multi(varan_input[0], "CNV", "sample_id.txt")
                clean_multi(varan_input[0], "SNV", "single_sample_vcf")
                clean_multi(varan_input[0], "SNV", "sample_id.txt")

        logger.success("The study is ready to be uploaded on cBioportal")
        return number_for_graph

    error_msg = "Validation Failed!"
    raise RuntimeError(error_msg)


def copy_maf(oldpath: str, output: str, copy_maf: bool, zip_maf: bool) -> None:
    """Copy MAF files from a previous path to the output folder.

    Args:
        oldpath (str): Path to the original MAF folder or ZIP archive.
        output (str): Output directory where MAF files will be copied.
        copy_maf (bool): Whether to copy the MAF files (from config, unused directly
        here).
        zip_maf (bool): Whether to compress MAF files after copying.

    Behavior:
        - Extracts ZIP if needed.
        - Matches sample IDs and copies the corresponding MAF files.
        - Recompresses the MAF folder if configured.

    Notes:
        - Reads SAMPLE_IDs from `data_clinical_sample.txt` in output folder.
        - Dynamically detects file suffix shared across MAF files.

    """
    copy_maf = config.get("Zip", "COPY_MAF")
    copy_maf = check_bool(copy_maf)
    if not copy_maf:
        return

    clin_sample = pd.read_csv(
        Path(output) / "data_clinical_sample.txt", sep="\t", header=4)
    sample_ids = clin_sample["SAMPLE_ID"]

    maf_dir = Path(oldpath) / "maf"
    maf_zip_path = Path(oldpath) / "maf.zip"
    output_maf_dir = Path(output) / "maf"
    final_zip = Path(output) / "maf.zip"

    extract_maf_zip_if_needed(maf_zip_path, maf_dir)


    if not maf_dir.exists():
        logger.warning(f"Unable to locate the MAF folder in {oldpath}.")
        logger.warning("MAF files will not be copied.")
        return

    zip_maf = config.get("Zip", "ZIP_MAF")
    zip_maf = check_bool(zip_maf)

    output_maf_dir.mkdir(parents=True, exist_ok=True)
    if zip_maf and final_zip.exists():
        with zipfile.ZipFile(final_zip, "r") as zip_existing:
            zip_existing.extractall(output_maf_dir)

    #TODO risolvere bug 
    # common_suffix = detect_common_suffix(maf_dir, sample_ids)
    common_suffix = ".hard-filtered.FILTERED.vcf.maf"

    for sample in sample_ids:
        file_name = f"{sample}{common_suffix}"
        old_file = maf_dir / file_name
        new_file = output_maf_dir / file_name
        if old_file.exists():
            shutil.copy2(old_file, new_file)

    if maf_zip_path.exists() and maf_dir.exists():
        shutil.rmtree(maf_dir)

    if zip_maf:
        logger.info(f"Zipping MAF files to {final_zip}...")
        zip_maf_files(output_maf_dir, final_zip)


def extract_maf_zip_if_needed(maf_zip_path: Path, maf_dir: Path) -> None:
    """Extract MAF ZIP archive if it exists.

    Args:
        maf_zip_path (Path): Path to the ZIP archive containing MAF files.
        maf_dir (Path): Directory where the MAF files will be extracted.

    Behavior:
        - Extracts the ZIP archive to the target directory if it exists.

    """
    if maf_zip_path.exists():
        with zipfile.ZipFile(maf_zip_path, "r") as zip_maf_file:
            zip_maf_file.extractall(maf_dir)


def detect_common_suffix(maf_dir: Path, sample_ids: list[str]) -> str:
    """Detect the most common suffix among MAF files based on sample IDs.

    Args:
        maf_dir (Path): Directory containing MAF files.
        sample_ids (list[str]): List of sample identifiers.

    Returns:
        str: The most common file suffix associated with the sample IDs.

    Behavior:
        - Iterates through MAF files and identifies the most frequent suffix
          following each sample ID.

    """
    sorted_samples = sorted(sample_ids, key=len, reverse=True)
    suffix_counter = collections.Counter()

    for file_path in maf_dir.iterdir():
        file_name = file_path.name
        for sample in sorted_samples:
            if file_name.startswith(sample):
                candidate = file_name[len(sample):]
                if candidate:
                    suffix_counter[candidate] += 1
                break

    common_suffix, _ = suffix_counter.most_common(1)[0]
    return common_suffix


def zip_maf_files(output_maf_dir: Path, zip_path: Path) -> None:
    """Compress the MAF directory into a ZIP archive.

    Args:
        output_maf_dir (Path): Path to the directory containing MAF files to compress.
        zip_path (Path): Destination path for the resulting ZIP file.

    Behavior:
        - Recursively compresses all files in the MAF directory.
        - Deletes the original uncompressed directory after zipping.

    """
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zip_ref:
        for root, _, files in os.walk(output_maf_dir):
            for file in files:
                file_path = Path(root) / file
                arcname = os.path.relpath(file_path, output_maf_dir)
                zip_ref.write(file_path, arcname)
    shutil.rmtree(output_maf_dir)


def remove_meta(output: str | Path) -> None:
    """Check existence of files and deletes metadatas if it doesn't.

    For each key-value pair in the internal mapping, the function checks whether
    the key file exists. If not, the corresponding value file is deleted.

    Parameters:
    ----------
    output : str or Path
        Path to the output directory where the files are located.
    """
    output = Path(output)

    file_map = {
        "data_cna_hg19.seg": "meta_cna_hg19_seg.txt",
        "data_mutations_extended.txt": "meta_mutations_extended.txt",
        "data_sv.txt": "meta_sv.txt",
        "data_cna.txt": "meta_cna.txt"
    }

    for data_file, meta_file in file_map.items():
        data_path = output / data_file
        meta_path = output / meta_file

        if not data_path.exists() and meta_path.exists():
                meta_path.unlink()


def check_all_data(output_folder: str | Path) -> None:
    """Check data files in the specified output folder.

    For each file in the list, if the file exists and contains only one data row,
    it is considered incomplete or invalid and is deleted.

    Parameters:
    ----------
    output_folder : str or Path
        Path to the directory containing the data files to be checked.

    """
    output_folder = Path(output_folder)

    file_to_delete = [
        "data_cna_hg19.seg",
        "data_cna_hg19.seg.fc.txt",
        "data_cna.txt",
        "data_mutations_extended.txt",
        "data_sv.txt"]

    for filename in file_to_delete:
        file_path = output_folder / filename
        if file_path.exists():
            try:
                file_df = pd.read_csv(file_path, sep="\t")
                if len(file_df) < 1:
                    file_path.unlink()
            except Exception as e:
                file_path.unlink()
