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

"""Module: update_functions.

Provides functions to update and synchronize study data files
across versioned folders in a cBioPortal-style workflow.

Key functionalities include:
  - Merging clinical sample and patient files.
  - Updating copy number alteration (CNA) and fold-change CNA data.
  - Synchronizing mutation and structural variation data.
  - Managing case list updates (CNA, sequenced, SV).

Each function reads the old and new versions of a file, merges or copies
as appropriate, and writes the result to the specified output folder.

"""

import re
import shutil
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from loguru import logger

from versioning import (
    create_newest_version_folder,
    get_version_list,
)


def update_clinical_samples(oldfile_path: str,
                            newfile_path: str,
                            output_folder: str) -> None:
    """Update data_clinical_sample.txt file.

    This function reads the original txt file from the given 'oldfile_path',
    insert new rows with the sample info founded inside the new txt file from
    the given 'newfile_path' and save the updated file.

    Args:
        oldfile_path (str): Path to the original data_clinical_sample.
        newfile_path (str): Path to the new data_clinical_sample.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    #header
    old_head = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=None, nrows=5)
    old_head.columns = old_head.iloc[4].tolist()
    old_head=old_head.set_index(old_head.columns[0])

    new_head = pd.read_csv(newfile_path, sep="\t", dtype = str, header=None, nrows=5)
    new_head.columns = new_head.iloc[4].tolist()
    new_head=new_head.set_index(new_head.columns[0])

    only_common_head = np.intersect1d(new_head.columns, old_head.columns)

    old_updated_head = old_head.drop(list(only_common_head), axis=1)
    old_updated_head.index = new_head.index
    new_head = pd.concat([new_head, old_updated_head], axis=1).fillna(value=np.nan)
    final_head = new_head.reset_index(drop=False)

    #body
    old_body = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=4)
    new_body = pd.read_csv(newfile_path, sep="\t", dtype = str, header=4)

    old_body_unique = old_body.drop(list(only_common_head), axis=1)
    final_body = old_body_unique.merge(new_body, how="outer", on="SAMPLE_ID")

    old_body=old_body.set_index(old_body.columns[0])
    new_body=new_body.set_index(new_body.columns[0])
    common_samples = old_body.index.intersection(new_body.index)
    old_body = old_body.drop(index=common_samples)
    final_body=final_body.set_index(final_body.columns[0])
    final_body.update(old_body, overwrite=True, filter_func=None, errors="ignore")
    final_body = final_body.reset_index(drop=False)

    #final
    final_patient = pd.concat([final_head, final_body], axis=0)
    outpath=Path(output_folder) / "data_clinical_sample.txt"
    final_patient.to_csv(outpath, header=False, index=False, sep="\t", na_rep="NaN")
    logger.info("data_clinical_sample.txt updated!")


def update_clinical_patient(oldfile_path: str,
                            newfile_path: str,
                            output_folder: str) -> None:
    """Update data_clinical_patient.txt file.

    This function reads the originaltxt file from the 'oldfile_path',
    insert new rows with the patients info founded inside the new txt file
    from the given 'newfile_path' and save the updated 'data_clinical_patient.txt'

    Args:
        oldfile_path (str): Path to the original data_clinical_patient.
        newfile_path (str): Path to the new data_clinical_patient.
        output_folder (str): Path to the output folder where the file will be saved.

    Returns:
        None

    """
    #header
    old_head = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=None, nrows=5)
    old_head.columns = old_head.iloc[4].tolist()
    old_head=old_head.set_index(old_head.columns[0])

    new_head = pd.read_csv(newfile_path, sep="\t", dtype = str, header=None, nrows=5)
    new_head.columns = new_head.iloc[4].tolist()
    new_head=new_head.set_index(new_head.columns[0])

    only_common_head = np.intersect1d(new_head.columns, old_head.columns)

    old_updated_head = old_head.drop(list(only_common_head), axis=1)
    old_updated_head.index = new_head.index
    new_head = pd.concat([new_head, old_updated_head], axis=1).fillna(value=np.nan)
    final_head = new_head.reset_index(drop=False)

    #body
    old_body = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=4)
    new_body = pd.read_csv(newfile_path, sep="\t", dtype = str, header=4)

    old_body_unique = old_body.drop(list(only_common_head), axis=1)
    final_body = old_body_unique.merge(new_body, how="outer", on="PATIENT_ID")

    old_body=old_body.set_index(old_body.columns[0])
    new_body=new_body.set_index(new_body.columns[0])
    common_samples = old_body.index.intersection(new_body.index)
    old_body = old_body.drop(index=common_samples)
    final_body=final_body.set_index(final_body.columns[0])
    final_body.update(old_body, overwrite=True, filter_func=None, errors="ignore")
    final_body = final_body.reset_index(drop=False)

    #final
    final_patient = pd.concat([final_head, final_body], axis=0)
    outpath=Path(output_folder) / "data_clinical_patient.txt"
    final_patient.to_csv(outpath, header=False, index=False, sep="\t", na_rep="NaN")
    logger.info("data_clinical_patient.txt updated!")


def update_cna_hg19(oldfile_path: str, newfile_path: str, output_folder: str) -> None:
    """Update data_cna_hg19.seg.

    This function reads the original CNA data from the given 'oldfile_path',
    insert new rows with the sample CNA data founded inside the new file from
    the given 'newfile_path' and save the updated file.

    Args:
        oldfile_path (str): Path to the original CNA data file.
        newfile_path (str): Path to the new CNA data file.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t")
    new = pd.read_csv(newfile_path, sep="\t")

    updated = pd.concat([old, new])
    updated = updated.drop_duplicates(
        subset=["ID","chrom","loc.start","loc.end"], keep="last")
    outpath=Path(output_folder) / "data_cna_hg19.seg"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_cna_hg19.seg updated!")


def update_cna_hg19_fc(oldfile_path: str,
                       newfile_path: str,
                       output_folder: str) -> None:
    """Update fold-change copy number alteration (CNA) data in hg19 format.

    Reads the original and new fold-change CNA segment files (TSV), merges rows,
    removes duplicates based on segment coordinates and gene, and writes the
    consolidated data to `data_cna_hg19.seg.fc.txt` in the output folder.

    Args:
        oldfile_path (str): Path to the original CNA fold-change file.
        newfile_path (str): Path to the new CNA fold-change file.
        output_folder (str): Directory where the updated file will be saved.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t")
    new = pd.read_csv(newfile_path, sep="\t")

    updated = pd.concat([old, new])
    updated = updated.drop_duplicates(
        subset=["ID","chrom","loc.start","loc.end", "gene"], keep="last")
    outpath=Path(output_folder) / "data_cna_hg19.seg.fc.txt"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_cna_hg19.seg.fc.txt updated!")


def update_cna(oldfile_path: str,
               newfile_path: str,
               output_folder: str) -> None:
    """Update sample inside a copy number alteration (CNA) data file.

    This function reads the original tab separated CNA data,
    insert new rows with the sample CNA data founded inside the new file
    from the given 'newfile_path' and save the updated file.

    Args:
        oldfile_path (str): Path to the original CNA data file.
        newfile_path (str): Path to the new CNA data file.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", index_col=0)
    new = pd.read_csv(newfile_path, sep="\t", index_col=0)

    sample_old = old.columns
    sample_new = new.columns
    to_remove = sample_old.intersection(sample_new)

    updated = pd.concat(
        [old.drop(columns=to_remove), new]
        , axis=1).replace(np.nan, 0).astype(int)
    outpath=Path(output_folder) / "data_cna.txt"
    updated.to_csv(outpath, index=True, sep="\t")
    logger.info("data_cna.txt updated!")


def update_mutations(oldfile_path: str,
                     newfile_path: str,
                     output_folder: str) -> None:
    """Update samples' mutation data inside data_mutations_extended.txt file.

    This function reads the original tab separated version txt file,
    insert new rows with the samples' mutation data founded inside the new txt file
    and save the updated file.

    Args:
        oldfile_path (str): Path to the original data_mutations_extended.
        newfile_path (str): Path to the new data_mutations_extended.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype=str)
    new = pd.read_csv(newfile_path, sep="\t", dtype=str)

    updated = pd.concat([old, new])
    updated=updated.drop_duplicates(
        subset=updated.loc[:, "Hugo_Symbol":"n_AF"].columns, keep="last")
    outpath=Path(output_folder) / "data_mutations_extended.txt"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_mutation_extended.txt updated!")


def update_sv(oldfile_path: str,
              newfile_path: str,
              output_folder: str) -> None:
    """Update samples' structural variation (SV) data inside data_sv.txt file.

    This function reads the original tab separated version txt file,
    insert new rows with the samples' SV data founded inside the new txt file
    and save the updated file named 'data_sv.txt' in the specified 'output_folder'.

    Args:
        oldfile_path (str): Path to the original data_sv.
        newfile_path (str): Path to the new data_sv.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    df_old = pd.read_csv(oldfile_path, sep="\t")
    df_new = pd.read_csv(newfile_path, sep="\t")
    merged_df = pd.concat(
        [df_old, df_new], axis=0, join="outer", ignore_index=True)
    merged_df=merged_df.drop_duplicates(
        subset=["Sample_Id",
                "Site1_Hugo_Symbol",
                "Site2_Hugo_Symbol",
                "SV_Status", "Class"], keep="last")
    output_file = Path(output_folder) / "data_sv.txt"
    merged_df.to_csv(output_file, sep="\t", index=False)

    logger.info("data_sv.txt updated!")


def update_caselist_cna(oldfile_path: str,
                        newfile_path: str,
                        output_folder: str) -> None:
    """Update cases' CNA data inside cases_cna.txt file.

    This function reads the original tab separated version txt file,
    insert new rows with the cases' CNA data founded inside the new txt
    file and save the updated file.

    Args:
        oldfile_path (str): Path to the original cases_cna.
        newfile_path (str): Path to the new cases_cna.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    with Path(oldfile_path).open() as old:
        with Path(newfile_path).open() as new:
            for raw_line in new:
                line = raw_line.strip()
                if line.startswith("case_list_ids"):
                    new_samples = line.split(":")[1]
                    len_new_sample = len(new_samples.split("\t"))

        outpath=Path(output_folder) / "cases_cna.txt"
        with outpath.open("w") as updated:
            for raw_line in old:
                if raw_line.startswith("case_list_description"):
                    n_old_samples = re.findall(r"\d+", raw_line)[0]
                    line = raw_line.replace(
                        n_old_samples, str(int(n_old_samples) + len_new_sample))
                if raw_line.startswith("case_list_ids"):
                    new_samples_filtered =[
                        sample for sample in new_samples.split("\t")
                        if sample not in raw_line]
                    line = "\t".join([raw_line, "\t".join(new_samples_filtered)])
                updated.write(line)


def update_caselist_sequenced(oldfile_path: str,
                              newfile_path: str,
                              output_folder: str) -> None:
    """Update cases' sequenced data inside cases_sequenced.txt file.

    This function reads the originalfile from the given 'oldfile_path',
    insert new rows founded inside the new txt file from the given 'newfile_path'
    and save the updated file named 'cases_sequenced.txt.

    Args:
        oldfile_path (str): Path to the original cases_sequenced.
        newfile_path (str): Path to the new cases_sequenced.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    with Path(oldfile_path).open() as old:
        with Path(newfile_path).open() as new:
            for raw_line in new:
                line = raw_line.strip()
                if line.startswith("case_list_ids"):
                    new_samples = line.split(":")[1]
                    len_new_sample = len(new_samples.split("\t"))
        outpath=Path(output_folder) / "cases_sequenced.txt"
        with outpath.open("w") as updated:
            for raw_line in old:
                if raw_line.startswith("case_list_description"):
                    n_old_samples = re.findall(r"\d+", raw_line)[0]
                    line = raw_line.replace(
                        n_old_samples, str(int(n_old_samples) + len_new_sample))
                if raw_line.startswith("case_list_ids"):
                    new_samples_filtered = [
                        sample for sample in new_samples.split("\t")
                        if sample not in raw_line]
                    line = "\t".join([raw_line, "\t".join(new_samples_filtered)])
                updated.write(line)


def update_caselist_sv(oldfile_path: str,
                       newfile_path: str,
                       output_folder: str) -> None:
    """Update cases' structural variation (SV) data inside cases_sv.txt file.

    This function reads the original file from the given 'oldfile_path',
    insert new rows with the cases' SV data founded inside 'newfile_path'
    and save the updated file named 'cases_sv.txt'.

    Args:
        oldfile_path (str): Path to the original cases_sv.
        newfile_path (str): Path to the new cases_sv.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    with Path(oldfile_path).open() as old:
        with Path(newfile_path).open() as new:
            for raw_line in new:
                line = raw_line.strip()
                if line.startswith("case_list_ids"):
                    new_samples = line.split(":")[1]
                    len_new_sample = len(new_samples.split("\t"))

        outpath=Path(output_folder) / "cases_sv.txt"
        with outpath.open("w") as updated:
            for raw_line in old:
                if raw_line.startswith("case_list_description"):
                    n_old_samples = re.findall(r"\d+", raw_line)[0]
                    line = raw_line.replace(
                        n_old_samples, str(int(n_old_samples) + len_new_sample))
                if raw_line.startswith("case_list_ids"):
                    new_samples_filtered = [
                        sample for sample in new_samples.split("\t")
                        if sample not in raw_line]
                    line = "\t".join([raw_line, "\t".join(new_samples_filtered)])
                updated.write(line)


def check_files(oldpath: str,
                newpath: str,
                output: str,
                file_name: str) -> None:
    """Dispatch and/or copy study data files between versions.

    Checks for the presence of a given file in both the old and new paths.
    - If present in both: calls the appropriate `update_*` function to merge.
    - If only in one: copies it to the output folder.
    - If in neither: logs a warning and skips.

    Args:
        oldpath (str): Path to the previous version folder.
        newpath (str): Path to the incoming data folder.
        output (str): Path to the target output folder.
        file_name (str): Name of the file to process, e.g.
            "data_clinical_sample.txt", "data_cna.txt", etc.

    Returns:
        None

    """
    file_updaters: dict[str, Callable[[Path, Path, Path], None]] = {
    "data_clinical_sample.txt": update_clinical_samples,
    "data_clinical_patient.txt": update_clinical_patient,
    "data_cna_hg19.seg": update_cna_hg19,
    "data_cna_hg19.seg.fc.txt": update_cna_hg19_fc,
    "data_cna.txt": update_cna,
    "data_mutations_extended.txt": update_mutations,
    "data_sv.txt": update_sv}

    o_data = Path(oldpath) / file_name
    n_data = Path(newpath) / file_name
    dest = Path(output) / file_name

    if o_data.exists() and n_data.exists():
        updater = file_updaters.get(file_name)
        if updater:
            updater(o_data, n_data, output)
        else:
            logger.warning(f"No updater function for '{file_name}'. Skipping merge.")
    elif o_data.exists():
        logger.warning(f"'{file_name}' only in old; copying forward.")
        shutil.copy(o_data, dest)
    elif n_data.exists():
        logger.warning(f"'{file_name}' only in new; copying forward.")
        shutil.copy(n_data, dest)
    else:
        logger.warning(f"'{file_name}' not found in either folder. Skipping.")


def check_files_cases(oldpath: str,
                      newpath: str,
                      output_caseslists: str,
                      file_name: str) -> None:
    """Dispatch and/or copy case list files between study versions.

    Checks for a given case list file in both old and new paths:
      - If present in both: calls `update_caselist_*` function to merge.
      - If only in one: copies it to the output 'case_lists' folder.
      - If in neither: logs a warning and skips.

    Args:
        oldpath (str): Path to the previous version folder.
        newpath (str): Path to the incoming data folder.
        output_caseslists (str): Path to the target 'case_lists' output folder.
        file_name (str): Name of the case list file to process,
            e.g. "cases_cna.txt", "cases_sequenced.txt", or "cases_sv.txt".

    Returns:
        None

    """
    o_data = Path(oldpath) / "case_lists" / file_name
    n_data = Path(newpath) / "case_lists" / file_name
    if o_data.exists() and n_data.exists():
        if file_name == "cases_cna.txt":
            update_caselist_cna(o_data,n_data,output_caseslists)
        elif file_name == "cases_sequenced.txt":
            update_caselist_sequenced(o_data,n_data,output_caseslists)
        elif file_name == "cases_sv.txt":
            update_caselist_sv(o_data,n_data,output_caseslists)
    elif o_data.exists() and not n_data.exists():
        logger.warning(
            f"{file_name} was not found in path 2 case_lists folder. "
                       "The file is being copied from path 1.")
        shutil.copy(o_data, output_caseslists)
    elif not o_data.exists() and n_data.exists():
        logger.warning(
            f"{file_name} was not found in path 1 case_lists folder. "
            "The file is being copied from path 2.")
        shutil.copy(n_data, output_caseslists)
    else:
        logger.warning(f"{file_name} not found in 'case_lists' folders. Skipping")

def copy_logo(oldpath: str, output: str) -> None:
    """Copy 'logo_VARAN.png' old folder to the new output folder.

    Args:
        oldpath (str): path to old folder containing the 'img/logo_VARAN.png' file.
        output (str): path to output folder where 'logo_VARAN.png' will be copied.

    Returns:
        None: This function performs file operations and doesn't return a value.

    Raises:
        FileNotFoundError: If the 'logo_VARAN.png' file is not found.

    """
    img_path = Path(oldpath) / "img" / "logo_VARAN.png"
    if Path(img_path).exists():
        img_output_dir = Path(output) / "img"
        img_output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(img_path, Path(img_output_dir) / "logo_VARAN.png")

def safe_check_file(oldpath: Path, newpath: Path, output: Path, file: str) -> None:
    """Safely checks a file for consistency between old and new study folders.

    This function wraps `check_files` in a try-except block.

    Args:
        oldpath (Path): Path to the old version of the study folder.
        newpath (Path): Path to the new version of the study folder.
        output (Path): Path to the output (updated) study folder.
        file (str): Name of the file to check (e.g., "data_clinical_sample.txt").

    Raises:
        IndexError: Raised when a parsing error is detected in the file structure.

    """
    try:
        check_files(oldpath, newpath, output, file)
    except pd.errors.ParserError as e:
        line_match = re.search(r"line (\d+)", str(e))
        line_number = line_match.group(1) if line_match else "unknown"
        logger.critical(
            f"Wrong column number in line {line_number} of {file} file")
        msg = "Exiting from Update script!"
        raise IndexError(msg) from e

def copy_metadata_files(oldpath: Path, output: Path) -> None:
    """Copy all metadata files from old study folder to new.

    This function is used to transfer configuration or metadata necessary
    for the updated study version.

    Args:
        oldpath (Path): Path to previous version of the study folder.
        output (Path): Path to output (new version) study folder.

    Returns:
        None

    """
    meta_files = list(oldpath.glob("*meta*"))
    if meta_files:
        for file in meta_files:
            shutil.copy(file, output)
    else:
        logger.warning("No meta files found!")

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
