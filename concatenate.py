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

"""Concatenate MAF files in a study folder and post-process them for cBioPortal use.

This script identifies mutation annotation files (e.g., `.maf`), merges them into
a single `data_mutations_extended.txt` file, and applies basic filtering logic based
on OncoKB or other flags. It handles folder detection, header checking, and cleanup.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
from loguru import logger


def concatenate_files(file_list: list[str], output_file: str) -> None:
    """Concatenate multiple MAF files into a single output text file.

    The function reads each file in `file_list`, removes headers except for the
    first file, strips ".bam" substrings from lines, and writes the combined
    content into `output_file`.

    Parameters
    ----------
    file_list : list of str
        List of file paths to be concatenated.
    output_file : str
        Path to the output file where concatenated data will be written.

    Returns
    -------
    None

    """
    with Path(output_file).open("w") as out_file:
        for i, file_name in enumerate(file_list):
            with Path(file_name).open() as in_file:
                lines = in_file.readlines()
                if "Hugo_Symbol" not in lines[0]:
                    del lines[0]
                lines = [x.replace(".bam", "") for x in lines]
                if i > 0:
                    lines = lines[1:]
                out_file.write("".join(lines))

    if not Path(output_file).exists():
        logger.critical(f"Something went wrong while writing {output_file}.")

    logger.info(f"{len(file_list)} maf file(s) concatenated")


def get_files_by_ext(folder: str, ext: str) -> list[str]:
    """Collect all files with a specific extension from the folder and subfolders.

    Recursively searches through `folder` and returns a list of file paths that
    end with the given extension `ext`. Logs the count of found files or warns if none.

    Parameters
    ----------
    folder : str
        Path to the directory to search for files.
    ext : str
        File extension to filter by (e.g., "maf", "txt").

    Returns
    -------
    list of str
        List of file paths matching the extension.

    """
    file_list = []
    for root, _, files in os.walk(folder):
        matched_files = [
            str(Path(root) / file)
            for file in files
            if file.endswith(ext)]
        file_list.extend(matched_files)

    if len(file_list)==0:
        logger.warning(f"No files found with .{ext} extension in {folder}")
    else:
        logger.info(f"{len(file_list)} {ext} file(s) found")
    return file_list


def extract_maf_folder(filters: str, oncokb: str | bool) -> str:
    """Determine the appropriate MAF folder name based on filters and OncoKB flag.

    Chooses between "MAF_Onco_filtered", "MAF_filtered", "MAF_OncoKB", or "maf"
    depending on the presence of filters and the OncoKB parameter.

    Parameters
    ----------
    filters : str
        Filter flags controlling which folder to select.
    oncokb : bool or str
        Indicates whether OncoKB filtering is active.

    Returns
    -------
    str
        The name of the folder containing the MAF files to process.

    """
    if oncokb:
        folder="MAF_OncoKB"
        if "o" in filters:
            folder="MAF_Onco_filtered"
    else:
        folder="maf"
        if filters not in {"d", ""}:
            folder="MAF_filtered"

    return folder


def concatenate_main(
    filters: str, output_folder: str, ext: str, oncokb: str | bool) -> None:
    """Run the full MAF concatenation workflow including file checks and cleanup.

    Executes folder selection, file retrieval, concatenation, duplicates removal,
    and moves the final output file to the specified output folder. Logs each step.

    Parameters
    ----------
    filters : str
        Filter flags to determine which MAF folder to use.
    output_folder : str
        Path to the base output directory.
    ext : str
        File extension of MAF files to concatenate.
    oncokb : bool or str
        OncoKB flag affecting folder selection and filtering.

    Returns
    -------
    None

    Raises
    ------
    Exception
        If the output file path is a directory or has incorrect extension.

    """
    logger.info("Starting concatenate_main script:")
    logger.info(
    f"concatenate_main args [filters:{filters}, folder:{output_folder}, "
    f"extension:{ext}, oncoKB:{oncokb}]")

    folder = extract_maf_folder(filters, oncokb)
    input_folder = Path(output_folder) / folder
    output_file = input_folder / "data_mutations_extended.txt"

    msg = "Exiting from filter_clinvar script!"
    if output_file.is_dir():
        logger.critical(
            f"It seems that the inserted output_file '{output_file}' is not a file, "
            "but a folder! Check your '-o/--output_file' field")
        raise RuntimeError(msg)

    if not output_file.name.endswith("txt"):
        logger.critical(
            f"It seems that the inserted output_file '{output_file}' has "
            "the wrong extension! Output file must be have a .txt extension.")
        raise RuntimeError(msg)

    file_list = get_files_by_ext(input_folder, ext)
    concatenate_files(file_list, output_file)

    logger.info("Checking data_mutations_extended...")
    with output_file.open() as data_mut:
        all_data_mut = data_mut.readlines()
        if (len(all_data_mut) == 1):
            output_file.unlink()
            logger.warning("data_mutations_extended is empty. File removed.")

    if output_file.exists():
        data_mut = pd.read_csv(output_file, sep="\t", dtype=str)
        data_mut = data_mut.drop_duplicates(keep="last")

    if Path(output_file).exists():
        logger.info(f"Extracting data_mutations_extended from {input_folder} folder")

        src = Path(input_folder) / "data_mutations_extended.txt"
        dst = Path(output_folder) / "data_mutations_extended.txt"
        src.rename(dst)

    logger.success("Concatenate script completed!\n")
