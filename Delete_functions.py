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

"""Clean clinical and genomic data by removing specific sample entries.

This module provides utility functions to remove samples from clinical,
genomic, and case list files. It ensures consistency by updating all relevant
files after removing given sample IDs. Intended for use in bioinformatics
pipelines for study customization or data sanitization.

"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger


def delete_clinical_samples(
    file_path: str,
    sample_ids: list[str],
    output_folder: str) -> None:
    """Remove specific samples from clinical sample data.

    Args:
        file_path (str): Path to input clinical sample file (TSV format).
        sample_ids (list[str]): List of sample IDs to remove.
        output_folder (str): Directory to save the updated clinical sample file.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t")
    idx_sample = np.argwhere(file.to_numpy() == "SAMPLE_ID")[0][1]
    filtered = file[~file.iloc[:, idx_sample].astype(str).isin(sample_ids)]
    output_path = Path(output_folder) / "data_clinical_sample.txt"
    filtered.to_csv(
        output_path,
        index=False,
        sep="\t",
    )


def delete_clinical_patient(
    oldpath: str,
    sample_ids: list[str],
    output_folder: str) -> None:
    """Remove patient entries associated with specific sample IDs.

    Args:
        oldpath (str): Directory containing input patient and sample files.
        sample_ids (list[str]): List of sample IDs used for patient filtering.
        output_folder (str): Directory to save the updated clinical patient file.

    Returns:
        None

    """
    oldpath = Path(oldpath)
    file = pd.read_csv(oldpath / "data_clinical_patient.txt", sep="\t")
    sample = pd.read_csv(oldpath / "data_clinical_sample.txt", sep="\t")

    idx_sample = np.argwhere(sample.to_numpy() == "SAMPLE_ID")[0][1]
    idx_patient = np.argwhere(sample.to_numpy() == "PATIENT_ID")[0][1]

    patient_ids = list(
        sample[
            sample.iloc[:, idx_sample].astype(str).isin(sample_ids)
        ].iloc[:, idx_patient])
    sample_ids = list(
        sample[
            sample.iloc[:, idx_sample].astype(str).isin(sample_ids)
            ].iloc[:, idx_sample])

    filtered_sample = sample.iloc[:, idx_sample].astype(str).isin(patient_ids)
    if len(sample[filtered_sample]) > len(sample_ids):
        pzt_list=sample[sample.iloc[:, idx_patient].astype(str).isin(patient_ids)]
        pzt_dup = [
            pzt_list[
                pzt_list.duplicated(subset=sample.iloc[0, idx_sample])
                ].iloc[0,1]]

        for p_dup in pzt_dup:
            df_dup = pzt_list[pzt_list.iloc[:,1]==p_dup]
            if len(
                df_dup[
                    df_dup.iloc[:, 0].astype(str).isin(sample_ids)
                    ])<df_dup.shape[0]:
                patient_ids.remove(p_dup)

    filtered = file[~file.iloc[:, 0].astype(str).isin(patient_ids)]
    output_path = Path(output_folder) / "data_clinical_patient.txt"
    filtered.to_csv(output_path, index=False, sep="\t", na_rep="NaN")


def delete_cna_hg19(file_path: str, sample_ids: list[str], output_folder: str) -> None:
    """Remove copy number alteration entries (hg19 format) by sample ID.

    Args:
        file_path (str): Path to the input CNA hg19 file (TSV format).
        sample_ids (list[str]): List of sample IDs to exclude.
        output_folder (str): Directory to save the filtered CNA file.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t")
    filtered = file[~file["ID"].astype(str).isin(sample_ids)]
    output_path = Path(output_folder) / "data_cna_hg19.seg"
    filtered.to_csv(output_path, index=False, sep="\t")


def delete_cna_hg19_fc(
    file_path: str,
    sample_ids: list[str],
    output_folder: str) -> None:
    """Remove CNA entries (hg19.fc format) associated with specific sample IDs.

    Args:
        file_path (str): Path to the input CNA hg19.fc file (TSV format).
        sample_ids (list[str]): List of sample IDs to exclude.
        output_folder (str): Directory to save the updated file.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t")
    filtered = file[~file["ID"].astype(str).isin(sample_ids)]

    output_path = Path(output_folder) / "data_cna_hg19.seg.fc.txt"
    filtered.to_csv(output_path, index=False, sep="\t")


def delete_cna(file_path: str, sample_ids: list[str], output_folder: str) -> None:
    """Remove copy number alteration columns for given sample IDs.

    Args:
        file_path (str): Path to the input CNA matrix file (TSV format).
        sample_ids (list[str]): List of sample IDs whose columns to remove.
        output_folder (str): Directory to save the updated CNA matrix.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t", index_col=0)
    filtered = file.drop(columns=sample_ids, axis=1, errors="ignore")
    filtered = filtered.loc[(filtered != 0).any(axis=1)]
    output_path = Path(output_folder) / "data_cna.txt"
    filtered.to_csv(output_path, sep="\t")


def delete_mutations(file_path: str, sample_ids: list[str], output_folder: str) -> None:
    """Remove mutation entries for specific tumor sample barcodes.

    Args:
        file_path (str): Path to the mutation file (TSV format).
        sample_ids (list[str]): Tumor sample barcodes to remove.
        output_folder (str): Directory to save the filtered mutation file.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t", dtype=str)
    filtered = file[~file["Tumor_Sample_Barcode"].astype(str).isin(sample_ids)]
    output_path = Path(output_folder) / "data_mutations_extended.txt"
    filtered.to_csv(output_path, index=False, sep="\t")


def delete_sv(file_path: str, sample_ids: list[str], output_folder: str) -> None:
    """Remove structural variation records containing specific sample IDs.

    Args:
        file_path (str): Path to the structural variation file (TXT format).
        sample_ids (list[str]): List of sample IDs to filter out.
        output_folder (str): Directory to save the updated SV file.

    Returns:
        None

    """
    file_path = Path(file_path)
    output_file = Path(output_folder) / "data_sv.txt"

    with file_path.open() as old_file, output_file.open("w") as of:
        for line in old_file:
            if not any(word in line for word in sample_ids):
                list_split = line.split("\t\t")
                list_strip = [elem.strip() for elem in list_split]
                new_row = "\t".join(list_strip) + "\n"
                of.write(new_row)


def delete_caselist_cna(
    file_path: str,
    sample_ids: list[str],
    output_folder: str) -> None:
    """Update CNA case list by removing specified sample IDs.

    Args:
        file_path (str): Path to the case list file for CNA.
        sample_ids (list[str]): Sample IDs to be excluded.
        output_folder (str): Directory to save the new case list file.

    Returns:
        None

    """
    file_path = Path(file_path)
    with file_path.open() as file:

        for raw_line in file:
            line = raw_line.strip()
            if line.startswith("case_list_ids"):
                samples = line.split(":")[1].split("\t")
                samples=list(map(str.strip, samples))
                updated = [elem for elem in samples if elem not in sample_ids]

        output_path = Path(output_folder) / "cases_cna.txt"
        with output_path.open("w") as filtered:
            file.seek(0)
            for raw_line in file:
                if raw_line.startswith("case_list_description"):
                    n_old_samples = re.findall(r"\d+", raw_line)[0]
                    line = raw_line.replace(n_old_samples, str(len(updated)))

                if raw_line.startswith("case_list_ids"):
                    line = "case_list_ids:" + "\t".join(updated)
                filtered.write(line)


def delete_caselist_sequenced(
    file_path: str,
    sample_ids: list[str],
    output_folder: str) -> None:
    """Update sequenced case list by removing specified sample IDs.

    Args:
        file_path (str): Path to the sequenced case list file.
        sample_ids (list[str]): Sample IDs to exclude.
        output_folder (str): Directory to save the updated case list.

    Returns:
        None

    """
    file_path = Path(file_path)
    with file_path.open() as file:
        for raw_line in file:
            line = raw_line.strip()
            if line.startswith("case_list_ids"):
                samples = line.split(":")[1].split("\t")
                samples=list(map(str.strip, samples))
                updated = [elem for elem in samples if elem not in sample_ids]

        output_path = Path(output_folder) / "cases_sequenced.txt"
        with output_path.open("w") as filtered:
            file.seek(0)
            for raw_line in file:
                if raw_line.startswith("case_list_description"):
                    n_old_samples = re.findall(r"\d+", raw_line)[0]
                    line = raw_line.replace(n_old_samples, str(len(updated)))

                if raw_line.startswith("case_list_ids"):
                    line = "case_list_ids:" + "\t".join(updated)
                filtered.write(line)


def delete_caselist_sv(
    file_path: str,
    sample_ids: list[str],
    output_folder: str) -> None:
    """Update structural variation case list by removing sample IDs.

    Args:
        file_path (str): Path to the SV case list file.
        sample_ids (list[str]): Sample IDs to be excluded.
        output_folder (str): Directory to save the new case list file.

    Returns:
        None

    """
    file_path = Path(file_path)
    with file_path.open() as file:
        for raw_line in file:
            line = raw_line.strip()
            if line.startswith("case_list_ids"):
                samples = line.split(":")[1].split("\t")
                samples = list(map(str.strip, samples))
                updated = [elem for elem in samples if elem not in sample_ids]

        output_path = Path(output_folder) / "cases_sv.txt"
        with output_path.open("w") as filtered:
            file.seek(0)
            for raw_line in file:
                if raw_line.startswith("case_list_description"):
                    n_old_samples = re.findall(r"\d+", raw_line)[0]
                    line = raw_line.replace(n_old_samples, str(len(updated)))

                if raw_line.startswith("case_list_ids"):
                    line = "case_list_ids:" + "\t".join(updated)
                filtered.write(line)


def check_sample_list(remove_path: str, oldpath: str) -> None:
    """Validate that sample list is correctly formatted and exists in the dataset.

    Args:
        remove_path (str): Path to file listing sample IDs to remove.
        oldpath (str): Directory containing the clinical sample file.

    Returns:
        None

    Raises:
        SystemExit: If the sample list format is invalid or results in an empty study.

    """
    with Path(remove_path).open() as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.critical(
                f"The file {remove_path} contains more than a column. "
                "It may not be in the correct format!")
            sys.exit(1)

    with Path(remove_path).open() as sample_list:
        all_samples_to_remove = {sample.strip() for sample in sample_list}

        sample_path = Path(oldpath) / "data_clinical_sample.txt"
        sample = pd.read_csv(sample_path, sep="\t", skiprows=4)
        old_samples = set(sample["SAMPLE_ID"])

        if not any(sample in old_samples for sample in all_samples_to_remove):
            logger.critical("The sample(s) you are trying to remove are not present "
            "in data_clinical_sample.txt file! Please check again!")
            sys.exit(1)

        missing_samples = set(all_samples_to_remove) - set(old_samples)
        if missing_samples:
            logger.warning("Some of the samples you are trying to remove may not be "
            "present in data_clinical_sample.txt file!")
            logger.warning(f"Missing sample(s): {', '.join(missing_samples)}")

        if len(set(old_samples) - set(all_samples_to_remove)) == 0:
            logger.critical("It looks like you are removing all the samples from "
            "original study! Cannot create an empty study!")
            sys.exit(1)


def delete_all_data(oldpath: str, sample_ids: list, output: str) -> None:
    """Delete data types for the provided sample IDs.

    Args:
        oldpath (str): Path to the original study folder.
        sample_ids (list): List of sample IDs to delete.
        output (str): Output folder path.

    Returns:
        None

    """
    file_to_delete = [
        ("data_clinical_patient.txt", delete_clinical_patient, False),
        ("data_clinical_sample.txt", delete_clinical_samples, True),
        ("data_cna_hg19.seg", delete_cna_hg19, True),
        ("data_cna_hg19.seg.fc.txt", delete_cna_hg19_fc, True),
        ("data_cna.txt", delete_cna, True),
        ("data_mutations_extended.txt", delete_mutations, True),
        ("data_sv.txt", delete_sv, True),
    ]

    for filename, deleter, use_path in file_to_delete:
        file_path = Path(oldpath) / filename
        if file_path.exists():
            if use_path:
                deleter(file_path, sample_ids, output)
            else:
                deleter(oldpath, sample_ids, output)
        else:
            logger.warning(f"{filename} not found in current folder. Skipping")
