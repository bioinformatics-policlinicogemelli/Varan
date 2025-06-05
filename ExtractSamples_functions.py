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

"""Module to extract genomic and clinical data from a study directory.

This script supports filtering and exporting:
- Clinical patient and sample information
- Copy number alterations (CNA) in different formats
- Mutation data
- Structural variation (SV) data

Intended for use in genomic data pre-processing pipelines.

"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger


def extract_clinical_samples(file_path: str,
                             sample_ids: list,
                             output_folder: str) -> None:
    """Extract specific clinical information for samples.

    This function reads a file from the given `file_path`, extracts rows
    corresponding to the provided list of `sample_ids`, and saves
    the extracted data, along with the header information, to a new file named
    'data_clinical_sample.txt' in the specified `output_folder`. If some sample
    identifiers are not found in the DataFrame, a warning message is printed.

    Args:
        file_path (str): Path to the original data_clinical_sample.
        sample_ids (list): List of sample identifiers to extract from the file.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t")
    idx_sample=np.argwhere(file.to_numpy() == "SAMPLE_ID")[0][1]
    header = file.loc[0:3, :]
    file = pd.read_csv(file_path, sep="\t")
    extracted = file[file.iloc[:, idx_sample].astype(str).isin(sample_ids)]

    extracted = pd.concat([header,extracted])
    outpath=Path(output_folder) / "data_clinical_sample.txt"
    extracted.to_csv(outpath, index=False, sep="\t")

def extract_clinical_patient(oldpath: str,
                            sample_ids: list,
                            output_folder: str) -> None:
    """Extract clinical patient data corresponding to specified sample identifiers.

    This function reads the 'data_clinical_patient.txt' file from the original
    data_clinical_patient, along with the 'data_clinical_sample.txt' file that
    should be present in the same directory. It extracts patient data rows based
    on the given `sample_ids`, determined by matching the '#Patient Identifier'
    from 'data_clinical_sample.txt'. The extracted patient data, including the header
    information, is saved to a new tab-separated file named 'data_clinical_patient.txt'
    in the specified `output_folder`. If some sample identifiers are not found in
    the DataFrame, a warning message is printed.

    Args:
        oldpath (str): Path to the old study folder.
        sample_ids (list): List of sample identifiers to extract.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    file = pd.read_csv(
        (Path(oldpath) / "data_clinical_patient.txt"), sep="\t", header=None)
    sample = pd.read_csv(
        (Path(oldpath) / "data_clinical_sample.txt"), sep="\t")

    idx_sample=np.argwhere(sample.to_numpy() == "SAMPLE_ID")[0][1]
    idx_patient=np.argwhere(sample.to_numpy() == "PATIENT_ID")[0][1]

    patient_ids = list(
        sample[sample.iloc[:, idx_sample].astype(str).isin(sample_ids)]
        [sample.iloc[0,idx_patient]])

    header = file.loc[0:4,:]
    data = file.loc[4:,:]

    idx_patient=np.argwhere(file.to_numpy() == "PATIENT_ID")[0][1]
    extracted = data[data.iloc[:, idx_patient].astype(str).isin(patient_ids)]

    extracted = pd.concat([header, extracted])
    outpath=Path(output_folder) / "data_clinical_patient.txt"
    extracted.to_csv(outpath, index=False, sep="\t", header=False, na_rep="NaN")


def extract_cna_hg19(file_path: str,
                     sample_ids: list,
                     output_folder: str) -> None:
    """Extract specific samples from CNA in hg19 genome format.

    This function reads a tab-separated CNA data file specified by 'file_path',
    extracts the rows corresponding to the provided 'sample_ids', and saves the
    extracted data into a new file named 'data_cna_hg19.seg' in the 'output_folder'.

    Args:
        file_path (str): Path to the input CNA data file in tab-separated format.
        sample_ids (list): List of sample IDs to be extracted from the input file.
        output_folder (str): Path to the output directory.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t")
    extracted = file[file["ID"].astype(str).isin(sample_ids)]
    outpath=Path(output_folder) / "data_cna_hg19.seg"
    extracted.to_csv(outpath, index=False, sep="\t")


def extract_cna_hg19_fc(file_path: str,
                        sample_ids: list,
                        output_folder: str) ->None:
    """Extract specific samples from CNA in hg19 genome format.

    Args:
        file_path (str): Path to the input CNA data file.
        sample_ids (list): List of sample IDs to be extracted from the input file.
        output_folder (str): Path to the output directory.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t")
    extracted = file[file["ID"].astype(str).isin(sample_ids)]
    outpath=Path(output_folder) / "data_cna_hg19.seg.fc.txt"
    extracted.to_csv(outpath, index=False, sep="\t")


def extract_cna(file_path: str,
                sample_ids: str,
                output_folder: str) -> None:
    """Extract specific samples' CNA data from the original CNA  file.

    This function reads the original CNA data file specified by 'file_path',
    extracts the columns corresponding to the provided 'sample_ids', and saves
    the extracted data into a new file named 'data_cna.txt' in the 'output_folder'.

    Args:
        file_path (str): Path to the input CNA data file in tab-separated format.
        sample_ids (list): List of sample IDs to be extracted from the input file.
        output_folder (str): Path to the output directory.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t", index_col=0)
    columns_to_keep = [sample for sample in sample_ids if sample in file.columns]
    extracted = file.loc[:, columns_to_keep]
    extracted = extracted.loc[(extracted != 0).any(axis=1)]
    outpath=Path(output_folder) / "data_cna.txt"
    extracted.to_csv(outpath, index=True, sep="\t")


def extract_mutations(file_path: str,
                      sample_ids: str,
                      output_folder: str) -> None:
    """Extract specific samples' mutation data from the `data_mutations_extended.txt`.

    This function reads a tab-separated mutation data file extracts the rows
    corresponding to the provided 'sample_ids', and saves the extracted data
    into a new file named 'data_mutations_extended.txt' in the 'output_folder'.

    Args:
        file_path (str): Path to the input mutation data file in tab-separated format.
        sample_ids (list): List of sample IDs to be extracted from the input file.
        output_folder (str): Path to the output directory.

    Returns:
        None

    """
    file = pd.read_csv(file_path, sep="\t", dtype=str)
    extracted = file[file["Tumor_Sample_Barcode"].astype(str).isin(sample_ids)]
    outpath=Path(output_folder) / "data_mutations_extended.txt"
    extracted.to_csv(outpath, index=False, sep="\t")


def extract_sv(file_path: str,
               sample_ids: str,
               output_folder: str) -> None:
    """Extract specific samples' SV data from the original `data_sv.txt`.

    This function reads a text file specified by 'file_path', extracts the lines
    corresponding to the provided 'sample_ids', and saves the extracted data into
    a new file named 'data_sv.txt' in the 'output_folder'.

    Args:
        file_path (str): Path to the input SV data file.
        sample_ids (list): List of sample IDs to be extracted from the input file.
        output_folder (str): Path to the output directory.

    Returns:
        None

    """
    old_file=pd.read_csv(file_path,sep="\t")
    new_file=old_file[old_file["Sample_Id"].isin(sample_ids)]
    outpath=Path(output_folder) / "data_sv.txt"
    new_file.to_csv(outpath,sep="\t",index=False)


def check_sample_list(extract_path: str, oldpath: str) -> None:
    """Validate the sample list file and check consistency with clinical sample data.

    This function performs several checks to ensure the list of samples to be extracted.
    It checks if this list contains only a single column; Includes at least one sample
    present in the 'data_clinical_sample.txt' file; Warns if some of the samples in the
    list are not found in the clinical sample data; Warns if the list includes all
    samples from the original study.

    Args:
        extract_path (str): Path to the text file containing sample ID to extract.
                            The file should contain one sample ID per line.
        oldpath (str): Path to the directory containing the original study data,
                            including the 'data_clinical_sample.txt' file.

    Returns:
        None

    """
    with Path(extract_path).open() as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.critical(
                f"The file {extract_path} contains more than a column."
                f"It may not be in the correct format!")
            sys.exit()

    with Path(extract_path).open() as sample_list:
        all_samples_to_extract = {sample.strip() for sample in sample_list}

        sample = pd.read_csv(
            (Path(oldpath) / "data_clinical_sample.txt"), sep="\t", skiprows=4)
        old_samples = set(sample["SAMPLE_ID"])

        if not any(sample in old_samples for sample in all_samples_to_extract):
            logger.critical(
                "The sample(s) you are trying to extract are not"
                "present in data_clinical_sample.txt file! Please check again!")
            sys.exit()

        missing_samples = set(all_samples_to_extract) - set(old_samples)
        if missing_samples:
            logger.warning(
                "Some of the samples you are trying to extract"
                "may not be present in data_clinical_sample.txt file!")
            logger.warning(f"Missing sample(s): {', '.join(missing_samples)}")

        if len(set(old_samples) - set(all_samples_to_extract)) == 0:
            logger.warning(
                "It looks like you are extracting all the samples from"
                "the original study, but this might be redundant, as it"
                "replicates the existing one.")

def extract_all_data(oldpath: str, sample_ids: list, output: str) -> None:
    """Extract available data types for the provided sample IDs.

    Args:
        oldpath (str): Path to the original study folder.
        sample_ids (list): List of sample IDs to extract.
        output (str): Output folder path.

    Returns:
        None

    """
    file_extractors = [
        ("data_clinical_patient.txt", extract_clinical_patient, False),
        ("data_clinical_sample.txt", extract_clinical_samples, True),
        ("data_cna_hg19.seg", extract_cna_hg19, True),
        ("data_cna_hg19.seg.fc.txt", extract_cna_hg19_fc, True),
        ("data_cna.txt", extract_cna, True),
        ("data_mutations_extended.txt", extract_mutations, True),
        ("data_sv.txt", extract_sv, True),
    ]

    for filename, extractor, use_path in file_extractors:
        file_path = Path(oldpath) / filename
        if file_path.exists():
            if use_path:
                extractor(file_path, sample_ids, output)
            else:
                extractor(oldpath, sample_ids, output)
        else:
            logger.warning(f"{filename} not found in current folder. Skipping")
