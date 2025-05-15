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

import os
import sys

import numpy as np
import pandas as pd
from loguru import logger


def extract_clinical_samples(file_path, sample_ids, output_folder):
    """Extracts specific clinical information for samples from the original data_clinical_samples.txt and saves the result to a new file.

    This function reads a tab-separated CSV file from the given `file_path`, extracts rows
    corresponding to the provided list of `sample_ids`, and saves the extracted data, along with
    the header information, to a new tab-separated file named 'data_clinical_sample.txt' in the
    specified `output_folder`. If some sample identifiers are not found in the DataFrame, a warning
    message is printed.

    Args:
        file_path (str): Path to the original data_clinical_sample.
        sample_ids (list): List of sample identifiers to extract from the file.
        output_folder (str): Path to the output folder where the extracted data will be saved.

    Example:
      >>>  extract_clinical_samples('input_data.csv', ['sample1', 'sample2'], 'output_folder/')

    """
    file = pd.read_csv(file_path, sep="\t")
    idx_sample=np.argwhere(file.values == "SAMPLE_ID")[0][1]
    header = file.loc[0:3, :]
    file = pd.read_csv(file_path, sep="\t")
    extracted = file[file.iloc[:, idx_sample].astype(str).isin(sample_ids)]

    extracted = pd.concat([header,extracted])
    extracted.to_csv(os.path.join(output_folder, "data_clinical_sample.txt"), index=False, sep="\t")

def extract_clinical_patient(oldpath, sample_ids, output_folder):
    """Extracts clinical patient data corresponding to specified sample identifiers and saves the result to a new file.
    
    This function reads the 'data_clinical_patient.txt' file from the original data_clinical_patient, along with the 'data_clinical_sample.txt'
    file that should be present in the same directory. It extracts patient data rows based on the given `sample_ids`, determined by
    matching the '#Patient Identifier' from 'data_clinical_sample.txt'. The extracted patient data, including the header information,
    is saved to a new tab-separated file named 'data_clinical_patient.txt' in the specified `output_folder`. If some sample identifiers
    are not found in the DataFrame, a warning message is printed.

    Args:
        oldpath (str): Path to the folder containing the input 'data_clinical_patient.txt' and 'data_clinical_sample.txt' files.
        sample_ids (list): List of sample identifiers for which patient data will be extracted.
        output_folder (str): Path to the output folder where the extracted patient data will be saved.

    Example:
      >>>  extract_clinical_patient('input_folder/', ['sample1', 'sample2'], 'output_folder/')

    """
    file = pd.read_csv(os.path.join(oldpath, "data_clinical_patient.txt"), sep="\t", header=None)
    sample = pd.read_csv(os.path.join(oldpath, "data_clinical_sample.txt"), sep="\t")

    idx_sample=np.argwhere(sample.values == "SAMPLE_ID")[0][1]
    idx_patient=np.argwhere(sample.values == "PATIENT_ID")[0][1]

    patient_ids = list(sample[sample.iloc[:, idx_sample].astype(str).isin(sample_ids)][sample.iloc[0,idx_patient]])

    header = file.loc[0:4,:]
    data = file.loc[4:,:]

    idx_patient=np.argwhere(file.values == "PATIENT_ID")[0][1]
    extracted = data[data.iloc[:, idx_patient].astype(str).isin(patient_ids)]

    extracted = pd.concat([header, extracted])
    extracted.to_csv(os.path.join(output_folder, "data_clinical_patient.txt"), index=False, sep="\t", header=False, na_rep="NaN")


def extract_cna_hg19(file_path, sample_ids, output_folder):
    """Extracts specific samples from a copy number alteration (CNA) data file in hg19 genome format.
    
    This function reads a tab-separated CNA data file specified by 'file_path', extracts the rows
    corresponding to the provided 'sample_ids', and saves the extracted data into a new file named
    'data_cna_hg19.seg' in the 'output_folder'.

    Args:
        file_path (str): Path to the input CNA data file in tab-separated format.
        sample_ids (list): List of sample IDs (as strings) to be extracted from the input file.
        output_folder (str): Path to the directory where the extracted data file will be saved.

    Notes:
        - The input file is expected to have a column named "ID" containing sample IDs.
        - If any sample ID provided in 'sample_ids' is not found in the input file, a warning
          will be printed indicating that some sample names are not present in the DataFrame.
        - The extracted data is saved in the 'data_cna_hg19.seg' file in tab-separated format.

    Example:
        extract_cna_hg19('input_data.tsv', ['sample1', 'sample2'], 'output_folder/')
      >>>  Extracts data for 'sample1' and 'sample2' from 'input_data.tsv' and saves it in 'output_folder/data_cna_hg19.seg'.

    """
    file = pd.read_csv(file_path, sep="\t")
    extracted = file[file["ID"].astype(str).isin(sample_ids)]
    extracted.to_csv(os.path.join(output_folder, "data_cna_hg19.seg"), index=False, sep="\t")


def extract_cna_hg19_fc(file_path, sample_ids, output_folder):
    """Args:
        file_path (str): Path to the input CNA data file in tab-separated format.
        sample_ids (list): List of sample IDs (as strings) to be extracted from the input file.
        output_folder (str): Path to the directory where the extracted data file will be saved.

    Notes:
        - The input file is expected to have a column named "ID" containing sample IDs.
        - If any sample ID provided in 'sample_ids' is not found in the input file, a warning
          will be printed indicating that some sample names are not present in the DataFrame.
        - The extracted data is saved in the 'data_cna_hg19.seg.fc.txt' file in tab-separated format.

    Example:
        extract_cna_hg19_fc('input_data.tsv', ['sample1', 'sample2'], 'output_folder/')
      >>>  Extracts data for 'sample1' and 'sample2' from 'input_data.tsv' and saves it in 'output_folder/data_cna_hg19.seg.fc.txt'.

    """
    file = pd.read_csv(file_path, sep="\t")
    extracted = file[file["ID"].astype(str).isin(sample_ids)]
    extracted.to_csv(os.path.join(output_folder, "data_cna_hg19.seg.fc.txt"), index=False, sep="\t")


def extract_cna(file_path, sample_ids, output_folder):
    """Extracts specific samples' copy number alteration (CNA) data from the original CNA  file.

    This function reads the original tab-separated CNA data file specified by 'file_path', extracts the columns
    corresponding to the provided 'sample_ids', and saves the extracted data into a new file named
    'data_cna.txt' in the 'output_folder'.

    Args:
        file_path (str): Path to the input CNA data file in tab-separated format.
        sample_ids (list): List of sample IDs (column names) to be extracted from the input file.
        output_folder (str): Path to the directory where the extracted data file will be saved.

    Notes:
        - The input file is expected to have a header row with column names.
        - If any sample ID provided in 'sample_ids' is not found in the input file's columns, a warning
          will be printed indicating that some sample names are not present in the DataFrame.
        - The extracted data is saved in the 'data_cna.txt' file in tab-separated format.

    Example:
      >>> extract_cna('input_data.tsv', ['sample1', 'sample2'], 'output_folder/')

    """
    file = pd.read_csv(file_path, sep="\t", index_col=0)
    columns_to_keep = [sample for sample in sample_ids if sample in file.columns]
    extracted = file.loc[:, columns_to_keep]
    extracted = extracted.loc[(extracted != 0).any(axis=1)]
    extracted.to_csv(os.path.join(output_folder, "data_cna.txt"), index=True, sep="\t")


def extract_mutations(file_path, sample_ids, output_folder):
    """Extracts specific samples' mutation data from the original `data_mutations_extended.txt`
    
    This function reads a tab-separated mutation data file specified by 'file_path', extracts the rows
    corresponding to the provided 'sample_ids', and saves the extracted data into a new file named
    'data_mutations_extended.txt' in the 'output_folder'.
    
    Args:
        file_path (str): Path to the input mutation data file in tab-separated format.
        sample_ids (list): List of sample IDs (Tumor_Sample_Barcode values) to be extracted from the input file.
        output_folder (str): Path to the directory where the extracted data file will be saved.

    Notes:
        - The input file is expected to have a column named "Tumor_Sample_Barcode" containing sample IDs.
        - If any sample ID provided in 'sample_ids' is not found in the input file, a warning
          will be printed indicating that some sample names are not present in the DataFrame.
        - The extracted data is saved in the 'data_mutations_extended.txt' file in tab-separated format.

    Example:
       >>> extract_mutations('input_mutations.tsv', ['sample1', 'sample2'], 'output_folder/')

    """
    file = pd.read_csv(file_path, sep="\t", dtype=str)
    extracted = file[file["Tumor_Sample_Barcode"].astype(str).isin(sample_ids)]
    extracted.to_csv(os.path.join(output_folder, "data_mutations_extended.txt"), index=False, sep="\t")


def extract_sv(file_path, sample_ids, output_folder):
    """Extracts specific samples' structural variation (SV) data from the original `data_sv.txt`.
    
    This function reads a text file specified by 'file_path', extracts the lines
    corresponding to the provided 'sample_ids', and saves the extracted data into a new file named
    'data_sv.txt' in the 'output_folder'.

    Args:
        file_path (str): Path to the input SV data file.
        sample_ids (list): List of sample IDs or keywords to be extracted from the input file.
        output_folder (str): Path to the directory where the extracted data file will be saved.

    Notes:
        - Each line in the input file is assumed to represent an SV record.
        - If any sample ID or keyword provided in 'sample_ids' is found in a line, that line is extracted.
        - The extracted data is saved in the 'data_sv.txt' file.

    Example:
      >>>  extract_sv('input_sv.txt', ['sample1', 'sample2'], 'output_folder/')

    """
    old_file=pd.read_csv(file_path,sep="\t")
    new_file=old_file[old_file["Sample_Id"].isin(sample_ids)]
    new_file.to_csv(os.path.join(output_folder,"data_sv.txt"),sep="\t",index=False)


def check_sample_list(extract_path, oldpath):
    with open(extract_path) as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.critical(f"The file {extract_path} contains more than a column. It may not be in the correct format!")
            sys.exit()

    with open(extract_path) as sample_list:
        all_samples_to_extract = set(sample.strip() for sample in sample_list.readlines())

        sample = pd.read_csv(os.path.join(oldpath, "data_clinical_sample.txt"), sep="\t", skiprows=4)
        old_samples = set(sample["SAMPLE_ID"])

        if not any(sample in old_samples for sample in all_samples_to_extract):
            logger.critical("The sample(s) you are trying to extract are not present in data_clinical_sample.txt file! Please check again!")
            sys.exit()

        missing_samples = set(all_samples_to_extract) - set(old_samples)
        if missing_samples:
            logger.warning("Some of the samples you are trying to extract may not be present in data_clinical_sample.txt file!")
            logger.warning(f"Missing sample(s): {", ".join(missing_samples)}")

        if len(set(old_samples) - set(all_samples_to_extract)) == 0:
            logger.warning("It looks like you are extracting all the samples from the original study, but this might be redundant, as it replicates the existing one.")
