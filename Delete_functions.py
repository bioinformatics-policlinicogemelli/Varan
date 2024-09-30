import os
import pandas as pd
import re
import numpy as np
import loguru
from loguru import logger 

def delete_clinical_samples(file_path, sample_ids, output_folder):
    """
    Delete specific clinical samples from a data file(data_clinical_sample).

    This function reads a data file in tab-separated format, filters out
    the rows with specified sample IDs, and saves the filtered data to
    a new file in the specified output folder.

    Args:
        file_path (str): Path to the input data file.
        sample_ids (list): List of sample IDs to be deleted.
        output_folder (str): Path to the folder where the output file will be saved.
        
    """    
    file = pd.read_csv(file_path, sep="\t")
    idx_sample=np.argwhere(file.values == "SAMPLE_ID")[0][1]
    filtered = file[~file.iloc[:, idx_sample].astype(str).isin(sample_ids)]
    filtered.to_csv(os.path.join(output_folder, "data_clinical_sample.txt"), index=False, sep="\t")


def delete_clinical_patient(oldpath, sample_ids, output_folder):
    """
    Delete clinical patients based on associated sample IDs.

    This function reads two data files: one for patients and one for samples,
    filters out patients based on sample IDs, and saves the filtered patient data
    to a new file in the specified output folder.

    Args:
        oldpath (str): Path to the folder containing the input data files.
        sample_ids (list): List of sample IDs used for patient filtering.
        output_folder (str): Path to the folder where the output patient file will be saved.

    """
    file = pd.read_csv(os.path.join(oldpath, "data_clinical_patient.txt"), sep="\t")
    sample = pd.read_csv(os.path.join(oldpath, "data_clinical_sample.txt"), sep="\t")
    
    idx_sample=np.argwhere(sample.values == "SAMPLE_ID")[0][1]
    idx_patient=np.argwhere(sample.values == "PATIENT_ID")[0][1]
    
    patient_ids = list(sample[sample.iloc[:, idx_sample].astype(str).isin(sample_ids)].iloc[:, idx_patient])
    
    # clean file from un-find samples
    sample_ids=list(sample[sample.iloc[:, idx_sample].astype(str).isin(sample_ids)].iloc[:, idx_sample])
    
    if len(sample[sample.iloc[:, idx_sample].astype(str).isin(patient_ids)]) > len(sample_ids):
        pzt_list=sample[sample.iloc[:, idx_patient].astype(str).isin(patient_ids)]
        pzt_dup=[pzt_list[pzt_list.duplicated(subset=sample.iloc[0, idx_sample])].iloc[0,1]]
        
        for p_dup in pzt_dup:
            df_dup=pzt_list[pzt_list.iloc[:,1]==p_dup]
            if len(df_dup[df_dup.iloc[:, 0].astype(str).isin(sample_ids)])<df_dup.shape[0]:
                patient_ids.remove(p_dup)
        
    filtered = file[~file.iloc[:, 0].astype(str).isin(patient_ids)]
    filtered.to_csv(os.path.join(output_folder, "data_clinical_patient.txt"), index=False, sep="\t")
   
def delete_cna_hg19(file_path, sample_ids, output_folder):
    """
    Delete specific copy number alteration data from a file in hg19 format.

    This function reads a data file in tab-separated format, filters out the rows
    with specified IDs, and saves the filtered data to a new file in the specified output folder.

    Args:
        file_path (str): Path to the input data file.
        sample_ids (list): List of IDs to be deleted.
        output_folder (str): Path to the folder where the output file will be saved.

    """
    file = pd.read_csv(file_path, sep="\t")
    filtered = file[~file["ID"].astype(str).isin(sample_ids)]
    filtered.to_csv(os.path.join(output_folder, "data_cna_hg19.seg"), index=False, sep="\t")


def delete_cna(file_path, sample_ids, output_folder):
    """
    Delete copy number alteration data associated with specific samples.

    This function reads a data file in tab-separated format, drops columns corresponding to
    the specified sample IDs, and saves the modified data to a new file in the specified output folder.

    Args:
        file_path (str): Path to the input data file.
        sample_ids (list): List of sample IDs for which data will be deleted.
        output_folder (str): Path to the folder where the output file will be saved.

    """
    file = pd.read_csv(file_path, sep="\t")
    filtered = file.drop(columns=sample_ids, axis=1, errors="ignore")
    filtered.to_csv(os.path.join(output_folder, "data_cna.txt"), index=False, sep="\t")
        
def delete_mutations(file_path, sample_ids, output_folder):
    """
    Delete specific mutation data associated with given sample IDs.

    This function reads a data file in tab-separated format, filters out rows
    with specified tumor sample barcodes, and saves the filtered data to a new file.

    Args:
        file_path (str): Path to the input data file.
        sample_ids (list): List of tumor sample barcodes to be deleted.
        output_folder (str): Path to the folder where the output file will be saved.

    """
    file = pd.read_csv(file_path, sep="\t")
    filtered = file[~file["Tumor_Sample_Barcode"].astype(str).isin(sample_ids)]
    filtered.to_csv(os.path.join(output_folder, "data_mutations_extended.txt"), index=False, sep="\t")


def delete_sv(file_path, sample_ids, output_folder):
    """
    Delete structural variation data associated with specified sample IDs.

    This function reads a data file line by line, filters out lines containing
    any of the specified sample IDs, and saves the filtered data to a new file.

    Args:
        file_path (str): Path to the input data file.
        sample_ids (list): List of sample IDs for which structural variation data will be deleted.
        output_folder (str): Path to the folder where the output file will be saved.

    """
    with open(file_path, "r") as old_file:
        with open(os.path.join(output_folder, "data_sv.txt"), "w") as of:
            for line in old_file:
                if not any(word in line for word in sample_ids):
                    list_split = line.split("\t\t")
                    list_strip = [elem.strip() for elem in list_split]
                    new_row = "\t".join(list_strip) + '\n'
                    of.write(new_row)


def delete_caselist_cna(file_path, sample_ids, output_folder):
    """
    Delete case list IDs associated with specific sample IDs for copy number alteration data.

    This function reads a case list file, identifies case list IDs, removes the specified sample IDs,
    updates the case list description, and saves the modified case list to a new file.

    Args:
        file_path (str): Path to the input case list file.
        sample_ids (list): List of sample IDs to be removed from the case list.
        output_folder (str): Path to the folder where the output case list file will be saved.

   
    """
    with open(file_path, "r") as file:
        
        for line in file:
            line = line.strip()
            if line.startswith("case_list_ids"):
                samples = line.split(":")[1].split("\t")
                samples=list(map(str.strip, samples))
                updated = [elem for elem in samples if elem not in sample_ids]

        with open(os.path.join(output_folder, "cases_cna.txt"), "w") as filtered:
            file.seek(0)
            for line in file:
                if line.startswith("case_list_description"):
                    n_old_samples = re.findall(r'\d+', line)[0]
                    line = line.replace(n_old_samples, str(len(updated)))

                if line.startswith("case_list_ids"):
                    line = "case_list_ids:" + "\t".join(updated)
                filtered.write(line)


def delete_caselist_sequenced(file_path, sample_ids, output_folder):
    """
    Delete case list IDs associated with specific sample IDs for sequenced data.

    This function reads a case list file, identifies case list IDs, removes the specified sample IDs,
    updates the case list description, and saves the modified case list to a new file.

    Args:
        file_path (str): Path to the input case list file.
        sample_ids (list): List of sample IDs to be removed from the case list.
        output_folder (str): Path to the folder where the output case list file will be saved.

    
    """
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("case_list_ids"):
                samples = line.split(":")[1].split("\t")
                samples=list(map(str.strip, samples))
                updated = [elem for elem in samples if elem not in sample_ids]

        with open(os.path.join(output_folder, "cases_sequenced.txt"), "w") as filtered:
            file.seek(0)
            for line in file:
                if line.startswith("case_list_description"):
                    n_old_samples = re.findall(r'\d+', line)[0]
                    line = line.replace(n_old_samples, str(len(updated)))

                if line.startswith("case_list_ids"):
                    line = "case_list_ids:" + "\t".join(updated)
                filtered.write(line)

    
    
def delete_caselist_sv(file_path, sample_ids, output_folder):
    """
    Delete case list IDs associated with specific sample IDs for structural variation data.

    This function reads a case list file, identifies case list IDs, removes the specified sample IDs,
    updates the case list description, and saves the modified case list to a new file.

    Args:
        file_path (str): Path to the input case list file.
        sample_ids (list): List of sample IDs to be removed from the case list.
        output_folder (str): Path to the folder where the output case list file will be saved.

    
    """
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("case_list_ids"):
                samples = line.split(":")[1].split("\t")
                samples=list(map(str.strip, samples))
                updated = [elem for elem in samples if elem not in sample_ids]

        with open(os.path.join(output_folder, "cases_sv.txt"), "w") as filtered:
            file.seek(0)
            for line in file:
                if line.startswith("case_list_description"):
                    n_old_samples = re.findall(r'\d+', line)[0]
                    line = line.replace(n_old_samples, str(len(updated)))

                if line.startswith("case_list_ids"):
                    line = "case_list_ids:" + "\t".join(updated)
                filtered.write(line)
