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
import pandas as pd
import argparse
from configparser import ConfigParser

def populate_cases_sv(project_id, folder, cases_list_dir, logger):
    """
        Function to populate cases_sv file
    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        cases_list_dir : path of case_list output dir
    """
    try:
        data_sv = pd.read_csv(os.path.join(folder, "data_sv.txt"), sep="\t")
    except pd.errors.EmptyDataError:
        logger.error("data_sv.txt is empty, skipping this step!")
        return
    nsamples = len(data_sv.Sample_Id.unique())
    sample_ids = list(data_sv.Sample_Id.unique())
    
    stable_id = project_id + "_sv"
    case_list_name = "Samples with SV data"
    case_list_category = "all_cases_with_sv_data"
    case_list_description = "Samples with SV data (" + str(nsamples) + " sample(s))"
    case_list_ids = "\t".join(sample_ids)

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "stable_id": stable_id,
        "case_list_category": case_list_category,
        "case_list_name": case_list_name,
        "case_list_description": case_list_description,
        "case_list_ids": case_list_ids,
    }
    
    case_sv_file = open(os.path.join(cases_list_dir, "cases_sv.txt"), "w")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=case_sv_file)
    case_sv_file.close()

def populate_cases_cna(project_id, folder, cases_list_dir, logger):
    """
        Function to populate cases_cna file
    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        cases_list_dir : path of case_list output dir
    """
    
    try:
        data_cna = pd.read_csv(os.path.join(folder, "data_cna.txt"), sep="\t")
    except pd.errors.EmptyDataError:
        logger.error("data_cna.txt is empty, skipping this step!")
        return
    
    nsamples = len(data_cna.columns)-1
    sample_ids = list(data_cna.columns)[1:]
    
    stable_id = project_id + "_cna"

    case_list_category = "all_cases_with_cna_data"
    case_list_name = "Samples with CNA data"
    case_list_description = "Samples with CNA data (" + str(nsamples) + " sample(s))"
    case_list_ids = "\t".join(sample_ids)

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "stable_id": stable_id,
        "case_list_category": case_list_category,
        "case_list_name": case_list_name,
        "case_list_description": case_list_description,
        "case_list_ids": case_list_ids,
    }

    case_cna_file = open(os.path.join(cases_list_dir, "cases_cna.txt"), "w")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=case_cna_file)
    case_cna_file.close()



def populate_cases_sequenced(project_id, folder, cases_list_dir, logger):
    """
        Function to populate cases_sequenced file
    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        cases_list_dir : path of case_list output dir
    """

    try:
        data_sequenced = pd.read_csv(os.path.join(folder, "data_mutations_extended.txt"), sep="\t", low_memory=False)
    except pd.errors.EmptyDataError:
        logger.error("data_mutations_extended.txt is empty, skipping this step!")
        return
    nsamples = len(data_sequenced["Tumor_Sample_Barcode"].unique())
    sample_ids = list(data_sequenced["Tumor_Sample_Barcode"].unique())

    stable_id = project_id + "_sequenced"

    case_list_category = "all_cases_with_mutation_data"
    case_list_name = "Sequenced Tumors"
    case_list_description = "All sequenced samples (" + str(nsamples) + " sample(s))"
    case_list_ids = "\t".join(sample_ids)

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "stable_id": stable_id,
        "case_list_category": case_list_category,
        "case_list_name": case_list_name,
        "case_list_description": case_list_description,
        "case_list_ids": case_list_ids,
    }

    case_sequenced_file = open(os.path.join(cases_list_dir, "cases_sequenced.txt"), "w")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=case_sequenced_file)
    case_sequenced_file.close()