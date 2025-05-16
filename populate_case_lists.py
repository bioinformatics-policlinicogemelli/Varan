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

"""Module for generating case lists from SV,CNA and mutation data.

Each function reads the input data file from a project directory,
extracts sample IDs, and writes a formatted `cases_*.txt` file into
the output directory for use in cBioPortal.

Functions:
- populate_cases_sv
- populate_cases_cna
- populate_cases_sequenced

"""

from pathlib import Path

import pandas as pd
from loguru import logger


def populate_cases_sv(project_id: str, folder: str, cases_list_dir: str,
                      logger: logger) -> None:
    """Populate the cases_sv.txt file based on SV data.

    Parameters
    ----------
    project_id : str
        Unique identifier for the cancer study.
    folder : str
        Path to the directory containing the `data_sv.txt` file.
    cases_list_dir : str
        Output directory where the `cases_sv.txt` will be saved.
    logger : logging.Logger
        Logger object for error reporting.

    """
    try:
        data_sv = pd.read_csv(Path(folder) / "data_sv.txt", sep="\t")
    except pd.errors.EmptyDataError:
        logger.exception("data_sv.txt is empty, skipping this step!")
        return()

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

    with (Path(cases_list_dir) / "cases_sv.txt").open("w") as case_sv_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=case_sv_file)

    return()


def populate_cases_cna(project_id: str, folder: str, cases_list_dir: str,
                       logger: logger) -> None:
    """Populate the cases_cna.txt file based on CNA data.

    Parameters
    ----------
    project_id : str
        Unique identifier for the cancer study.
    folder : str
        Path to the directory containing the `data_cna.txt` file.
    cases_list_dir : str
        Output directory where the `cases_cna.txt` will be saved.
    logger : logging.Logger
        Logger object for error reporting.

    """
    try:
        data_cna = pd.read_csv(Path(folder) / "data_cna.txt", sep="\t")
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

    with (Path(cases_list_dir) / "cases_cna.txt").open("w") as case_cna_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=case_cna_file)


def populate_cases_sequenced(project_id: str, folder: str, cases_list_dir: str,
                             logger: logger) -> None:
    """Populate the cases_sequenced.txt file based on mutation data.

    Parameters
    ----------
    project_id : str
        Unique identifier for the cancer study.
    folder : str
        Path to the directory containing the `data_mutations_extended.txt` file.
    cases_list_dir : str
        Output directory where the `cases_sequenced.txt` will be saved.
    logger : logging.Logger
        Logger object for error reporting.

    """
    try:
        data_sequenced = pd.read_csv(Path(folder) / "data_mutations_extended.txt",
                                    sep="\t", low_memory=False)
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

    with (Path(cases_list_dir) / "cases_sequenced.txt").open("w") as case_sequence_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=case_sequence_file)
