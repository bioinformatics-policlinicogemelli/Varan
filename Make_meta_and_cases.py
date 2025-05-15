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
from configparser import ConfigParser
from datetime import datetime

from loguru import logger

from populate_case_lists import (
    populate_cases_cna,
    populate_cases_sequenced,
    populate_cases_sv,
)
from versioning import extract_version_str


def create_meta_study(cancer, project_name, project_id, description, output_dir, version):
    """Function to create meta_study_file
    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        description : optional description to overwrite default description
        output_dir : path of output dir

    """
    date = datetime.today().strftime("%d-%m-%Y")
    name = "project"

    if project_name == "":
        project_name = os.path.basename(os.path.normpath(output_dir))
    project_name = project_name.upper() + " (" + date + ")"

    if project_id == "":
        project_id = cancer + "_" + name + version

    add_global_case_list = "true"

    if description == "":
        description = "Comprehensive profiling of " + cancer.capitalize() + " cancer samples."

    dictionary_file = {
        "type_of_cancer": cancer,
        "cancer_study_identifier": project_id,
        "name": project_name,
        "add_global_case_list": add_global_case_list,
        "description": description,
        }

    meta_file = open(os.path.join(output_dir, "meta_study.txt"), "w")
    logger.info("Writing meta_study.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()
    logger.info("meta_study.txt created!")

    return project_id


def create_meta_clinical_patient(project_id, output_dir):
    """Function to create meta_clinical_patient

    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        output_dir : path of output dir

    """
    alteration_type = "CLINICAL"
    datatype = "PATIENT_ATTRIBUTES"
    data_filename = "data_clinical_patient.txt"
    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "data_filename": data_filename,
    }
    meta_file = open(os.path.join(output_dir, "meta_clinical_patient.txt"), "w")
    logger.info("Writing meta_clinical_patient.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_clinical_sample(project_id, output_dir):
    """Function to create meta_clinical_sample

    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        output_dir : path of output dir

    """
    alteration_type = "CLINICAL"
    datatype = "SAMPLE_ATTRIBUTES"
    data_filename = "data_clinical_sample.txt"
    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "data_filename": data_filename,
    }
    meta_file = open(os.path.join(output_dir, "meta_clinical_sample.txt"), "w")
    logger.info("Writing meta_clinical_sample.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_mutations(cancer, project_id, profile, output_dir):
    """Function to create meta_mutations

    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        profile: Description to overwrite default description
        output_dir : path of output dir

    """
    alteration_type = "MUTATION_EXTENDED"
    datatype = "MAF"
    stable_id = "mutations"
    profile_in_analysis_tab = "TRUE"
    profile_name = "Mutations"

    if profile == " ":
        profile = "Sequencing of " + cancer.capitalize() + " tumor."

    data_filename = "data_mutations_extended.txt"

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "stable_id": stable_id,
        "show_profile_in_analysis_tab": profile_in_analysis_tab,
        "profile_name": profile_name,
        "profile_description": profile,
        "data_filename": data_filename,
    }
    meta_file = open(os.path.join(output_dir, "meta_mutations_extended.txt"), "w")
    logger.info("Writing meta_mutations_extended.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_sv(project_id, profile, output_dir):
    """Function to create meta_sv

    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        output_dir : path of output dir

    """
    alteration_type = "STRUCTURAL_VARIANT"
    datatype = "SV"
    stable_id = "structural_variants"
    profile_in_analysis_tab = "true"
    profile_name = "Structural variants"

    if profile == "":
        profile = "Structural Variant Data for " + project_id

    data_filename = "data_sv.txt"

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "stable_id": stable_id,
        "show_profile_in_analysis_tab": profile_in_analysis_tab,
        "profile_name": profile_name,
        "profile_description": profile,
        "data_filename": data_filename,
    }

    meta_file = open(os.path.join(output_dir, "meta_sv.txt"), "w")
    logger.info("Writing meta_sv.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_cna(cancer, project_id, profile, output_dir):
    """Function to create meta_cna

    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        output_dir : path of output dir

    """
    alteration_type = "COPY_NUMBER_ALTERATION"
    datatype = "DISCRETE"
    stable_id = "cna"
    profile_in_analysis_tab = "true"
    profile_name = "Putative copy-number alterations from GISTIC of " + cancer + " tumor"

    if profile == "":
        profile = "Putative copy-number from GISTIC 2.0. Values: -2 = hemizygous deletion; 0 = neutral / no change; 2 = gain."

    data_filename = "data_cna.txt"

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "stable_id": stable_id,
        "show_profile_in_analysis_tab": profile_in_analysis_tab,
        "profile_name": profile_name,
        "profile_description": profile,
        "data_filename": data_filename,
    }

    meta_file = open(os.path.join(output_dir, "meta_cna.txt"), "w")
    logger.info("Writing meta_cna.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_cna_hg19(project_id, profile, output_dir):
    """Function to create meta_cna_hg19

    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        output_dir : path of output dir

    """
    alteration_type = "COPY_NUMBER_ALTERATION"
    datatype = "SEG"
    genome_id = "hg19"
    if profile == "":
        profile = "Somatic CNA data (copy number ratio from tumor samples minus ratio from matched normals)."
    data_filename = "data_cna_hg19.seg"

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "reference_genome_id": genome_id,
        "description": profile,
        "data_filename": data_filename,
    }

    meta_file = open(os.path.join(output_dir, "meta_cna_hg19_seg.txt"), "w")
    logger.info("Writing meta_cna_hg19_seg.txt file...")
    for key, value in dictionary_file.items():
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def meta_case_main(cancer, output_folder, old_study_info=[], rename=""):

    logger.info("Starting meta_case_main script:")
    logger.info(f"meta_case_main args [cancer:{cancer}, output_file:{output_folder}]")

    version = extract_version_str(output_folder)

    config = ConfigParser()
    config.read("conf.ini")

    if len(old_study_info)==0:
        project_id = config.get("Project","PROJECT_ID")
        project_name = config.get("Project","PROJECT_NAME")

    elif len(old_study_info)>0 and not old_study_info[-1]:
        if rename!="":
            logger.info("Study will be renamed in meta")
            project_id = rename
            project_name = rename
        else:
            project_id =  os.path.basename(output_folder)
            project_name = os.path.basename(output_folder).replace("_"," ")

    elif len(old_study_info)>0 and old_study_info[-1]:
        version = extract_version_str(output_folder)
        project_id = old_study_info[0]+version
        project_name = old_study_info[1]+version.replace("_"," ")

    description=config.get("Project","DESCRIPTION")
    profile_mut=config.get("Project","PROFILE_MUT")
    profile_cna=config.get("Project","PROFILE_CNA")
    profile_cna19=config.get("Project","PROFILE_CNA_HG19")
    profile_sv=config.get("Project","PROFILE_SV")

    logger.info("Creating case list folder...")
    cases_list_dir = os.path.join(output_folder, "case_lists")
    if os.path.exists(cases_list_dir):
        pass
    else:
        os.mkdir(cases_list_dir)


    ###########  METAFILE FUNCTIONS ###########
    project_id = create_meta_study(cancer, project_name, project_id, description, output_folder, version)

    if os.path.exists(os.path.join(output_folder,"data_clinical_patient.txt")):
        create_meta_clinical_patient(project_id, output_folder)
        logger.info("meta_clinical_patient.txt created!")
    if os.path.exists(os.path.join(output_folder,"data_clinical_sample.txt")):
        create_meta_clinical_sample(project_id, output_folder)
        logger.info("meta_clinical_sample.txt created!")
    if os.path.exists(os.path.join(output_folder,"data_mutations_extended.txt")):
        create_meta_mutations(cancer, project_id, profile_mut, output_folder)
        logger.info("meta_mutations_extended.txt created!")
    if os.path.exists(os.path.join(output_folder,"data_sv.txt")):
        create_meta_sv(project_id, profile_sv, output_folder)
        logger.info("meta_sv.txt created!")
    if os.path.exists(os.path.join(output_folder,"data_cna.txt")):
        create_meta_cna(cancer, project_id, profile_cna, output_folder)
        logger.info("meta_cna.txt created!")
    if os.path.exists(os.path.join(output_folder,"data_cna_hg19.seg")):
        create_meta_cna_hg19(project_id, profile_cna19, output_folder)
        logger.info("meta_cna_hg19.txt created!")

    ########### CASE LIST FUNCTION ###########
    if os.path.exists(os.path.join(output_folder,"data_mutations_extended.txt")):
        populate_cases_sequenced(project_id, output_folder, cases_list_dir, logger)
    else: logger.warning("data_mutations_extended.txt file not found, not writing cases_sequenced.txt!")

    if os.path.exists(os.path.join(output_folder,"data_cna.txt")):
        populate_cases_cna(project_id, output_folder, cases_list_dir, logger)
    else: logger.warning("data_cna.txt file not found, not writing cases_cna.txt!")

    if os.path.exists(os.path.join(output_folder,"data_sv.txt")):
        populate_cases_sv(project_id, output_folder, cases_list_dir, logger)
    else: logger.warning("data_sv.txt file not found, not writing cases_sv.txt!")

    check_cases(output_folder)
    logger.success("Make_meta_and_cases script completed!\n")


def check_cases(output_folder):
    cases_path = os.path.join(output_folder,"case_lists")
    for case_file in os.listdir(cases_path):
        with open(os.path.join(cases_path, case_file)) as file:
            lines = file.readlines()
            if lines[-1].strip() == "case_list_ids:":
                os.remove(os.path.join(cases_path, case_file))
                logger.info(f"No samples found in {case_file}. This file will not be created")
