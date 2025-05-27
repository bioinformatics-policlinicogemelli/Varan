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

"""Meta file and case list generator for cancer genomics projects.

This script generates meta files and case lists for studies using the cBioPortal
format. It supports mutation, CNA, SV, and clinical data, based on the contents of
the output directory. Also handles study renaming and versioning from configuration.

"""
from __future__ import annotations

from configparser import ConfigParser
from datetime import datetime, timezone
from pathlib import Path

from loguru import logger

from populate_case_lists import (populate_cases_cna, populate_cases_sequenced,
                                 populate_cases_sv)
from versioning import extract_version_str


def create_meta_study(cancer: str, project_name: str, project_id: str,
                      description: str, output_dir: str, version: str) -> str:
    """Create the meta_study.txt file for a cBioPortal study.

    Args:
        cancer: Cancer type (e.g., 'lung').
        project_name: Display name of the project.
        project_id: Study identifier (used by cBioPortal).
        description: Study description.
        output_dir: Path to the output directory.
        version: Version string to append.

    Returns:
        A string containing the resolved project ID.

    """
    date = datetime.now(timezone.utc).strftime("%d-%m-%Y")
    name = "project"

    if project_name == "":
        project_name = Path(output_dir).name
    project_name = project_name.upper() + " (" + date + ")"

    if project_id == "":
        project_id = cancer + "_" + name + version

    add_global_case_list = "true"

    if description == "":
        description = (
            f"Comprehensive profiling of {cancer.capitalize()} cancer samples.")

    dictionary_file = {
        "type_of_cancer": cancer,
        "cancer_study_identifier": project_id,
        "name": project_name,
        "add_global_case_list": add_global_case_list,
        "description": description,
        }

    meta_file_path = Path(output_dir) / "meta_study.txt"
    logger.info("Writing meta_study.txt file...")

    with meta_file_path.open("w") as meta_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=meta_file)

    logger.info("meta_study.txt created!")

    return project_id


def create_meta_clinical_patient(project_id: str, output_dir: str) -> None:
    """Create the meta_clinical_patient.txt file.

    Args:
        project_id: Study identifier.
        output_dir: Path to the output directory.

    Returns:
        None.

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

    file_path = Path(output_dir) / "meta_clinical_patient.txt"
    logger.info("Writing meta_clinical_patient.txt file...")
    with file_path.open("w") as meta_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=meta_file)


def create_meta_clinical_sample(project_id: str, output_dir: str) -> None:
    """Create the meta_clinical_sample.txt file.

    Args:
        project_id: Study identifier.
        output_dir: Path to the output directory.

    Returns:
        None.

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

    file_path = Path(output_dir) / "meta_clinical_sample.txt"
    logger.info("Writing meta_clinical_sample.txt file...")
    with file_path.open("w") as meta_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=meta_file)


def create_meta_mutations(cancer: str, project_id: str, profile: str,
                          output_dir: str) -> None:
    """Create the meta_mutations_extended.txt file.

    Args:
        cancer: Cancer type.
        project_id: Study identifier.
        profile: Description of sequencing profile.
        output_dir: Path to the output directory.

    Returns:
        None.

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

    file_path = Path(output_dir) / "meta_mutations_extended.txt"
    logger.info("Writing meta_mutations_extended.txt file...")
    with file_path.open("w") as meta_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=meta_file)


def create_meta_sv(project_id: str, profile: str, output_dir: str) -> None:
    """Create the meta_sv.txt file for structural variants.

    Args:
        project_id: Study identifier.
        profile: Description of SV profile.
        output_dir: Path to the output directory.

    Returns:
        None.

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

    file_path = Path(output_dir) / "meta_sv.txt"
    logger.info("Writing meta_sv.txt file...")
    with file_path.open("w") as meta_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=meta_file)


def create_meta_cna(cancer: str, project_id: str, profile: str,
                    output_dir: str) -> None:
    """Create the meta_cna.txt file for discrete CNA.

    Args:
        cancer: Cancer type.
        project_id: Study identifier.
        profile: CNA profile description.
        output_dir: Path to the output directory.

    Returns:
        None.

    """
    alteration_type = "COPY_NUMBER_ALTERATION"
    datatype = "DISCRETE"
    stable_id = "cna"
    profile_in_analysis_tab = "true"
    profile_name = (f"Putative copy-number alterations from GISTIC of {cancer} tumor")

    if profile == "":
        profile = (
            "Putative copy-number from GISTIC 2.0. Values: -2 = "
            "hemizygous deletion; 0 = neutral / no change; 2 = gain.")

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

    meta_file_path = output_dir / "meta_cna.txt"
    with meta_file_path.open("w") as meta_file:
        logger.info("Writing meta_cna.txt file...")
        for key, value in dictionary_file.items():
            meta_file.write(f"{key}: {value}\n")


def create_meta_cna_hg19(project_id: str, profile: str, output_dir: str) -> None:
    """Create the meta_cna_hg19_seg.txt file for hg19 SEG data.

    Args:
        project_id: Study identifier.
        profile: SEG file description.
        output_dir: Path to the output directory.

    Returns:
        None.

    """
    alteration_type = "COPY_NUMBER_ALTERATION"
    datatype = "SEG"
    genome_id = "hg19"
    if profile == "":
        profile = (
            "Somatic CNA data (copy number ratio from tumor "
            "samples minus ratio from matched normals).")
    data_filename = "data_cna_hg19.seg"

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "reference_genome_id": genome_id,
        "description": profile,
        "data_filename": data_filename,
    }

    meta_file_path = output_dir / "meta_cna_hg19_seg.txt"
    with meta_file_path.open("w") as meta_file:
        logger.info("Writing meta_cna_hg19_seg.txt file...")
        for key, value in dictionary_file.items():
            meta_file.write(f"{key}: {value}\n")


def meta_case_main(
    cancer: str,
    output_folder: str,
    old_study_info: list[str] | None = None,
    rename: str = "",
) -> None:
    """Generate all metadata files and case lists for a cancer study.

    This function reads the project configuration, manages versioning,
    initializes project details, and triggers the creation of necessary
    meta files and case list files based on available input data.

    Args:
        cancer: Cancer type identifier.
        output_folder: Directory path where output files will be written.
        old_study_info: Optional list containing old study ID and name;
            if None or empty, uses config defaults.
        rename: Optional string to rename the study project.

    Returns:
        None

    """
    logger.info("Starting meta_case_main script:")
    logger.info(f"meta_case_main args [cancer:{cancer}, output_file:{output_folder}]")

    version = extract_version_str(output_folder)

    config = ConfigParser()
    config.read("conf.ini")

    if old_study_info is None:
        old_study_info = []

    if len(old_study_info)==0:
        project_id = config.get("Project","PROJECT_ID")
        project_name = config.get("Project","PROJECT_NAME")

    elif len(old_study_info)>0 and not old_study_info[-1]:
        if rename!="":
            logger.info("Study will be renamed in meta")
            project_id = rename
            project_name = rename
        else:
            project_id = Path(output_folder).name
            project_name = project_id.replace("_", " ")

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
    cases_list_dir = Path(output_folder) / "case_lists"
    cases_list_dir.mkdir(exist_ok=True)


    ###########  METAFILE FUNCTIONS  ###########
    project_id = create_meta_study(
        cancer,
        project_name,
        project_id,
        description,
        output_folder,
        version)

    output_folder_path = Path(output_folder)

    if (output_folder_path / "data_clinical_patient.txt").exists():
        create_meta_clinical_patient(project_id, output_folder)
        logger.info("meta_clinical_patient.txt created!")

    if (output_folder_path / "data_clinical_sample.txt").exists():
        create_meta_clinical_sample(project_id, output_folder)
        logger.info("meta_clinical_sample.txt created!")

    if (output_folder_path / "data_mutations_extended.txt").exists():
        create_meta_mutations(cancer, project_id, profile_mut, output_folder)
        logger.info("meta_mutations_extended.txt created!")

    if (output_folder_path / "data_sv.txt").exists():
        create_meta_sv(project_id, profile_sv, output_folder)
        logger.info("meta_sv.txt created!")

    if (output_folder_path / "data_cna.txt").exists():
        create_meta_cna(cancer, project_id, profile_cna, output_folder)
        logger.info("meta_cna.txt created!")

    if (output_folder_path / "data_cna_hg19.seg").exists():
        create_meta_cna_hg19(project_id, profile_cna19, output_folder)
        logger.info("meta_cna_hg19.txt created!")

    ########### CASE LIST FUNCTION ###########

    if (Path(output_folder) / "data_mutations_extended.txt").exists():
        populate_cases_sequenced(project_id, output_folder, cases_list_dir, logger)
    else:
        logger.warning("data_mutations_extended.txt file not found, "
        "not writing cases_sequenced.txt!")

    if (Path(output_folder) / "data_cna.txt").exists():
        populate_cases_cna(project_id, output_folder, cases_list_dir, logger)
    else:
        logger.warning("data_cna.txt file not found, not writing cases_cna.txt!")

    if(Path(output_folder) / "data_sv.txt").exists():
        populate_cases_sv(project_id, output_folder, cases_list_dir, logger)
    else:
        logger.warning("data_sv.txt file not found, not writing cases_sv.txt!")

    check_cases(output_folder)
    logger.success("Make_meta_and_cases script completed!\n")


def check_cases(output_folder: str) -> None:
    """Remove empty case list files from the case_lists folder.

    Args:
        output_folder: Path to the output directory.

    Returns:
        None.

    """
    cases_path = Path(output_folder) / "case_lists"
    for case_file in cases_path.iterdir():
        if case_file.is_file():
            with case_file.open() as file:
                lines = file.readlines()
                if lines and lines[-1].strip() == "case_list_ids:":
                    case_file.unlink()
                    logger.info(f"No samples found in {case_file.name}. "
                                "This file will not be created")
