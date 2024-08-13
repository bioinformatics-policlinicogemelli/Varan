# cBioportal metafile
import argparse
import os
from configparser import ConfigParser
from populate_case_lists import populate_cases_cna,populate_cases_sequenced,populate_cases_sv
from loguru import logger
import sys
from versioning import extract_version_str
from datetime import datetime


def create_meta_study(cancer, project_name, project_id, description, output_dir, version):
    """
        Function to create meta_study_file
    Args:
        cancer : cancer type
        vus : Flag to select Vus inclusion
        description : optional description to overwrite default description
        output_dir : path of output dir
    """
    data = datetime.today().strftime('%d-%m-%Y')
    name = "project"

    if project_name == "":
        project_name = cancer.capitalize() + " Cancer " + name.upper() + " " + version.upper().replace("_", "")
    project_name = project_name + " (" + data + ")"

    if project_id == "":
        project_id = cancer + "_" + name + version 
    
    add_global_case_list = "true"

    if description == " ":
        description = "Comprehensive profiling of " + cancer.capitalize() + " cancer samples."
                
    dictionary_file = {
        "type_of_cancer": cancer, 
        "cancer_study_identifier": project_id,
        "name": project_name, 
        "add_global_case_list": add_global_case_list,
        "description": description
        }

    meta_file = open(os.path.join(output_dir, "meta_study.txt"), "w")
    logger.info("Writing meta_study.txt file...")
    for key, value in dictionary_file.items():
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()

    return project_id


def create_meta_clinical_patient(project_id, output_dir):
    """
        Function to create meta_clinical_patient
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
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_clinical_sample(project_id, output_dir):
    """
        Function to create meta_clinical_sample
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
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_mutations(cancer, project_id, profile, output_dir):
    """
        Function to create meta_mutations
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
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_sv(project_id, profile, output_dir):
    """
        Function to create meta_sv
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
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_cna(cancer, project_id, profile, output_dir):
    """
        Function to create meta_cna
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
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


def create_meta_cna_hg19(project_id, profile, output_dir):
    """
        Function to create meta_cna_hg19
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
        logger.info(f"{key}: {value}", file=meta_file)
        print(f"{key}: {value}", file=meta_file)
    meta_file.close()


# def create_cases_sequenced(cases_list_dir, project_id):
#     """
#         Function to cases_sequenced
#     Args:
#         cancer : cancer type
#         vus : Flag to select Vus inclusion
#         cases_list_dir : path of case_list output dir
#     """
    

#     stable_id = project_id+"_sequenced"

#     case_list_category = "all_cases_with_mutation_data"
#     case_list_name = "Sequenced Tumors"
#     case_list_description = "All sequenced samples (   samples)"
#     case_list_ids = " "

#     dictionary_file = {
#         "cancer_study_identifier": study_id,
#         "stable_id": stable_id,
#         "case_list_category": case_list_category,
#         "case_list_name": case_list_name,
#         "case_list_description": case_list_description,
#         "case_list_ids": case_list_ids,
#     }

#     meta_file = open(os.path.join(cases_list_dir, "cases_sequenced.txt"), "w")
#     for key, value in dictionary_file.items():
#         logger.info(f"{key}: {value}", file=meta_file)
#         print(f"{key}: {value}", file=meta_file)
#     meta_file.close()


# def create_cases_cna(cancer, cases_list_dir, version):
#     """
#         Function to cases_cna
#     Args:
#         cancer : cancer type
#         vus : Flag to select Vus inclusion
#         cases_list_dir : path of case_list output dir
#     """

#     study_id = cancer+project_name + version
    
#     stable_id = study_id+"_cna"

#     case_list_category = "all_cases_with_cna_data"
#     case_list_name = "Samples with CNA data"
#     case_list_description = "Samples with CNA data (   samples)"
#     case_list_ids = " "

#     dictionary_file = {
#         "cancer_study_identifier": study_id,
#         "stable_id": stable_id,
#         "case_list_category": case_list_category,
#         "case_list_name": case_list_name,
#         "case_list_description": case_list_description,
#         "case_list_ids": case_list_ids,
#     }

#     meta_file = open(os.path.join(cases_list_dir, "cases_cna.txt"), "w")
#     for key, value in dictionary_file.items():
#         logger.info(f"{key}: {value}", file=meta_file)
#         print(f"{key}: {value}", file=meta_file)
#     meta_file.close()
    


# def create_cases_sv(cancer, cases_list_dir, version):
#     """
#         Function to cases_sv
#     Args:
#         cancer : cancer type
#         vus : Flag to select Vus inclusion
#         cases_list_dir : path of case_list output dir
#     """
    
#     study_id = cancer + project_name + version

#     stable_id = study_id + "_sv"
#     case_list_name = "Samples with SV data"
#     case_list_description = "All samples ( ) samples"
#     case_list_category = "all_cases_with_sv_data"
#     case_list_ids = " "

#     dictionary_file = {
#         "cancer_study_identifier": study_id,
#         "stable_id": stable_id,
#         "case_list_name": case_list_name,
#         "case_list_description": case_list_description,
#         "case_list_category": case_list_category,
#         "case_list_ids": case_list_ids,
#     }

#     meta_file = open(os.path.join(cases_list_dir, "cases_sv.txt"), "w")
#     for key, value in dictionary_file.items():
#         logger.info(f"{key}: {value}", file=meta_file)
#         print(f"{key}: {value}", file=meta_file)
#     meta_file.close()



def meta_case_main(cancer, output_folder):
    logger.info("Starting meta_case_main script:")
    logger.info(f"meta_case_main args [cancer:{cancer}, output_file:{output_folder}]")
    
    config = ConfigParser()
    config.read("conf.ini")
    
    project_name = config.get("Project","PROJECT_NAME")
    project_id = config.get("Project","PROJECT_ID")
    description=config.get("Project","DESCRIPTION")
    profile_mut=config.get("Project","PROFILE_MUT")
    profile_cna=config.get("Project","PROFILE_CNA")
    profile_cna19=config.get("Project","PROFILE_CNA_HG19")
    profile_sv=config.get("Project","PROFILE_SV")

    logger.info("Creating case list folder...")
    cases_list_dir = os.path.join(output_folder, "case_lists")
    if os.path.exists(cases_list_dir):
        logger.info("It seemps that this folder already exists")
        pass
    else:
        os.mkdir(cases_list_dir)

    version = extract_version_str(output_folder)
    # output_folder = output_folder.split("_v")[0]

    ###########  METAFILE FUNCTIONS ###########
    project_id = create_meta_study(cancer, project_name, project_id, description, output_folder, version)
    
    if os.path.exists(os.path.join(output_folder,"data_clinical_patient.txt")):
        create_meta_clinical_patient(project_id, output_folder)
        logger.info("meta_clinical_patient.txt created!")
    if os.path.exists(os.path.join(output_folder,"data_clinical_sample.txt")):
        create_meta_clinical_sample(project_id, output_folder)
        logger.info("meta_clinical_sample.txt created!")
    if os.path.exists(os.path.join(os.path.join(output_folder,"data_mutations_extended.txt"))):
        create_meta_mutations(cancer, project_id, profile_mut, output_folder)
        logger.info("meta_mutations_extended.txt created!")
    if os.path.exists(os.path.join(os.path.join(output_folder,"data_sv.txt"))):
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
    else: logger.warning("data_mutations_extended.txt file not found!")
    
    if os.path.exists(os.path.join(output_folder,"data_cna.txt")):
        populate_cases_cna(project_id, output_folder, cases_list_dir, logger)
    else: logger.warning("data_cna.txt file not found!")
    
    if os.path.exists(os.path.join(output_folder,"data_sv.txt")):
        populate_cases_sv(project_id, output_folder, cases_list_dir, logger)
    else: logger.warning("data_sv.txt file not found!")

    logger.success("Make_meta_and_cases script completed!")

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__=="__main__":
    
    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='cBioportal arguments')

    parser.add_argument('-c', '--Cancer', required=False,
                        help='Cancer Name')
    parser.add_argument('-o', '--Output', required=False,
                        help='Output path')

    try:
        args = parser.parse_args()
    except Exception as err:
        logger.remove()
        logfile = "make_meta_and_cases_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}", colorize=True, catch=True)
        logger.add(os.path.join('Logs', logfile), format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}", mode="w")
        logger.critical(f"error: {err}", file=sys.stderr)
  
    cancer = args.Cancer
    output = args.Output
    
    # config = ConfigParser()

    # parse existing file
    # read config file
    # configFile = config.read("conf.ini")
    # project = config.get("Project","PROJECT_NAME")
    # project_name = "_" + project
    # description = config.get("Project","DESCRIPTION")
    # profile = config.get("Project","PROFILE")


    meta_case_main(cancer,output)