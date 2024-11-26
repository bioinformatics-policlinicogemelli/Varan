import os
from ExtractSamples_functions import *
from ValidateFolder import validateOutput
from Make_meta_and_cases import meta_case_main
from versioning import *
from loguru import logger
import shutil
import sys

def extract_main(oldpath, extract_path, output, study_id, overwrite):
    
    logger.info(f"extract_main args [old_path:{oldpath}, extract_path:{extract_path}, output_folder:{output}]")	
    
    logger.info("Checking input...")
    if not os.path.isdir(oldpath):
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit()
     
    if output!="":
        no_out=False
        if os.path.exists(oldpath):
            logger.info("Original folder found")
        if os.path.exists(extract_path):
            logger.info("Sample list to extract found")
    else:
        no_out=True
        output=re.split(r'_v[0-9]$',oldpath)[0]
    
    check_sample_list(extract_path, oldpath)
    old_versions=get_version_list(output)

    if len(old_versions)>0 and os.path.exists(old_versions[-1]):
        if overwrite:
            logger.info(f"Overwrite option set. Start removing folder")
            shutil.rmtree(old_versions[-1])

    output=create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")     
    # os.mkdir(output)
    output_caseslists=os.path.join(output, "case_lists")
    os.mkdir(output_caseslists)

    logger.info("Great! Everything is ready to start")
    os.system("cp " + oldpath + "/*meta* " + output)

    sampleIds = open(extract_path, "r").readlines()
    sampleIds = [sample.strip() for sample in sampleIds]

    o_clinical_patient = os.path.join(oldpath,"data_clinical_patient.txt")
    if os.path.exists(o_clinical_patient):
        extract_clinical_patient(oldpath,sampleIds,output)
    else:
        logger.warning("data_clinical_patient.txt not found in current folder. Skipping")
    
    o_clinical_sample = os.path.join(oldpath,"data_clinical_sample.txt")
    if os.path.exists(o_clinical_sample) :
        extract_clinical_samples(o_clinical_sample,sampleIds,output)
    else:
        logger.warning("data_clinical_sample.txt not found in current folder. Skipping")
    
    o_cna_hg19 = os.path.join(oldpath,"data_cna_hg19.seg")
    if os.path.exists(o_cna_hg19):
        extract_cna_hg19(o_cna_hg19,sampleIds,output)
    else:
        logger.warning("data_cna_hg19.seg not found in current folder. Skipping")
    
    o_cna_hg19_fc = os.path.join(oldpath,"data_cna_hg19.seg.fc.txt")
    if os.path.exists(o_cna_hg19_fc):
        extract_cna_hg19_fc(o_cna_hg19_fc,sampleIds,output)
    else:
        logger.warning("data_cna_hg19.seg.fc.txt not found in current folder. Skipping")

    o_cna = os.path.join(oldpath, "data_cna.txt")
    if os.path.exists(o_cna):
        extract_cna(o_cna, sampleIds, output)
    else:
        logger.warning("data_cna.txt not found in current folder. Skipping")
    
    o_mutations = os.path.join(oldpath, "data_mutations_extended.txt")
    if os.path.exists(o_mutations):
        extract_mutations(o_mutations, sampleIds, output)
    else:
        logger.warning("data_mutations_extended.txt not found in current folder. Skipping")
    
    o_sv = os.path.join(oldpath, "data_sv.txt")
    if os.path.exists(o_sv):
        extract_sv(o_sv,sampleIds, output)
    else:
        logger.warning("data_sv.txt not found in current folder. Skipping")
    
    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)
    meta_case_main(cancer, output, study_info, study_id)
    
    # if len(old_versions)>=1:
    #     old_version=old_versions[-1]
    #     compare_version(output, old_version, "extract", output)

    compare_version(oldpath, output, "extract")
    
    logger.info("Starting Validation Folder...")
    validateOutput(output, None, False, True)

    report_file_path = os.path.join(output, "report.txt")  

    with open(report_file_path, "r") as file:
        val_report = file.readlines()

    with open(report_file_path, "w") as file:
        file.write("This is the report from cBioPortal Validator. The numbers indicated are the rows where the error occurred.\n")
        file.writelines(val_report)
    
    logger.success("The process ended without errors")
    logger.success("Successfully extracted sample(s)!")
    

    
    