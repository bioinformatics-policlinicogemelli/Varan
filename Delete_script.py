import os
import sys
import argparse
import shutil
from Make_meta_and_cases import meta_case_main
from Delete_functions import *
from ValidateFolder import validateOutput, copy_maf
from versioning import *
from loguru import logger
from write_report import *


config = ConfigParser()
configFile = config.read("conf.ini")

ZIP_MAF = config.getboolean('Zip', 'ZIP_MAF')
COPY_MAF = config.getboolean('Zip', 'COPY_MAF')


def delete_main(oldpath, removepath, output, study_id, overwrite):
    
    logger.info(f"delete_main args [old_path:{oldpath}, remove_path:{removepath}, destination_folder:{output}]")
    logger.info("Checking input...")

    oldpath = oldpath.rstrip("/")

    if not os.path.isdir(oldpath):
        logger.critical(f"{oldpath} is not a valid folder!")
        sys.exit()	

    if output!="":
        no_out=False
        if os.path.exists(oldpath):
            logger.info("Original folder found")
        if os.path.exists(removepath):
            logger.info("Sample list to remove found")
    else:
        no_out=True
        output=re.split(r'_v[0-9]+$',oldpath)[0]

    check_sample_list(removepath, oldpath)
    old_versions = get_version_list(output)

    if len(old_versions)>0 and os.path.exists(old_versions[-1]):
        if overwrite:
            logger.info(f"Overwrite option set. Start removing folder")
            shutil.rmtree(old_versions[-1])
    
    output = create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")
    output_caseslists = os.path.join(output, "case_lists")
    os.mkdir(output_caseslists)
    
    logger.info("Great! Everything is ready to start")

    os.system("cp " + oldpath + "/*meta* " + output)

    with open(removepath) as sample_list:
        first_line = sample_list.readline()
        if len(first_line.split("\t")) > 1:
            logger.warning(f"The file {removepath} contains more than a column. It may not be in the correct format!")

    sampleIds = open(removepath,"r").readlines()
    sampleIds = [sample.strip() for sample in sampleIds]

    o_clinical_patient=os.path.join(oldpath,"data_clinical_patient.txt")
    if os.path.exists(o_clinical_patient):
        delete_clinical_patient(oldpath,sampleIds,output)
    else:
        logger.warning("data_clinical_patient.txt not found in current folder. Skipping")
    #
    o_clinical_sample=os.path.join(oldpath,"data_clinical_sample.txt")
    if os.path.exists(o_clinical_sample) :
        delete_clinical_samples(o_clinical_sample,sampleIds,output)
    else:
        logger.warning("data_clinical_sample.txt not found in current folder. Skipping")
    
    #
    o_cna_hg19=os.path.join(oldpath,"data_cna_hg19.seg")
    if os.path.exists(o_cna_hg19):
        delete_cna_hg19(o_cna_hg19,sampleIds,output)
    else:
        logger.warning("data_cna_hg19.seg not found in current folder. Skipping")

    #
    o_cna_hg19_fc=os.path.join(oldpath,"data_cna_hg19.seg.fc.txt")
    if os.path.exists(o_cna_hg19_fc):
        delete_cna_hg19_fc(o_cna_hg19_fc,sampleIds,output)
    else:
        logger.warning("data_cna_hg19.seg.fc.txt not found in current folder. Skipping")

    #
    o_cna=os.path.join(oldpath,"data_cna.txt")
    if os.path.exists(o_cna):
        delete_cna(o_cna,sampleIds,output)
    else:
        logger.warning("data_cna.txt not found in current folder. Skipping")
    
    #
    o_mutations=os.path.join(oldpath,"data_mutations_extended.txt")
    if os.path.exists(o_mutations):
        delete_mutations(o_mutations,sampleIds,output)
    else:
        logger.warning("data_mutations_extended.txt not found in current folder. Skipping")
    #
    o_sv=os.path.join(oldpath,"data_sv.txt")
    if os.path.exists(o_sv):
        delete_sv(o_sv,sampleIds,output)
    else:
        logger.warning("data_sv.txt not found in current folder. Skipping")
    
    cancer, study_info = extract_info_from_meta(oldpath)
    study_info.append(oldpath)
    study_info.append(no_out)
    
    meta_case_main(cancer, output, study_info, study_id)
    
    # if len(old_versions)>=1:
    #     old_version=old_versions[-1]
    #     compare_version(output, old_version, "delete", output)

    copy_maf(oldpath, output, COPY_MAF, ZIP_MAF)

    logger.info("Starting Validation Folder...")
    number_for_graph = validateOutput(output, None, False, True, None, None, None)

    logger.info("Starting writing report_VARAN.html...")
    write_report_remove(oldpath, output, number_for_graph)

    logger.success("The process ended without errors")
    logger.success("Successfully removed sample(s)!")