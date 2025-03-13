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
from ExtractSamples_functions import *
from ValidateFolder import validateOutput, copy_maf
from Make_meta_and_cases import meta_case_main
from versioning import *
from loguru import logger
import shutil
import sys
from write_report import *
from filter_clinvar import check_bool

config = ConfigParser()
configFile = config.read("conf.ini")



def extract_main(oldpath, extract_path, output, study_id, overwrite):
    
    logger.info(f"extract_main args [old_path:{oldpath}, extract_path:{extract_path}, output_folder:{output}]")	
    logger.info("Checking input...")
    oldpath = oldpath.rstrip("/")

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
        output=re.split(r'_v[0-9]+$',oldpath)[0]
    
    check_sample_list(extract_path, oldpath)
    old_versions=get_version_list(output)

    if len(old_versions)>0 and os.path.exists(old_versions[-1]):
        if overwrite:
            logger.info(f"Overwrite option set. Start removing folder")
            shutil.rmtree(old_versions[-1])

    output=create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")     

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

    ZIP_MAF = config.get('Zip', 'ZIP_MAF')
    ZIP_MAF = check_bool(ZIP_MAF)
    COPY_MAF = config.get('Zip', 'COPY_MAF')
    COPY_MAF = check_bool(COPY_MAF)
    copy_maf(oldpath, output, COPY_MAF, ZIP_MAF)
    
    logger.info("Starting Validation Folder...")
    number_for_graph = validateOutput(output, None, False, True, None, None, None)

    logger.info("Starting writing report_VARAN.html...")
    write_report_extract(oldpath, output, number_for_graph)
    
    logger.success("The process ended without errors")
    logger.success("Successfully extracted sample(s)!")