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
import re 
from loguru import logger
import pandas as pd
from datetime import datetime
from collections import defaultdict


def extract_version_str(foldername):
    version = extract_version_int(foldername)
    versionname  = "_v" + str(version)
    return versionname

def extract_version_int(foldername):
    try:
        version = int(re.search(r'_v(\d+)$', foldername).group(1))
    except:
        version = foldername
    return version


def get_version_list(output_folder):
    foldername = re.split(r'_v[0-9]+$',os.path.basename(output_folder))[0]
    outputfolderpath = os.path.dirname(output_folder)
    if outputfolderpath == "":
        outputfolderpath = os.getcwd()
    old_versions = [file for file in os.listdir(os.path.realpath(outputfolderpath)) if re.split(fr'{output_folder}_v[0-9]+$', os.path.basename(output_folder))[0] + "_v" in file]

    global old_version_exists
    if old_versions == []:
        old_version_exists = False
    else:
        old_version_exists = True

    version_n = [extract_version_int(version) for version in old_versions]
    version_n=[elem for elem in version_n if type(elem) == int]
    sorted_version = sorted(list(set(version_n)), key=int)
    version_name_ordered = list(map(lambda x: foldername + "_v" + str(x), sorted_version))
    return version_name_ordered

def get_newest_version(output_folder):
    foldername = re.split(r'_v[0-9]+$',os.path.basename(output_folder))[0]
    outputfolderpath = os.path.dirname(output_folder)
    if outputfolderpath == "":
        outputfolderpath = os.getcwd()

    old_versions = [file for file in os.listdir(os.path.realpath(outputfolderpath)) if re.match(rf"^{re.escape(foldername)}_v[0-9]+$", file)]

    logger.info(f"{len(old_versions)} version of the selected output folder found: {old_versions}")
    old_versions_number = list(map(extract_version_int, old_versions))
    if old_versions_number == []:
        v = "_v0"
        version = "_v1"
    else:
        v = max(old_versions_number)
        version = "_v" + str(v + 1)
    output_folder_version = foldername + version
    return os.path.join(os.path.dirname(output_folder),output_folder_version), "_v" + str(v)

def create_newest_version_folder(outputfolder):
    
    if len(get_version_list(outputfolder)) == 0:
        version = "_v1"
        outputfolder_newest_version = outputfolder + version
        os.mkdir(outputfolder_newest_version)
    else:
        outputfolder_newest_version, _= get_newest_version(outputfolder)
        os.mkdir(outputfolder_newest_version)
    
    return outputfolder_newest_version


def extract_info_from_meta(folder):
    file_meta = os.path.join(folder, "meta_study.txt")
    with open(file_meta, 'r') as meta:
        for line in meta:
            if line.startswith("type_of_cancer"):
                cancer = line.split(" ")[1].strip()
            if line.startswith("cancer_study_identifier"):
                study_id = re.split(r'_v[0-9]+$',line.split(" ")[1])[0].strip()
            if line.startswith("name"):
                study_name = re.split(r'V[0-9]+$',line.split(":")[1].split("(")[0].strip())[0].strip()
 
    return cancer, [study_id, study_name]