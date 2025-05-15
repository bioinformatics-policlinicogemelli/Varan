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
from loguru import logger


def concatenate_files(file_list, output_file):
    with open(output_file, "w") as out_file:
        for i, file_name in enumerate(file_list):
            with open(file_name) as in_file:
                lines = in_file.readlines()
                if "Hugo_Symbol" not in lines[0]:
                    del lines[0]
                lines=list(map(lambda x: x.replace(".bam",""),lines))
                if i > 0:
                    lines = lines[1:]
                out_file.write("".join(lines))

    if not os.path.exists(output_file):
        logger.critical(f"Something went wrong while writing {output_file}.")

    logger.info(f"{len(file_list)} maf file(s) concatenated")

def get_files_by_ext(folder, ext):
    file_list = []
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith(ext):
                file_list.append(os.path.join(root, file))
    if len(file_list)==0:
        logger.warning(f"No files found with .{ext} extension in {folder}")
    else:
        logger.info(f"{len(file_list)} {ext} file(s) found")
    return file_list

def extract_maf_folder(filters, oncoKB):
    if oncoKB and "o" in filters:
        folder="MAF_Onco_filtered"
    elif filters!="d" and filters!="":
        folder="MAF_filtered"
    elif oncoKB:
        folder="MAF_OncoKB"
    else:
        folder="maf"
    return folder

def concatenate_main(filters, output_folder, ext, oncoKB):

    logger.info("Starting concatenate_main script:")
    logger.info(f"concatenate_main args [filters:{filters}, folder:{output_folder}, extension:{ext}, oncoKB:{oncoKB}]")

    folder = extract_maf_folder(filters, oncoKB)
    input_folder=os.path.join(output_folder,folder)
    output_file=os.path.join(input_folder,"data_mutations_extended.txt")

    if os.path.isdir(output_file):
        logger.critical(f"It seems that the inserted output_file '{output_file}' is not a file, but a folder! Check your '-o/--output_file' field")
        raise Exception("Exiting from filter_clinvar script!")
    if not output_file.endswith("txt"):
        logger.critical(f"It seems that the inserted output_file '{output_file}' has the wrong extension! Output file must be have a .txt extension.")
        raise Exception("Exiting from filter_clinvar script!")

    file_list = get_files_by_ext(input_folder, ext)
    concatenate_files(file_list, output_file)

    logger.info("Checking data_mutations_extended...")
    with open(output_file) as data_mut:
        all_data_mut = data_mut.readlines()
        if (len(all_data_mut) == 1):
            os.remove(output_file)
            logger.warning("data_mutations_extended is empty. File removed.")

    if os.path.exists(output_file):
        data_mut = pd.read_csv(output_file, sep="\t", dtype=str)
        data_mut.drop_duplicates(keep="last", inplace=True)

    if os.path.exists(output_file):
        logger.info(f"Extracting data_mutations_extended from {input_folder} folder")
        os.system("mv "+os.path.join(input_folder,"data_mutations_extended.txt")+" "+ output_folder )

    logger.success("Concatenate script completed!\n")
