import os 
import re 
from loguru import logger
import pandas as pd


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
    outputfolderpath = os.getcwd()
    if outputfolderpath == "":
        outputfolderpath = os.getcwd()
        
    old_versions = [file for file in os.listdir(os.path.realpath(outputfolderpath)) \
                    if re.split(r'_v[0-9]+$',os.path.basename(output_folder))[0] in file]
    version_n = [extract_version_int(version) for version in old_versions]
    sorted_version = sorted(version_n, key=int)
    version_name_ordered = list(map(lambda x: foldername + "_v" + str(x), sorted_version))
    return version_name_ordered

def get_newest_version(output_folder):
    
    foldername = re.split(r'_v[0-9]+$',os.path.basename(output_folder))[0]
    outputfolderpath = os.path.dirname(output_folder)
    if outputfolderpath == "":
        outputfolderpath = os.getcwd()
    old_versions = [file for file in os.listdir(os.path.realpath(outputfolderpath)) if foldername == re.split(r'_v[0-9]+$', file)[0]]

    logger.info(f"{len(old_versions)} output folder versions found: {old_versions}")

    old_versions_number = list(map(extract_version_int, old_versions))
    if old_versions_number == []:
        v = "_v0"
        version = "_v1"
    else:
        v = max(old_versions_number)
        version = "_v" + str(v + 1)
    output_folder_version = foldername + version
    return output_folder_version, "_v" + str(v)

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
    #vus = False
    with open(file_meta, 'r') as meta:
        for line in meta:
            if line.startswith("type_of_cancer"):
                cancer = line.split(" ")[1].strip()
            if line.startswith("cancer_study_identifier"):
                study_id = re.split(r'_v[0-9]+$',line.split(" ")[1])[0].strip()
            if line.startswith("name"):
                study_name = re.split(r'V[0-9]$',line.split(":")[1].split("(")[0].strip())[0].strip()  
 
    return cancer, [study_id, study_name]#, vus 
        
        
def extract_sample_list(filecase):
    with open(filecase, 'r') as meta:
        for line in meta:
            if line.startswith("case_list_ids:"):
               sample_part = line.split(":")[1]
               samples = sample_part.split("\t")
               sample_list = [sample.strip() for sample in samples]
    return sample_list       
        
        

def compare_sample_file(input_folder, file1, file2, filename, action, outputfolder, summary_file):
    
    if os.path.exists(file1) and os.path.exists(file2):
        samples_file1 = extract_sample_list(file1)  
        samples_file2 = extract_sample_list(file2)
        
        clin_sam_old_df = pd.read_csv(os.path.join(input_folder, "data_clinical_sample.txt"), sep="\t")
        clin_sam_new_df = pd.read_csv(os.path.join(outputfolder, "data_clinical_sample.txt"), sep="\t")
         
        if action == "delete":
            clin_pat = set(clin_sam_old_df.iloc[4:, 1]).difference(clin_sam_new_df.iloc[4:, 1])
            clin_sample = set(clin_sam_old_df.iloc[4:, 0]).difference(clin_sam_new_df.iloc[4:, 0])
            if filename == "cases_cna":
                print(f"Patients: {clin_pat} were {action}d", file=summary_file)       
                print(f"Samples: {clin_sample} were {action}d\n", file=summary_file)

            removed_samples = [sample for sample in samples_file2 if not sample in samples_file1 and sample != ""]
            if len(removed_samples) > 0:
                print(f"{len(removed_samples)} sample(s) removed from {filename}: {removed_samples}", file=summary_file)
                print(f"There are now {len(set(samples_file1)-set(removed_samples))} sample(s) in {filename}\n", 
                      file=summary_file)
            
        elif action == "extract":
            clin_pat = set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
            clin_sample = set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])
            if filename == "cases_cna":
                print(f"Patients: {clin_pat} were {action}ed", file=summary_file)       
                print(f"Samples: {clin_sample} were {action}ed\n", file=summary_file)
            extracted_samples = [sample for sample in samples_file2 if sample in samples_file1 and sample != ""]
            if len(extracted_samples) > 0:
                print(f"{len(extracted_samples)} sample(s) extracted from {filename}: {extracted_samples}")
                print(f"{len(extracted_samples)} sample(s) extracted from {filename}: {extracted_samples}\n", file=summary_file)
        
        elif action == "update":
            updated_clin_pat = set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
            updated_clin_sample = set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])
            added_clin_pat = set(clin_sam_new_df.iloc[4:, 1]).difference(clin_sam_old_df.iloc[4:, 1])
            added_clin_sample = set(clin_sam_new_df.iloc[4:, 0]).difference(clin_sam_old_df.iloc[4:, 0])

            if filename == "cases_cna":
                print(f"Patients: {list(updated_clin_pat)} were updated", file=summary_file)       
                print(f"Samples: {list(updated_clin_sample)} were updated\n", file=summary_file)
                print(f"Patients: {list(added_clin_pat)} were added", file=summary_file)       
                print(f"Samples: {list(added_clin_sample)} were added\n", file=summary_file)

            new_samples = [sample for sample in samples_file1 if not sample in samples_file2 and sample != ""]
            removed_samples = [sample for sample in samples_file2 if not sample in samples_file1 and sample != ""]
            if len(new_samples) > 0:
                print(f"{len(new_samples)} sample(s) added in {filename}: {new_samples}")
                print(f"{len(new_samples)} sample(s) added in {filename}: {new_samples}", file=summary_file)
            if len(removed_samples) > 0:
                print(f"{len(removed_samples)} sample(s) removed from {filename}: {removed_samples}")
                print(f"{len(removed_samples)} sample(s) removed from {filename}: {removed_samples}", file=summary_file)
        
        if filename == "cases_sv":
            return len(set(clin_sam_new_df.iloc[4:, 1])), len(set(clin_sam_new_df.iloc[4:, 0]))
   
    else:   
        if not os.path.exists(file1):
            print(f"{file1} does not exist")
            print(f"{file1} does not exist", file=summary_file)
        if not os.path.exists(file2):
            print(f"{file2} does not exist")
            print(f"{file2} does not exist", file=summary_file)

    
        
def compare_version(folder1, folder2, action):
    with open(os.path.join(folder1, "summary.txt"), "a") as summary_file:
        summary_file.write(f"Summary for {action.upper()} action\n\n")
        summary_file.write(f"ORIGINAL STUDY: {folder2}\n")
        summary_file.write(f"NEW STUDY: {folder1}\n\n")

        case_list1 = os.path.join(folder1, "case_lists")
        case_list2 = os.path.join(folder2, "case_lists")

        # Compare case_list_cna
        cna_1 = os.path.join(case_list1, "cases_cna.txt")
        cna_2 = os.path.join(case_list2, "cases_cna.txt")
        compare_sample_file(folder2, cna_1, cna_2, "cases_cna", action, folder1, summary_file)
        
        # Compare case_list_sequenced
        sequenced_1 = os.path.join(case_list1, "cases_sequenced.txt")
        sequenced_2 = os.path.join(case_list2, "cases_sequenced.txt")
        compare_sample_file(folder2, sequenced_1, sequenced_2, "cases_sequenced", action, folder1, summary_file)
        
        # Compare case_list_sv
        sv_1 = os.path.join(case_list1, "cases_sv.txt")
        sv_2 = os.path.join(case_list2, "cases_sv.txt")
        np, ns = compare_sample_file(folder2, sv_1, sv_2, "cases_sv", action, folder1, summary_file)
        
        summary_file.write(f"TOTAL PATIENT: {ns}\n")
        summary_file.write(f"TOTAL SAMPLE: {np}")
        