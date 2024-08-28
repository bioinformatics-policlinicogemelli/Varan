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
        

def compare_version(folder1, folder2, action):
    with open(os.path.join(folder2, "summary.txt"), "a") as summary_file:
        summary_file.write(f"Summary for {action.upper()} action\n\n")
        summary_file.write(f"ORIGINAL STUDY: {folder1}\n")
        summary_file.write(f"NEW STUDY: {folder2}\n\n")

        case_list1 = os.path.join(folder1, "case_lists")
        case_list2 = os.path.join(folder2, "case_lists")

        # Compare case_list_cna
        cna_1 = os.path.join(case_list1, "cases_cna.txt")
        cna_2 = os.path.join(case_list2, "cases_cna.txt")
        compare_sample_file(cna_1, cna_2, folder1, folder2, "cases_cna", action, summary_file)
        
        # Compare case_list_sequenced
        sequenced_1 = os.path.join(case_list1, "cases_sequenced.txt")
        sequenced_2 = os.path.join(case_list2, "cases_sequenced.txt")
        compare_sample_file(sequenced_1, sequenced_2, folder1, folder2, "cases_sequenced", action, summary_file)
        
        # Compare case_list_sv
        sv_1 = os.path.join(case_list1, "cases_sv.txt")
        sv_2 = os.path.join(case_list2, "cases_sv.txt")
        np, ns = compare_sample_file(sv_1, sv_2, folder1, folder2, "cases_sv", action, summary_file)
        
        summary_file.write(f"\nTOTAL PATIENT: {ns}\n")
        summary_file.write(f"TOTAL SAMPLE: {np}")



def compare_sample_file(file1, file2, input_folder, outputfolder, filename, action, summary_file):

    clin_sam_old_df = pd.read_csv(os.path.join(input_folder, "data_clinical_sample.txt"), sep="\t")
    clin_sam_new_df = pd.read_csv(os.path.join(outputfolder, "data_clinical_sample.txt"), sep="\t")
    
    if os.path.exists(file1) and os.path.exists(file2):
        samples_file1 = extract_sample_list(file1)  
        samples_file2 = extract_sample_list(file2)
         
        if action == "delete":
            clin_pat = set(clin_sam_old_df.iloc[4:, 1]).difference(clin_sam_new_df.iloc[4:, 1])
            clin_sample = set(clin_sam_old_df.iloc[4:, 0]).difference(clin_sam_new_df.iloc[4:, 0])
            if filename == "cases_cna":
                summary_file.write(f"Patients: {list(clin_pat)} was/were removed\n")
                summary_file.write(f"Samples: {list(clin_sample)} was/were removed\n\n")

            removed_samples = [sample for sample in samples_file1 if not sample in samples_file2 and sample != ""]
            summary_file.write(f"{len(removed_samples)} sample(s) removed from {filename}: {removed_samples}\n")
            summary_file.write(f"There is/are now {len(set(samples_file1)-set(removed_samples))} sample(s) in {filename}\n")
            
        elif action == "extract":
            clin_pat = set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
            clin_sample = set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])
            if filename == "cases_cna":
                summary_file.write(f"Patients: {list(clin_pat)} was/were extracted\n")
                summary_file.write(f"Samples: {list(clin_sample)} was/were extracted\n\n")     

            extracted_samples = [sample for sample in samples_file2 if sample in samples_file1 and sample != ""]
            if len(extracted_samples) > 0:
                summary_file.write(f"{len(extracted_samples)} sample(s) extracted from {filename}: {extracted_samples}\n")
    
    else:
        if not os.path.exists(file1):
            summary_file.write(f"{file1} does not exist\n")

    if filename == "cases_sv":
        return len(set(clin_sam_new_df.iloc[4:, 1])), len(set(clin_sam_new_df.iloc[4:, 0]))


def compare_version_update(folder1, folder2, output, action):
    
    with open(os.path.join(output, "summary.txt"), "a") as summary_file:
        summary_file.write(f"Summary for {action.upper()} action\n\n")
        summary_file.write(f"ORIGINAL STUDY: {folder1}\n")
        summary_file.write(f"UPDATING WITH: {folder2}\n")
        summary_file.write(f"NEW STUDY: {output}\n\n")

        case_list1 = os.path.join(folder1, "case_lists")
        case_list2 = os.path.join(folder2, "case_lists")

        # Compare case_list_cna
        cna_1 = os.path.join(case_list1, "cases_cna.txt")
        cna_2 = os.path.join(case_list2, "cases_cna.txt")
        compare_sample_file_update(cna_1, cna_2, folder1, folder2, output, "cases_cna", summary_file)
        
        # Compare case_list_sequenced
        sequenced_1 = os.path.join(case_list1, "cases_sequenced.txt")
        sequenced_2 = os.path.join(case_list2, "cases_sequenced.txt")
        compare_sample_file_update(sequenced_1, sequenced_2, folder1, folder2, output, "cases_sequenced", summary_file)
        
        # Compare case_list_sv
        sv_1 = os.path.join(case_list1, "cases_sv.txt")
        sv_2 = os.path.join(case_list2, "cases_sv.txt")
        np, ns = compare_sample_file_update(sv_1, sv_2, folder1, folder2, output, "cases_sv", summary_file)
        
        summary_file.write(f"\nTOTAL PATIENT: {ns}\n")
        summary_file.write(f"TOTAL SAMPLE: {np}")


def compare_sample_file_update(file1, file2, input_folder1, input_folder2, outputfolder, filename, summary_file):

    if os.path.exists(file1) and os.path.exists(file2):
        samples_file1 = extract_sample_list(file1)  
        samples_file2 = extract_sample_list(file2)
        clin_sam_old_df = pd.read_csv(os.path.join(input_folder1, "data_clinical_sample.txt"), sep="\t")
        clin_sam_new_df = pd.read_csv(os.path.join(input_folder2, "data_clinical_sample.txt"), sep="\t")

        updated_clin_pat = set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
        updated_clin_sample = set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])
        added_clin_pat = set(clin_sam_new_df.iloc[4:, 1]).difference(clin_sam_old_df.iloc[4:, 1])
        added_clin_sample = set(clin_sam_new_df.iloc[4:, 0]).difference(clin_sam_old_df.iloc[4:, 0])

        if filename == "cases_cna":
            summary_file.write(f"Updated patients: {list(updated_clin_pat)}\n")
            summary_file.write(f"Updated samples: {list(updated_clin_sample)}\n")
            summary_file.write(f"Added patients: {list(added_clin_pat)}\n")
            summary_file.write(f"Added samples: {list(added_clin_sample)}\n\n")

        updated_samples = [sample for sample in samples_file1 if sample in samples_file2 and sample != ""]
        new_samples = [sample for sample in samples_file2 if not sample in samples_file1 and sample != ""]
        summary_file.write(f"{len(updated_samples)} sample(s) updated from {filename}: {updated_samples}\n")
        summary_file.write(f"{len(new_samples)} sample(s) added in {filename}: {new_samples}\n")

    else:   
        if not os.path.exists(file1):
            summary_file.write(f"{file1} does not exist\n")
            
        if not os.path.exists(file2):
            summary_file.write(f"{file2} does not exist\n")
    
    if filename == "cases_sv":
        clin_sam_output_df = pd.read_csv(os.path.join(outputfolder, "data_clinical_sample.txt"), sep="\t")
        return len(set(clin_sam_output_df.iloc[4:, 1])), len(set(clin_sam_output_df.iloc[4:, 0]))