import os
from loguru import logger

def concatenate_files(file_list, output_file):
    with open(output_file, 'w') as out_file:
        for i, file_name in enumerate(file_list):
            with open(file_name) as in_file:
                lines = in_file.readlines()
                lines=list(map(lambda x: x.replace(".bam",""),lines))
                if i > 0:
                    lines = lines[1:]     
                out_file.write("".join(lines))
                
    if not os.path.exists(output_file):
        logger.critical(f"Something went wrong while writing {output_file}.")
    
    logger.info(f"#{len(file_list)} maf file(s) concatenated")
    
def get_files_by_ext(folder, ext):
    file_list = []
    for root, _, files in os.walk(folder):
        for file in files:
            if file.endswith(ext):
                file_list.append(os.path.join(root, file))
    if len(file_list)==0:
        logger.warning(f"No files found with .{ext} extension in {folder}")
    else:
        logger.info(f"#{len(file_list)} {ext} file(s) found: {file_list}")
    return file_list

def extract_maf_folder(filters, oncoKB):
    if oncoKB and "o" in filters:
        folder="MAF_Onco_filtered"
    elif "o" not in filters and (filters!="d" and filters!=""):
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
    if not output_file.endswith('txt'):
        logger.critical(f"It seems that the inserted output_file '{output_file}' has the wrong extension! Output file must be have a .txt extension.")
        raise Exception("Exiting from filter_clinvar script!")
        
    file_list = get_files_by_ext(input_folder, ext)
    concatenate_files(file_list, output_file)
    
    logger.info(f"Extracting data_mutations_extended from {input_folder} folder") 
    os.system("mv "+os.path.join(input_folder,"data_mutations_extended.txt")+" "+ output_folder )
    
    logger.success("Concatenate script completed!\n")