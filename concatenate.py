#####################################
# NAME: concatenate.py
# AUTHOR: Luciano Giaco'
# Date: 18/01/2023
version = "1.0"
# ===================================

import os
import argparse
from loguru import logger
import sys


def concatenate_files(file_list, output_file):
    with open(output_file, 'w') as out_file:
        for i, file_name in enumerate(file_list):
            with open(file_name) as in_file:
                lines = in_file.readlines()
                if i > 0:
                    lines = lines[1:]
                out_file.write("".join(lines))
                
    if not os.path.exists(output_file):
        logger.critical(f"Something went wrong while writing {output_file}.")
    
    logger.info(f"#{len(file_list)} maf file(s) concatenated")
    
def get_files_by_ext(folder, ext):
    file_list = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(ext):
                file_list.append(os.path.join(root, file))
    if len(file_list)==0:
        logger.warning(f"No files found with .{ext} extension in {folder}")
    else:
        logger.info(f"#{len(file_list)} {ext} file(s) found: {file_list}")
    return file_list



def concatenate_main(folder, ext, output_file, log=False):
    
    if not log:
        logger.remove()
        logfile="concatenate_main_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}")
    
    logger.info("Starting concatenate_main script:")
    logger.info(f"concatenate_main args [folder:{folder}, extension:{ext}, output_file:{output_file}]") 
    
    if os.path.isdir(output_file):
        logger.critical(f"It seems that the inserted output_file '{output_file}' is not a file, but a folder! Check your '-o/--output_file' field")
        raise Exception("Exiting from filter_clinvar script!")
    if not output_file.endswith('txt'):
        logger.critical(f"It seems that the inserted output_file '{output_file}' has the wrong extension! Output file must be have a .txt extension.")
        raise Exception("Exiting from filter_clinvar script!")
        
    file_list = get_files_by_ext(folder, ext)

    concatenate_files(file_list, output_file)
    
    logger.success("Concatenate script completed!\n")

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits"""
  def error(self, message):
    raise ValueError(message)

if __name__ == '__main__':

    parser = MyArgumentParser(add_help=False, exit_on_error=False, usage=None, description='Concatenate several files maf froma a given folder')
    
    parser.add_argument('-f', '--folder', required=True,
                                            help='Path folder containing the maf files')
    parser.add_argument('-e', '--extension', required=True,
                                            help='Extension of the files to concatenate (eg. maf)')
    parser.add_argument('-o', '--output_file', required=True,
                                            help='Output txt file (eg data_mutations_extended.txt)')

    try:
        args = parser.parse_args()
    except Exception as err:
        logger.remove()
        logfile="concatenate_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
        logger.level("INFO", color="<green>")
        logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True,catch=True)
        logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w")
        logger.critical(f"error: {err}", file=sys.stderr)

    folder = args.folder
    ext = args.extension
    output_file = args.output_file

    concatenate_main(folder, ext, output_file, log=False)