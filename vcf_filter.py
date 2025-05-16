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

"""Module for parsing and processing VCF files.

This module provides functionality to extract and filter data from VCF files,
write header lines to an output file, and save filtered results in tab-delimited format.

"""

from pathlib import Path

import pandas as pd


def parsing_vcf(file_input: str, file_output: str) -> None:
    """Parse a VCF file, filters variants, and writes filtered data to an output file.

    Args:
        file_input (str): Path to the input VCF file.
        file_output (str): Path to the output file where filtered data will be saved.

    This function processes the input VCF file by:
    - Removing metadata lines starting with "##".
    - Filtering 'ALT' column when not equal to "." and 'FILTER' column is "PASS".
    - Saving the filtered data into a tab-delimited file.

    Returns:
        None: Writes the filtered data to the specified output file.

    """
    with Path(file_input).open() as f:
        correct_data=[line.strip().split("\t")
                      for line in f if not line.startswith("##")]

    data=pd.DataFrame(correct_data[1:],columns=correct_data[0])
    data=data[data["ALT"]!="."]
    data=data[data["FILTER"]=="PASS"]
    data.to_csv(file_output,sep="\t", mode="a", index=False)


def write_header_lines(input_vcf: str, output_vcf: str) -> None:
    """Write the header lines from an input VCF file to an output file.

    Args:
        input_vcf (str): Path to the input VCF file.
        output_vcf (str): Path to the output file where header lines will be written.

    Returns:
        None: Writes the header lines to the specified output file.

    """
    with Path(input_vcf).open() as f_in, Path(output_vcf).open("w") as f_out:
        for line in f_in:
            if line.startswith("##"):
                f_out.write(line)


def main(input_path: str, output_path: str) -> None:
    """Write header then parse and save VCF file.

    Args:
        input_path (str): Path to the input VCF file.
        output_path (str): Path to the output folder.

    Returns:
        None: Calls the helper functions to process and write VCF data to output file.

    """
    write_header_lines(input_path, output_path)
    parsing_vcf(input_path, output_path)
