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

"""Convert VCF CNV data into tab-delimited tables for analysis.

This script parses CNV (Copy Number Variation) data from VCF files and
converts them into flat table formats suitable for downstream processing.
It supports different VCF versions and handles edge cases such as negative
fold changes or missing values. It also logs unusual values for review.

Functions:
- is_positive: Check if a fold change is positive, with logging for negatives.
- vcf_to_table: Convert CNV data from a VCF to a simplified table format.
- vcf_to_table_fc: Enhanced version of vcf_to_table including gene and discrete info.
- load_table: Load a tab-delimited table into a Pandas DataFrame.
- main: Run the complete VCF to table conversion pipeline.

"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd
from loguru import logger


def is_positive(number: float, sample: str) -> bool:
    """Check if a number is positive, otherwise log a warning for the given sample.

    Parameters
    ----------
    number : float
        The value to check for positivity (e.g. a fold change).
    sample : str
        The sample ID related to the value.

    Returns
    -------
    bool
        True if the number is positive, False otherwise.

    """
    if number >= 0:
        return True
    with Path("negative_FC.log").open("a") as n:
        n.write(f"[WARNING] the sample {sample} has a Fold change in CNV with "
        "negative value\n")
    return False


def extract_fc_from_fields(
    version: str,
    file_format: str,
    info: str,
    fields: list[str],
    sample: str,
    ) -> float | None:
    """Extract the fc value from a VCF record based on the VCF version and format.

    Parameters
    ----------
    version : str
        Version of the VCF file (e.g., "VCFv4.1", "VCFv4.2").
    file_format : str
        FORMAT field string from the VCF record, indicating the structure of
        sample-specific data.
    info : str
        Sample-specific data string corresponding to the FORMAT field.
    fields : list of str
        The full list of tab-separated values from the VCF line.
    sample : str
        Sample identifier, used for logging in case of errors.

    Returns
    -------
    float
        The extracted fold change (fc) value as a string (to be converted to float by
        the caller).

    Raises
    ------
    SystemExit
        If the VCF version is unsupported or expected keys are not found in the format
        strings.

    """
    if version == "VCFv4.1":
        if file_format == "FC":
            return fields[-1]
        return None
    if version == "VCFv4.2":
        format_infos = file_format.split(":")
        fc_position = format_infos.index("SM")
        sample_info = info.split(":")
        return sample_info[fc_position]

    version_number = version.split("v")[-1]
    logger.critical(
        f"The VCF file for sample {sample} is in an unsupported "
        f"version (v{version_number}). Supported versions: [4.1, 4.2]")
    sys.exit()


def extract_header_positions(header_line: str, sample: str) -> dict:
    """Extract column positions from the VCF header line.

    Parameters
    ----------
    header_line : str
        The line starting with '#' from the VCF file.
    sample : str
        The sample name to locate in the header.

    Returns
    -------
    dict
        Dictionary with column names as keys and their indices as values.

    """
    fields_names = [x.strip() for x in header_line.split("\t")]
    return {
        "chrom": fields_names.index("#CHROM"),
        "info": fields_names.index("INFO"),
        "start": fields_names.index("POS"),
        "qual": fields_names.index("QUAL"),
        "format_infos": fields_names.index("FORMAT"),
        "sample_format": fields_names.index(sample),
    }


def vcf_to_table(vcf_file: str, table_file: str, sample: str, mode: str) -> None:
    """Convert a VCF CNV file into a tab-delimited file with log2 fc values.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file.
    table_file : str
        Path to the output table file.
    sample : str
        Sample ID to extract from the VCF.
    mode : str
        File mode, "a" for append or "w" for write.

    Returns
    -------
    None

    """
    table_path = Path(table_file)
    vcf_path = Path(vcf_file)
    mode = "a" if table_path.exists() else "w"

    with vcf_path.open() as vcf, table_path.open(mode) as table:
        if mode != "a":
            table.write("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n")

        sample = sample.split(".")[0]
        header_positions = {}

        for line in vcf:
            if line.startswith("##fileformat"):
                version = (line.split("=")[-1]).strip()

            if line.startswith("##"):
                continue

            if line.startswith("#"):
                header_positions = extract_header_positions(line, sample)
                continue

            fields = line.strip().split("\t")
            chrom = fields[header_positions["chrom"]].strip("chr")
            start = fields[header_positions["start"]]
            infos = fields[header_positions["info"]].split(";")
            end = next(info.split("=")[-1] for info in infos if "END" in info)
            qual = fields[header_positions["qual"]]
            info = fields[header_positions["sample_format"]]
            file_format = fields[header_positions["format_infos"]]

            fc = extract_fc_from_fields(version, file_format, info, fields, sample)

            if fc != ".":
                fc = float(fc)
                if is_positive(fc, sample):
                    log2fc = math.log2(fc)
                else:
                    fc = 0.0001
                    log2fc = math.log2(fc)
            else:
                continue

            table.write(f"{sample}\t{chrom}\t{start}\t{end}\t{qual}\t{log2fc}\n")


def parse_fc_and_gene(
    version: str,
    infos: list[str],
    sample_infos: str,
    file_format: str,
    fields: list[str],
    sample: str) -> tuple[str, str]:
    """Extract the fold change (FC) and gene name from VCF data fields.

    Parameters
    ----------
    version : str
        VCF file version.
    infos : list[str]
        List of semicolon-separated INFO field entries.
    sample_infos : str
        Sample-specific field content.
    file_format : str
        The FORMAT field value.
    fields : list[str]
        Entire parsed VCF line.
    sample : str
        Sample ID.

    Returns
    -------
    tuple[str, str]
        FC value and gene name.

    """
    if version == "VCFv4.1":
        gene = next(info.split("=")[-1] for info in infos if "ANT" in info)
        fc = fields[-1] if file_format == "FC" else "."
    elif version == "VCFv4.2":
        gene = next(info.split("=")[-1] for info in infos if "SEGID" in info)
        format_infos = file_format.split(":")
        sample_infos_split = sample_infos.split(":")
        fc_position = format_infos.index("SM")
        fc = sample_infos_split[fc_position]
    else:
        version_number = version.split("v")[-1]
        logger.critical(
            f"The VCF file for sample {sample} is in an unsupported version "
            f"(v{version_number}).")
        sys.exit()
    return fc, gene


def vcf_to_table_fc(vcf_file: str, table_file: str, sample: str, mode: str) -> None:
    """Convert VCF CNV to table format with FC, gene, and discrete CNV type values.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file.
    table_file : str
        Path to the output table file.
    sample : str
        Sample ID to extract from the VCF.
    mode : str
        File mode, "a" for append or "w" for write.

    Returns
    -------
    None

    """
    table_path = Path(table_file)
    vcf_path = Path(vcf_file)

    mode = "a" if table_path.exists() else "w"
    with vcf_path.open() as vcf, table_path.open(mode) as table:
        if mode != "a":
            table.write("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\tgene\tdiscrete\n")
        sample=sample.split(".")[0]

        for line in vcf:
            if line.startswith("##fileformat"):
                version = (line.split("=")[-1]).strip()

            if line.startswith("##"):
                continue

            if line.startswith("#"):
                fields_names = line.split("\t")
                fields_names = [x.strip() for x in fields_names]
                chrom_position = fields_names.index("#CHROM")
                info_position = fields_names.index("INFO")
                start_position = fields_names.index("POS")
                qual_position = fields_names.index("QUAL")
                alt_position = fields_names.index("ALT")
                format_position = fields_names.index("FORMAT")
                sample_position = fields_names.index(sample)
                continue

            fields = line.strip().split("\t")

            chrom = fields[chrom_position].strip("chr")
            start = fields[start_position]
            infos = fields[info_position].split(";")
            end = next(info.split("=")[-1] for info in infos if "END" in info)
            qual = fields[qual_position]
            alt = fields[alt_position]
            file_format = fields[format_position]
            sample_infos = fields[sample_position]

            fc, gene = parse_fc_and_gene(
                version,
                infos,
                sample_infos,
                file_format,
                fields,
                sample)

            if fc != ".":
                fc = float(fc)
                if not is_positive(fc, sample):
                    fc = 0.0001
            else:
                continue

            ########################
            # manage discrete data #
            ########################

            if alt == "<DUP>":
                discr = "2"
            elif alt == "<DEL>":
                discr = "-2"
            else:
                discr = "0"

            table.write(f"{sample}\t{chrom}\t{start}\t{end}\t{qual}\t{fc}\t{gene}\t{discr}\n")


def load_table(file_path: str) -> pd.DataFrame:
    """Load a tab-delimited table file into a Pandas DataFrame.

    Parameters
    ----------
    file_path : str
        Path to the file to be read.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the table contents.

    """
    return pd.read_csv(file_path, sep="\t", header=0)


def main(my_input: str, my_output: str, sample: str, mode: str) -> None:
    """Run the full VCF to table conversion pipeline for a sample.

    Parameters
    ----------
    my_input : str
        Path to the input VCF file.
    my_output : str
        Path to the output table file.
    sample : str
        Sample ID to process.
    mode : str
        Mode to open the output file: 'a' for append, 'w' for write.

    Returns
    -------
    None

    """
    vcf_to_table(my_input, my_output, sample, mode)
    vcf_to_table_fc(my_input, my_output, sample, mode)
