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

"""Provide the main entry point for the Varan bioinformatics pipeline.

This script parses command-line arguments, configures logging, and runs
various stages of the Varan pipeline, such as data walking, filtering,
concatenation, table creation, validation, updating, deleting, and extracting
samples from studies.

It supports flexible options for analysis types, filters, and sample handling.

"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, NoReturn

if TYPE_CHECKING:
    from collections.abc import Sequence

from loguru import logger

from concatenate import concatenate_main
from Delete_script import delete_main
from ExtractSamples_script import extract_main
from filter_clinvar import filter_main
from Make_meta_and_cases import meta_case_main
from Update_script import update_main
from ValidateFolder import validate_output
from walk import walk_folder


def logo() -> None:
    """Print the ASCII art logo for the Varan pipeline."""
    logo_text = r"""
__| |_______________________________________________________________________| |__
__   _______________________________________________________________________   __
  | |                                                                       | |
  | | █████   █████   █████████   ███████████     █████████   ██████   █████| |
  | |░░███   ░░███   ███░░░░░███ ░░███░░░░░███   ███░░░░░███ ░░██████ ░░███ | |
  | | ░███    ░███  ░███    ░███  ░███    ░███  ░███    ░███  ░███░███ ░███ | |
  | | ░███    ░███  ░███████████  ░██████████   ░███████████  ░███░░███░███ | |
  | | ░░███   ███   ░███░░░░░███  ░███░░░░░███  ░███░░░░░███  ░███ ░░██████ | |
  | |  ░░░█████░    ░███    ░███  ░███    ░███  ░███    ░███  ░███  ░░█████ | |
  | |    ░░███      █████   █████ █████   █████ █████   █████ █████  ░░█████| |
  | |     ░░░      ░░░░░   ░░░░░ ░░░░░   ░░░░░ ░░░░░   ░░░░░ ░░░░░    ░░░░░ | |
__| |_______________________________________________________________________| |__
__   _______________________________________________________________________   __
  | |                                                                       | |
"""
    logger.info(logo_text)


def varan(
    varan_input: Sequence[str] | None,
    cancer: str | None,
    output_folder: str,
    oncokb: bool,
    filters: str,
    analysis_type: str | None = None,
    overwrite_output: bool = False,
    resume: bool = False,
    multiple: bool = False,
    update: bool = False,
    extract: bool = False,
    remove: bool = False,
) -> None:
    """Run the full Varan pipeline workflow based on provided arguments.

    Execute various steps such as preparation, filtering, concatenation,
    table creation, validation, updating, deleting, and extracting samples.

    Parameters
    ----------
    varan_input : Optional[Sequence[str]]
        List of input paths (folders or files), or None if not applicable.
    cancer : Optional[str]
        Name of the cancer type or None if not provided.
    output_folder : str
        Path to the output directory.
    oncokb : bool
        Flag indicating whether to apply OncoKB annotation.
    filters : str
        String of filter options to apply during processing.
    analysis_type : Optional[str], optional
        Type of analysis to perform, by default None.
        Valid options include "snv", "cnv", "fus", or "tab".
    overwrite_output : bool, optional
        Whether to overwrite the output folder if it exists, by default False.
    resume : bool, optional
        Whether to resume an existing analysis, by default False.
    multiple : bool, optional
        Whether multiple sample VCF files are expected, by default False.
    update : bool, optional
        Whether to run the update study process, by default False.
    extract : bool, optional
        Whether to run the extract samples process, by default False.
    remove : bool, optional
        Whether to run the remove samples process, by default False.

    Returns
    -------
    None

    """
    if not any([update, extract, remove]):

        logger.info(
            f"Varan args [input:{varan_input}, output_folder:{output_folder}, "
            f"filters:{filters}, cancer:{cancer}, oncoKB:{oncokb}, "
            f"analysis_type:{analysis_type}, overwrite_output:{overwrite_output}, "
            f"resume:{resume}, multiple:{multiple}, update:{update}, "
            f"extract:{extract}, remove:{remove}]")

        ###########################
        #        1.  WALK         #
        ###########################

        logger.info("Starting preparation study folder")
        output_folder, varan_input, _ = walk_folder(
            varan_input, multiple, output_folder, oncokb, cancer,
            overwrite_output, resume, analysis_type, filters,
            )


        ###########################
        #       2. FILTER         #
        ###########################

        logger.info("Starting MAF filtering")
        if analysis_type not in ["cnv", "fus", "tab"]:
            filter_main(
                varan_input,
                output_folder,
                output_folder,
                oncokb,
                filters,
                cancer,
                resume)


        ############################
        #      3. CONCATENATE      #
        ############################

        maf_path = Path(output_folder) / "maf"

        if maf_path.exists() and analysis_type not in ["cnv", "fus", "tab"]:
            logger.info("Concatenating mutation file")
            concatenate_main(filters, output_folder, "maf", oncokb)


        ###########################################
        #      4. MAKE AND POPULATE TABLES        #
        ###########################################

        logger.info("It's time to create tables!")
        meta_case_main(cancer, output_folder)


        ############################
        #      5. VALIDATION       #
        ############################

        logger.info("Starting validation...")
        validate_output(
            output_folder,
            varan_input,
            multiple,
            False,
            cancer,
            oncokb,
            filters,
            )


    ############################
    #         UPDATE           #
    ############################

    if update:
        logger.info("Starting Update study")
        oldpath=args.Path
        new=args.NewPath
        new_name = args.Name
        output_folder=args.output_folder
        update_main(oldpath, new, output_folder, new_name, overwrite_output)


    ############################
    #         DELETE           #
    ############################

    if remove:
        logger.info("Starting Delete sample(s) from study")
        oldpath=args.Path
        removepath=args.SampleList
        new_name = args.Name
        output_folder=args.output_folder
        delete_main(oldpath, removepath, output_folder, new_name, overwrite_output)


    ############################
    #         EXTRACT          #
    ############################

    if extract:
        logger.info("Starting Extract sample(s) from study")
        oldpath = args.Path
        removepath = args.SampleList
        new_name = args.Name
        output_folder = args.output_folder
        extract_main(oldpath, removepath, output_folder, new_name, overwrite_output)


#################################################################################################################

class MyArgumentParser(argparse.ArgumentParser):
  """An argument parser that raises an error, instead of quits."""

def error(self, message: str) -> NoReturn:
    """Raise a ValueError with the given error message."""
    raise ValueError(message)

if __name__ == "__main__":

    logger.remove()
    logfile="Varan_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
    logger.level("INFO", color="<green>")

    logger.add(
        sys.stderr,
        format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",
        colorize=True,
        catch=True)
    logger.add(
        Path("Logs") / logfile,
        format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",
        mode="w")

    logo()

    logger.info("Welcome to VARAN 🦎\n")

    parser = MyArgumentParser(
        add_help=True,
        exit_on_error=False,
        usage=None,
        description="Argument of Varan script")

    # WALK BLOCK
    parser.add_argument("-c", "--Cancer", required=False,
                        help="Cancer Name")
    parser.add_argument("-i", "--varan_input", nargs="+", required=False, type=str,
    help=("list with 1) input folder/sample file tsv (required) "
    "2) patient tsv 3) fusion file"))

    parser.add_argument("-t", "--analysis_type", required=False,
    choices=["snv", "cnv", "fus", "tab"],
    help=("Select the analysis to follow (snv -> snv analysis; "
    "cnv -> cnv analysis; fus  -> fusions analysis; tab  -> table creation)"))

    parser.add_argument("-w", "--overWrite", required=False, action="store_true",
                        help="Overwrite output folder if it exists")
    parser.add_argument("-R", "--resume", required=False, action="store_true",
                        help="Resume an already started analysis")

    # ANNOTATION BLOCK
    parser.add_argument("-k", "--oncokb", required=False, action="store_true",
                        help="OncoKB annotation")
    parser.add_argument("-m", "--multiple", required=False, action="store_true",
                        help="Multiple sample VCF?")

    # FILTER BLOCK
    parser.add_argument("-f", "--Filter", required=False, default="",
                        help=("Select filter for SNV [d -> filter, p -> filter==PASS, "
                        "v-> vaf, o-> Oncokb, a -> AF, q -> Consequence, y-> polyphens,"
                        " c -> clin_sig, n -> novel, i -> Impact]"))

    # UPDATE BLOCK
    parser.add_argument("-u", "--Update", required=False,action="store_true",
                        help="Add this argument if you want to concatenate two studies")
    parser.add_argument("-n", "--NewPath", required=False,
                        help="Path of new study folder to add")

    # DELETE BLOCK
    parser.add_argument(
        "-r", "--Remove", required=False,action="store_true",
        help="Add this argument if you want to remove samples from a study")

    # EXTRACT BLOCK
    parser.add_argument(
        "-e", "--Extract", required=False, action="store_true",
        help="Add this argument if you want to extract samples from a study")

    # COMMON BLOCK
    parser.add_argument("-o", "--output_folder", required=False, default="",
                        help="Output folder")
    parser.add_argument("-s", "--SampleList", required=False,
                        help="Path of file with list of SampleIDs")
    parser.add_argument("-p", "--Path", required=False,
                        help="Path of original study folder")
    parser.add_argument(
        "-N", "--Name", required=False, default="",
        help=(
            "Add this argument if you want to give a custom name to the extract study"))

    try:
        args = parser.parse_args()

        cancer = args.Cancer
        varan_input = args.varan_input
        filters=args.Filter
        output_folder = args.output_folder
        analysis_type=args.analysis_type
        overwrite_output=args.overWrite
        resume=args.resume
        oncokb=args.oncokb
        multiple=args.multiple

        update=args.Update
        extract=args.Extract
        remove=args.Remove

        if sum([args.Update, args.Extract, args.Remove]) > 1:
            logger.critical(
                "Please select only one option between Update, Extract and "
                "Remove")
            sys.exit(1)

        if not any([args.Update, args.Extract, args.Remove]) and (
            args.varan_input is None or args.varan_input[0].strip() == ""):
            logger.critical("Error Argument: Valid Input is required")
            sys.exit(1)

        if not any([args.Update, args.Extract, args.Remove]) and args.output_folder=="":
            logger.critical("Error Argument: Output is required")
            sys.exit(1)

        if not any([args.Update, args.Extract, args.Remove]) and args.Cancer is None:
            logger.critical("Error Argument: Cancer name is required")
            sys.exit(1)

        if args.Update and (args.Path is None or args.NewPath is None):
            logger.critical(
                "To update a study, you need to specify both original "
                "and new folder paths")
            sys.exit(1)

        if (any([args.Remove, args.Extract]) and args.Path is None) or \
        (any([args.Remove, args.Extract]) and args.SampleList is None):
            logger.critical(
                "To remove/extract samples from a study, you need to specify both "
                "original folder path and samples' list")
            sys.exit(1)

        if (args.output_folder=="" and args.Name!=""):
            logger.critical("To use -N option it's required to set also -o")
            sys.exit(1)

        if "n" in filters and "v" not in filters:
            logger.critical(
                'To use the "n" option in filters it\'s required to set also the "v"')
            sys.exit(1)

        if "o" in filters and not oncokb:
            logger.critical(
                'To use the "o" option in filters it\'s required to set also -k')
            sys.exit(1)

        if resume and overwrite_output:
            logger.critical(
                "Both resume and overwrite options are selected. "
                "Please select only one!")
            sys.exit(1)


        varan(
            varan_input,
            cancer,
            output_folder,
            oncokb,
            filters,
            analysis_type,
            overwrite_output,
            resume,
            multiple,
            update,
            extract,
            remove)

    except ValueError as err:
        logger.critical(f"ValueError: {err}", file=sys.stderr)
    except FileNotFoundError as err:
        logger.critical(f"File not found: {err}", file=sys.stderr)
