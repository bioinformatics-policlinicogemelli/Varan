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
from datetime import datetime
from pathlib import Path
import subprocess
from typing import TYPE_CHECKING, NoReturn

if TYPE_CHECKING:
    from collections.abc import Sequence

from loguru import logger

from config_loader import get_config, set_config_path
from vendor_adapters import ADAPTERS
from versioning import get_git_version

# NOTE: the rest of Varan's own modules (concatenate, Delete_script,
# ExtractSamples_script, filter_clinvar, Make_meta_and_cases, Update_script,
# ValidateFolder, walk, write_report) are intentionally *not* imported here.
# They each read conf.ini at their own module's import time, so importing
# them before --config has been parsed would make that flag a no-op. They're
# imported further down, inside `if __name__ == "__main__":`, right after
# set_config_path() is called. get_git_version() lives in versioning.py
# specifically so the --version banner (built before argument parsing) can
# be printed without pulling in a conf.ini-dependent module early.
#
# vendor_adapters (ADAPTERS) is safe to import up here despite that rule:
# neither vendor_adapters/__init__.py nor its adapter modules read conf.ini
# at import time - only their run() functions take conf.ini-derived
# overrides as plain arguments, resolved lazily inside run_vendor_adapter()
# below, well after set_config_path() has run. It needs to be available
# before argument parsing anyway, to build --pipeline's `choices` list.

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
    path: str | None = None,
    new_path: str | None = None,
    name: str = "",
    sample_list: str | None = None,
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
    path : Optional[str], optional
        Original study folder path, required when update/extract/remove is True.
    new_path : Optional[str], optional
        Incoming data folder path, required when update is True.
    name : str, optional
        New study name/rename, used by update/extract/remove.
    sample_list : Optional[str], optional
        Path to the sample ID list, required when extract/remove is True.

    Returns
    -------
    None

    """
    if not any([update, extract, remove]):
        start_time = datetime.now().astimezone().strftime("%d/%m/%Y, %H:%M:%S")

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
            start_time,
            analysis_type,
            )


    ############################
    #         UPDATE           #
    ############################

    if update:
        logger.info("Starting Update study")
        update_main(path, new_path, output_folder, name, overwrite_output)


    ############################
    #         DELETE           #
    ############################

    if remove:
        logger.info("Starting Delete sample(s) from study")
        delete_main(path, sample_list, output_folder, name, overwrite_output)


    ############################
    #         EXTRACT          #
    ############################

    if extract:
        logger.info("Starting Extract sample(s) from study")
        extract_main(path, sample_list, output_folder, name, overwrite_output)


#################################################################################################################

def resolve_pipeline(cli_pipeline: str | None) -> str:
    """Decide which vendor pipeline this run uses.

    Precedence, matching how -C/--config already overrides the default
    conf.ini path elsewhere in this codebase: an explicit `--pipeline` CLI
    value always wins; otherwise fall back to conf.ini's `[Vendor]
    PIPELINE` (meant for a future GUI to render as a dropdown of known,
    registered options); otherwise default to "native" - today's existing,
    unchanged behavior (-i already points at Varan's own sample.tsv/
    patient.tsv/[fusion.tsv] files, or an already-shaped SNV/CNV/
    CombinedOutput input folder).

    Deliberately NOT auto-detected from the input file's shape/content:
    which vendor produced a given raw input is a consequential
    classification (it picks an entire parsing code path), so it's made an
    explicit choice here rather than an implicit guess - the same
    reasoning already applied to medullo/ewing's do_mva handling elsewhere
    in this codebase. See run_vendor_adapter()'s handling of an adapter
    finding no data for the supplementary sanity-check heuristic (catches
    a human picking the wrong vendor, without making that guess the actual
    dispatch mechanism).

    Args:
        cli_pipeline (str | None): The parsed `--pipeline` CLI value, or
            None if the flag wasn't given.

    Returns:
        str: "native", or one of vendor_adapters.ADAPTERS's registered
            vendor names.

    """
    if cli_pipeline:
        return cli_pipeline

    conf_pipeline = get_config().get("Vendor", "PIPELINE", fallback="").strip()
    if conf_pipeline:
        return conf_pipeline

    return "native"


def run_vendor_adapter(
    pipeline: str, varan_input: Sequence[str], output_folder: str,
) -> tuple[list[str], str]:
    """Convert a vendor's raw run input into Varan's own sample.tsv/
    patient.tsv/fusion.tsv shape, then return a `varan_input`-shaped list
    so the rest of `varan.py` can proceed exactly as if `-i` had been given
    those files directly - see MULTIVENDOR_INTEGRATION_NOTES.md for the
    full design rationale.

    `varan_input[0]` is reinterpreted as the vendor's own raw input (e.g.
    an S3 run folder for the Guardant adapter) rather than an already-built
    sample.tsv, and is passed straight through to the adapter's own
    `run(folder=...)` - unchanged from vendor_adapters' existing interface,
    no new input shape invented here. `varan_input[1]` (patient.tsv), if
    given, is passed through as-is (no vendor adapter currently generates
    its own patient.tsv). A third element (a fusion file) is not supported
    in vendor mode - the adapter always generates its own - and is ignored
    with a warning if present, rather than silently misused.

    Intermediate files (the generated sample.tsv/fusions.tsv, and for
    Guardant, the converted per-sample VCFs) are written under this run's
    own `<output_folder>/scratch/<random>/` folder - the same namespaced
    scratch convention `create_random_name_folder`/`clear_scratch` already
    use elsewhere in walk.py for VEP's temp files - and removed again once
    `varan()` has finished reading them (that copy-out-of-scratch step
    happens inside `transform_input()`/`check_folders()`, called from
    `_walk_setup()`, so nothing downstream depends on scratch surviving
    past that point).

    Args:
        pipeline (str): A registered vendor name (never "native" - callers
            should only call this when a vendor pipeline was selected).
        varan_input (Sequence[str]): The raw `-i` argument list.
        output_folder (str): This run's (not-yet-versioned) output folder.

    Returns:
        tuple[list[str], str]: A 3-element [sample_tsv, patient_tsv,
            fusion_tsv] list, ready to hand to `varan()` exactly like a
            native `-i` input (empty strings mean "not applicable",
            matching the existing file-mode convention - see
            walk.input_extraction_file) - together with the scratch
            folder path, so the caller can remove it once `varan()` has
            finished with it.

    """
    raw_input = varan_input[0]
    patient_passthrough = varan_input[1].strip() if len(varan_input) > 1 else ""
    if len(varan_input) > 2 and varan_input[2].strip():
        logger.warning(
            "A third -i argument (fusion file) was given together with "
            f"--pipeline {pipeline}, but vendor mode always generates its "
            "own fusion file from the raw input - the one you passed is "
            "being ignored.")

    scratch_dir = create_random_name_folder(output_folder)

    adapter = ADAPTERS[pipeline]
    config = get_config()

    # Defaults route every intermediate file into this run's own scratch
    # folder, so a normal run leaves nothing behind in a shared/production
    # location. A deployment that wants persistent paths instead (e.g. to
    # keep the converted VCFs around, or to point at a non-default
    # dict.csv) can override any of these - or any other keyword
    # run() accepts - via a [Vendor.<pipeline>] conf.ini section, e.g.:
    #   [Vendor.guardant]
    #   dict_path = /path/to/dict.csv
    scratch_defaults = {
        "report_base_dir": str(Path(scratch_dir) / "sample_tsv"),
        "vcf_base_dir": str(Path(scratch_dir) / "vcf"),
        "temp_local_dir": str(Path(scratch_dir) / "tmp"),
    }
    section = f"Vendor.{pipeline}"
    overrides = dict(config.items(section)) if config.has_section(section) else {}
    run_kwargs = {**scratch_defaults, **overrides}

    logger.info(f"Running '{pipeline}' vendor adapter against '{raw_input}'...")
    result = adapter.run(folder=raw_input, **run_kwargs)

    if result is None:
        clear_scratch(scratch_dir)
        msg = (
            f"The '{pipeline}' vendor adapter found no usable data in "
            f"'{raw_input}'. Check that this is really a {pipeline} run "
            "folder/selection file, and that --pipeline (or conf.ini's "
            "[Vendor] PIPELINE) matches the vendor that actually produced "
            "this input - this is not auto-detected.")
        raise ValueError(msg)

    sample_tsv = str(result["report_path"])
    fusion_tsv_raw = result.get("fusion_table_path")
    fusion_tsv = (
        str(fusion_tsv_raw)
        if fusion_tsv_raw and Path(fusion_tsv_raw).exists() else "")

    logger.success(
        f"'{pipeline}' adapter produced {sample_tsv}"
        + (f" and {fusion_tsv}" if fusion_tsv else " (no fusions file)"))

    return [sample_tsv, patient_passthrough, fusion_tsv], scratch_dir

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

    __version__ = get_git_version()

    parser = MyArgumentParser(
        add_help=True,
        exit_on_error=False,
        usage=None,
        description="Argument of Varan script")

    # VERSION BLOCK
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"Varan {__version__}"
    )

    # WALK BLOCK
    parser.add_argument("-c", "--Cancer", required=False,
                        help="Cancer Name")

    parser.add_argument("-i", "--varan_input", nargs="+", required=False, type=str,
    help=("list with 1) input folder/sample file tsv (required) "
    "2) patient tsv 3) fusion file. If --pipeline selects a vendor "
    "(not 'native'), 1) is instead that vendor's raw run input (e.g. an "
    "S3 folder) and 3) is not used - the vendor adapter generates its "
    "own fusion file."))

    parser.add_argument(
        "--pipeline", required=False, default=None,
        choices=["native", *sorted(ADAPTERS)],
        help=(
            "Which vendor produced -i's input. 'native' (the default if "
            "neither this flag nor conf.ini's [Vendor] PIPELINE is set) "
            "means -i already points at Varan's own sample.tsv/patient.tsv/"
            "[fusion.tsv] or an already-shaped SNV/CNV/CombinedOutput "
            "folder - today's unchanged behavior. Any other value runs "
            "that vendor's preprocessing adapter (vendor_adapters/<name>.py) "
            "against -i's raw input automatically before the rest of the "
            "pipeline, so you no longer need to run create_Varan_input.py "
            "yourself first. This is an explicit choice, not "
            "auto-detected from the input's content - a CLI value here "
            "always overrides conf.ini's [Vendor] PIPELINE."))

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
                        help=("Select filter for SNV [d -> drop rows with ALT=='.' or "
                        "FILTER!='PASS' (on the raw VCF), p -> filter==PASS (on the MAF), "
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

    # CONFIG BLOCK
    parser.add_argument(
        "-C", "--config", required=False, default="conf.ini",
        help="Path to the conf.ini file to use for this run (default: ./conf.ini)")

    try:
        args = parser.parse_args()

        # Set the conf.ini path *before* importing any other Varan module -
        # every one of them reads conf.ini at import time, so this ordering
        # is what makes --config actually take effect instead of being
        # silently ignored. See config_loader.py for the full explanation.
        set_config_path(args.config)

        from concatenate import concatenate_main
        from Delete_script import delete_main
        from ExtractSamples_script import extract_main
        from filter_clinvar import filter_main
        from Make_meta_and_cases import meta_case_main
        from Update_script import update_main
        from ValidateFolder import validate_output
        from walk import clear_scratch, create_random_name_folder, walk_folder

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

        # VENDOR PIPELINE DISPATCH
        #
        # Explicit selection only (--pipeline, falling back to conf.ini's
        # [Vendor] PIPELINE, falling back to "native") - never auto-detected
        # from -i's content. See resolve_pipeline()'s docstring for why.
        pipeline = resolve_pipeline(args.pipeline)
        if pipeline not in {"native", *ADAPTERS}:
            logger.critical(
                f"conf.ini's [Vendor] PIPELINE is '{pipeline}', which isn't "
                f"a registered vendor (known: native, "
                f"{', '.join(sorted(ADAPTERS))}). Fix conf.ini, or pass "
                "--pipeline explicitly.")
            sys.exit(1)

        scratch_dir = None
        if pipeline != "native" and not any([update, extract, remove]):
            varan_input, scratch_dir = run_vendor_adapter(
                pipeline, varan_input, output_folder)

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
            remove,
            args.Path,
            args.NewPath,
            args.Name,
            args.SampleList)

        # Only clean up on success, matching the existing scratch-folder
        # convention elsewhere (walk._walk_process_snv's own VEP scratch
        # folder is likewise left in place if an exception interrupts the
        # run, so a failed run can be inspected before retrying).
        if scratch_dir:
            clear_scratch(scratch_dir)

    except ValueError as err:
        logger.critical(f"ValueError: {err}", file=sys.stderr)
    except FileNotFoundError as err:
        logger.critical(f"File not found: {err}", file=sys.stderr)
