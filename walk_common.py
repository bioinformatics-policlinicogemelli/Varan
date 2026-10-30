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

"""Shared pieces of the walk stage: WalkContext (the state setup hands to
every stage), conf.ini setting accessors, scratch folders and small helpers
used by several stages.
"""

from __future__ import annotations

import secrets
import shutil
import string
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from loguru import logger

# conf.ini settings, read when used (not at import - see config_loader.py).
from config_loader import config


def _path_setting(key: str) -> str:
    return config.get("Paths", key)


def sample_type() -> str:
    """[Sample_Type] TYPE, upper-cased: "SOLID" or "LIQUID"."""
    return config.get("Sample_Type", "TYPE").strip().strip('"').strip("'").upper()


def ploidy() -> int:
    return int(config.get("Cna", "PLOIDY"))


def threshold_msi_liquid() -> float:
    return float(config.get("MSI", "THRESHOLD_MSI_LIQUID"))


output_filtered = "snv_filtered"


def create_random_name_folder(output_folder: str) -> str:
    """Create a temporary scratch folder with a random name inside output_folder.

    Scratch lives at `output_folder/scratch/<random>` instead of a shared,
    cwd-relative `scratch/` folder, so each run's temp files are namespaced
    under its own output folder - no more collision risk when several
    varan.py processes run in parallel from the same working directory.

    Args:
        output_folder (str): The study's output folder for this run.

    Returns:
        str: Path to the created temporary folder.

    """
    folder_name = "".join(
        secrets.choice(string.ascii_lowercase + string.digits) for _ in range(10))
    temporary = Path(output_folder) / "scratch" / folder_name

    try:
        temporary.mkdir(parents=True)
    except Exception as err:
        logger.critical("Something went wrong while creating the vep tmp folder")
        msg = "Error in create_random_name_folder: exiting from walk script!"
        raise(Exception(msg)) from err
    return(str(temporary))


def clear_scratch(folder: str | None = None) -> None:
    """Remove this run's scratch subfolder, then the (now-empty) scratch/
    parent directory itself, so nothing scratch-related is left in the
    study's output folder.

    Only ever removes `<output_folder>/scratch/<this run's random name>`
    and then `<output_folder>/scratch` itself - never scans/clears a shared,
    cwd-relative `scratch/` directory, since that would risk deleting other
    runs' still-in-use temporary files when multiple varan.py processes run
    in parallel from the same working directory. Removing the parent here is
    still safe: it's namespaced under this run's own output_folder, and this
    is the only place a scratch subfolder is ever created for that folder.

    Args:
        folder (str | None): Path to this run's own scratch subfolder, or None
            if this run never created one (e.g. resume=True).

    Returns:
        None

    """
    if folder is None:
        return
    to_rem = Path(folder)
    if to_rem.exists():
        shutil.rmtree(to_rem)
    scratch_parent = to_rem.parent
    if scratch_parent.name == "scratch" and scratch_parent.exists():
        shutil.rmtree(scratch_parent, ignore_errors=True)
        # create_random_name_folder()'s mkdir(parents=True) may have had to
        # create the study's output folder itself just to hold this
        # scratch/ subfolder - e.g. vendor-adapter preprocessing
        # (run_vendor_adapter() in varan.py) runs against the not-yet-
        # versioned output folder name, before create_output_folder() ever
        # creates the real "<name>_v1" folder. Removing scratch/ above can
        # then leave that now-empty output folder behind. Only ever remove
        # it when it's completely empty, so a real output folder that
        # already has pipeline output in it (e.g. walk.py's own mid-run
        # scratch cleanup) is never touched.
        output_folder = scratch_parent.parent
        if output_folder.exists() and not any(output_folder.iterdir()):
            output_folder.rmdir()


def flatten(nested_list: list[list]) -> list:
    """Flatten a nested list into a single flat list.

    Args:
        nested_list (list[list]): A list of lists to flatten.

    Returns:
        list: A single flattened list.

    """
    return [item for sublist in nested_list for item in sublist]


def check_field_tsv(row: pd.Series, name: str) -> str:
    """Return a field from a TSV row. Exit if field is missing.

    Args:
        row (pd.Series): A row from a DataFrame.
        name (str): Name of the column to extract.

    Returns:
        str: The value in the specified column.

    """
    try:
        field=str(row[name])
    except KeyError as e:
        logger.critical(f"KeyError: {e} not found! Check if column name is correctly "
        "spelled or if there are tabs/spaces before or after the coloumn key: "
        f"\n{row.index}. \nThis error may also occur if the table columns have "
        "not been separated by tabs!")
        sys.exit(1)
    return field


def get_combined_variant_output_from_folder(
    input_folder: str, file: pd.DataFrame, isinputfile: bool) -> dict:
    """Create a dict mapping SAMPLE_ID to CombinedVariantOutput file paths.

    Args:
        input_folder (str): Base folder for default paths.
        file (pd.DataFrame): Input sample information.
        isinputfile (bool): Whether to use the path from the file or default.

    Returns:
        dict: Mapping from SAMPLE_ID to CombinedOutput file paths.

    """
    combined_dict = {}
    for _,row in file.iterrows():
        sample_id = check_field_tsv(row, "SAMPLE_ID")
        patient_id = check_field_tsv(row, "PATIENT_ID")
        if isinputfile:
            combined_path = Path(check_field_tsv(row, "comb_path"))
        else:
            combined_path = (
                Path(input_folder)
                / "CombinedOutput"
                / f"{patient_id}_CombinedVariantOutput.tsv")

        if combined_path.exists():
            pass
        else:
            logger.warning("comb_path in conf.ini does not exists")
        combined_dict[sample_id] = combined_path
    return combined_dict


@dataclass
class WalkContext:
    """Shared state produced once by _walk_setup() and read by the four
    independent stages below (_walk_process_cnv/_snv/_fusion and
    _walk_write_clinical_tables).

    None of those four stages depends on another's *output* - CNV, SNV,
    fusion/splice and the clinical tables are each derived only from this
    setup context - they are simply called in sequence by walk_folder() for
    simplicity and readable logs. That independence is exactly what lets
    them be exposed as separate Snakemake rules (see Snakefile) with real
    per-stage resume/parallelism, instead of the previous single
    ~400-line walk_folder() that Snakemake could only shell out to as one
    opaque unit.

    isinputfile: whether -i was a sample.tsv file (True) or an input folder.

    start_time, vendor_scratch_dir and old_version_exists are never read by
    walk.py itself: walk_stage.py fills them in after setup so
    pipeline_stage.py's later processes can pass the run's start time to
    the report, remove the vendor adapter's scratch folder once the run
    succeeds, and restore versioning.py's `old_version_exists` global -
    which is only set inside setup's process (by get_version_list()) but
    read by write_report.py during validation.
    analysis_type is the -t the user actually passed (None if none):
    vcf_type above may have been auto-detected by setup, but varan()
    gates filter/concatenate and the report's sections on the CLI value.
    command_line is the varan.py command equivalent to this run, shown as
    the report's "Command Line" (see walk_stage.py's setup).
    snv_sample_ids / cnv_sample_ids map each SNV / CNV input file name to
    its sample.tsv SAMPLE_ID (see walk.sample_ids_by_file()); files not in
    them fall back to the id derived from the file name.
    """

    output_folder: str
    input_folder: Path
    input_path: object
    patient_tsv: str
    fusion_tsv: str
    input_folder_snv: Path
    input_folder_cnv: Path
    case_folder_arr: dict
    case_folder_arr_cnv: dict | None
    clin_file: pd.DataFrame
    clin_sample_path: Path
    multiple: bool
    oncokb: bool
    cancer: str
    filters: str
    vcf_type: str | None
    resume: bool
    isinputfile: bool
    sigma: bool = False
    start_time: str = ""
    vendor_scratch_dir: str = ""
    old_version_exists: bool = False
    analysis_type: str | None = None
    command_line: list | None = None
    snv_sample_ids: dict = field(default_factory=dict)
    cnv_sample_ids: dict = field(default_factory=dict)


def get_table_from_folder(tsvpath: str) -> dict[str, list[str]]:
    """Extract sample-to-patient mapping from a TSV file.

    Args:
        tsvpath (str): Path to the input TSV file.

    Returns:
        dict[str, list[str]]: Mapping from SAMPLE_ID to list with
        corresponding PATIENT_ID.

    """
    table_dict = {}
    file = pd.read_csv(tsvpath, sep="\t", index_col=False, dtype=str)
    for _, row in file.iterrows():
        sample_id = str(row["SAMPLE_ID"])
        if ".bam" in sample_id:
           sample_id = sample_id.replace(".bam", "")
        if sample_id not in table_dict:
            table_dict[sample_id] = [str(row["PATIENT_ID"])]
    return table_dict
