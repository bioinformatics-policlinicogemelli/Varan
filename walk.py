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

"""Walk stage entry points: input handling, setup and the whole-walk driver.

_walk_setup() resolves the input (sample.tsv or input folder, vendor
adapter output included), creates/resumes the versioned study folder and
returns the WalkContext every stage runs from; walk_folder() then runs the
four stages, each in its own module and independent of the others:

- walk_cnv.py: CNV calls -> data_cna*.txt / .seg, BRCA exon-level CNVs
- walk_snv.py: SNV calls -> per-sample vcf2maf/VEP MAFs (+ SigMA)
- walk_fusion.py: fusions and splice variants -> data_sv.txt
- walk_clinical.py: data_clinical_patient.txt / data_clinical_sample.txt

Shared pieces (WalkContext, conf.ini setting accessors, scratch folders)
live in walk_common.py. Everything is re-exported here, so callers keep
using `walk.<name>` (walk_stage.py runs each stage as its own Snakemake
job through it).
"""
from __future__ import annotations

import os
import shutil
import sys
import zipfile
from pathlib import Path

import pandas as pd
from loguru import logger

from config_loader import config
from filter_clinvar import check_bool
from populate_case_lists import write_profiled_samples
from versioning import get_newest_version, get_version_list

# Re-exported so `import walk` / `walk.<name>` keep working for every
# caller (walk_stage.py, pipeline_stage.py, varan.py) after the split.
from walk_common import (  # noqa: F401
    describe_unknown_oncotree_codes,
    unknown_oncotree_codes,
    get_table_from_folder,
    WalkContext,
    _path_setting,
    check_field_tsv,
    clear_scratch,
    create_random_name_folder,
    flatten,
    get_combined_variant_output_from_folder,
    output_filtered,
    ploidy,
    sample_type,
    threshold_msi_liquid,
)
from walk_cnv import (  # noqa: F401
    _walk_process_cnv,
    build_exon_df,
    check_data_cna,
    cnv_type_from_folder,
    extract_multiple_cnv,
    get_cnv_from_folder,
    get_sample_id_from_cnv,
    write_exon_brca,
    write_exon_brca_generic_assay,
)
from walk_snv import (  # noqa: F401
    _read_vcf_header_and_data,
    _run_sigma_for_snv_sample,
    _sigma_sample_rows_dir,
    _walk_process_snv,
    _walk_process_snv_sample,
    _walk_snv_finish,
    _walk_snv_samples,
    _write_sigma_intermediate,
    extract_multiple_snv,
    get_sample_id_from_snv,
    get_snv_from_folder,
    recover_colocated_variants,
    run_vcf2maf,
    snv_job_keys,
    snv_type_from_folder,
    split_colocated_variants,
    vcf2maf_constructor,
    vcf_filtering,
)
from walk_fusion import (  # noqa: F401
    _walk_process_fusion,
    annotate_fusion,
    fill_fusion_from_combined,
    fill_fusion_from_temp,
    fill_splice_from_combined,
)
from walk_clinical import (  # noqa: F401
    _is_blank_value,
    _resolve_biomarker_threshold,
    _walk_write_clinical_tables,
    add_header_patient_long,
    add_header_patient_short,
    add_header_patient_type,
    fill_from_combined,
    fill_from_file,
    update_data_clinical_with_exon_info,
    update_data_clinical_with_hrd_info,
    write_clinical_sample,
    write_default_clinical_patient,
)








































def create_folder(output_folder: str,
                  overwrite_output: bool, resume: bool) -> Path:
    """Create versioned output folder for storing results.

    Args:
        output_folder (str): Base output path.
        overwrite_output (bool): If True, delete existing version and create new.
        resume (bool): If True, resume from last existing version if available.

    Returns:
        Path: Path to the versioned output folder.

    """
    output_list = [
        Path(output_folder).parent / x
        for x in get_version_list(output_folder)
    ]

    if output_list and output_list[-1].exists():
        output = output_list[-1]
        logger.warning(
            f"It seems that a version of the folder '{output_folder}' already exists.")
        if overwrite_output:
            logger.info("Overwrite option set. Start removing folder")
            shutil.rmtree(output)
        elif resume:
            _,current_version = get_newest_version(output_folder)
            return Path(output_folder + current_version)

    if not output_list:
        version = "_v1"
        output_folder_version = Path(f"{output_folder}{version}")


    else:
        output_folder_version, _ = get_newest_version(output_folder)

    logger.info(f"Creating the output folder '{output_folder_version}'...")

    output_folder_version.mkdir(parents=True, exist_ok=True)

    maf_path = output_folder_version / "maf"
    maf_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"The folder '{output_folder_version}' was correctly created!")

    return output_folder_version






def check_multiple_file(input_file: str, multiple: bool) -> None:
    """Check sample.tsv and conf.ini for multiple sample input.

    Args:
        input_file (str): Path to the sample.tsv file.
        multiple (bool): Whether the input is in multiple-sample mode.

    Raises:
        Exception: If validation of configuration and sample.tsv fails.

    """
    conf_snv = config.get("Multiple", "SNV")
    conf_cnv = config.get("Multiple", "CNV")
    file_paths = pd.read_csv(input_file, sep="\t", header=0)
    snv_file_path = file_paths["snv_path"].isna().all()
    cnv_file_path = file_paths["cnv_path"].isna().all()

    #CASE 1: multiple = True; neither snv or cnv are filled in sample.tsv;
    #multiple snv or cnv are filled in conf.ini
    if (multiple
    and not (snv_file_path or cnv_file_path)
    and not (conf_snv == "" or conf_cnv == "")):
        logger.critical(
            "-m was selected and Muliple section in conf.ini was filled but the "
            "file doesn't looks like a multiVCF.")
        msg = "Input error"
        raise Exception(msg)
    #CASE 2: multiple = True; neither snv or cnv are filled in sample.tsv;
    #neither multiple snv or cnv are filled in conf.ini
    if (multiple
    and not (snv_file_path or cnv_file_path)
    and (conf_snv == "" or conf_cnv == "")):
        logger.critical("-m was selected but Muliple section in conf.ini wasn't "
        "filled and the file doesn't looks like a multiVCF.")
        msg = "Input error"
        raise Exception (msg)
    #CASE 3: multiple = False; snv or cnv are filled in sample.tsv;
    #multiple snv or cnv are filled in conf.ini
    if (not multiple
    and (snv_file_path or cnv_file_path)
    and not (conf_snv == "" or conf_cnv == "")):
        logger.critical("-m was not selected but both sample.tsv and Muliple "
        "section in conf.ini were filled.")
        msg = "Input error"
        raise Exception (msg)
    #CASE 4: multiple = False; snv or cnv are filled in sample.tsv;
    #neither multiple snv or cnv are filled in conf.ini
    if (not multiple
    and (snv_file_path or cnv_file_path)
    and (conf_snv == "" or conf_cnv == "")):
        logger.warning("SNV and/or CNV columns in sample.tsv were not filled.")


def check_multiple_folder(input_dir: str, multiple: bool) -> None:
    """Validate whether SNV/CNV VCFs are single or multi-sample based on -m option.

    Args:
        input_dir (str): Path to the base input directory.
        multiple (bool): Whether the input should be multi-sample.

    Raises:
        Exception: If actual file type does not match `multiple` expectation.

    """
    snv_mulitple, snv_single, cnv_mulitple, cnv_single = False, False, False, False
    snv_folder = Path(input_dir) / "SNV"
    try:
        vcf_files = (f.name for f in snv_folder.iterdir() if f.suffix == ".vcf")
        multiple_vcf_snv = next(vcf_files)
    except StopIteration:
        multiple_vcf_snv = ""
        snv_folder.mkdir(parents=True, exist_ok=True)

    snv_file = Path(input_dir) / "sample_id_snv.txt"
    cmd_snv = f"vcf-query -l {(snv_folder / multiple_vcf_snv)} > {snv_file}"
    os.system(cmd_snv)

    with snv_file.open() as file:
        lines = file.readlines()

    minimum_lines_for_multiple = 2

    if len(lines) >= minimum_lines_for_multiple and not multiple:
        snv_mulitple = True
        logger.error("-m option was not selected but the SNV file is multiple!")
    elif len(lines) < minimum_lines_for_multiple and multiple:
        snv_single = True
        logger.error("-m option was selected but the SNV file is not multiple!")

    snv_file.unlink()

    cnv_folder = Path(input_dir) / "CNV"
    try:
        vcf_files = (f.name for f in cnv_folder.iterdir() if f.suffix == ".vcf")
        multiple_vcf_cnv = next(vcf_files)
    except StopIteration:
        multiple_vcf_cnv = ""
        cnv_folder.mkdir(parents=True, exist_ok=True)

    cnv_file = Path(input_dir) / "sample_id_cnv.txt"
    cmd_cnv = f"vcf-query -l {(Path(cnv_folder) / multiple_vcf_cnv)} > {cnv_file}"
    os.system(cmd_cnv)

    with cnv_file.open() as file:
        lines = file.readlines()

    if len(lines) >= minimum_lines_for_multiple and not multiple:
        cnv_mulitple = True
        logger.error("-m option was not selected but the CNV file is multiple!")
    elif len(lines) < minimum_lines_for_multiple and multiple:
        cnv_single = True
        logger.error("-m option was selected but the CNV file is not multiple!")

    cnv_file.unlink()

    msg = "Input error"
    if cnv_mulitple or snv_mulitple:
        raise Exception(msg)
    if cnv_single or snv_single:
        raise Exception(msg)





















def check_input_file(
    output_folder: str, file: str, copy_to: str, sample_id: str,
    relevant: bool = True) -> bool:
    """Copy an input file to a temporary folder, or report it's missing.

    Args:
        output_folder (str): Base output directory.
        file (str): File path to copy.
        copy_to (str): Subfolder name.
        sample_id (str): ID of the sample (for logging).
        relevant (bool): Whether this file type is actually part of the
            selected analysis type (-t). When False, a missing/unset path
            is expected and not worth a warning - only an unexpected state
            (path given but the file doesn't exist) still gets logged.

    Returns:
        bool: True when `file` is unset for this sample and `relevant` is
            True - the caller aggregates this flag across the whole sample
            batch (see check_folders/transform_input) so a config-wide gap
            (e.g. every Guardant sample's blank comb_path) logs ONE warning
            instead of repeating the identical line once per sample. A
            path that was given but doesn't exist on disk is still warned
            about immediately here, since that message differs per sample.

    """
    file_path = Path(file)
    destination = Path(output_folder) / "temp" / copy_to

    # `file` empty (e.g. Guardant's always-blank comb_path) must be checked
    # BEFORE file_path.exists(): Path("") == Path(".") and the current
    # directory always exists, so file_path.exists() would be True and this
    # would try to `cp .` the whole cwd - producing a confusing "cp: -r not
    # specified" shell error instead of the intended "no path set" warning,
    # without actually copying anything.
    if not file:
        return relevant

    if file_path.exists():
        # Copied under the source file's own original name - no longer
        # forced into a fixed SAMPLE_ID-based filename. That rename (once
        # here, for get_sample_id_from_snv()/get_sample_id_from_cnv() to
        # recover the sample id from a fixed suffix downstream) hid the
        # file's real origin and wasn't needed for this pipeline's actual
        # inputs, which already carry the sample id in their own name.
        shutil.copy(file_path, destination)
    else:
        logger.warning(f"{file_path} not found")
    return False


def check_folders(
    output_folder: str,
    snv_path: str,
    cnv_path: str,
    combout: str,
    sample_id: str,
    vcf_type: str | None = None) -> dict[str, bool]:
    """Check presence of SNV, CNV, and CombinedOutput files for a sample.

    Only flags a file type that's actually missing or unset for a sample
    where it's needed - a missing/unset path for a file type excluded by
    the selected analysis type (-t) is expected, not a problem (see
    check_input_file's `relevant` flag).

    Args:
        output_folder (str): Base output directory.
        snv_path (str): Path to SNV file.
        cnv_path (str): Path to CNV file.
        combout (str): Path to CombinedOutput file.
        sample_id (str): ID of the sample.
        vcf_type (str | None): Selected analysis type restriction (-t), or
            None if no restriction was given.

    Returns:
        dict[str, bool]: For each of "SNV"/"CNV"/"CombinedOutput", True if
            this sample has no path set for it (and it's relevant) - the
            caller aggregates these across the batch (see transform_input)
            to log one warning per file type instead of one per sample.

    """
    return {
        "SNV": check_input_file(
            output_folder, snv_path, "SNV", sample_id,
            relevant=vcf_type not in ["cnv", "fus", "tab"]),
        "CNV": check_input_file(
            output_folder, cnv_path, "CNV", sample_id,
            relevant=vcf_type not in ["snv", "fus", "tab"]),
        "CombinedOutput": check_input_file(
            output_folder, combout, "CombinedOutput", sample_id),
    }


def transform_input(
    tsv: str,
    clin_pzt: str,
    fusion_tsv: str,
    output_folder: str,
    multiple: bool,
    vcf_type: str | None = None) -> str:
    """Prepare temp folder structure and copy necessary files.

    Args:
        tsv (str): Sample TSV file path.
        clin_pzt (str): Clinical patient TSV file path.
        fusion_tsv (str): Fusion TSV file path.
        output_folder (str): Output directory.
        multiple (bool): Whether input is from multi-sample VCF.
        vcf_type (str | None): Selected analysis type restriction (-t), used
            to silence expected-missing-file warnings for file types the
            current run doesn't need (see check_folders).

    Returns:
        str: Path to the temporary folder created.

    """
    base_temp = Path(output_folder) / "temp"
    (base_temp).mkdir(parents=True, exist_ok=True)
    (base_temp / "SNV").mkdir(parents=True, exist_ok=True)
    (base_temp / "CNV").mkdir(parents=True, exist_ok=True)
    (base_temp / "CombinedOutput").mkdir(parents=True, exist_ok=True)
    (base_temp / "FUSIONS").mkdir(parents=True, exist_ok=True)

    sample_path = Path(output_folder) / "temp" / "sample.tsv"
    patient_path = Path(output_folder) / "temp" / "patient.tsv"
    fusion_path = Path(output_folder) / "temp" / "FUSIONS" / "fusion.tsv"

    shutil.copy(tsv, sample_path)

    if clin_pzt != "":
        shutil.copy(clin_pzt, patient_path)

    if fusion_tsv != "":
        shutil.copy(fusion_tsv, fusion_path)

    if multiple:
        snv_path = config.get("Multiple", "SNV")
        cnv_path = config.get("Multiple", "CNV")
        combout = config.get("Multiple", "COMBOUT")

        missing = check_folders(
            output_folder, snv_path, cnv_path, combout, "multiple", vcf_type)
        for copy_to, was_missing in missing.items():
            if was_missing:
                logger.warning(
                    f"No final_path set in conf.ini for {copy_to}!")

    else:
        tsv_file = pd.read_csv(tsv, sep="\t", dtype="string", keep_default_na=False)

        # Collected per file type instead of warning inline per sample: a
        # config-wide gap (e.g. every Guardant sample's blank comb_path,
        # since Guardant never has a CombinedOutput) would otherwise log
        # the identical line once per sample - see check_input_file.
        missing_by_type: dict[str, list[str]] = {
            "SNV": [], "CNV": [], "CombinedOutput": []}
        for _,row in tsv_file.iterrows():
            sample_id = row["SAMPLE_ID"]
            snv_path = row["snv_path"]
            cnv_path = row["cnv_path"]
            combout = row["comb_path"]

            missing = check_folders(
                output_folder, snv_path, cnv_path, combout, sample_id, vcf_type)
            for copy_to, was_missing in missing.items():
                if was_missing:
                    missing_by_type[copy_to].append(sample_id)

        total_samples = len(tsv_file)
        for copy_to, missing_samples in missing_by_type.items():
            if not missing_samples:
                continue
            if len(missing_samples) == total_samples:
                logger.warning(
                    f"No final_path set in conf.ini for {copy_to} on any "
                    "of the samples!")
            else:
                for sample_id in missing_samples:
                    logger.warning(
                        f"No final_path set in conf.ini for sample "
                        f"{sample_id}'s {copy_to}!")

    return str(Path(output_folder) / "temp")




















def input_extraction_file(input_f: list) -> tuple:
    """Extract sample, patient, and fusion TSV file paths from input list.

    Args:
        input_f (list): List of file path strings (up to 3 elements).

    Returns:
        tuple: Paths to sample, patient, and fusion TSV files.

    """
    sample_tsv = input_f[0]
    patient_tsv, fusion_tsv = "",""

    limit = 2
    if len(input_f) > 1:
        patient_tsv = input_f[1].strip()
        if len(input_f) > limit:
            fusion_tsv=input_f[limit].strip()

    return sample_tsv, patient_tsv, fusion_tsv


def input_extraction_folder(input_path: str) -> tuple:
    """Extract sample, patient, and fusion TSV file paths from input folder.

    Args:
        input_path (str): Path to input folder.

    Returns:
        tuple: Paths to sample.tsv, patient.tsv, and fusion TSV files.

    """
    input_path = Path(input_path)
    patient_tsv, fusion_tsv = "",""
    sample_tsv = input_path / "sample.tsv"
    if (input_path / "patient.tsv").exists():
        patient_tsv = input_path / "patient.tsv"
    if (input_path / "FUSIONS" / "Fusions.tsv").exists():
        fusion_tsv = input_path / "FUSIONS" / "Fusions.tsv"
    return sample_tsv, patient_tsv, fusion_tsv


def validate_input(
    oncokb: bool,
    vcf_type: str | None,
    filters: str,
    cancer: str,
    input_path: str) -> None:
    """Validate the input configuration and required files for the pipeline.

    Args:
        oncokb (bool): Whether OncoKB annotation is enabled.
        vcf_type (str | None): VCF type or None.
        filters (str): Filter options as string.
        cancer (str): Cancer ID to validate against cBioPortal.
        input_path (str): Path to input file or folder.

    Raises:
        Exception: If required columns are missing or config fields are empty.
        SystemExit: If cancer ID is not recognized by cBioPortal.

    """
    isinputfile = Path(input_path).is_file()
    if not isinputfile:
        sample_path = Path(input_path) / "sample.tsv"
        file_tsv = pd.read_csv(sample_path, sep="\t", dtype=str)
    elif isinputfile:
        file_tsv = pd.read_csv(input_path, sep="\t", dtype=str)

    column_list = [
        "SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path",
        "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR"]
    if not all(name in file_tsv.columns for name in column_list):
        logger.critical('Required columns: "SAMPLE_ID", "PATIENT_ID",'
        '"ONCOTREE_CODE", "snv_path", "cnv_path", "comb_path", "MSI", '
        '"TMB", "MSI_THR", "TMB_THR"')
        msg = "The input file is missing some important columns!"
        raise(ValueError(msg))

    if oncokb and config.get("OncoKB", "ONCOKB") == "":
        msg = "oncokb option was set but ONCOKB field in conf.ini is empty!"
        raise ValueError(msg)

    unknown_codes = unknown_oncotree_codes(file_tsv)
    if unknown_codes:
        logger.warning(describe_unknown_oncotree_codes(unknown_codes))

    # Fail here, before the expensive VEP/vcf2maf step, rather than deep
    # inside filter_oncokb() after that work is already done: a conf.ini
    # that predates the SNV/CNV/FUSION split (a single ONCOKB_FILTER used
    # to cover all three) is otherwise only caught with a raw
    # configparser.NoOptionError once filtering actually runs.
    oncokb_filter_checks = []
    if oncokb and "o" in filters:
        oncokb_filter_checks.append(("Filters", "ONCOKB_FILTER_SNV"))
    if oncokb and "o" in filters and vcf_type not in ["snv", "fus", "tab"]:
        oncokb_filter_checks.append(("Cna", "ONCOKB_FILTER_CNV"))
    if oncokb and "o" in filters and vcf_type not in ["cnv", "snv", "tab"]:
        oncokb_filter_checks.append(("FUSION", "ONCOKB_FILTER_FUSION"))

    for section, key in oncokb_filter_checks:
        if not config.has_option(section, key):
            logger.critical(
                f"conf.ini is missing '{key}' under [{section}] - required "
                "for this run's OncoKB filtering. If this conf.ini predates "
                "the SNV/CNV/FUSION split (a single ONCOKB_FILTER used to "
                "cover all three), add ONCOKB_FILTER_SNV under [Filters], "
                "ONCOKB_FILTER_CNV under [Cna] and ONCOKB_FILTER_FUSION "
                "under [FUSION] - see the conf.ini template.")
            msg = f"Missing conf.ini option: [{section}] {key}"
            raise ValueError(msg)

    if sample_type() not in {"SOLID", "LIQUID"}:
        raise ValueError('Please select a sample type between "Solid" and "Liquid" in conf.ini.')

    if (vcf_type is None or "snv" in vcf_type) and (not _path_setting("VEP_PATH") or not _path_setting("VEP_DATA")):
        msg = "VEP_PATH and/or VEP_DATA field in conf.ini is empty!"
        raise ValueError(msg)

    if not _path_setting("REF_FASTA"):
        msg = "REF_FASTA field in conf.ini is empty!"
        raise ValueError(msg)

    if (vcf_type is None or "snv" in vcf_type) and not _path_setting("CLINV").strip():
        logger.warning(
            "[Paths] CLINV in conf.ini is empty: VEP won't annotate variants "
            "with ClinVar, so CLIN_SIG stays empty and the 'c' (clinical "
            "significance) filter has nothing to filter on.")

    if "o" in filters and not oncokb:
        logger.warning("OncoKB filter was selected in filters options but -k option "
        "was not set. This filtering will be ignored.")

    cancer_cbio = pd.read_csv("cancer_list.txt", sep="\t")
    cancer_cbio = cancer_cbio["TYPE_OF_CANCER_ID"].to_numpy().tolist()

    if cancer not in cancer_cbio:
        logger.critical(f"The cancer_id '{cancer}' is not recognize by cBioPortal. "
        "Check the cancer_list.txt to find the correct cancer id")
        sys.exit(1)














def _walk_setup(
    input_path: list,
    multiple: bool,
    output_folder: str,
    oncokb: bool,
    cancer: str,
    overwrite_output: bool,
    resume: bool,
    vcf_type: str | None,
    filters: str,
    sigma: bool = False) -> WalkContext:
    """Stage 0: resolve input, create/resume the output folder, and check
    which of SNV/CNV/CombinedOutput are actually present.

    Every other stage (_walk_process_cnv/_snv/_fusion,
    _walk_write_clinical_tables) depends only on this stage's output, not
    on each other's - see WalkContext and walk_folder().
    """
    logger.info("Starting walk_folder script:")
    logger.info(
        f"walk_folder args [input:{input_path}, output_folder:{output_folder}, "
        f"overwrite:{overwrite_output}, resume:{resume}, vcf_type:{vcf_type}, "
        f"filters:{filters}, multiple:{multiple}]")

    if not Path(input_path[0]).exists():
        msg = f"No valid file/folder {input_path} found. Check your input path"
        raise FileNotFoundError(msg)


    path = Path(input_path[0])
    isinputfile = path.is_file()
    validate_input(oncokb, vcf_type, filters, cancer, input_path[0])

    if not resume or not (Path(output_folder) / "temp").exists():
        output_folder = create_folder(output_folder, overwrite_output, resume)
    else:
        get_version_list(output_folder)

    if not isinputfile:
        input_folder = input_path[0]
        input_path, patient_tsv, fusion_tsv = input_extraction_folder(input_folder)
        check_multiple_folder(input_folder, multiple)

    elif isinputfile:
        input_path, patient_tsv, fusion_tsv = input_extraction_file(input_path)
        check_multiple_file(input_path, multiple)
        input_folder = transform_input(
            input_path, patient_tsv, fusion_tsv, output_folder, multiple, vcf_type)

    else:
        logger.critical(f"The input {input_path} isn't a file nor a folder")
        msg = "Exiting from walk script!"
        raise(FileNotFoundError(msg))

    img_path = Path("docs") / "img" / "logo_VARAN.png"
    if img_path.exists():
        img_output_dir = Path(output_folder) / "img"
        img_output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(img_path, img_output_dir / "logo_VARAN.png")

    input_folder = Path(input_folder)
    input_folder_snv = (input_folder / "SNV").resolve()
    input_folder_cnv = (input_folder / "CNV").resolve()
    input_folder_comb_out = (input_folder / "CombinedOutput").resolve()

    if not (
        any(input_folder_snv.iterdir()) or
        any(input_folder_cnv.iterdir()) or
        any(input_folder_comb_out.iterdir())):
        msg = ("Empty input folder(s)! Either SNV, CNV or CombinedOutput folder should"
        "contain input files.")
        raise ValueError(msg)

    maf_path = Path(output_folder) / "maf"
    maf_zip_path = Path(output_folder) / "maf.zip"
    clin_sample_path = Path(input_folder) / "sample.tsv"

    try:
        clin_file = pd.read_csv(clin_sample_path, sep="\t", dtype=str)
    except Exception as err:
        logger.critical(f"Something went wrong while reading {clin_sample_path}!")
        msg = "Error in reading the input file! Please check again."
        raise OSError(msg) from err

    # MAF/clinical consistency only makes sense when this run actually
    # processes SNV data (vcf_type is None, or explicitly "snv") - the same
    # gate _walk_process_snv itself uses. A CNV-only/fusion-only/tab-only
    # resume (-t cnv/fus/tab -R) never produces a MAF in the first place, so
    # checking for one here was both pointless and a crash risk (maf_samples
    # was left undefined - see below - if neither maf/ nor maf.zip existed).
    if resume and vcf_type not in ["cnv", "fus", "tab"]:
        maf_samples = None
        if maf_path.exists():
            try:
                maf_samples = {
                    entry.name
                    for entry in maf_path.iterdir()
                    if entry.is_file()}
            except Exception as err:
                logger.critical(
                    "Error accessing the maf folder. Please check its integrity.")
                msg = "Error reading the maf folder."
                raise Exception(msg) from err

        elif maf_zip_path.exists():
            try:
                with zipfile.ZipFile(maf_zip_path, "r") as zipped_maf:
                    zipped_maf.extractall(maf_path)
                maf_samples = {
                    entry.name
                    for entry in maf_path.iterdir()
                    if entry.is_file()}
            except zipfile.BadZipFile as err:
                logger.critical("Corrupted maf.zip file.")
                msg = ("Unable to extract the ZIP file. Try decompressing it manually "
                "and run Varan again.")
                raise Exception(msg) from err

        if maf_samples is None:
            logger.warning(
                "Resuming an SNV-inclusive run, but no maf/ folder or maf.zip "
                "was found to resume from - proceeding as if this were the "
                "first run for the MAF step.")
        else:
            clin_samples = set(clin_file["SAMPLE_ID"])
            clin_in_maf = all(
                any(clin_sample in maf_sample for maf_sample in maf_samples)
                for clin_sample in clin_samples)
            maf_in_clin = all(
                any(clin_sample in maf_sample for clin_sample in clin_samples)
                for maf_sample in maf_samples)

            if not (clin_in_maf and maf_in_clin) and len(maf_samples) != 0:
                logger.critical("It seems you are resuming an existing study with a "
                "different set of input samples. Please verify the sample consistency!")
                msg = "Sample mismatch detected."
                raise FileNotFoundError(msg)

        zip_maf = config.get("Zip", "ZIP_MAF")
        zip_maf = check_bool(zip_maf)
        if maf_zip_path.exists() and not zip_maf:
            maf_zip_path.unlink()

    input_folder_snv = Path(input_folder_snv)
    input_folder_cnv = Path(input_folder_cnv)

    snv_folder_empty = len(list(input_folder_snv.iterdir())) == 0
    cnv_folder_empty = len(list(input_folder_cnv.iterdir())) == 0

    if vcf_type is None:
        if snv_folder_empty and not cnv_folder_empty:
            vcf_type = "cnv"
            logger.info("SNV path was empty, the analysis will exclude SNV")
        elif cnv_folder_empty and not snv_folder_empty:
            vcf_type = "snv"
            logger.info("CNV path was empty, the analysis will exclude CNV")
        elif snv_folder_empty and cnv_folder_empty:
            # Both empty is a valid CombinedOutput-only run (fusions, splice
            # variants, MSI/TMB all come from CombinedOutput, not SNV/CNV
            # VCFs). Leaving vcf_type as None here is required: setting it to
            # "cnv" or "snv" (as a naive single `if` chain would, since the
            # SNV check would win) makes _walk_process_fusion's
            # `vcf_type in ["cnv", "snv", "tab"]` gate skip fusion/splice
            # processing entirely, silently dropping data_sv.txt.
            logger.info(
                "SNV and CNV paths were both empty; proceeding with "
                "CombinedOutput-only processing.")

    case_folder_arr_cnv = None
    if input_folder_cnv.exists() and vcf_type not in ["snv", "fus", "tab"]:
        if multiple:
            multivcf = next(
                (i for i in input_folder_cnv.iterdir() if i.name.endswith(".vcf")),
                None)
            if multivcf:
                # Split inside this study's temp/ (removed at the end of the
                # run), never inside the user's own input folder.
                work_dir = Path(output_folder) / "temp" / "CNV"
                extract_multiple_cnv(multivcf, work_dir)
                input_folder_cnv = work_dir / "single_sample_vcf"
        logger.info("Checking CNV files...")
        case_folder_arr_cnv = get_cnv_from_folder(input_folder_cnv)
        logger.info("Everything ok!")

    input_folder_snv = Path(input_folder_snv)

    if input_folder_snv.exists() and vcf_type not in ["cnv", "fus", "tab"]:
        if multiple:
            multivcf = next(f for f in input_folder_snv.iterdir() if f.suffix == ".vcf")
            work_dir = Path(output_folder) / "temp" / "SNV"
            extract_multiple_snv(multivcf, work_dir)
            input_folder_snv = work_dir / "single_sample_vcf"

        logger.info("Checking SNV files...")

    case_folder_arr = get_snv_from_folder(input_folder_snv)
    logger.info("Everything ok!")

    ctx = WalkContext(
        output_folder=output_folder,
        input_folder=Path(input_folder),
        input_path=input_path,
        patient_tsv=patient_tsv,
        fusion_tsv=fusion_tsv,
        input_folder_snv=Path(input_folder_snv),
        input_folder_cnv=Path(input_folder_cnv),
        case_folder_arr=case_folder_arr,
        case_folder_arr_cnv=case_folder_arr_cnv,
        clin_file=clin_file,
        clin_sample_path=clin_sample_path,
        multiple=multiple,
        oncokb=oncokb,
        cancer=cancer,
        filters=filters,
        vcf_type=vcf_type,
        resume=resume,
        isinputfile=isinputfile,
        sigma=sigma,
    )
    ctx.snv_sample_ids = sample_ids_by_file(ctx, "snv_path", "SNV")
    ctx.cnv_sample_ids = sample_ids_by_file(ctx, "cnv_path", "CNV")
    write_profiled_samples(ctx.output_folder, _profiled_samples(ctx))
    return ctx


def _profiled_samples(ctx: WalkContext) -> dict[str, list[str]]:
    """Which of sample.tsv's SAMPLE_IDs this run tests for each data type,
    whether or not they turn out to have any event - cBioPortal's
    "profiled" samples, listed in the case lists (see
    populate_case_lists.py). A sample counts as tested for a type only
    when its input for that type is really there:

    - sequenced / cna: its SNV / CNV file (snv_path / cnv_path, or - in
      folder input with those columns blank - a file in SNV/ / CNV/ named
      after the SAMPLE_ID); every sample with -m/--multiple's multi-sample
      VCF.
    - sv: its CombinedVariantOutput has a [Fusions]/[Data Fusions]/
      [Splice Variants] section that isn't "NA" (TSO500 writes NA when RNA
      wasn't run). Samples only listed in a fusion.tsv are left to
      data_sv.txt itself: that file lists positives only.

    sample.tsv's optional cnv_tested / fusion_tested columns ("yes"/"no",
    written by vendor adapters - see vendor_adapters/common.py) override
    the cna / sv rules above for that sample when filled: e.g. Guardant
    always writes a CNV VCF, even without a CNV report, and its fusions
    file only lists the positives.

    A data type this run doesn't process at all (-t, or auto-detected from
    empty input folders) gets no entry.
    """
    clin = ctx.clin_file.fillna("")
    sample_ids = [str(s).strip() for s in clin["SAMPLE_ID"] if str(s).strip()]
    rows = list(clin.to_dict("records"))
    profiled: dict[str, list[str]] = {}

    def has_input(column: str, folder: str, flag_column: str = "") -> list[str]:
        type_dir = Path(ctx.input_folder) / folder
        if ctx.multiple and type_dir.is_dir() and any(f.is_file() for f in type_dir.iterdir()):
            return sample_ids
        by_sample = input_files_by_sample(ctx, column, folder)
        tested = []
        for row in rows:
            sample_id = str(row.get("SAMPLE_ID", "")).strip()
            if not sample_id:
                continue
            flag = _tested_flag(row, flag_column) if flag_column else None
            if flag is not None:
                if flag:
                    tested.append(sample_id)
                continue
            if sample_id in by_sample:
                tested.append(sample_id)
        return tested

    if ctx.vcf_type not in ["cnv", "fus", "tab"]:
        profiled["sequenced"] = has_input("snv_path", "SNV")
    if ctx.vcf_type not in ["snv", "fus", "tab"]:
        profiled["cna"] = _with_cna_evaluable(
            ctx, has_input("cnv_path", "CNV", "cnv_tested"))
    if ctx.vcf_type not in ["cnv", "snv", "tab"]:
        comb_dir = Path(ctx.input_folder) / "CombinedOutput"
        files = [f for f in comb_dir.iterdir() if f.is_file()] if comb_dir.is_dir() else []
        tested = []
        for row in rows:
            sample_id = str(row.get("SAMPLE_ID", "")).strip()
            if not sample_id:
                continue
            flag = _tested_flag(row, "fusion_tested")
            if flag is not None:
                if flag:
                    tested.append(sample_id)
                continue
            value = str(row.get("comb_path", "")).strip()
            path = _resolve_input_file(value, comb_dir)
            if path is None and not value:
                # Folder input's documented naming is per patient
                # (<PATIENT_ID>_CombinedVariantOutput.tsv - the one
                # get_combined_variant_output_from_folder() reads MSI/TMB
                # from); a file named after the sample is accepted too.
                by_patient = comb_dir / (
                    f"{str(row.get('PATIENT_ID', '')).strip()}_CombinedVariantOutput.tsv")
                path = by_patient if by_patient.is_file() else next(
                    (f for f in files if _named_after(sample_id, f.name)), None)
            if path is not None and _tests_fusions(path):
                tested.append(sample_id)
        profiled["sv"] = tested
    return profiled


def input_files_by_sample(ctx: WalkContext, column: str, folder: str) -> dict[str, Path]:
    """SAMPLE_ID -> that sample's own input file of one type: sample.tsv's
    path column (column, e.g. "snv_path"), or - folder input with that
    column blank - the file in <input>/<folder>/ named after the sample.
    The one place a sample's identity is tied to its input files."""
    type_dir = Path(ctx.input_folder) / folder
    files = [f for f in type_dir.iterdir() if f.is_file()] if type_dir.is_dir() else []
    by_sample = {}
    for row in ctx.clin_file.fillna("").to_dict("records"):
        sample_id = str(row.get("SAMPLE_ID", "")).strip()
        if not sample_id:
            continue
        value = str(row.get(column, "")).strip()
        path = _resolve_input_file(value, type_dir)
        if path is None and not value:
            path = next((f for f in files if _named_after(sample_id, f.name)), None)
        if path is not None:
            by_sample[sample_id] = path
    return by_sample


def sample_ids_by_file(ctx: WalkContext, column: str, folder: str) -> dict[str, str]:
    """Input file name -> SAMPLE_ID for one data type (see
    input_files_by_sample()), so the SNV/CNV stages label each file with
    the sample.tsv SAMPLE_ID it belongs to instead of guessing it from the
    file name (e.g. a DRAGEN library id). Empty with -m/--multiple (one
    multi-sample VCF); a file name shared by several samples is left out."""
    if ctx.multiple:
        return {}
    by_sample = input_files_by_sample(ctx, column, folder)
    names = [p.name for p in by_sample.values()]
    shared = {n for n in names if names.count(n) > 1}
    if shared:
        logger.warning(
            f"{folder} input file name(s) shared by several samples: "
            f"{', '.join(sorted(shared))} - their sample ids are taken from "
            "the file names instead.")
    return {p.name: s for s, p in by_sample.items() if p.name not in shared}


def _with_cna_evaluable(ctx: WalkContext, sample_ids: list[str]) -> list[str]:
    """Drop the CNV samples whose discrete CNA calls won't be computed: with
    [Cna] CNVKIT_algorithm = True the calls need the sample's tumor content
    (sample.tsv TC) and a sample without one is skipped - listing it as
    tested would tell cBioPortal "no alteration" instead of "not evaluated"."""
    if not check_bool(config.get("Cna", "CNVKIT_algorithm", fallback="False")):
        return sample_ids
    clin = ctx.clin_file.fillna("")
    tc = {str(r.get("SAMPLE_ID", "")).strip(): str(r.get("TC", "")).strip()
          for r in clin.to_dict("records")}
    missing = [s for s in sample_ids if not tc.get(s)]
    if missing:
        logger.warning(
            f"No TC (tumor content) in sample.tsv for {', '.join(missing)}: with "
            "[Cna] CNVKIT_algorithm = True their discrete CNA calls can't be "
            "computed, so they are not listed as tested for CNAs.")
    return [s for s in sample_ids if tc.get(s)]


def _tested_flag(row: dict, column: str) -> bool | None:
    """sample.tsv's cnv_tested / fusion_tested value for one row: True
    (yes), False (no), or None when blank/absent - then the file-based
    rules in _profiled_samples() decide."""
    value = str(row.get(column, "")).strip().lower()
    if value in ("yes", "true", "1"):
        return True
    if value in ("no", "false", "0"):
        return False
    return None


def _resolve_input_file(value: str, type_dir: Path) -> Path | None:
    """sample.tsv's path for one input file, as an existing file: as given,
    or - folder input - by its file name inside type_dir."""
    if not value:
        return None
    path = Path(value)
    if path.is_file():
        return path
    if (type_dir / path.name).is_file():
        return type_dir / path.name
    return None


def _named_after(sample_id: str, file_name: str) -> bool:
    """Whether file_name is sample_id's own file (S1_..., S1.vcf), not
    another sample's that merely starts the same way (S10_...)."""
    return file_name == sample_id or (
        file_name.startswith(sample_id) and file_name[len(sample_id)] in "._-")


def _tests_fusions(combined_output: Path) -> bool:
    """Whether a CombinedVariantOutput file reports a fusion/splice-variant
    analysis: a [Fusions]/[Data Fusions]/[Splice Variants] section whose
    content isn't "NA" - read the same way tsv.get_fusions() and
    tsv.get_splice_variants() read them."""
    try:
        lines = combined_output.read_text(errors="replace").splitlines()
    except OSError:
        return False
    for i, line in enumerate(lines):
        if any(h in line for h in ("[Fusions]", "[Data Fusions]", "[Splice Variants]")):
            if i + 2 < len(lines) and lines[i + 2].strip() == "NA":
                continue
            if i + 1 < len(lines):
                return True
    return False
























def walk_folder(
    input_path: list,
    multiple: bool,
    output_folder: str,
    oncokb: bool,
    cancer: str,
    overwrite_output: bool = False,
    resume: bool = False,
    vcf_type: str | None = None,
    filters: str = "",
    sigma: bool = False) -> tuple:
    """Process input files/folders for SNV, CNV, fusions, and prepare output.

    Thin orchestrator over four independent stages - see WalkContext,
    _walk_setup, _walk_process_cnv, _walk_process_snv, _walk_process_fusion
    and _walk_write_clinical_tables. CNV, SNV, fusion/splice and the
    clinical tables each depend only on _walk_setup()'s output, never on
    each other, which is what makes them safe to expose as independent
    Snakemake rules (see Snakefile) instead of one opaque shell-out.

    Args:
        input_path (list): List with input path(s), either folder(s) or file(s).
        multiple (bool): Whether multiple samples per file are expected.
        output_folder (str): Path to the output directory.
        oncokb (bool): Whether to enable OncoKB annotation.
        cancer (str): Cancer ID used for validations and annotations.
        overwrite_output (bool, optional): Overwrite output folder if exists.
        Defaults to False.
        resume (bool, optional): Resume from previous run if True. Defaults to False.
        vcf_type (str | None, optional): Type of VCF to process (snv, cnv, fus, tab).
        Defaults to None.
        filters (str, optional): Filters to apply on VCF data. Defaults to "".
        sigma (bool, optional): Whether to run SigMA mutational-signature
        analysis per sample (see sigma_runner.run_sigma_for_sample). Defaults
        to False. All of SigMA's own parameters live in conf.ini's [SigMA]
        section - this is only the on/off switch.

    Returns:
        tuple: Returns output folder path, input path or file, and fusion TSV path.

    """
    ctx = _walk_setup(
        input_path, multiple, output_folder, oncokb, cancer,
        overwrite_output, resume, vcf_type, filters, sigma)

    _walk_process_cnv(ctx)
    _walk_process_snv(ctx)
    _walk_process_fusion(ctx)
    _walk_write_clinical_tables(ctx)

    logger.success("Walk script completed!\n")

    return ctx.output_folder, ctx.input_path, ctx.fusion_tsv
