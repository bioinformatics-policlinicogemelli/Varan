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

"""Walk stage: data_clinical_patient.txt / data_clinical_sample.txt (MSI,
TMB, exon-level and HRD/GIS information).
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

import specimen
import tsv
from clinical_file import set_sample_column
from config_loader import config
from thresholds import meets_threshold
from walk_cnv import build_exon_df
from walk_common import (
    flatten,
    get_combined_variant_output_from_folder,
    get_table_from_folder,
    WalkContext,
)


def write_clinical_sample(
    clin_samp_path: str,
    output_folder: str,
    table_dict: dict,
    combined_output_used: bool = False,
) -> None:
    """Write the `data_clinical_sample.txt` file by merging sample and metrics data.

    The function reads the clinical sample file, optionally merges it with MSI/TMB
    metrics from a provided table, validates headers, and writes the output file
    with appropriate headers as required by cBioPortal format.

    Args:
        clin_samp_path (str): Path to the input clinical sample TSV file.
        output_folder (str): Directory where the final file will be saved.
        table_dict (dict): Dictionary containing metrics (MSI, TMB) by sample ID.
        combined_output_used (bool): True if `table_dict` was built from a real
            CombinedOutput folder (fill_from_combined); False if it was built as a
            fallback straight from sample.tsv (fill_from_file). Used to log which
            source MSI/TMB actually came from, and to avoid comparing sample.tsv
            against itself when no CombinedOutput was ever provided.

    Raises:
        NameError: If expected columns or header definitions are missing or incorrect.
        KeyError: If input template does not contain required fields.

    """
    logger.info("Writing data_clinical_sample.txt file...")
    conf_header_short = config.get("ClinicalSample", "HEADER_SAMPLE_SHORT")
    conf_header_long = config.get("ClinicalSample", "HEADER_SAMPLE_LONG")
    conf_header_type = config.get("ClinicalSample", "HEADER_SAMPLE_TYPE")

    data_clin_samp = pd.read_csv(clin_samp_path, sep="\t", header=0, dtype=str)
    data_clin_samp["ONCOTREE_CODE"] = data_clin_samp["ONCOTREE_CODE"].str.upper()

    try:
        data_clin_samp = data_clin_samp.drop(
            columns=["snv_path", "cnv_path", "comb_path"])
    except KeyError as err:
        logger.critical("snv_path, cnv_path or comb_path columns were removed or "
        "modified from template. Please use the correct template!")
        msg = "Exiting from script!"
        raise NameError(msg) from err
    # Optional input-control column (see _profiled_samples()), not a
    # clinical attribute.
    data_clin_samp = data_clin_samp.drop(
        columns=["cnv_tested", "fusion_tested"], errors="ignore")

    data_clin_samp.columns = data_clin_samp.columns.str.upper()
    # Optional input column (matched white blood cells' VCF, see
    # clonal_hematopoiesis.py), not a clinical attribute.
    data_clin_samp = data_clin_samp.drop(columns=["WBC_PATH"], errors="ignore")

    combout_df = pd.DataFrame.from_dict(table_dict).transpose().reset_index()
    combout_df = combout_df.rename(columns={"index": "SAMPLE_ID", 0: "PATIENT_ID"})

    final_data_sample = data_clin_samp
    limit = 2

    if len(combout_df.columns) > limit:
        combout_df = combout_df.rename(columns={
            1: "MSI",
            2: "TMB",
            3: "MSI_THR",
            4:"TMB_THR"})
        source = "CombinedOutput" if combined_output_used else "sample.tsv"
        logger.info(
            f"MSI/TMB values written to data_clinical_sample.txt were taken from {source}.")

        try:
            if combined_output_used:
                msi_sample = pd.to_numeric(data_clin_samp["MSI"], errors="coerce")
                msi_combined = pd.to_numeric(combout_df["MSI"], errors="coerce")
                tmb_sample = pd.to_numeric(data_clin_samp["TMB"], errors="coerce")
                tmb_combined = pd.to_numeric(combout_df["TMB"], errors="coerce")

                msi_mismatch = ((msi_sample != msi_combined)
                                 & msi_sample.notna() & msi_combined.notna()).any()
                tmb_mismatch = ((tmb_sample != tmb_combined)
                                 & tmb_sample.notna() & tmb_combined.notna()).any()

                if msi_mismatch or tmb_mismatch:
                    logger.warning(
                        "MSI and/or TMB values reported in sample.tsv differ from "
                        "the ones computed from CombinedOutput! CombinedOutput "
                        "values were used, sample.tsv values were discarded.")
            try:
                data_clin_samp = data_clin_samp.drop(
                    columns=["MSI", "TMB", "MSI_THR", "TMB_THR"])
            except KeyError as err:
                logger.critical(
                    "MSI_THR or TMB_THR columns were removed or modified from template."
                    " Please use the correct template!")
                msg = "Exiting from script!"
                raise NameError(msg) from err
        except KeyError as err:
            logger.warning("No MSI or TMB columns found in template.")
            msg = "Exiting from script!"
            raise(KeyError(msg)) from err

        final_data_sample = data_clin_samp.merge(
            combout_df, on=["PATIENT_ID", "SAMPLE_ID"])

    # SigMA columns only ever appear when -g/--sigma was actually passed for
    # this run: sigma_runner/_walk_process_snv only ever writes this file
    # when ctx.sigma is True and at least one sample went through the SigMA
    # pipeline (see _write_sigma_intermediate). No file -> no columns added,
    # matching this codebase's existing convention of only showing/recording
    # a filter's info when that filter was actually applied (see
    # ONCOKB_FILTER_CNV/ONCOKB_FILTER_FUSION in write_report.py).
    sigma_path = Path(output_folder) / "intermediate" / "sigma" / "data_sigma.txt"
    if sigma_path.exists():
        sigma_df = pd.read_csv(sigma_path, sep="\t", dtype=str)
        if "SAMPLE_ID" in sigma_df.columns:
            final_data_sample = final_data_sample.merge(
                sigma_df, on="SAMPLE_ID", how="left")
            logger.info(
                f"Merged SigMA results for {len(sigma_df)} sample(s) into "
                "data_clinical_sample.txt.")
        else:
            logger.warning(
                f"{sigma_path} exists but has no SAMPLE_ID column - "
                "skipping SigMA column merge into data_clinical_sample.txt.")

    basic_columns = ["SAMPLE_ID", "PATIENT_ID", "MSI", "TMB", "MSI_THR", "TMB_THR"]

    # A pipeline that has no native TMB (or MSI) of its own - e.g. Guardant,
    # which reports MSI as a text call but never reports TMB at all - must
    # not be forced to ship a clinical column full of "NA" placeholders.
    # Only applies to the sample.tsv-driven path (no CombinedOutput source
    # to fall back on): native TSO500 CombinedOutput ingestion always has
    # real MSI/TMB values.
    if not combined_output_used:
        for marker in ("MSI", "TMB"):
            thr_col = f"{marker}_THR"
            if marker in final_data_sample.columns and (
                final_data_sample[marker].apply(_is_blank_value).all()
                and final_data_sample[thr_col].apply(_is_blank_value).all()):
                logger.info(
                    f"No {marker} value or {marker}_THR was provided for any "
                    f"sample in this run, and there is no CombinedOutput to "
                    f"compute it from - dropping {marker}/{thr_col} from "
                    "data_clinical_sample.txt entirely instead of writing NA.")
                basic_columns.remove(marker)
                basic_columns.remove(thr_col)
                final_data_sample = final_data_sample.drop(columns=[marker, thr_col])

    other_columns = [c for c in final_data_sample.columns
                      if c not in ("SAMPLE_ID", "PATIENT_ID", "MSI", "TMB",
                                   "MSI_THR", "TMB_THR")]
    new_cols = [*basic_columns, *other_columns]

    final_data_sample = final_data_sample[new_cols]
    dataclin_columns = list(final_data_sample.columns)

    # Add header's fifth row
    default_row = pd.DataFrame([dataclin_columns], columns=dataclin_columns)
    final_data_sample = pd.concat([default_row, final_data_sample], ignore_index=True)

    # Add header's fourth row (SERIES OF 1s)
    header_numbers = pd.DataFrame(
        [[1] * len(final_data_sample.columns)], columns=dataclin_columns)
    final_data_sample = pd.concat(
        [header_numbers, final_data_sample], ignore_index=True)

    # Add header's third row (HEADER_SAMPLE_TYPE)
    if not conf_header_type:
        header_row = ["STRING"] * len(dataclin_columns)
        for marker in ("MSI", "TMB", "SIGMA_TOTAL_SNVS", "SIGMA_SIGNATURE3_MVA"):
            if marker in dataclin_columns:
                header_row[dataclin_columns.index(marker)] = "NUMBER"
        for marker in ("SIGMA_DO_MVA",):
            if marker in dataclin_columns:
                header_row[dataclin_columns.index(marker)] = "BOOLEAN"
        sample_header_type = pd.DataFrame([header_row], columns=dataclin_columns)
    else:
        types_list = conf_header_type.split(",")
        types_list = [x.strip() for x in types_list]
        for types in types_list:
            if types.upper() not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{types} is not a valid type. Please check the given "
                "input in conf.ini. Valid types: STRING, NUMBER, BOOLEAN")
                msg = "The type is not valid: exiting from walk script!"
                raise(NameError(msg))
        try:
            types_list = [x.upper() for x in types_list]
            sample_header_type = pd.DataFrame([types_list], columns=dataclin_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(types_list)}) in "
            "HEADER_SAMPLE_TYPE is different from the effective number of columns "
            f"({len(final_data_sample.columns)}).")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    final_data_sample = pd.concat(
        [sample_header_type, final_data_sample], ignore_index=True)

    # Add header's second row (HEADER_SAMPLE_LONG)
    if not conf_header_long:
        sample_header_long = default_row
    else:
        try:
            combined_headers = conf_header_long.split(",")
            sample_header_long = pd.DataFrame(
                [combined_headers], columns=dataclin_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(combined_headers)}) in "
            "HEADER_SAMPLE_LONG in conf.ini is different from the effective number of "
            f"columns ({len(final_data_sample.columns)}).")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err
    final_data_sample = pd.concat(
        [sample_header_long, final_data_sample], ignore_index=True)

    # Add header's first row (HEADER_SAMPLE_SHORT)
    if not conf_header_short:
        sample_header_short = default_row
    else:
        try:
            combined_headers = conf_header_short.split(",")
            sample_header_short = pd.DataFrame(
                [combined_headers], columns=dataclin_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(combined_headers)}) in "
            "HEADER_SAMPLE_SHORT in conf.ini is different from the effective number of "
            f"columns ({len(final_data_sample.columns)}).")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err
    final_data_sample = pd.concat(
        [sample_header_short, final_data_sample], ignore_index=True)
    final_data_sample.loc[0:3, "SAMPLE_ID"] = (
        final_data_sample.loc[0:3, "SAMPLE_ID"].apply(lambda x: f"#{x}"))

    output_folder = Path(output_folder)
    data_clin_txt = output_folder / "data_clinical_sample.txt"
    final_data_sample.to_csv(data_clin_txt, sep="\t", index=False, header=False)


def write_specimen_type(output_folder: str, types: dict[str, str]) -> None:
    """Write every sample's specimen type (Tissue/Plasma, see specimen.py)
    to data_clinical_sample.txt's SPECIMEN_TYPE column - sample.tsv's own
    SPECIMEN_TYPE column, if it had one, comes through with its raw values
    (blank, "Solid"...), so it's rewritten here too."""
    default = specimen.study_specimen_type()
    set_sample_column(
        output_folder, specimen.COLUMN, types, "Specimen Type",
        "Specimen type of the sequenced sample (Tissue or Plasma)",
        default=default)


def write_default_clinical_patient(output_folder: str, table_dict: dict) -> None:
    """Write a default `data_clinical_patient.txt` file.

    Args:
        output_folder (str): Path to the output directory.
        table_dict (dict): Dictionary containing sample info mapped to patient IDs.

    """
    logger.info("Writing data_clinical_patient.txt file...")

    output_folder = Path(output_folder)
    data_clin_samp = output_folder / "data_clinical_patient.txt"

    with data_clin_samp.open("w", encoding="utf-8") as cil_sample:
        cil_sample.write("#Patient Identifier\tAge\tGender\n")
        cil_sample.write("#Patient identifier\tAge\tGender\n")
        cil_sample.write("#STRING\tNUMBER\tSTRING\n")
        cil_sample.write("#1\t1\t1\n")
        cil_sample.write("PATIENT_ID\tAGE\tGENDER\n")

        nested_list = list(table_dict.values())
        list_patients = sorted(set(flatten(nested_list)))

        for v in list_patients:
            cil_sample.write(f"{v}\tNaN\tNaN\n")


def add_header_patient_type(
    patient_tsv: str,
    datapat_columns: list[str],
    conf_header_type: str,
    final_data_pat: pd.DataFrame) -> pd.DataFrame:
    """Add HEADER_PATIENT_TYPE row to the patient DataFrame.

    Args:
        patient_tsv (str): File name used in log messages for context.
        datapat_columns (list[str]): List of column names.
        conf_header_type (str): Comma-separated types from config.
        final_data_pat (pd.DataFrame): The data to which headers are appended.

    Returns:
        pd.DataFrame: DataFrame with the header row prepended.

    Raises:
        NameError: If types are invalid or column count mismatches.

    """
    if not conf_header_type:
        def_type = ["STRING", "NUMBER", "STRING"]
        header_type = def_type + ["STRING"] * (len(datapat_columns)-3)
        header_type_df = pd.DataFrame([header_type], columns=datapat_columns)
    else:
        types_list = conf_header_type.split(",")
        types_list = [x.strip().upper() for x in types_list]
        for types in types_list:
            if types not in ["STRING", "BOOLEAN", "NUMBER"]:
                logger.critical(f"{types} is not a valid type. Please check the given "
                "input in conf.ini. Valid types: STRING, NUMBER, BOOLEAN")
                msg = "The type is not valid: exiting from walk script!"
                raise(NameError(msg))
        try:
            header_type_df = pd.DataFrame([types_list], columns=datapat_columns)
        except ValueError as err:
            logger.critical(f"The number of column names ({len(types_list)}) in "
            "HEADER_PATIENT_TYPE is different from the effective number of "
            f"columns ({len(datapat_columns)}) in {patient_tsv}.")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    return pd.concat([header_type_df, final_data_pat], ignore_index=True)


def add_header_patient_short(
    patient_tsv: str,
    datapat_columns: list[str],
    conf_header_short: str,
    default_row: pd.DataFrame,
    final_data_pat: pd.DataFrame) -> pd.DataFrame:
    """Add HEADER_PATIENT_SHORT row to the patient DataFrame.

    Args:
        patient_tsv (str): File name for context in logging.
        datapat_columns (list[str]): List of column names.
        conf_header_short (str): Header from config, comma-separated.
        default_row (pd.DataFrame): Default row to use if config not provided.
        final_data_pat (pd.DataFrame): The target DataFrame.

    Returns:
        pd.DataFrame: Updated DataFrame.

    """
    if not conf_header_short:
        pat_header_short = default_row
    else:
        try:
            pat_header_short = pd.DataFrame(
                [conf_header_short.split(", ")], columns=datapat_columns)
        except ValueError as err:
            logger.critical(f"The number of column names "
            f"({len(conf_header_short.split(', '))}) in HEADER_PATIENT_SHORT "
            "in conf.ini is different from the effective number of columns "
            f"({len(datapat_columns)}) in {patient_tsv}.")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    return pd.concat([pat_header_short, final_data_pat], ignore_index=True)


def add_header_patient_long(
    patient_tsv: str,
    datapat_columns: list[str],
    conf_header_long: str,
    default_row: pd.DataFrame,
    final_data_pat: pd.DataFrame) -> pd.DataFrame:
    """Add HEADER_PATIENT_LONG row to the patient DataFrame.

    Args:
        patient_tsv (str): File name for context.
        datapat_columns (list[str]): List of column names.
        conf_header_long (str): Header from config, comma-separated.
        default_row (pd.DataFrame): Default row if config missing.
        final_data_pat (pd.DataFrame): DataFrame to prepend header to.

    Returns:
        pd.DataFrame: DataFrame with header row added.

    """
    if not conf_header_long:
        pat_header_long = default_row
    else:
        try:
            pat_header_long = pd.DataFrame(
                [conf_header_long.split(", ")], columns=datapat_columns)
        except ValueError as err:
            logger.critical(f"The number of column names "
            f"({len(conf_header_long.split(', '))}) in HEADER_PATIENT_LONG "
            "in conf.ini is different from the effective number of columns "
            f"({len(datapat_columns)}) in {patient_tsv}.")
            msg = "Different number of columns: exiting from walk script!"
            raise(NameError(msg)) from err

    return pd.concat([pat_header_long, final_data_pat], ignore_index=True)


def _is_blank_value(value: object) -> bool:
    """Return True if a MSI/TMB VALUE or THR cell should be treated as empty."""
    if value is None:
        return True
    if isinstance(value, float) and np.isnan(value):
        return True
    return str(value).strip().upper() in ("", "NAN", "NA")


def _resolve_biomarker_threshold(
    sample_id: str,
    marker_name: str,
    value: object,
    thr_given: object,
    compute_thr) -> str:
    """Resolve a MSI/TMB threshold label for one sample.

    Precedence (identical for MSI and TMB, per explicit user design so the
    input can carry a pipeline's own pre-computed threshold - e.g. Guardant
    - instead of forcing everything through conf.ini's thresholds):

    - THR pre-filled in the input -> keep it as-is, never overwritten by
      conf.ini. Logged as a warning (conf.ini threshold ignored) if VALUE is
      also populated, or as info (THR trusted alone) if VALUE is empty.
    - THR empty, VALUE populated -> compute from conf.ini's own threshold,
      exactly as before.
    - Neither populated -> "NA".

    Args:
        sample_id (str): Sample ID, for logging.
        marker_name (str): "MSI" or "TMB", for logging.
        value (object): The raw MSI/TMB numeric value cell (may be blank).
        thr_given (object): The raw MSI_THR/TMB_THR cell (may be blank).
        compute_thr (Callable[[str], str]): conf.ini-threshold-based
            classifier, called only when thr_given is blank.

    Returns:
        str: The resolved threshold label to write to data_clinical_sample.txt.

    """
    value_present = not _is_blank_value(value)
    thr_present = not _is_blank_value(thr_given)

    if thr_present:
        thr_given = str(thr_given).strip()
        if value_present:
            logger.warning(
                f"Sample {sample_id}: {marker_name}_THR is pre-filled "
                f"('{thr_given}') and {marker_name} value is also populated "
                f"- keeping the provided {marker_name}_THR and ignoring "
                f"conf.ini's {marker_name} threshold.")
        else:
            logger.info(
                f"Sample {sample_id}: {marker_name} value is empty but "
                f"{marker_name}_THR is pre-filled ('{thr_given}') - keeping "
                "it as-is.")
        return thr_given

    if value_present:
        return compute_thr(value)

    return "NA"


def fill_from_file(
    table_dict_patient: dict[str, list],
    file_input_clinical: pd.DataFrame,
    types: dict[str, str] | None = None) -> dict[str, list]:
    """Populate clinical table dictionary with MSI and TMB values and statuses.

    Honors a pre-filled MSI_THR/TMB_THR column in the input (if present) over
    conf.ini's own thresholds - see _resolve_biomarker_threshold(). This
    replaces the previous behavior of always recomputing Stable/Unstable and
    the TMB category from conf.ini, which forced non-Illumina pipelines
    (e.g. Guardant) that already know their own threshold to encode it via
    a fake VALUE just to get the right label out.

    Args:
        table_dict_patient (dict): Patient clinical data dictionary.
        file_input_clinical (pd.DataFrame): Clinical input with SAMPLE_ID,
            MSI, TMB and, if present, MSI_THR/TMB_THR.
        types (dict | None): {SAMPLE_ID: Tissue/Plasma}; each sample's
            thresholds are [MSI]/[TMB]'s, or its [MSI.<type>]/[TMB.<type>]
            override (see specimen.setting()).

    Returns:
        dict: Updated clinical data dictionary.

    """
    has_msi_thr_col = "MSI_THR" in file_input_clinical.columns
    has_tmb_thr_col = "TMB_THR" in file_input_clinical.columns

    types = types or {}

    for _, row in file_input_clinical.iterrows():
        k = row["SAMPLE_ID"]
        kind = types.get(str(k).strip())
        msi_thr = specimen.setting("MSI", "THRESHOLD_MSI", kind)
        tmb_thr = ast.literal_eval(specimen.setting("TMB", "THRESHOLD_TMB", kind))

        def compute_msi(m: object, msi_thr: str = msi_thr) -> str:
            return "Stable" if meets_threshold(m, msi_thr) else "Unstable"

        def compute_tmb(t: object, tmb_thr: dict = tmb_thr, k: str = k) -> str:
            for _k, _v in tmb_thr.items():
                if meets_threshold(t, _v):
                    return _k
            logger.warning(f"TMB {t} out of range for {k}")
            return "Out of threshold ranges"

        m = row["MSI"]
        t = row["TMB"]
        table_dict_patient[k].append(m)
        table_dict_patient[k].append(t)

        msi_thr_given = row["MSI_THR"] if has_msi_thr_col else None
        table_dict_patient[k].append(
            _resolve_biomarker_threshold(k, "MSI", m, msi_thr_given, compute_msi))

        tmb_thr_given = row["TMB_THR"] if has_tmb_thr_col else None
        table_dict_patient[k].append(
            _resolve_biomarker_threshold(k, "TMB", t, tmb_thr_given, compute_tmb))

    return table_dict_patient


def fill_from_combined(
    combined_dict: dict[str, str],
    table_dict_patient: dict[str, list],
    types: dict[str, str] | None = None) -> dict[str, list]:
    """Extract MSI and TMB from combined output files and update clinical dict.

    Each sample's thresholds are [MSI]/[TMB]'s, or its [MSI.<type>]/
    [TMB.<type>] override (see specimen.setting()). Plasma samples are
    called with [MSI.Plasma] THRESHOLD_MSI when set (same syntax as
    [MSI] THRESHOLD_MSI: the condition for "Stable"), otherwise with the
    older [MSI] THRESHOLD_MSI_LIQUID (Unstable at or above it) when set,
    otherwise with [MSI] THRESHOLD_MSI.

    Args:
        combined_dict (dict): Map sample ID to CombinedVariantOutput path.
        table_dict_patient (dict): Clinical data dictionary to update.
        types (dict | None): {SAMPLE_ID: Tissue/Plasma}; samples not in it
            take the study's [Sample_Type] TYPE.

    Returns:
        dict: Updated clinical data dictionary.

    """
    types = types or {}
    default_type = specimen.study_specimen_type()

    for k, v in combined_dict.items():
        kind = types.get(str(k).strip(), default_type)
        msi_sites_thr = specimen.setting("MSI", "THRESHOLD_SITES", kind)
        msi_thr = specimen.setting("MSI", "THRESHOLD_MSI", kind)
        tmb = ast.literal_eval(specimen.setting("TMB", "THRESHOLD_TMB", kind))
        try:
            tmv_msi = tsv.get_msi_tmb(Path(v), kind)
        except Exception as err:
            logger.error(f"Something went wrong reading CombinedOutput for sample "
                         f"{k} ({err})! MSI/TMB set to NA for this sample.")
            # Keep table_dict_patient[k] the same length as every successfully
            # processed sample (4 appends below) - otherwise
            # pd.DataFrame.from_dict(table_dict_patient) in write_clinical_sample
            # raises "All arrays must be of the same length" for the whole batch.
            table_dict_patient[k].extend(["NA", "NA", "NA", "NA"])
            continue

        if (
            tmv_msi["MSI"][0][1] != "NA" and
            meets_threshold(tmv_msi["MSI"][0][1], msi_sites_thr)):
            table_dict_patient[k].append(tmv_msi["MSI"][1][1])
        else:
            table_dict_patient[k].append("NA")

        table_dict_patient[k].append(tmv_msi["TMB_Total"])

        if (tmv_msi["MSI"][1][1] != "NA" and table_dict_patient[k][-2] != "NA"):
                msi_value = float(tmv_msi["MSI"][1][1])

                if kind == specimen.PLASMA and not config.has_option(
                        f"MSI.{specimen.PLASMA}", "THRESHOLD_MSI") and config.get(
                        "MSI", "THRESHOLD_MSI_LIQUID", fallback="").strip():
                    if msi_value >= float(config.get("MSI", "THRESHOLD_MSI_LIQUID")):
                        table_dict_patient[k].append("Unstable")
                    else:
                        table_dict_patient[k].append("Stable")

                else:
                    if meets_threshold(msi_value, msi_thr):
                        table_dict_patient[k].append("Stable")
                    else:
                        table_dict_patient[k].append("Unstable")

        else:
            table_dict_patient[k].append("NA")

        if tmv_msi["TMB_Total"] != "NA":
            found = False
            for _k, _v in tmb.items():
                if meets_threshold(tmv_msi["TMB_Total"], _v):
                    table_dict_patient[k].append(_k)
                    found = True
                    break
            if not found:
                logger.warning(
                    f"TMB {tmv_msi['TMB_Total']} out of range for {k}")
                table_dict_patient[k].append("Out of threshold ranges")
        else:
            table_dict_patient[k].append("NA")

    return table_dict_patient


def update_data_clinical_with_exon_info(
    combined_output: Path,
    combined_dict: dict,
    output_folder: Path):
    """Updates clinical_sample WITH exon-level CNV information for BRCA1 and BRCA2

    Args:
        combined_output (Path): Directory containing combined output files.
        combined_dict (dict): Mapping of {sample_id: path_to_tsv}.
        output_folder (Path): Output directory where data_clinical_sample.txt is located.

    """
    data_clin_path = Path(output_folder) / "data_clinical_sample.txt"

    if combined_output.is_dir() and any(combined_output.iterdir()):
        all_exon_dfs = []

        for sample_id, tsv_path in combined_dict.items():
            try:
                exonic = tsv.get_exons(Path(tsv_path))
            except Exception as e:
                logger.error(f"Error while reading exon-level CNVs for sample {sample_id}: {e}")
                exonic = []

            exon_df = build_exon_df(output_folder, exonic, sample_id=sample_id)
            all_exon_dfs.append(exon_df)

        if all_exon_dfs:
            all_exon_dfs = pd.concat(all_exon_dfs, ignore_index=True)

            with open(data_clin_path, "r", encoding="utf-8") as f:
                header_lines = [next(f) for _ in range(4)]

            data_clin_df = pd.read_csv(
                data_clin_path, 
                sep="\t", 
                header=4, 
                dtype={"SAMPLE_ID": str}
            )
            merged_data_clin = pd.merge(data_clin_df, all_exon_dfs, on="SAMPLE_ID", how="left")

            new_columns = ["Exonic_BRCA1", "Exonic_BRCA2", "BRCA1_details", "BRCA2_details"]

            updated_headers = []
            updated_headers.append(header_lines[0].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n")
            updated_headers.append(header_lines[1].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n")
            updated_headers.append(header_lines[2].rstrip("\n") + "\t" + "\t".join(["STRING"] * len(new_columns)) + "\n")
            updated_headers.append(header_lines[3].rstrip("\n") + "\t" + "\t".join(["1"] * len(new_columns)) + "\n")

            with open(data_clin_path, "w", encoding="utf-8") as f:
                f.writelines(updated_headers)
                # lineterminator="\n": f is a text-mode handle that already
                # turns "\n" into the platform's line ending - pandas' own
                # default (os.linesep) would double it into "\r\r\n" on
                # Windows, which cBioPortal's validator rejects.
                merged_data_clin.to_csv(f, sep="\t", index=False, lineterminator="\n")

        else:
            logger.warning("No exon CNV data available to update the clinical file.")


def update_data_clinical_with_hrd_info(
    combined_output: Path,
    combined_dict: dict,
    output_folder: Path) -> None:
    """Update data_clinical_sample.txt with HRD/GIS biomarker info.

    Genomic Instability Score, Tumor Fraction and Ploidy are only reported by
    DRAGEN for samples run with the TSO500 HRD feature enabled - for every
    other sample these three new columns are just "NA", the same way
    update_data_clinical_with_exon_info() layers BRCA exon-level info onto
    this same file without touching the rest of it.

    Args:
        combined_output (Path): Directory containing combined output files.
        combined_dict (dict): Mapping of {sample_id: path_to_tsv}.
        output_folder (Path): Output directory where data_clinical_sample.txt is located.

    """
    data_clin_path = Path(output_folder) / "data_clinical_sample.txt"

    if combined_output.is_dir() and any(combined_output.iterdir()):
        all_gis_rows = []

        for sample_id, tsv_path in combined_dict.items():
            try:
                gis = tsv.get_gis(Path(tsv_path))
            except Exception as e:
                logger.error(f"Error while reading GIS/HRD info for sample {sample_id}: {e}")
                gis = {"GIS": "NA", "Tumor_Fraction": "NA", "Ploidy": "NA"}

            all_gis_rows.append({
                "SAMPLE_ID": sample_id,
                "GENOMIC_INSTABILITY_SCORE": gis["GIS"],
                "TUMOR_FRACTION": gis["Tumor_Fraction"],
                "PLOIDY": gis["Ploidy"],
            })

        if all_gis_rows:
            all_gis_df = pd.DataFrame(all_gis_rows)

            with open(data_clin_path, "r", encoding="utf-8") as f:
                header_lines = [next(f) for _ in range(4)]

            data_clin_df = pd.read_csv(
                data_clin_path,
                sep="\t",
                header=4,
                dtype={"SAMPLE_ID": str}
            )
            merged_data_clin = pd.merge(data_clin_df, all_gis_df, on="SAMPLE_ID", how="left")

            new_columns = ["GENOMIC_INSTABILITY_SCORE", "TUMOR_FRACTION", "PLOIDY"]

            updated_headers = [
                header_lines[0].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n",
                header_lines[1].rstrip("\n") + "\t" + "\t".join(new_columns) + "\n",
                header_lines[2].rstrip("\n") + "\t" + "\t".join(["NUMBER"] * len(new_columns)) + "\n",
                header_lines[3].rstrip("\n") + "\t" + "\t".join(["1"] * len(new_columns)) + "\n",
            ]

            with open(data_clin_path, "w", encoding="utf-8") as f:
                f.writelines(updated_headers)
                # lineterminator="\n": f is a text-mode handle that already
                # turns "\n" into the platform's line ending - pandas' own
                # default (os.linesep) would double it into "\r\r\n" on
                # Windows, which cBioPortal's validator rejects.
                merged_data_clin.to_csv(f, sep="\t", index=False, lineterminator="\n")
        else:
            logger.warning("No sample data available to update HRD/GIS clinical info.")


def _walk_write_clinical_tables(ctx: WalkContext) -> None:
    """Stage: data_clinical_patient.txt + data_clinical_sample.txt (MSI/TMB,
    exon, HRD info).

    Depends only on WalkContext, independent of _walk_process_cnv/_snv/_fusion
    - it reads sample.tsv/patient.tsv and CombinedOutput directly, not their
    output files.
    """
    table_dict_patient = get_table_from_folder(ctx.clin_sample_path)
    logger.info("Writing clinical files...")

    if Path(ctx.patient_tsv).name.strip() != "":
        logger.info("Writing data_clinical_patient.txt file...")

        input_file_path = Path(ctx.input_folder) / "patient.tsv"
        data_clin_pat = pd.read_csv(input_file_path, sep="\t", header=0, dtype=str)

        data_clin_pat.columns = data_clin_pat.columns.str.upper()
        datapat_columns = list(data_clin_pat.columns)

        # Get headers from conf.ini if they're present
        conf_header_short = config.get("ClinicalPatient", "HEADER_PATIENT_SHORT")
        conf_header_long = config.get("ClinicalPatient", "HEADER_PATIENT_LONG")
        conf_header_type = config.get("ClinicalPatient", "HEADER_PATIENT_TYPE")

        # Add header's fifth row
        default_row = pd.DataFrame([datapat_columns], columns=datapat_columns)
        final_data_pat = pd.concat([default_row, data_clin_pat], ignore_index=True)

        # Add header's fourth row (1s)
        header_numbers = pd.DataFrame(
            [[1] * len(datapat_columns)], columns=datapat_columns)
        final_data_pat = pd.concat([header_numbers, final_data_pat], ignore_index=True)

        # Add header's third row (HEADER_PATIENT_TYPE)
        final_data_pat = add_header_patient_type(
            ctx.patient_tsv, datapat_columns, conf_header_type,
            final_data_pat)

        # Add header's second row (HEADER_PATIENT_LONG)
        final_data_pat = add_header_patient_long(
            ctx.patient_tsv, datapat_columns, conf_header_long,
            default_row, final_data_pat)

        # Add header's first row (HEADER_PATIENT_SHORT)
        final_data_pat = add_header_patient_short(
            ctx.patient_tsv, datapat_columns, conf_header_short,
            default_row, final_data_pat)

        final_data_pat.loc[0:3, "PATIENT_ID"] = final_data_pat.loc[
            0:3, "PATIENT_ID"].apply(lambda x: f"#{x}")

        data_clin_txt = Path(ctx.output_folder) / "data_clinical_patient.txt"
        final_data_pat.to_csv(data_clin_txt, sep="\t", index=False, header=False)

    else:
        write_default_clinical_patient(ctx.output_folder, table_dict_patient)

    file_input_sample = pd.read_csv(
        ctx.clin_sample_path, sep="\t", index_col=False, dtype=str)

    types = specimen.specimen_types(ctx.clin_file)

    combined_output = Path(ctx.input_folder) / "CombinedOutput"
    if combined_output.exists() and len(list(combined_output.iterdir())) > 0:
        combined_dict = get_combined_variant_output_from_folder(
            ctx.input_folder, ctx.clin_file, ctx.isinputfile)
        new_table_dict_patient = fill_from_combined(
            combined_dict, table_dict_patient, types)
    else:
        new_table_dict_patient = fill_from_file(
            table_dict_patient, file_input_sample, types)
        combined_dict = {}

    write_clinical_sample(ctx.clin_sample_path, ctx.output_folder, new_table_dict_patient,
                           combined_output_used=bool(combined_dict))

    write_specimen_type(ctx.output_folder, types)

    update_data_clinical_with_exon_info(combined_output, combined_dict, ctx.output_folder)

    update_data_clinical_with_hrd_info(combined_output, combined_dict, ctx.output_folder)
