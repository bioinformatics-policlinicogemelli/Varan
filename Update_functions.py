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

"""Module: update_functions.

Provides functions to update and synchronize study data files
across versioned folders in a cBioPortal-style workflow.

Key functionalities include:
  - Merging clinical sample and patient files.
  - Updating copy number alteration (CNA) and fold-change CNA data.
  - Synchronizing mutation and structural variation data.
  - Managing case list updates (CNA, sequenced, SV).

Each function reads the old and new versions of a file, merges or copies
as appropriate, and writes the result to the specified output folder.

"""

import re
import shutil
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from loguru import logger

from versioning import (
    create_newest_version_folder,
    get_version_list,
)


def update_clinical_samples(oldfile_path: str,
                            newfile_path: str,
                            output_folder: str) -> None:
    """Update data_clinical_sample.txt file.

    This function reads the original txt file from the given 'oldfile_path',
    insert new rows with the sample info founded inside the new txt file from
    the given 'newfile_path' and save the updated file.

    Args:
        oldfile_path (str): Path to the original data_clinical_sample.
        newfile_path (str): Path to the new data_clinical_sample.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    #header
    # Indexed by column *name* (SAMPLE_ID), not position: this file's usual
    # column order happens to put SAMPLE_ID first, but assuming that instead
    # of naming it explicitly meant a differently-ordered (e.g. hand-edited,
    # or PATIENT_ID-first) clinical file would have SAMPLE_ID silently
    # dropped from old_body_unique below, crashing the later merge on
    # "SAMPLE_ID" with a KeyError.
    old_head = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=None, nrows=5)
    old_head.columns = old_head.iloc[4].tolist()
    old_head=old_head.set_index("SAMPLE_ID")

    new_head = pd.read_csv(newfile_path, sep="\t", dtype = str, header=None, nrows=5)
    new_head.columns = new_head.iloc[4].tolist()
    new_head=new_head.set_index("SAMPLE_ID")

    only_common_head = np.intersect1d(new_head.columns, old_head.columns)

    old_updated_head = old_head.drop(list(only_common_head), axis=1)
    old_updated_head.index = new_head.index
    new_head = pd.concat([new_head, old_updated_head], axis=1).fillna(value=np.nan)
    final_head = new_head.reset_index(drop=False)

    #body
    old_body = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=4)
    new_body = pd.read_csv(newfile_path, sep="\t", dtype = str, header=4)

    old_body_unique = old_body.drop(list(only_common_head), axis=1)
    final_body = old_body_unique.merge(new_body, how="outer", on="SAMPLE_ID")

    old_body=old_body.set_index("SAMPLE_ID")
    new_body=new_body.set_index("SAMPLE_ID")
    common_samples = old_body.index.intersection(new_body.index)
    old_body = old_body.drop(index=common_samples)
    final_body=final_body.set_index("SAMPLE_ID")
    final_body.update(old_body, overwrite=True, filter_func=None, errors="ignore")
    final_body = final_body.reset_index(drop=False)

    #final
    final_patient = pd.concat([final_head, final_body], axis=0)
    outpath=Path(output_folder) / "data_clinical_sample.txt"
    final_patient.to_csv(outpath, header=False, index=False, sep="\t", na_rep="NaN")
    logger.info("data_clinical_sample.txt updated!")


def update_clinical_patient(oldfile_path: str,
                            newfile_path: str,
                            output_folder: str) -> None:
    """Update data_clinical_patient.txt file.

    This function reads the originaltxt file from the 'oldfile_path',
    insert new rows with the patients info founded inside the new txt file
    from the given 'newfile_path' and save the updated 'data_clinical_patient.txt'

    Args:
        oldfile_path (str): Path to the original data_clinical_patient.
        newfile_path (str): Path to the new data_clinical_patient.
        output_folder (str): Path to the output folder where the file will be saved.

    Returns:
        None

    """
    #header
    # Indexed by column *name* (PATIENT_ID), not position - see the matching
    # comment in update_clinical_samples for why position-based indexing is
    # fragile here.
    old_head = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=None, nrows=5)
    old_head.columns = old_head.iloc[4].tolist()
    old_head=old_head.set_index("PATIENT_ID")

    new_head = pd.read_csv(newfile_path, sep="\t", dtype = str, header=None, nrows=5)
    new_head.columns = new_head.iloc[4].tolist()
    new_head=new_head.set_index("PATIENT_ID")

    only_common_head = np.intersect1d(new_head.columns, old_head.columns)

    old_updated_head = old_head.drop(list(only_common_head), axis=1)
    old_updated_head.index = new_head.index
    new_head = pd.concat([new_head, old_updated_head], axis=1).fillna(value=np.nan)
    final_head = new_head.reset_index(drop=False)

    #body
    old_body = pd.read_csv(oldfile_path, sep="\t", dtype = str, header=4)
    new_body = pd.read_csv(newfile_path, sep="\t", dtype = str, header=4)

    old_body_unique = old_body.drop(list(only_common_head), axis=1)
    final_body = old_body_unique.merge(new_body, how="outer", on="PATIENT_ID")

    old_body=old_body.set_index("PATIENT_ID")
    new_body=new_body.set_index("PATIENT_ID")
    common_samples = old_body.index.intersection(new_body.index)
    old_body = old_body.drop(index=common_samples)
    final_body=final_body.set_index("PATIENT_ID")
    final_body.update(old_body, overwrite=True, filter_func=None, errors="ignore")
    final_body = final_body.reset_index(drop=False)

    #final
    final_patient = pd.concat([final_head, final_body], axis=0)
    outpath=Path(output_folder) / "data_clinical_patient.txt"
    final_patient.to_csv(outpath, header=False, index=False, sep="\t", na_rep="NaN")
    logger.info("data_clinical_patient.txt updated!")


def update_cna_hg19(oldfile_path: str, newfile_path: str, output_folder: str) -> None:
    """Update data_cna_hg19.seg.

    This function reads the original CNA data from the given 'oldfile_path',
    insert new rows with the sample CNA data founded inside the new file from
    the given 'newfile_path' and save the updated file.

    Any sample (ID) present in both files is treated as fully reprocessed:
    all of its old segments are dropped before merging, rather than relying
    on old and new segments to be identical on (ID, chrom, loc.start,
    loc.end) to collapse via drop_duplicates. If a sample was re-annotated
    with different results, its segmentation can shift slightly - matching
    only on coordinates left stale old segments in place whenever a
    boundary changed, alongside the new ones.

    Args:
        oldfile_path (str): Path to the original CNA data file.
        newfile_path (str): Path to the new CNA data file.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype={"ID": str})
    new = pd.read_csv(newfile_path, sep="\t", dtype={"ID": str})

    updated_samples = set(old["ID"]) & set(new["ID"])
    old = old[~old["ID"].isin(updated_samples)]

    updated = pd.concat([old, new], ignore_index=True)
    outpath=Path(output_folder) / "data_cna_hg19.seg"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_cna_hg19.seg updated!")


def update_cna_hg19_fc(oldfile_path: str,
                       newfile_path: str,
                       output_folder: str) -> None:
    """Update fold-change copy number alteration (CNA) data in hg19 format.

    Reads the original and new fold-change CNA segment files (TSV) and
    writes the consolidated data to `data_cna_hg19.seg.fc.txt`.

    Any sample (ID) present in both files is treated as fully reprocessed:
    all of its old segments are dropped before merging - see update_cna_hg19
    for why matching only on segment coordinates isn't enough.

    Args:
        oldfile_path (str): Path to the original CNA fold-change file.
        newfile_path (str): Path to the new CNA fold-change file.
        output_folder (str): Directory where the updated file will be saved.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype={"ID": str})
    new = pd.read_csv(newfile_path, sep="\t", dtype={"ID": str})

    updated_samples = set(old["ID"]) & set(new["ID"])
    old = old[~old["ID"].isin(updated_samples)]

    updated = pd.concat([old, new], ignore_index=True)
    outpath=Path(output_folder) / "data_cna_hg19.seg.fc.txt"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_cna_hg19.seg.fc.txt updated!")


def update_cna(oldfile_path: str,
               newfile_path: str,
               output_folder: str) -> None:
    """Update sample inside a copy number alteration (CNA) data file.

    This function reads the original tab separated CNA data,
    insert new rows with the sample CNA data founded inside the new file
    from the given 'newfile_path' and save the updated file.

    Args:
        oldfile_path (str): Path to the original CNA data file.
        newfile_path (str): Path to the new CNA data file.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", index_col=0)
    new = pd.read_csv(newfile_path, sep="\t", index_col=0)

    sample_old = old.columns
    sample_new = new.columns
    to_remove = sample_old.intersection(sample_new)

    updated = pd.concat(
        [old.drop(columns=to_remove), new]
        , axis=1).replace(np.nan, 0).astype(int)
    outpath=Path(output_folder) / "data_cna.txt"
    updated.to_csv(outpath, index=True, sep="\t")
    logger.info("data_cna.txt updated!")


def update_exon_brca_cna(oldfile_path: str,
                         newfile_path: str,
                         output_folder: str) -> None:
    """Update sample columns in the BRCA exon-level CNV Generic Assay matrix.

    Mirrors update_cna's column-intersection replace logic, but for
    data_exon_brca_cna.txt's shape: two leading id columns
    (ENTITY_STABLE_ID, NAME) that must be kept as-is rather than treated as
    a single index column, and CATEGORICAL string values (LOSS/GAIN/NEUTRAL)
    rather than data_cna.txt's discrete integers - so it can't reuse
    update_cna directly.

    Args:
        oldfile_path (str): Path to the original data_exon_brca_cna.txt.
        newfile_path (str): Path to the new data_exon_brca_cna.txt.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t")
    new = pd.read_csv(newfile_path, sep="\t")

    # Must be present in *both* files, not just old: merge(on=id_cols)
    # requires every column in id_cols to exist on both sides, so an
    # id_cols entry old has but new doesn't (e.g. a schema change that
    # drops/renames NAME) would otherwise raise KeyError on the missing
    # column in new instead of just merging on whichever id columns both
    # files still agree on.
    id_cols = [c for c in ["ENTITY_STABLE_ID", "NAME"]
               if c in old.columns and c in new.columns]
    sample_old = [c for c in old.columns if c not in id_cols]
    sample_new = [c for c in new.columns if c not in id_cols]
    to_remove = set(sample_old) & set(sample_new)

    old = old.drop(columns=list(to_remove))
    updated = old.merge(new, on=id_cols, how="outer") if id_cols else pd.concat(
        [old, new], axis=1)
    outpath=Path(output_folder) / "data_exon_brca_cna.txt"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_exon_brca_cna.txt updated!")


def update_mutations(oldfile_path: str,
                     newfile_path: str,
                     output_folder: str) -> None:
    """Update samples' mutation data inside data_mutations_extended.txt file.

    This function reads the original tab separated version txt file,
    insert new rows with the samples' mutation data founded inside the new txt file
    and save the updated file.

    Any sample (Tumor_Sample_Barcode) present in both files is treated as
    fully reprocessed: all of its old rows are dropped before merging.
    Relying on old and new rows to be byte-identical on the Hugo_Symbol:n_AF
    column range to collapse via drop_duplicates is fragile across two
    independent annotation runs (a changed ONCOTREE_CODE, a ClinVar/VEP/
    OncoKB version bump between runs, ...) even when the underlying variant
    calls are the same - and silently leaves the old, possibly stale, rows
    in place for any sample whose re-annotated output differs even
    slightly, alongside the new ones (duplicates, or stale rows for calls
    that no longer exist).

    Args:
        oldfile_path (str): Path to the original data_mutations_extended.
        newfile_path (str): Path to the new data_mutations_extended.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype=str)
    new = pd.read_csv(newfile_path, sep="\t", dtype=str)

    updated_samples = set(old["Tumor_Sample_Barcode"]) & set(new["Tumor_Sample_Barcode"])
    old = old[~old["Tumor_Sample_Barcode"].isin(updated_samples)]

    updated = pd.concat([old, new], ignore_index=True)
    outpath=Path(output_folder) / "data_mutations_extended.txt"
    updated.to_csv(outpath, index=False, sep="\t")
    logger.info("data_mutation_extended.txt updated!")


def update_sv(oldfile_path: str,
              newfile_path: str,
              output_folder: str) -> None:
    """Update samples' structural variation (SV) data inside data_sv.txt file.

    This function reads the original tab separated version txt file,
    insert new rows with the samples' SV data founded inside the new txt file
    and save the updated file named 'data_sv.txt' in the specified 'output_folder'.

    Any sample (Sample_Id) present in both files is treated as fully
    reprocessed: all of its old rows are dropped before merging, rather
    than relying on the row-content dedup below (which matches cBioPortal's
    own uniqueness rules, not sample identity) to catch a re-annotated
    sample whose fusion/splice calls changed - a call that disappeared
    between runs would otherwise never collide with anything and would be
    left behind as a stale row.

    Args:
        oldfile_path (str): Path to the original data_sv.
        newfile_path (str): Path to the new data_sv.
        output_folder (str): Path to the output folder.

    Returns:
        None

    """
    df_old = pd.read_csv(oldfile_path, sep="\t", dtype=str)
    df_new = pd.read_csv(newfile_path, sep="\t", dtype=str)

    updated_samples = set(df_old["Sample_Id"]) & set(df_new["Sample_Id"])
    df_old = df_old[~df_old["Sample_Id"].isin(updated_samples)]

    merged_df = pd.concat(
        [df_old, df_new], axis=0, join="outer", ignore_index=True)

    # Dedup exactly like the fresh-creation path in walk.py: match cBioPortal's own
    # StructuralVariantValidator.UNIQUENESS_COLUMNS (which includes breakpoint
    # position/chromosome, not just the gene pair), so we neither merge two
    # genuinely distinct fusion calls between the same genes, nor leave in place
    # rows cBioPortal itself would reject as duplicates. Among true duplicates,
    # keep the one with the highest Normal_Paired_End_Read_Count.
    if "Normal_Paired_End_Read_Count" in merged_df.columns:
        merged_df["Normal_Paired_End_Read_Count"] = pd.to_numeric(
            merged_df["Normal_Paired_End_Read_Count"], errors="coerce")
        merged_df = merged_df.sort_values(
            by="Normal_Paired_End_Read_Count", ascending=False)
        col_subset = [col for col in merged_df.columns
                      if col != "Normal_Paired_End_Read_Count"]
        merged_df = merged_df.drop_duplicates(subset=col_subset, keep="first")
        merged_df["Normal_Paired_End_Read_Count"] = (
            merged_df["Normal_Paired_End_Read_Count"]
            .astype(str).str.replace(r"\.0$", "", regex=True)
            .replace("nan", ""))
    else:
        merged_df = merged_df.drop_duplicates(keep="first")

    output_file = Path(output_folder) / "data_sv.txt"
    merged_df.to_csv(output_file, sep="\t", index=False)

    if "Class" in merged_df.columns:
        n_splice = (merged_df["Class"] == "SPLICE").sum()
        n_fusion = (merged_df["Class"] == "FUSION").sum()
        logger.info(
            f"data_sv.txt updated! ({n_fusion} fusion, {n_splice} splice "
            "variant row(s) after merge/dedup)")
    else:
        logger.info("data_sv.txt updated!")


def update_generic_by_sample_id(
    oldfile_path: Path, newfile_path: Path, output_folder: Path,
    id_column: str = "Sample_Id") -> None:
    """Update a generic Sample_Id-keyed file by replacing reprocessed samples.

    Works regardless of how many rows a sample has (e.g. exon-level CNA data can have
    one row per gene per sample) - reusable for any future per-sample multi-row file
    without writing a new update_* function each time.

    Any sample (`id_column`) present in both files is treated as fully
    reprocessed: all of its old rows are dropped before merging, instead of
    relying on old and new rows being byte-identical to collapse via
    drop_duplicates - which left stale rows behind for any sample whose
    re-annotated output changed even slightly, exactly like the
    data_mutations_extended.txt/data_sv.txt bug this mirrors.

    Args:
        oldfile_path (Path): Path to the file in the original study folder.
        newfile_path (Path): Path to the same-named file in the incoming study folder.
        output_folder (Path): Path to the output folder.
        id_column (str): Name of the sample-id column in the file.

    Returns:
        None

    """
    old = pd.read_csv(oldfile_path, sep="\t", dtype=str)
    new = pd.read_csv(newfile_path, sep="\t", dtype=str)

    if id_column in old.columns and id_column in new.columns:
        updated_samples = set(old[id_column]) & set(new[id_column])
        old = old[~old[id_column].isin(updated_samples)]

    merged = pd.concat([old, new], ignore_index=True)
    output_file = Path(output_folder) / Path(oldfile_path).name
    merged.to_csv(output_file, sep="\t", index=False)
    logger.info(f"{Path(oldfile_path).name} updated!")


def check_files(oldpath: str,
                newpath: str,
                output: str,
                file_name: str) -> None:
    """Dispatch and/or copy study data files between versions.

    Checks for the presence of a given file in both the old and new paths.
    - If present in both: calls the appropriate `update_*` function to merge.
    - If only in one: copies it to the output folder.
    - If in neither: logs a warning and skips.

    Args:
        oldpath (str): Path to the previous version folder.
        newpath (str): Path to the incoming data folder.
        output (str): Path to the target output folder.
        file_name (str): Name of the file to process, e.g.
            "data_clinical_sample.txt", "data_cna.txt", etc.

    Returns:
        None

    """
    file_updaters: dict[str, Callable[[Path, Path, Path], None]] = {
    "data_clinical_sample.txt": update_clinical_samples,
    "data_clinical_patient.txt": update_clinical_patient,
    "data_cna_hg19.seg": update_cna_hg19,
    "data_cna_hg19.seg.fc.txt": update_cna_hg19_fc,
    "data_cna.txt": update_cna,
    "data_mutations_extended.txt": update_mutations,
    "data_sv.txt": update_sv,
    "exon_CNA_data.txt": update_generic_by_sample_id,
    "data_exon_brca_cna.txt": update_exon_brca_cna}

    o_data = Path(oldpath) / file_name
    n_data = Path(newpath) / file_name
    dest = Path(output) / file_name

    if o_data.exists() and n_data.exists():
        updater = file_updaters.get(file_name)
        if updater:
            updater(o_data, n_data, output)
        else:
            logger.warning(f"No updater function for '{file_name}'. Skipping merge.")
    elif o_data.exists():
        logger.warning(f"'{file_name}' only in old; copying forward.")
        shutil.copy(o_data, dest)
    elif n_data.exists():
        logger.warning(f"'{file_name}' only in new; copying forward.")
        shutil.copy(n_data, dest)
    else:
        logger.warning(f"'{file_name}' not found in either folder. Skipping.")


def safe_check_file(oldpath: Path, newpath: Path, output: Path, file: str) -> None:
    """Safely checks a file for consistency between old and new study folders.

    This function wraps `check_files` in a try-except block.

    Args:
        oldpath (Path): Path to the old version of the study folder.
        newpath (Path): Path to the new version of the study folder.
        output (Path): Path to the output (updated) study folder.
        file (str): Name of the file to check (e.g., "data_clinical_sample.txt").

    Raises:
        IndexError: Raised when a parsing error is detected in the file structure.

    """
    try:
        check_files(oldpath, newpath, output, file)
    except pd.errors.ParserError as e:
        line_match = re.search(r"line (\d+)", str(e))
        line_number = line_match.group(1) if line_match else "unknown"
        logger.critical(
            f"Wrong column number in line {line_number} of {file} file")
        msg = "Exiting from Update script!"
        raise IndexError(msg) from e

def copy_metadata_files(oldpath: Path, newpath: Path, output: Path) -> None:
    """Copy meta_*.txt files forward from both source studies.

    Copies from oldpath, then adds any meta file that only exists in
    newpath - e.g. meta_exon_brca_cna.txt when oldpath predates that data
    type but newpath (a freshly walked batch being merged in) already has
    it. Only scanning oldpath silently dropped any such new-in-this-batch
    meta file. Which source wins on a same-named file doesn't matter: every
    file copied here is just a placeholder for remove_meta's cleanup pass
    and meta_case_main's proper regeneration afterward, both of which run
    later in update_main based on the actual merged data.

    Args:
        oldpath (Path): Path to previous version of the study folder.
        newpath (Path): Path to the incoming data folder being merged in.
        output (Path): Path to output (new version) study folder.

    Returns:
        None

    """
    meta_files = {f.name: f for f in newpath.glob("*meta*")}
    meta_files.update({f.name: f for f in oldpath.glob("*meta*")})
    if meta_files:
        for file in meta_files.values():
            shutil.copy(file, output)
    else:
        logger.warning("No meta files found in either study folder!")

def prepare_output_folder(oldpath: str,
                          output: str,
                          overwrite: bool) -> tuple[str, bool, str]:
    """Prepare output folder for the updated study version.

    This function handles:
    - Deriving the base output path if not provided.
    - Removing the last versioned folder if overwrite is enabled.
    - Creating a new versioned output folder.
    - Creating a 'case_lists' subfolder inside the output directory.

    Args:
        oldpath (str): Path to the previous version of the study folder.
        output (str): Base output path for the new version.
        overwrite (bool): Whether to overwrite the latest version folder.

    Returns:
        tuple[str, bool, str]: A tuple containing:
            - The path to the newly created output folder.
            - Boolean flag indicating whether output was auto-inferred.
            - Path to the 'case_lists' subfolder inside the output folder.

    """
    if output != "":
        no_out = False
    else:
        no_out = True
        output = re.split(r"_v[0-9]+$", oldpath)[0]

    old_versions = get_version_list(output)
    if old_versions and Path(old_versions[-1]).exists() and overwrite:
        logger.info("Overwrite option set. Start removing folder")
        shutil.rmtree(old_versions[-1])

    output = create_newest_version_folder(output)
    logger.info(f"Creating a new folder: {output}")

    output_caseslists = Path(output) / "case_lists"
    Path(output_caseslists).mkdir()

    return output, no_out, output_caseslists
