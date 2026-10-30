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

"""Per-sample study data files, and the three operations on them that
update / extract / remove need: keep some samples, drop some samples, or
merge two studies.

STUDY_FILES is the one registry of those files - add a new per-sample
file here and all three flows handle it. Most files are generic:

- "rows": one row per event, a column holding the sample id (MAF, SV,
  .seg, gene panel matrix, ...);
- "columns": one column per sample next to entity columns (data_cna.txt,
  data_exon_brca_cna.txt).

The clinical files have their own layout (cBioPortal's 4 commented header
lines) and rules - a patient is only removed once none of its samples is
left - and data_sv.txt deduplicates on merge; those are the only special
cases.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger


@dataclass(frozen=True)
class StudyFile:
    """One per-sample study data file.

    layout: "rows", "columns", "clinical_sample" or "clinical_patient".
    id_column: rows layout - the column holding the sample id.
    entity_columns: columns layout - the leading non-sample columns.
    integer_values: columns layout - values are integers (0 for a sample
        with no call), and rows left without any non-zero value are
        dropped when samples are removed (data_cna.txt).
    merge: optional replacement for the generic merge (old, new, output).
    """

    name: str
    layout: str
    id_column: str = ""
    entity_columns: tuple[str, ...] = ()
    integer_values: bool = False
    merge: Callable[[Path, Path, Path], None] | None = None


# ---------------------------------------------------------------- subsetting

def _subset_rows(spec: StudyFile, src: Path, out: Path, keep: Callable) -> None:
    df = pd.read_csv(src, sep="\t", dtype=str)
    if spec.id_column not in df.columns:
        logger.warning(f"{spec.name}: column '{spec.id_column}' not found, skipping.")
        return
    kept = df[df[spec.id_column].astype(str).map(keep)]
    if "Class" in df.columns:
        dropped = df.drop(kept.index)
        for label, part in (("kept", kept), ("removed", dropped)):
            counts = part["Class"].value_counts()
            if counts.get("FUSION", 0) or counts.get("SPLICE", 0):
                logger.info(
                    f"{spec.name}: {label} {counts.get('FUSION', 0)} fusion and "
                    f"{counts.get('SPLICE', 0)} splice variant row(s).")
    kept.to_csv(out / spec.name, sep="\t", index=False)


def _subset_columns(spec: StudyFile, src: Path, out: Path, keep: Callable) -> None:
    if spec.integer_values:
        df = pd.read_csv(src, sep="\t", index_col=0)
        kept = df[[c for c in df.columns if keep(str(c))]]
        kept = kept.loc[(kept != 0).any(axis=1)]
        # Integers, as cBioPortal requires: also turns an older study's
        # "2.0"-style values into "2".
        kept.fillna(0).astype(int).to_csv(out / spec.name, sep="\t")
        return
    df = pd.read_csv(src, sep="\t")
    entity = [c for c in spec.entity_columns if c in df.columns]
    df[entity + [c for c in df.columns if c not in entity and keep(str(c))]].to_csv(
        out / spec.name, sep="\t", index=False)


def _subset_clinical_sample(spec: StudyFile, src: Path, out: Path, keep: Callable) -> None:
    df = pd.read_csv(src, sep="\t")
    idx_sample = np.argwhere(df.to_numpy() == "SAMPLE_ID")[0][1]
    is_header = df.index <= 3  # the other 3 "#" lines + the SAMPLE_ID header
    kept = df[is_header | df.iloc[:, idx_sample].astype(str).map(keep)]
    kept.to_csv(out / spec.name, sep="\t", index=False)


def _subset_clinical_patient(old_folder: Path, out: Path, keep: Callable,
                             extracting: bool) -> None:
    """Extract keeps the patients of the kept samples; remove drops a
    patient only once all of its samples are removed (patients without any
    sample are never dropped by a remove)."""
    patients = pd.read_csv(old_folder / "data_clinical_patient.txt", sep="\t", header=None)
    samples = pd.read_csv(old_folder / "data_clinical_sample.txt", sep="\t")
    idx_sample = np.argwhere(samples.to_numpy() == "SAMPLE_ID")[0][1]
    idx_patient = np.argwhere(samples.to_numpy() == "PATIENT_ID")[0][1]
    body = samples.iloc[4:]
    sample_kept = body.iloc[:, idx_sample].astype(str).map(keep)
    with_kept = set(body.iloc[:, idx_patient][sample_kept].astype(str))
    with_dropped = set(body.iloc[:, idx_patient][~sample_kept].astype(str))
    idx_pt = np.argwhere(patients.to_numpy() == "PATIENT_ID")[0][1]
    pt = patients.iloc[:, idx_pt].astype(str)
    stays = pt.isin(with_kept) if extracting else pt.isin(with_kept) | ~pt.isin(with_dropped)
    stays.iloc[:5] = True  # 4 "#" lines + the PATIENT_ID header
    patients[stays].to_csv(
        out / "data_clinical_patient.txt", sep="\t", index=False, header=False,
        na_rep="NaN")


def subset_study(old_folder: str, sample_ids: list[str], output: str, keep: bool) -> None:
    """Write `output`'s copy of every per-sample file of `old_folder`,
    keeping only `sample_ids` (keep=True, extract) or every sample but them
    (keep=False, remove)."""
    wanted = {str(s) for s in sample_ids}
    def stays(sample_id: str) -> bool:
        return (sample_id in wanted) == keep
    old, out = Path(old_folder), Path(output)
    for spec in STUDY_FILES:
        src = old / spec.name
        if not src.exists():
            logger.warning(f"{spec.name} not found in current folder. Skipping")
            continue
        if spec.layout == "clinical_patient":
            _subset_clinical_patient(old, out, stays, extracting=keep)
        elif spec.layout == "clinical_sample":
            _subset_clinical_sample(spec, src, out, stays)
        elif spec.layout == "columns":
            _subset_columns(spec, src, out, stays)
        else:
            _subset_rows(spec, src, out, stays)


# ------------------------------------------------------------------- merging

def _merge_rows(spec: StudyFile, old_path: Path, new_path: Path, out: Path) -> None:
    """A sample in both studies is fully replaced by its new rows."""
    old = pd.read_csv(old_path, sep="\t", dtype=str)
    new = pd.read_csv(new_path, sep="\t", dtype=str)
    if spec.id_column in old.columns and spec.id_column in new.columns:
        old = old[~old[spec.id_column].isin(set(new[spec.id_column]))]
    pd.concat([old, new], ignore_index=True).to_csv(out / spec.name, sep="\t", index=False)
    logger.info(f"{spec.name} updated!")


def _merge_columns(spec: StudyFile, old_path: Path, new_path: Path, out: Path) -> None:
    """A sample in both studies is fully replaced by its new column."""
    if spec.integer_values:
        old = pd.read_csv(old_path, sep="\t", index_col=0)
        new = pd.read_csv(new_path, sep="\t", index_col=0)
        merged = pd.concat([old.drop(columns=old.columns.intersection(new.columns)), new],
                           axis=1).replace(np.nan, 0).astype(int)
        merged.to_csv(out / spec.name, index=True, sep="\t")
    else:
        old = pd.read_csv(old_path, sep="\t")
        new = pd.read_csv(new_path, sep="\t")
        entity = [c for c in spec.entity_columns if c in old.columns and c in new.columns]
        old = old.drop(columns=[c for c in old.columns if c not in entity and c in new.columns])
        merged = old.merge(new, on=entity, how="outer") if entity else pd.concat([old, new], axis=1)
        merged.to_csv(out / spec.name, index=False, sep="\t")
    logger.info(f"{spec.name} updated!")


def _merge_clinical(old_path: Path, new_path: Path, out: Path, key: str) -> None:
    """Merge two cBioPortal clinical files (4 "#" header lines + header
    row) on `key`: the new study's header wins for shared attributes, a
    record in both studies takes its new values, attributes only one study
    has are kept."""
    def head(path: Path) -> pd.DataFrame:
        h = pd.read_csv(path, sep="\t", dtype=str, header=None, nrows=5)
        h.columns = h.iloc[4].tolist()
        return h.set_index(key)
    old_head, new_head = head(old_path), head(new_path)
    common = np.intersect1d(new_head.columns, old_head.columns)
    old_only_head = old_head.drop(list(common), axis=1)
    old_only_head.index = new_head.index
    final_head = pd.concat([new_head, old_only_head], axis=1).fillna(value=np.nan).reset_index()

    old_body = pd.read_csv(old_path, sep="\t", dtype=str, header=4)
    new_body = pd.read_csv(new_path, sep="\t", dtype=str, header=4)
    final_body = old_body.drop(list(common), axis=1).merge(new_body, how="outer", on=key)
    old_body = old_body.set_index(key)
    old_body = old_body.drop(index=old_body.index.intersection(new_body.set_index(key).index))
    final_body = final_body.set_index(key)
    final_body.update(old_body, overwrite=True, filter_func=None, errors="ignore")
    pd.concat([final_head, final_body.reset_index()], axis=0).to_csv(
        out / Path(old_path).name, header=False, index=False, sep="\t", na_rep="NaN")
    logger.info(f"{Path(old_path).name} updated!")


def _merge_sv(old_path: Path, new_path: Path, out: Path) -> None:
    """_merge_rows(), then drop duplicate events - keeping the copy with the
    highest supporting read count when only that differs."""
    old = pd.read_csv(old_path, sep="\t", dtype=str)
    new = pd.read_csv(new_path, sep="\t", dtype=str)
    old = old[~old["Sample_Id"].isin(set(old["Sample_Id"]) & set(new["Sample_Id"]))]
    merged = pd.concat([old, new], axis=0, join="outer", ignore_index=True)
    reads = "Normal_Paired_End_Read_Count"
    if reads in merged.columns:
        merged[reads] = pd.to_numeric(merged[reads], errors="coerce")
        merged = merged.sort_values(by=reads, ascending=False)
        merged = merged.drop_duplicates(
            subset=[c for c in merged.columns if c != reads], keep="first")
        merged[reads] = (merged[reads].astype(str).str.replace(r"\.0$", "", regex=True)
                         .replace("nan", ""))
    else:
        merged = merged.drop_duplicates(keep="first")
    merged.to_csv(out / "data_sv.txt", sep="\t", index=False)
    if "Class" in merged.columns:
        logger.info(
            f"data_sv.txt updated! ({(merged['Class'] == 'FUSION').sum()} fusion, "
            f"{(merged['Class'] == 'SPLICE').sum()} splice variant row(s) after "
            "merge/dedup)")
    else:
        logger.info("data_sv.txt updated!")


def merge_study_file(old_folder: str, new_folder: str, output: str, name: str) -> None:
    """Write `output`'s `name` from the two studies: merged when both have
    it, copied forward when only one does."""
    spec = next(s for s in STUDY_FILES if s.name == name)
    old_path, new_path = Path(old_folder) / name, Path(new_folder) / name
    out = Path(output)
    if old_path.exists() and new_path.exists():
        if spec.merge:
            spec.merge(old_path, new_path, out)
        elif spec.layout == "columns":
            _merge_columns(spec, old_path, new_path, out)
        else:
            _merge_rows(spec, old_path, new_path, out)
    elif old_path.exists():
        logger.warning(f"'{name}' only in old; copying forward.")
        shutil.copy(old_path, out / name)
    elif new_path.exists():
        logger.warning(f"'{name}' only in new; copying forward.")
        shutil.copy(new_path, out / name)
    else:
        logger.warning(f"'{name}' not found in either folder. Skipping.")


STUDY_FILES = [
    StudyFile("data_clinical_sample.txt", "clinical_sample",
              merge=lambda o, n, out: _merge_clinical(o, n, out, "SAMPLE_ID")),
    StudyFile("data_clinical_patient.txt", "clinical_patient",
              merge=lambda o, n, out: _merge_clinical(o, n, out, "PATIENT_ID")),
    StudyFile("data_cna_hg19.seg", "rows", id_column="ID"),
    StudyFile("data_cna_hg19.seg.fc.txt", "rows", id_column="ID"),
    StudyFile("data_cna.txt", "columns", entity_columns=("Hugo_Symbol",), integer_values=True),
    StudyFile("data_mutations_extended.txt", "rows", id_column="Tumor_Sample_Barcode"),
    StudyFile("data_sv.txt", "rows", id_column="Sample_Id", merge=_merge_sv),
    StudyFile("exon_CNA_data.txt", "rows", id_column="Sample_Id"),
    StudyFile("data_exon_brca_cna.txt", "columns", entity_columns=("ENTITY_STABLE_ID", "NAME")),
    StudyFile("data_gene_panel_matrix.txt", "rows", id_column="SAMPLE_ID"),
]
