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

"""Module for generating case lists from SV,CNA and mutation data.

Each function reads the input data file from a project directory,
extracts sample IDs, and writes a formatted `cases_*.txt` file into
the output directory for use in cBioPortal.

Functions:
- populate_cases_sv
- populate_cases_cna
- populate_cases_sequenced

cBioPortal treats these case lists as the samples *profiled* for each data
type - the denominators of every alteration frequency - so a sample tested
for a data type but without any event of that type must still be listed.
The data files alone can't say that (data_mutations_extended.txt and
data_sv.txt only have rows for events, data_cna.txt only columns for
samples with calls), so each populate_cases_*() also takes the samples
known to be profiled - see PROFILED_SAMPLES_FILE and
profiled_from_case_lists() - and lists the union of both.

"""

from pathlib import Path

import pandas as pd
from loguru import logger

# Written by walk.py at study creation: one "<profile>\t<SAMPLE_ID>" line
# per sample tested for that data type, whether or not it had any event -
# profile being one of PROFILE_CASE_LISTS' keys. Read by meta_case_main().
PROFILED_SAMPLES_FILE = Path("service") / "profiled_samples.tsv"

# profile -> the case list file it feeds.
PROFILE_CASE_LISTS = {
    "sequenced": "cases_sequenced.txt",
    "cna": "cases_cna.txt",
    "sv": "cases_sv.txt",
}


def write_profiled_samples(folder: str, profiled: dict[str, list[str]]) -> None:
    """Write PROFILED_SAMPLES_FILE under the study folder `folder`."""
    path = Path(folder) / PROFILED_SAMPLES_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for profile, sample_ids in profiled.items():
            for sample_id in sample_ids:
                f.write(f"{profile}\t{sample_id}\n")


def read_profiled_samples(folder: str) -> dict[str, list[str]]:
    """Read PROFILED_SAMPLES_FILE from the study folder `folder` ({} if the
    study has none, e.g. one made by update/extract/remove)."""
    path = Path(folder) / PROFILED_SAMPLES_FILE
    profiled: dict[str, list[str]] = {}
    if not path.exists():
        return profiled
    for line in path.read_text().splitlines():
        if "\t" in line:
            profile, sample_id = line.split("\t", 1)
            profiled.setdefault(profile, []).append(sample_id)
    return profiled


def profiled_from_case_lists(study_folders: list[str]) -> dict[str, list[str]]:
    """Profiled samples per data type, as listed by the existing case lists
    of `study_folders` - how update/extract/remove carry forward the
    samples tested without events, which their merged/filtered data files
    no longer show."""
    profiled: dict[str, list[str]] = {}
    for folder in study_folders:
        for profile, case_file in PROFILE_CASE_LISTS.items():
            path = Path(folder) / "case_lists" / case_file
            if not path.exists():
                continue
            for line in path.read_text().splitlines():
                if line.startswith("case_list_ids:"):
                    ids = profiled.setdefault(profile, [])
                    ids.extend(s.strip() for s in line.split(":", 1)[1].split("\t")
                               if s.strip() and s.strip() not in ids)
    return profiled


def _study_sample_ids(folder: str) -> set[str] | None:
    """SAMPLE_IDs in the study's data_clinical_sample.txt, or None if it
    has none to check against."""
    path = Path(folder) / "data_clinical_sample.txt"
    if not path.exists():
        return None
    try:
        clin = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    except pd.errors.EmptyDataError:
        return None
    if "SAMPLE_ID" not in clin.columns:
        return None
    return set(clin["SAMPLE_ID"].dropna().str.strip())


# data_sv.txt's header as Varan writes it (see walk.fill_fusion_from_combined).
SV_HEADER = ("Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\tSite2_Hugo_Symbol\t"
             "Normal_Paired_End_Read_Count\tEvent_Info\tRNA_Support\n")


def profiled_in_study(folder: str, profiled: list[str] | None) -> list[str]:
    """The samples of `profiled` that the study still contains."""
    in_study = _study_sample_ids(folder)
    return [s for s in dict.fromkeys(profiled or [])
            if in_study is None or s in in_study]


def ensure_sv_file(folder: str, profiled: list[str] | None) -> None:
    """Write a header-only data_sv.txt when the study has samples tested for
    fusions/splice variants but none had any: without the file there's no
    structural-variant profile at all, so cBioPortal can't show "0 of N
    tested" - an empty file plus cases_sv.txt is how it says that (checked
    with cBioPortal's own validator). No-op when data_sv.txt exists."""
    path = Path(folder) / "data_sv.txt"
    if path.exists():
        return
    tested = profiled_in_study(folder, profiled)
    if not tested:
        return
    path.write_text(SV_HEADER)
    logger.info(
        f"No fusions or splice variants in this study, but {len(tested)} "
        "sample(s) were tested for them: wrote a header-only data_sv.txt so "
        "they're shown as tested without events.")


# Data files carrying a cBioPortal profile, by PROFILE_CASE_LISTS key.
PROFILE_DATA_FILES = {
    "sequenced": "data_mutations_extended.txt",
    "cna": "data_cna.txt",
    "sv": "data_sv.txt",
}
OTHER_DATA_FILES = ["data_cna_hg19.seg", "data_cna_hg19.seg.fc.txt",
                    "data_exon_brca_cna.txt", "exon_CNA_data.txt"]


def drop_unusable_data_files(folder: str, profiled: dict[str, list[str]]) -> None:
    """Remove only the data files cBioPortal couldn't load - empty files
    without events are otherwise kept (see ensure_sv_file()):

    - a 0-byte file (not even a header line: an interrupted write);
    - a mutation/CNA/SV file without any event when no sample of the study
      was tested for that data type: its profile would need a case list
      (cases_sequenced/cna/sv), and that list would be empty.
    """
    for name in [*PROFILE_DATA_FILES.values(), *OTHER_DATA_FILES]:
        path = Path(folder) / name
        if path.exists() and path.stat().st_size == 0:
            path.unlink()
            logger.warning(f"{name} was an empty (0-byte) file - removed.")
    for profile, name in PROFILE_DATA_FILES.items():
        path = Path(folder) / name
        if not path.exists() or profiled_in_study(folder, profiled.get(profile)):
            continue
        if profile == "sequenced":
            from filter_clinvar import read_maf
            no_events = read_maf(path).empty
        else:
            no_events = pd.read_csv(path, sep="\t", dtype=str).empty
        if no_events:
            path.unlink()
            logger.info(
                f"{name} has no events and no sample of this study was tested "
                "for that data type - removed (cBioPortal would need a case "
                "list of tested samples for it).")


def ensure_cna_file(folder: str, profiled: list[str] | None) -> None:
    """Same as ensure_sv_file(), for CNAs: when the study has samples tested
    for CNVs but data_cna.txt is gone - no sample had any call, or
    extract/remove only kept samples without one - write it with just its
    header (Hugo_Symbol + one column per tested sample) and no gene row,
    so the CNA profile still exists with 0 alterations over them (checked
    with cBioPortal's own validator). No-op when data_cna.txt exists."""
    path = Path(folder) / "data_cna.txt"
    if path.exists():
        return
    tested = profiled_in_study(folder, profiled)
    if not tested:
        return
    path.write_text("\t".join(["Hugo_Symbol", *tested]) + "\n")
    logger.info(
        f"No CNV call in this study, but {len(tested)} sample(s) were tested "
        "for CNVs: wrote a data_cna.txt with only its header so they're "
        "shown as tested without alterations.")


def pad_unaltered_cna_samples(folder: str, profiled: list[str] | None) -> None:
    """Add an all-0 column to data_cna.txt for every CNV-profiled sample
    of the study that has no call at all: cBioPortal expects the discrete
    CNA file to contain every profiled sample, and data_cna.txt is only
    built from the calls. Leaves the file untouched when nothing's missing.
    """
    path = Path(folder) / "data_cna.txt"
    if not profiled or not path.exists():
        return
    try:
        data_cna = pd.read_csv(path, sep="\t", dtype=str)
    except pd.errors.EmptyDataError:
        return
    in_study = _study_sample_ids(folder)
    missing = [s for s in dict.fromkeys(profiled)
               if s not in data_cna.columns and (in_study is None or s in in_study)]
    if not missing:
        return
    for sample_id in missing:
        data_cna[sample_id] = "0"
    data_cna.to_csv(path, sep="\t", index=False)
    logger.info(
        f"Added {len(missing)} CNV-profiled sample(s) without any call to "
        f"data_cna.txt as all-0 columns: {', '.join(missing)}")


def with_profiled(sample_ids: list[str], profiled: list[str] | None,
                  folder: str) -> list[str]:
    """`sample_ids` (the samples with events, from the data file) followed
    by every other profiled sample still in the study. Samples with events
    are never dropped; profiled ones are only added if the study's
    data_clinical_sample.txt defines them (e.g. not after `remove`)."""
    if not profiled:
        return sample_ids
    in_study = _study_sample_ids(folder)
    result = list(sample_ids)
    seen = {str(s).strip() for s in sample_ids}
    for sample_id in profiled:
        if sample_id not in seen and (in_study is None or sample_id in in_study):
            result.append(sample_id)
            seen.add(sample_id)
    return result


def populate_cases_sv(project_id: str, folder: str, cases_list_dir: str,
                      logger: logger, profiled: list[str] | None = None) -> None:
    """Populate the cases_sv.txt file based on SV data.

    Parameters
    ----------
    project_id : str
        Unique identifier for the cancer study.
    folder : str
        Path to the directory containing the `data_sv.txt` file.
    cases_list_dir : str
        Output directory where the `cases_sv.txt` will be saved.
    logger : logging.Logger
        Logger object for error reporting.
    profiled : list[str], optional
        Samples tested for fusions/splice variants, listed even without
        any event - see with_profiled().

    """
    try:
        data_sv = pd.read_csv(Path(folder) / "data_sv.txt", sep="\t", dtype={"Sample_Id": str})
    except pd.errors.EmptyDataError:
        logger.exception("data_sv.txt is empty, skipping this step!")
        return()

    # Filter out NaN before counting, not just when building case_list_ids
    # below - otherwise a NaN Sample_Id (e.g. from a malformed splice-variant
    # row) counts as "1 sample" in nsamples while being excluded from the
    # actual listed IDs, so the reported count and the list length disagree.
    sample_ids = [s for s in data_sv.Sample_Id.unique() if pd.notna(s)]
    sample_ids = with_profiled(sample_ids, profiled, folder)
    nsamples = len(sample_ids)

    stable_id = project_id + "_sv"
    case_list_name = "Samples with SV data"
    case_list_category = "all_cases_with_sv_data"
    case_list_description = "Samples with SV data (" + str(nsamples) + " sample(s))"
    case_list_ids = "\t".join([str(s).strip() for s in sample_ids])

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "stable_id": stable_id,
        "case_list_category": case_list_category,
        "case_list_name": case_list_name,
        "case_list_description": case_list_description,
        "case_list_ids": case_list_ids,
    }

    with (Path(cases_list_dir) / "cases_sv.txt").open("w") as case_sv_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=case_sv_file)

    return()


def populate_cases_cna(project_id: str, folder: str, cases_list_dir: str,
                       logger: logger, profiled: list[str] | None = None) -> None:
    """Populate the cases_cna.txt file based on CNA data.

    Parameters
    ----------
    project_id : str
        Unique identifier for the cancer study.
    folder : str
        Path to the directory containing the `data_cna.txt` file.
    cases_list_dir : str
        Output directory where the `cases_cna.txt` will be saved.
    logger : logging.Logger
        Logger object for error reporting.
    profiled : list[str], optional
        Samples tested for CNV, listed even without any call - see
        with_profiled() and pad_unaltered_cna_samples().

    """
    try:
        data_cna = pd.read_csv(Path(folder) / "data_cna.txt", sep="\t", dtype=str)
    except pd.errors.EmptyDataError:
        logger.error("data_cna.txt is empty, skipping this step!")
        return

    sample_ids = with_profiled(list(data_cna.columns)[1:], profiled, folder)
    nsamples = len(sample_ids)

    stable_id = project_id + "_cna"

    case_list_category = "all_cases_with_cna_data"
    case_list_name = "Samples with CNA data"
    case_list_description = "Samples with CNA data (" + str(nsamples) + " sample(s))"
    case_list_ids = "\t".join([str(s) for s in sample_ids])

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "stable_id": stable_id,
        "case_list_category": case_list_category,
        "case_list_name": case_list_name,
        "case_list_description": case_list_description,
        "case_list_ids": case_list_ids,
    }

    with (Path(cases_list_dir) / "cases_cna.txt").open("w") as case_cna_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=case_cna_file)


def populate_cases_sequenced(project_id: str, folder: str, cases_list_dir: str,
                             logger: logger, profiled: list[str] | None = None) -> None:
    """Populate the cases_sequenced.txt file based on mutation data.

    Parameters
    ----------
    project_id : str
        Unique identifier for the cancer study.
    folder : str
        Path to the directory containing the `data_mutations_extended.txt` file.
    cases_list_dir : str
        Output directory where the `cases_sequenced.txt` will be saved.
    logger : logging.Logger
        Logger object for error reporting.
    profiled : list[str], optional
        Samples sequenced for SNVs/indels, listed even without any
        mutation left after filtering - see with_profiled().

    """
    try:
        data_sequenced = pd.read_csv(Path(folder) / "data_mutations_extended.txt",
                                    sep="\t", low_memory=False, dtype={"Tumor_Sample_Barcode": str})
    except pd.errors.EmptyDataError:
        logger.error("data_mutations_extended.txt is empty, skipping this step!")
        return

    sample_ids = with_profiled(
        list(data_sequenced["Tumor_Sample_Barcode"].unique()), profiled, folder)
    nsamples = len(sample_ids)

    stable_id = project_id + "_sequenced"

    case_list_category = "all_cases_with_mutation_data"
    case_list_name = "Sequenced Tumors"
    case_list_description = "All sequenced samples (" + str(nsamples) + " sample(s))"
    case_list_ids = "\t".join([str(s) for s in sample_ids])

    dictionary_file = {
        "cancer_study_identifier": project_id,
        "stable_id": stable_id,
        "case_list_category": case_list_category,
        "case_list_name": case_list_name,
        "case_list_description": case_list_description,
        "case_list_ids": case_list_ids,
    }

    with (Path(cases_list_dir) / "cases_sequenced.txt").open("w") as case_sequence_file:
        for key, value in dictionary_file.items():
            print(f"{key}: {value}", file=case_sequence_file)
