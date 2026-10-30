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

"""Optional gene panel support: tell cBioPortal which genes each sample was
tested for.

Without it, cBioPortal assumes every sample profiled for a data type (see
the case lists, populate_case_lists.py) was tested for every gene - fine
for a single-assay study, wrong once samples from different assays share a
study (e.g. TSO500 and Guardant360 merged with -u): a gene outside a
sample's assay counts as "not altered" in that gene's frequency instead of
"not profiled".

Configured in conf.ini's [GenePanel] section - all optional; with nothing
set Varan works exactly as before, only logging what it can't do:

- MUTATIONS / CNA / SV: the gene panel id (cBioPortal stable_id) this run's
  samples were tested with, per data type. Writes
  data_gene_panel_matrix.txt + meta_gene_panel_matrix.txt: one row per
  sample, one column per profile (mutations / cna / structural_variants),
  holding that panel id - or NA where the sample wasn't tested for that
  data type (see populate_case_lists.PROFILED_SAMPLES_FILE). Every panel id
  used must already be loaded in the target cBioPortal instance, or the
  study import fails.
- PANEL_DIR: a folder with one plain gene list per panel, named
  <panel id>.txt (one HUGO symbol per line; "#" lines ignored). Optional:
  with it, Varan also writes each panel in cBioPortal's own gene panel
  file format (service/gene_panels/, for the portal admin to load) and
  warns about events in genes outside the sample's panel.

update/extract/remove carry data_gene_panel_matrix.txt forward like any
other per-sample file; finalize_gene_panel_matrix() then drops it with a
warning when it no longer covers every sample of the study (e.g. merging
a study with panels and one without) - cBioPortal needs every sample in it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from loguru import logger

from config_loader import get_config

MATRIX_FILE = "data_gene_panel_matrix.txt"
META_FILE = "meta_gene_panel_matrix.txt"

# Matrix column (the profile's stable_id in its meta file) -> conf.ini
# [GenePanel] key, profile name in PROFILED_SAMPLES_FILE, and the data
# file holding that profile.
PROFILES = {
    "mutations": ("MUTATIONS", "sequenced", "data_mutations_extended.txt"),
    "cna": ("CNA", "cna", "data_cna.txt"),
    "structural_variants": ("SV", "sv", "data_sv.txt"),
}

NO_PANEL_WARNING = (
    "No gene panel configured ([GenePanel] in conf.ini): cBioPortal will "
    "assume every sample was tested for every gene. Per-gene alteration "
    "frequencies, and the OncoPrint's 'not profiled' marks, will be wrong "
    "for genes outside a sample's assay - notably once studies from "
    "different assays are merged (e.g. Illumina TSO500 + Guardant360). To "
    "enable it, set [GenePanel] MUTATIONS / CNA / SV to the panel ids "
    "loaded in your cBioPortal, and optionally PANEL_DIR to a folder of "
    "<panel id>.txt gene lists.")


def configured_panels() -> dict[str, str]:
    """{matrix column: panel id} for every data type [GenePanel] sets."""
    config = get_config()
    panels = {}
    for column, (key, _, _) in PROFILES.items():
        value = config.get("GenePanel", key, fallback="").strip().strip('"').strip("'")
        if value:
            panels[column] = value
    return panels


def read_gene_list(panel_id: str) -> list[str] | None:
    """The genes of `panel_id` from [GenePanel] PANEL_DIR/<panel_id>.txt,
    or None if there's no such file."""
    panel_dir = get_config().get("GenePanel", "PANEL_DIR", fallback="").strip()
    if not panel_dir:
        return None
    path = Path(panel_dir) / f"{panel_id}.txt"
    if not path.is_file():
        return None
    genes = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.extend(g for g in line.replace(",", "\t").split("\t") if g.strip())
    return list(dict.fromkeys(g.strip() for g in genes))


def build_gene_panel_matrix(folder: str, profiled: dict[str, list[str]]) -> None:
    """At study creation: write data_gene_panel_matrix.txt from [GenePanel]
    and the samples tested per data type, export the panels' gene lists,
    and check events against them. Only logs NO_PANEL_WARNING when no
    panel is configured."""
    panels = configured_panels()
    if not panels:
        logger.warning(NO_PANEL_WARNING)
        return

    sample_ids = _study_sample_ids(folder)
    if not sample_ids:
        logger.warning(
            f"{MATRIX_FILE} not written: no SAMPLE_ID found in "
            "data_clinical_sample.txt.")
        return

    matrix = pd.DataFrame({"SAMPLE_ID": sample_ids})
    for column, panel_id in panels.items():
        tested = set(profiled.get(PROFILES[column][1], []))
        matrix[column] = [panel_id if s in tested else "NA" for s in sample_ids]
    matrix.to_csv(Path(folder) / MATRIX_FILE, sep="\t", index=False)
    logger.info(
        f"Wrote {MATRIX_FILE}: "
        + ", ".join(f"{c} = {p}" for c, p in panels.items()))

    for panel_id in dict.fromkeys(panels.values()):
        genes = read_gene_list(panel_id)
        if genes is None:
            logger.warning(
                f"No gene list for panel '{panel_id}' (expected "
                f"<[GenePanel] PANEL_DIR>/{panel_id}.txt): the study still "
                "references it - it must already be loaded in cBioPortal - "
                "but Varan can't export it or check events against it.")
            continue
        _export_panel(folder, panel_id, genes)
    check_events_in_panels(folder)


def _export_panel(folder: str, panel_id: str, genes: list[str]) -> None:
    """Write `panel_id` in cBioPortal's gene panel file format, under the
    study's service/gene_panels/ (not part of the study import itself:
    panels are loaded into the portal once, by its admin)."""
    out_dir = Path(folder) / "service" / "gene_panels"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"data_gene_panel_{panel_id}.txt"
    out.write_text(
        f"stable_id: {panel_id}\n"
        f"description: {panel_id} ({len(genes)} genes, exported by Varan)\n"
        f"gene_list: {chr(9).join(genes)}\n", encoding="utf-8")
    logger.info(f"Gene panel '{panel_id}' ({len(genes)} genes) exported to {out}")


def check_events_in_panels(folder: str) -> None:
    """Warn about events in genes outside the sample's panel for that data
    type - usually a wrong panel assignment or a gene-name mismatch."""
    matrix_path = Path(folder) / MATRIX_FILE
    if not matrix_path.exists():
        return
    matrix = _read_matrix(matrix_path).set_index("SAMPLE_ID")
    gene_lists: dict[str, set[str] | None] = {}

    for column, (_, _, data_file) in PROFILES.items():
        path = Path(folder) / data_file
        if column not in matrix.columns or not path.exists():
            continue
        for sample_id, genes in _event_genes(path, column).items():
            panel_id = matrix[column].get(sample_id)
            if not panel_id or panel_id == "NA":
                continue
            if panel_id not in gene_lists:
                listed = read_gene_list(panel_id)
                gene_lists[panel_id] = set(listed) if listed is not None else None
            panel_genes = gene_lists[panel_id]
            if panel_genes is None:
                continue
            outside = sorted(genes - panel_genes)
            if outside:
                logger.warning(
                    f"{sample_id}: {column} event(s) in gene(s) outside panel "
                    f"'{panel_id}': {', '.join(outside)}. Check the panel "
                    "assigned to this sample, or the gene names.")


def _event_genes(path: Path, column: str) -> dict[str, set[str]]:
    """{sample: genes with at least one event} for one data file."""
    events: dict[str, set[str]] = {}
    try:
        if column == "mutations":
            from filter_clinvar import read_maf
            df = read_maf(path)
            for sample_id, gene in zip(df["Tumor_Sample_Barcode"], df["Hugo_Symbol"]):
                events.setdefault(str(sample_id), set()).add(str(gene))
        elif column == "cna":
            df = pd.read_csv(path, sep="\t", dtype=str).set_index("Hugo_Symbol")
            for sample_id in df.columns:
                altered = df.index[pd.to_numeric(df[sample_id], errors="coerce").fillna(0) != 0]
                events[str(sample_id)] = {str(g) for g in altered}
        else:
            df = pd.read_csv(path, sep="\t", dtype=str)
            for _, row in df.iterrows():
                genes = {str(row.get(c)) for c in ("Site1_Hugo_Symbol", "Site2_Hugo_Symbol")
                         if pd.notna(row.get(c)) and str(row.get(c)).strip()}
                events.setdefault(str(row["Sample_Id"]), set()).update(genes)
    except (KeyError, pd.errors.EmptyDataError) as err:
        logger.warning(f"Could not check {path.name} against gene panels: {err}")
    return events


def finalize_gene_panel_matrix(folder: str) -> None:
    """Keep data_gene_panel_matrix.txt consistent with the study it ends up
    in (any flow: create/update/extract/remove): drop columns whose profile
    has no data file anymore, and drop the whole matrix - with a warning -
    when it doesn't list every sample of the study, since cBioPortal needs
    each one in it (e.g. after merging a study with panels and one
    without: the panels of the latter's samples aren't known)."""
    path = Path(folder) / MATRIX_FILE
    if not path.exists():
        return
    matrix = _read_matrix(path)
    keep = ["SAMPLE_ID"] + [
        c for c in matrix.columns
        if c in PROFILES and (Path(folder) / PROFILES[c][2]).exists()]
    matrix = matrix[keep]
    # update/extract/remove's generic per-sample functions read "NA" as a
    # missing value and write it back blank; cBioPortal wants "NA".
    matrix[keep[1:]] = matrix[keep[1:]].replace("", "NA")

    missing = sorted(set(_study_sample_ids(folder)) - set(matrix["SAMPLE_ID"]))
    if missing or len(keep) == 1:
        path.unlink()
        reason = (f"no gene panel known for sample(s) {', '.join(missing)}"
                  if missing else "none of its data types is in the study")
        logger.warning(
            f"{MATRIX_FILE} dropped: {reason}. cBioPortal will assume every "
            "sample was tested for every gene. To keep gene panels, every "
            "sample must come from a study created with [GenePanel] set.")
        return
    matrix.to_csv(path, sep="\t", index=False)


def create_meta_gene_panel_matrix(project_id: str, folder: str) -> None:
    """Write meta_gene_panel_matrix.txt for an existing matrix file."""
    if not (Path(folder) / MATRIX_FILE).exists():
        return
    (Path(folder) / META_FILE).write_text(
        f"cancer_study_identifier: {project_id}\n"
        "genetic_alteration_type: GENE_PANEL_MATRIX\n"
        "datatype: GENE_PANEL_MATRIX\n"
        f"data_filename: {MATRIX_FILE}\n", encoding="utf-8")
    logger.info(f"{META_FILE} created!")


def _read_matrix(path: Path) -> pd.DataFrame:
    """Read data_gene_panel_matrix.txt keeping "NA" as text - cBioPortal's
    value for "not profiled" - instead of turning it into a blank cell."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _study_sample_ids(folder: str) -> list[str]:
    path = Path(folder) / "data_clinical_sample.txt"
    if not path.exists():
        return []
    try:
        clin = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    except pd.errors.EmptyDataError:
        return []
    if "SAMPLE_ID" not in clin.columns:
        return []
    return list(dict.fromkeys(clin["SAMPLE_ID"].dropna().str.strip()))
