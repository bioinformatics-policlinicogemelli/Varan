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

"""OncoKB actionability per sample, and the derived case lists.

With OncoKB annotation (-k), each sample gets (data_clinical_sample.txt):
- ONCOKB_HIGHEST_LEVEL: highest therapeutic sensitivity level among its
  mutations, fusions and amplifications/deep deletions (1, 2, 3A, 3B, 4;
  "None" if annotated without a level, NA if not annotated).
- ONCOKB_RESISTANCE_LEVEL: highest resistance level (R1, R2).
- ONCOKB_N_ACTIONABLE / ONCOKB_ACTIONABLE: the alterations at an
  [Actionability] ONCOKB_LEVELS level (default 1 and 2, standard care).
Levels are those OncoKB gave for the sample's ONCOTREE_CODE. The OncoKB data
version is fetched when annotating and kept in service/oncokb_info.json (a
study merged from batches annotated with different versions lists them all).
CNV levels are kept in service/oncokb_cna_levels.tsv, since data_cna.txt
can't hold them; a level counts only if data_cna.txt still has that call.

Case lists (when they have samples): OncoKB actionable, CIViC level A-B
(civic.py), tumour-agnostic biomarker (tumor_agnostic.py), MSI-H, TMB-high
(the lab's TMB_THR category "High"), Tissue, Plasma.

Reference: Chakravarty et al., JCO Precis Oncol 2017, doi:10.1200/PO.17.00011.
"""

from __future__ import annotations

import json
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd
from loguru import logger

import study_alterations as sa
from clinical_file import drop_sample_columns, set_sample_column
from config_loader import config

INFO_FILE = Path("service") / "oncokb_info.json"
CNA_LEVELS_FILE = Path("service") / "oncokb_cna_levels.tsv"
SUMMARY_FILE = Path("service") / "actionability.json"
CNA_ANNOTATED = Path("intermediate") / "02_CNV_annotated_merged.txt"
SENSITIVE = ["1", "2", "3A", "3B", "4"]
RESISTANCE = ["R1", "R2"]
COLUMNS = {
    "ONCOKB_HIGHEST_LEVEL": ("OncoKB Highest Level",
                             "Highest OncoKB therapeutic level (sensitivity) of the sample's alterations",
                             "STRING"),
    "ONCOKB_RESISTANCE_LEVEL": ("OncoKB Resistance Level",
                                "Highest OncoKB resistance level of the sample's alterations", "STRING"),
    "ONCOKB_N_ACTIONABLE": ("OncoKB Actionable Alterations",
                            "Alterations at the OncoKB levels set in [Actionability] ONCOKB_LEVELS",
                            "NUMBER"),
    "ONCOKB_ACTIONABLE": ("OncoKB Actionable", "Those alterations, with their level", "STRING"),
}
CASE_LIST_PREFIX = "cases_varan_"


def actionable_levels() -> list[str]:
    value = config.get("Actionability", "ONCOKB_LEVELS", fallback="").strip() or "1,2"
    return [v.strip().upper().replace("LEVEL_", "") for v in value.split(",") if v.strip()]


def _level(value: object) -> str:
    text = str(value or "").strip().upper()
    return text.replace("LEVEL_", "") if text.startswith("LEVEL_") else ""


def _best(levels: list[str], order: list[str]) -> str:
    present = [lv for lv in order if lv in levels]
    return present[0] if present else ""


# ----------------------------------------------------- at annotation

def record_oncokb_info(folder: str | Path) -> None:
    """Fetch OncoKB's data version (needs the token and the network) into
    service/oncokb_info.json. A failure is a warning, not an error."""
    token = config.get("OncoKB", "ONCOKB", fallback="").strip()
    entry = {"data_version": "unknown", "data_date": "", "retrieved": date.today().isoformat()}
    try:
        request = urllib.request.Request(
            "https://www.oncokb.org/api/v1/info", headers={"Authorization": f"Bearer {token}"})
        with urllib.request.urlopen(request, timeout=30) as response:
            info = json.loads(response.read().decode("utf-8"))
        entry["data_version"] = info.get("dataVersion", {}).get("version", "unknown")
        entry["data_date"] = info.get("dataVersion", {}).get("date", "")
    except Exception as err:  # noqa: BLE001 - recorded as unknown, never fatal
        logger.warning(f"Couldn't read OncoKB's data version ({err}): recorded as unknown.")
    path = Path(folder) / INFO_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"versions": [entry]}, indent=2), encoding="utf-8")
    logger.info(f"OncoKB data version: {entry['data_version']} {entry['data_date']}")


def save_cna_levels(folder: str | Path) -> None:
    """At creation: the OncoKB levels of the CNV calls (intermediate/) into
    service/, where update/extract/remove can still find them."""
    source = Path(folder) / CNA_ANNOTATED
    if not source.exists():
        return
    table = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    keep = [c for c in ("HIGHEST_SENSITIVE_LEVEL", "HIGHEST_RESISTANCE_LEVEL", "ONCOGENIC") if c in table]
    out = pd.DataFrame({
        "SAMPLE_ID": table["Tumor_Sample_Barcode"].str.replace(".cnv.bam", "", regex=False),
        "GENE": table["Hugo_Symbol"],
        "VALUE": table["Copy_Number_Alteration"],
        **{c: table[c] for c in keep}})
    out = out[out["VALUE"].isin(["2", "-2"])]
    path = Path(folder) / CNA_LEVELS_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)


def carry_service(sources: list[str | Path], output: str | Path) -> None:
    """update/extract/remove: OncoKB versions and CNV levels of the sources."""
    versions, frames = [], []
    for source in sources:
        info = Path(source) / INFO_FILE
        if info.exists():
            for entry in json.loads(info.read_text(encoding="utf-8")).get("versions", []):
                if entry not in versions:
                    versions.append(entry)
        levels = Path(source) / CNA_LEVELS_FILE
        if levels.exists():
            frames.append(pd.read_csv(levels, sep="\t", dtype=str, keep_default_na=False))
    out = Path(output)
    if versions:
        (out / INFO_FILE).parent.mkdir(parents=True, exist_ok=True)
        (out / INFO_FILE).write_text(json.dumps({"versions": versions}, indent=2), encoding="utf-8")
    if frames:
        merged = pd.concat(frames).drop_duplicates(["SAMPLE_ID", "GENE", "VALUE"], keep="last")
        (out / CNA_LEVELS_FILE).parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out / CNA_LEVELS_FILE, sep="\t", index=False)


def oncokb_versions(folder: str | Path) -> list[dict]:
    path = Path(folder) / INFO_FILE
    if not path.exists():
        return []
    return json.loads(path.read_text(encoding="utf-8")).get("versions", [])


# ------------------------------------------------------------- study

def annotate_study(folder: str | Path) -> dict | None:
    """Write the sample attributes described in the module docstring; None
    (and the attributes removed) when the study has no OncoKB annotation."""
    folder = Path(folder)
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    samples = _samples(folder)
    found: dict[str, list[tuple[str, str, str]]] = {}   # sample -> [(label, sens, res)]
    annotated: set[str] = set()
    has_oncokb = False

    maf = sa.read_maf(folder)
    if maf is not None and "HIGHEST_SENSITIVE_LEVEL" in maf.columns:
        has_oncokb = True
        for r in maf.to_dict("records"):
            sample = str(r["Tumor_Sample_Barcode"]).strip()
            if _was_annotated(r):
                annotated.add(sample)
                found.setdefault(sample, []).append((
                    sa.mutation_label(r), _level(r.get("HIGHEST_SENSITIVE_LEVEL")),
                    _level(r.get("HIGHEST_RESISTANCE_LEVEL"))))
    sv = sa.fusions(sa.read_sv(folder))
    if "HIGHEST_SENSITIVE_LEVEL" in sv.columns:
        has_oncokb = True
        for r in sv.to_dict("records"):
            sample = str(r["Sample_Id"]).strip()
            annotated.add(sample)
            found.setdefault(sample, []).append((
                sa.fusion_label(r), _level(r.get("HIGHEST_SENSITIVE_LEVEL")),
                _level(r.get("HIGHEST_RESISTANCE_LEVEL"))))
    levels_path = folder / CNA_LEVELS_FILE
    if levels_path.exists():
        has_oncokb = True
        levels = pd.read_csv(levels_path, sep="\t", dtype=str, keep_default_na=False)
        calls = {(c["sample"], c["gene"], str(c["value"])) for c in sa.read_cna(folder)}
        for r in levels.to_dict("records"):
            annotated.add(r["SAMPLE_ID"])
            if (r["SAMPLE_ID"], r["GENE"], r["VALUE"]) in calls:
                found.setdefault(r["SAMPLE_ID"], []).append((
                    sa.cna_label({"gene": r["GENE"], "value": int(r["VALUE"])}),
                    _level(r.get("HIGHEST_SENSITIVE_LEVEL")),
                    _level(r.get("HIGHEST_RESISTANCE_LEVEL"))))

    if not has_oncokb or samples is None:
        drop_sample_columns(folder, list(COLUMNS))
        return None

    wanted = actionable_levels()
    attrs: dict[str, dict[str, object]] = {c: {} for c in COLUMNS}
    per_sample = {}
    maf_samples = set() if maf is None else set(maf["Tumor_Sample_Barcode"].astype(str).str.strip())
    for sample in samples["SAMPLE_ID"]:
        alterations = found.get(sample, [])
        if sample not in annotated and sample in maf_samples:
            # Variants without OncoKB annotation: merged from a batch run
            # without -k. A sample without any alteration is "None".
            for c in COLUMNS:
                attrs[c][sample] = "NA"
            continue
        best = _best([s for _, s, _ in alterations], SENSITIVE)
        resistance = _best([r for _, _, r in alterations], RESISTANCE)
        actionable = sorted({f"{label} ({s})" for label, s, _ in alterations if s in wanted})
        attrs["ONCOKB_HIGHEST_LEVEL"][sample] = best or "None"
        attrs["ONCOKB_RESISTANCE_LEVEL"][sample] = resistance or "None"
        attrs["ONCOKB_N_ACTIONABLE"][sample] = len(actionable)
        attrs["ONCOKB_ACTIONABLE"][sample] = ", ".join(actionable) or "None"
        per_sample[sample] = {"highest": best or "None", "resistance": resistance or "None",
                              "actionable": actionable}

    for column, (short, long, datatype) in COLUMNS.items():
        set_sample_column(folder, column, attrs[column], short, long, datatype)

    summary = cohort_summary(samples, per_sample, wanted, SENSITIVE)
    summary["versions"] = oncokb_versions(folder)
    (folder / SUMMARY_FILE).parent.mkdir(parents=True, exist_ok=True)
    (folder / SUMMARY_FILE).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info(f"OncoKB actionability: {summary['n_actionable']}/{summary['n_samples']} sample(s) "
                f"with a level {'/'.join(wanted)} alteration.")
    return summary


def _was_annotated(row: dict) -> bool:
    flag = str(row.get("ANNOTATED", "")).strip().lower()
    if flag:
        return flag == "true"
    return any(str(row.get(c, "")).strip() not in ("", "nan") for c in
               ("HIGHEST_SENSITIVE_LEVEL", "ONCOGENIC", "MUTATION_EFFECT"))


def _samples(folder: Path) -> pd.DataFrame | None:
    import specimen
    table = specimen.read_study_samples(folder)
    if table is None:
        return None
    table = table.copy()
    table["SAMPLE_ID"] = table["SAMPLE_ID"].astype(str).str.strip()
    return table


def cohort_summary(samples: pd.DataFrame, per_sample: dict[str, dict], wanted: list[str],
                   order: list[str]) -> dict:
    """Counts for the report: samples with an actionable alteration, by
    highest level and by tumour type (ONCOTREE_CODE)."""
    evaluated = [s for s in samples["SAMPLE_ID"] if s in per_sample]
    with_hit = [s for s in evaluated if per_sample[s]["actionable"]]
    by_level = {lv: sum(1 for s in evaluated if per_sample[s]["highest"] == lv) for lv in order}
    by_level["None"] = sum(1 for s in evaluated if per_sample[s]["highest"] == "None")
    tumour = dict(zip(samples["SAMPLE_ID"],
                      samples.get("ONCOTREE_CODE", pd.Series("NA", index=samples.index)).fillna("NA")))
    by_tumour: dict[str, list[int]] = {}
    for s in evaluated:
        counts = by_tumour.setdefault(str(tumour.get(s, "NA")), [0, 0])
        counts[0] += 1
        counts[1] += bool(per_sample[s]["actionable"])
    return {"levels": wanted, "n_samples": len(evaluated), "n_total": len(samples),
            "n_actionable": len(with_hit),
            "by_level": by_level, "by_tumour": by_tumour, "samples": per_sample}


def read_summary(folder: str | Path) -> dict | None:
    path = Path(folder) / SUMMARY_FILE
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


# -------------------------------------------------------- case lists

def write_case_lists(project_id: str, folder: str | Path) -> None:
    """The derived case lists (module docstring), from
    data_clinical_sample.txt; one without samples isn't written."""
    folder = Path(folder)
    cases = folder / "case_lists"
    for old in cases.glob(f"{CASE_LIST_PREFIX}*.txt"):
        old.unlink()
    path = folder / "data_clinical_sample.txt"
    if not path.exists():
        return
    table = pd.read_csv(path, sep="\t", header=4, dtype=str, keep_default_na=False)

    def where(column: str, test) -> list[str]:
        if column not in table.columns:
            return []
        return [s for s, v in zip(table["SAMPLE_ID"], table[column]) if test(str(v).strip())]

    def number(value: str) -> float:
        try:
            return float(value)
        except ValueError:
            return 0.0

    civic_levels = [lv.strip().upper() for lv in config.get(
        "Actionability", "CIVIC_LEVELS", fallback="").split(",") if lv.strip()] or ["A", "B"]
    lists = [
        ("oncokb_actionable", "OncoKB actionable",
         f"Samples with an alteration at OncoKB level {'/'.join(actionable_levels())}",
         where("ONCOKB_N_ACTIONABLE", lambda v: number(v) > 0)),
        ("civic_actionable", "CIViC level " + "/".join(civic_levels),
         "Samples with an alteration with CIViC predictive evidence at level "
         f"{'/'.join(civic_levels)} (any tumour type)",
         where("CIVIC_N_ACTIONABLE", lambda v: number(v) > 0)),
        ("tumor_agnostic", "Tumour-agnostic biomarker",
         "Samples with a biomarker of a tumour-agnostic approval (see report)",
         where("TUMOR_AGNOSTIC_BIOMARKERS", lambda v: v not in ("", "NA", "None detected"))),
        ("msi_h", "MSI-H", "Samples with MSI status Unstable",
         where("MSI_THR", lambda v: v.lower() == "unstable")),
        ("tmb_high", "TMB-high", "Samples in the High TMB category (lab thresholds, TMB_THR)",
         where("TMB_THR", lambda v: v.lower() == "high")),
        ("tissue", "Tissue samples", "Samples with SPECIMEN_TYPE Tissue",
         where("SPECIMEN_TYPE", lambda v: v == "Tissue")),
        ("plasma", "Plasma samples", "Samples with SPECIMEN_TYPE Plasma",
         where("SPECIMEN_TYPE", lambda v: v == "Plasma")),
    ]
    for key, name, description, sample_ids in lists:
        if not sample_ids:
            continue
        cases.mkdir(exist_ok=True)
        (cases / f"{CASE_LIST_PREFIX}{key}.txt").write_text(
            f"cancer_study_identifier: {project_id}\n"
            f"stable_id: {project_id}_{key}\n"
            "case_list_category: other\n"
            f"case_list_name: {name}\n"
            f"case_list_description: {description} ({len(sample_ids)} sample(s))\n"
            f"case_list_ids: {chr(9).join(sample_ids)}\n", encoding="utf-8")
