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

"""Tumour-agnostic biomarkers per sample.

tumor_agnostic.tsv lists the biomarkers of tumour-agnostic approvals
("agnostic": NTRK and RET fusions, BRAF V600E, MSI-H, TMB-H) and biomarkers
relevant across tumour types without one ("pan-tumour": ALK, ROS1, FGFR,
ERBB2, HRD), with the FDA and EMA approvals and the date they were checked.
The two groups are never mixed.

Each biomarker is, for each sample:
- Detected (with the alteration found);
- Not detected - not the same as negative: coverage and detection limits
  aren't checked;
- Not tested - the sample has no data of that kind (no fusion/CNA data,
  no MSI/TMB/GIS value) or the gene isn't in its panel (when known).

data_clinical_sample.txt gets TUMOR_AGNOSTIC_BIOMARKERS (detected,
"agnostic" group), TUMOR_AGNOSTIC_NOT_TESTED and PAN_TUMOR_BIOMARKERS; the
report has the full table. [TumorAgnostic] TMB_HIGH / GIS_HIGH override the
table's thresholds. Research use only.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from loguru import logger

import gene_panels
import specimen
import study_alterations as sa
from clinical_file import drop_sample_columns, set_sample_column
from config_loader import VARAN_DIR, config

TABLE = VARAN_DIR / "tumor_agnostic.tsv"
SUMMARY_FILE = Path("service") / "tumor_agnostic.json"
DETECTED, NOT_DETECTED, NOT_TESTED = "Detected", "Not detected", "Not tested"
COLUMNS = {
    "TUMOR_AGNOSTIC_BIOMARKERS": ("Tumor-Agnostic Biomarkers",
                                  "Biomarkers of tumour-agnostic approvals detected (see report)", "STRING"),
    "TUMOR_AGNOSTIC_NOT_TESTED": ("Tumor-Agnostic Not Tested",
                                  "Tumour-agnostic biomarkers this sample couldn't be tested for", "STRING"),
    "PAN_TUMOR_BIOMARKERS": ("Pan-Tumor Biomarkers",
                             "Biomarkers relevant across tumour types (no tumour-agnostic approval) detected",
                             "STRING"),
}
DISCLAIMER = ("\"Not detected\" isn't a negative result: coverage and detection limits aren't "
              "checked. Approvals change: check the current labels. Research use only.")
# data_gene_panel_matrix.txt column and profiled-samples key per data kind.
DATA = {sa.MUTATION: ("mutations", "sequenced"), sa.FUSION: ("structural_variants", "sv"),
        sa.AMPLIFICATION: ("cna", "cna")}


def load_table(path: Path = TABLE) -> list[dict]:
    rows = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False).to_dict("records")
    for r in rows:
        r["genes"] = [g.strip() for g in r["genes"].split(",") if g.strip()]
        r["specs"] = sa.parse(r["alteration"]) if r["alteration"] else []
        override = {"tmb": "TMB_HIGH", "gis": "GIS_HIGH"}.get(r["kind"])
        if override:
            r["threshold"] = config.get("TumorAgnostic", override, fallback="").strip() or r["threshold"]
    return rows


def _number(value: object) -> float | None:
    try:
        number = float(str(value).strip())
    except ValueError:
        return None
    return None if number != number else number


def annotate_study(folder: str | Path, profiled: dict[str, list[str]] | None = None) -> dict | None:
    """Write the sample attributes and service/tumor_agnostic.json (for the
    report). `profiled`: the samples with each kind of data (read from the
    case lists / service/ by Make_meta_and_cases)."""
    folder = Path(folder)
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    samples = specimen.read_study_samples(folder)
    if samples is None or not config.get("TumorAgnostic", "ENABLED", fallback="True").strip().lower() \
            in ("true", "t", "yes", "1"):
        drop_sample_columns(folder, list(COLUMNS))
        return None
    clinical = samples.fillna("").astype(str).set_index("SAMPLE_ID")
    types = specimen.study_types(folder)
    maf = sa.read_maf(folder)
    fusions = sa.fusions(sa.read_sv(folder))
    cna = sa.read_cna(folder)
    rules = load_table()
    profiled = profiled or {}

    def tested(sample: str, kind: str, gene: str) -> str:
        """"" if `sample` was tested for `gene`, else why not."""
        column, key = DATA[kind]
        if profiled.get(key) is not None and sample not in set(profiled[key]):
            return {"sequenced": "no mutation data", "sv": "no fusion data", "cna": "no CNA data"}[key]
        if kind == sa.MUTATION and maf is None:
            return "no mutation data"
        if kind == sa.FUSION and not (folder / sa.SV_FILE).exists():
            return "no fusion data"
        if kind == sa.AMPLIFICATION and not (folder / sa.CNA_FILE).exists():
            return "no CNA data"
        genes = gene_panels.sample_genes(folder, sample, types.get(sample), column)
        return f"{gene} not in panel" if genes is not None and gene not in genes else ""

    results: dict[str, dict[str, dict]] = {}
    for sample in clinical.index:
        row = clinical.loc[sample]
        results[sample] = {}
        for rule in rules:
            kind = rule["kind"]
            if kind in ("msi", "tmb", "gis"):
                status, detail = _clinical_status(rule, row)
            else:
                status, detail = _gene_status(rule, sample, maf, fusions, cna, tested)
            results[sample][rule["biomarker"]] = {"status": status, "detail": detail}

    attrs: dict[str, dict[str, object]] = {c: {} for c in COLUMNS}
    for sample, found in results.items():
        groups = {g: [r["biomarker"] for r in rules if r["group"] == g] for g in ("agnostic", "pan-tumour")}
        hits = [b for b in groups["agnostic"] if found[b]["status"] == DETECTED]
        untested = [b for b in groups["agnostic"] if found[b]["status"] == NOT_TESTED]
        pan = [b for b in groups["pan-tumour"] if found[b]["status"] == DETECTED]
        attrs["TUMOR_AGNOSTIC_BIOMARKERS"][sample] = "; ".join(hits) or "None detected"
        attrs["TUMOR_AGNOSTIC_NOT_TESTED"][sample] = "; ".join(untested) or "None"
        attrs["PAN_TUMOR_BIOMARKERS"][sample] = "; ".join(pan) or "None detected"
    for column, (short, long, datatype) in COLUMNS.items():
        set_sample_column(folder, column, attrs[column], short, long, datatype)

    summary = {
        "rules": [{k: r[k] for k in ("biomarker", "group", "kind", "threshold", "fda", "ema", "note")}
                  for r in rules],
        "checked": _checked_date(),
        "samples": results,
    }
    (folder / SUMMARY_FILE).parent.mkdir(parents=True, exist_ok=True)
    (folder / SUMMARY_FILE).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    n = sum(1 for s in attrs["TUMOR_AGNOSTIC_BIOMARKERS"].values() if s != "None detected")
    logger.info(f"Tumour-agnostic biomarkers detected in {n}/{len(results)} sample(s).")
    return summary


def _clinical_status(rule: dict, row: pd.Series) -> tuple[str, str]:
    if rule["kind"] == "msi":
        value = str(row.get("MSI_THR", "")).strip().lower()
        if value == "unstable":
            return DETECTED, f"MSI {row.get('MSI', '')}".strip()
        if value == "stable":
            return NOT_DETECTED, ""
        return NOT_TESTED, "no MSI status"
    column = "TMB" if rule["kind"] == "tmb" else "GENOMIC_INSTABILITY_SCORE"
    value = _number(row.get(column, ""))
    if value is None:
        return NOT_TESTED, f"no {column} value"
    threshold = float(rule["threshold"])
    if value >= threshold:
        return DETECTED, f"{column} {value:g} (>= {threshold:g})"
    return NOT_DETECTED, f"{column} {value:g}"


def _gene_status(rule: dict, sample: str, maf: pd.DataFrame | None, fusions: pd.DataFrame,
                 cna: list[dict], tested) -> tuple[str, str]:
    kind = rule["kind"]
    found = []
    if kind == sa.MUTATION and maf is not None:
        rows = maf[maf["Tumor_Sample_Barcode"].astype(str).str.strip() == sample]
        for r in rows.to_dict("records"):
            if r.get("Hugo_Symbol") in rule["genes"] and any(s.matches_mutation(r) for s in rule["specs"]):
                found.append(sa.mutation_label(r))
    elif kind == sa.FUSION:
        rows = fusions[fusions["Sample_Id"].astype(str).str.strip() == sample]
        for r in rows.to_dict("records"):
            if any(sa.fusion_matches(g, r) for g in rule["genes"]):
                found.append(sa.fusion_label(r))
    elif kind == sa.AMPLIFICATION:
        found = [sa.cna_label(c) for c in cna
                 if c["sample"] == sample and c["gene"] in rule["genes"] and c["value"] == 2]
    if found:
        return DETECTED, ", ".join(sorted(set(found)))
    reasons = [tested(sample, kind, g) for g in rule["genes"]]
    if all(reasons):
        return NOT_TESTED, "; ".join(sorted(set(reasons)))
    missing = [r for r in reasons if r]
    return NOT_DETECTED, ("partly tested: " + "; ".join(missing)) if missing else ""


def _checked_date() -> str:
    for line in TABLE.read_text(encoding="utf-8").splitlines():
        if "Last checked:" in line:
            return line.split("Last checked:", 1)[1].split(".")[0].strip()
    return ""


def read_summary(folder: str | Path) -> dict | None:
    path = Path(folder) / SUMMARY_FILE
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
