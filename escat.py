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

"""ESCAT tiers from a table the lab curates ([ESCAT] FILE).

Varan ships no ESCAT content - only this mechanism and the format, in
escat_template.tsv. The table (TSV, "#" lines are comments):
  GENE        HUGO symbol
  ALTERATION  the syntax of study_alterations.py (V600E, G12, exon 19
              deletion, fusion, amplification, ...; several with ",")
  TUMOR       OncoTree codes (",") - a row matches samples coded with that
              code or a descendant (NSCLC matches LUAD, see
              oncotree_parents.tsv) - or * for any tumour
  LEVEL       I-A, I-B, I-C, II-A, II-B, III-A, III-B, IV, V or X
  NOTE, SOURCE  optional, shown in the report
"# version: ...", "# date: ..." and "# source: ..." comment lines are shown
in the report. The table is checked before the run starts: an unknown
level, alteration or OncoTree code stops it.

Output: MAF namespace ESCAT (ESCAT.LEVEL, best level of the variant for the
sample's tumour type); data_clinical_sample.txt ESCAT_HIGHEST_LEVEL and
ESCAT_ALTERATIONS (mutations, fusions, amplifications/deletions).

Reference: Mateo et al., Ann Oncol 2018, doi:10.1093/annonc/mdy263.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from loguru import logger

import study_alterations as sa
from clinical_file import drop_sample_columns, set_sample_column
from config_loader import VARAN_DIR, config

LEVELS = ["I-A", "I-B", "I-C", "II-A", "II-B", "III-A", "III-B", "IV", "V", "X"]
ONCOTREE = VARAN_DIR / "oncotree_parents.tsv"
SUMMARY_FILE = Path("service") / "escat.json"
NAMESPACE = "ESCAT"
MAF_COLUMN = "ESCAT.LEVEL"
COLUMNS = {
    "ESCAT_HIGHEST_LEVEL": ("ESCAT Highest Level",
                            "Best ESCAT tier of the sample's alterations for its tumour type (lab table)",
                            "STRING"),
    "ESCAT_ALTERATIONS": ("ESCAT Alterations", "The sample's alterations with an ESCAT tier", "STRING"),
}


def table_path() -> Path | None:
    value = config.get("ESCAT", "FILE", fallback="").strip()
    if not value:
        return None
    path = Path(value)
    return path if path.is_absolute() or path.exists() else VARAN_DIR / path


def parents() -> dict[str, str]:
    table = pd.read_csv(ONCOTREE, sep="\t", comment="#", dtype=str, keep_default_na=False)
    return dict(zip(table["code"].str.upper(), table["parent"].str.upper()))


def lineage(code: str, tree: dict[str, str]) -> set[str]:
    """`code` and its ancestors."""
    out, code = set(), code.strip().upper()
    while code and code not in out:
        out.add(code)
        code = tree.get(code, "")
    return out


def _level(value: str) -> str:
    text = value.strip().upper().replace(" ", "")
    if "-" not in text and len(text) > 1 and text[-1] in "ABC" and text[:-1] in ("I", "II", "III"):
        text = f"{text[:-1]}-{text[-1]}"
    return text


def load(path: Path | None = None) -> tuple[list[dict], dict[str, str]]:
    """(rules, header info) from the table; ValueError for a wrong row."""
    path = path or table_path()
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    info = {}
    for line in lines:
        if line.startswith("#") and ":" in line:
            key, value = line[1:].split(":", 1)
            if key.strip().lower() in ("version", "date", "source"):
                info[key.strip().lower()] = value.strip()
    table = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    table.columns = [c.strip().upper() for c in table.columns]
    missing = {"GENE", "ALTERATION", "TUMOR", "LEVEL"} - set(table.columns)
    if missing:
        msg = f"{Path(path).name}: missing column(s) {', '.join(sorted(missing))} (see escat_template.tsv)."
        raise ValueError(msg)
    tree = parents()
    rules = []
    for n, row in enumerate(table.to_dict("records"), start=2):
        where = f"{Path(path).name} row {n}"
        level = _level(row["LEVEL"])
        if level not in LEVELS:
            msg = f"{where}: LEVEL '{row['LEVEL']}' - use one of {', '.join(LEVELS)}."
            raise ValueError(msg)
        tumours = {t.strip().upper() for t in row["TUMOR"].split(",") if t.strip()}
        unknown = sorted(t for t in tumours if t != "*" and t not in tree)
        if unknown or not tumours:
            msg = f"{where}: TUMOR {', '.join(unknown) or '(empty)'} isn't an OncoTree code (or *)."
            raise ValueError(msg)
        try:
            specs = sa.parse(row["ALTERATION"])
        except ValueError as err:
            raise ValueError(f"{where}: {err}") from err
        rules.append({"gene": row["GENE"].strip(), "specs": specs, "alteration": row["ALTERATION"],
                      "tumours": tumours, "level": level,
                      "note": row.get("NOTE", ""), "source": row.get("SOURCE", "")})
    return rules, info


def check_table() -> None:
    """Before the run: stop on a wrong table (ValueError)."""
    path = table_path()
    if path is None:
        return
    if not path.exists():
        msg = f"[ESCAT] FILE not found: {path}"
        raise ValueError(msg)
    rules, info = load(path)
    logger.info(f"ESCAT table {path.name}: {len(rules)} row(s) {info.get('version', '')}".rstrip())


def _best(levels: list[str]) -> str:
    return next((lv for lv in LEVELS if lv in levels), "")


def annotate_study(folder: str | Path) -> dict | None:
    folder = Path(folder)
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    maf_path = folder / sa.MAF_FILE
    maf = sa.read_maf(folder)
    if maf is not None and MAF_COLUMN in maf.columns:
        maf = maf.drop(columns=[MAF_COLUMN])
        maf.to_csv(maf_path, sep="\t", index=False)
    path = table_path()
    from actionability import _samples
    samples = _samples(folder)
    if path is None or samples is None:
        drop_sample_columns(folder, list(COLUMNS))
        return None
    rules, info = load(path)
    tree = parents()
    tumour = {s: lineage(str(c), tree) for s, c in zip(
        samples["SAMPLE_ID"], samples.get("ONCOTREE_CODE", pd.Series("", index=samples.index)).fillna(""))}

    def matching(gene: str, sample: str, test) -> list[dict]:
        codes = tumour.get(sample, set())
        return [r for r in rules if r["gene"] == gene
                and ("*" in r["tumours"] or r["tumours"] & codes) and any(test(s) for s in r["specs"])]

    found: dict[str, list[dict]] = {}

    def add(sample: str, label: str, hits: list[dict]) -> str:
        level = _best([h["level"] for h in hits])
        if level:
            found.setdefault(sample, []).append({
                "alteration": label, "level": level,
                "rows": [f"{h['gene']} {h['alteration']} ({', '.join(sorted(h['tumours']))}): {h['level']}"
                         + (f" - {h['note']}" if h["note"] else "") for h in hits]})
        return level

    if maf is not None and not maf.empty:
        levels = []
        for r in maf.to_dict("records"):
            sample = str(r["Tumor_Sample_Barcode"]).strip()
            hits = matching(str(r.get("Hugo_Symbol", "")), sample, lambda s, r=r: s.matches_mutation(r))
            levels.append(add(sample, sa.mutation_label(r), hits))
        maf[MAF_COLUMN] = levels
        maf.to_csv(maf_path, sep="\t", index=False)
    for r in sa.fusions(sa.read_sv(folder)).to_dict("records"):
        sample = str(r["Sample_Id"]).strip()
        hits = [h for g in {r.get("Site1_Hugo_Symbol", ""), r.get("Site2_Hugo_Symbol", "")}
                for h in matching(g, sample, lambda s: s.kind == sa.FUSION)]
        add(sample, sa.fusion_label(r), hits)
    for call in sa.read_cna(folder):
        kind = sa.AMPLIFICATION if call["value"] > 0 else sa.DELETION
        add(call["sample"], sa.cna_label(call),
            matching(call["gene"], call["sample"], lambda s, kind=kind: s.kind == kind))

    attrs: dict[str, dict[str, object]] = {c: {} for c in COLUMNS}
    for sample in samples["SAMPLE_ID"]:
        hits = found.get(sample, [])
        attrs["ESCAT_HIGHEST_LEVEL"][sample] = _best([h["level"] for h in hits]) or "None"
        attrs["ESCAT_ALTERATIONS"][sample] = ", ".join(
            sorted({f"{h['alteration']} ({h['level']})" for h in hits})) or "None"
    for column, (short, long, datatype) in COLUMNS.items():
        set_sample_column(folder, column, attrs[column], short, long, datatype)

    summary = {"file": path.name, "info": info, "rows": len(rules), "samples": found}
    (folder / SUMMARY_FILE).parent.mkdir(parents=True, exist_ok=True)
    (folder / SUMMARY_FILE).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info(f"ESCAT: {len(found)} sample(s) with a tiered alteration ({path.name}).")
    return summary


def read_summary(folder: str | Path) -> dict | None:
    path = Path(folder) / SUMMARY_FILE
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
