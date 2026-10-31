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

"""Serial ctDNA: compare the plasma samples of a patient over time.

Input: sample.tsv's optional COLLECTION_DATE column (YYYY-MM-DD or
DD/MM/YYYY) and optional BASELINE column (yes/true/1 marks a patient's
baseline plasma sample - by default the earliest one). The dates are kept
in the study's service/collection_dates.tsv, never in the files cBioPortal
imports, so that a later batch merged with -u is ordered correctly; what
cBioPortal gets is days from the patient's first sample.

Output, only when at least one sample has a COLLECTION_DATE:
- data_timeline_specimen.txt: one SPECIMEN event per dated sample, so
  cBioPortal's patient view ("Genomic Evolution") orders the samples in time.
- data_clinical_sample.txt: DAYS_FROM_FIRST_SAMPLE; for plasma samples
  CTDNA_MAX_VAF_PCT, CTDNA_MEAN_VAF_PCT, CTDNA_N_VARIANTS; for patients with
  two or more dated plasma samples CTDNA_CHANGE_PCT (vs the baseline),
  MOLECULAR_RESPONSE, CTDNA_N_NEW, CTDNA_N_LOST and RESISTANCE_VARIANTS.
- data_mutations_extended.txt ("SERIAL" namespace): SERIAL.STATUS (Baseline /
  Persistent / New / New (not covered before)) and SERIAL.RESISTANCE.

Metrics leave out variants flagged as clonal hematopoiesis (see
clonal_hematopoiesis.py) and germline-like VAFs. MOLECULAR_RESPONSE follows
the ctMoniTR aggregate analysis of TKI-treated NSCLC (Clin Cancer Res 2025;
31(11):2162): Decrease = the metric drops by more than [Serial]
DECREASE_PCT (50%), Increase = rises by more than INCREASE_PCT (20%),
Intermediate otherwise. There is no consensus definition of molecular
response; the thresholds are configurable. "Not detected" (no variant left)
is kept apart from Decrease. Research use only.

Treatments (treatments.tsv, see clinical_timeline.py), when the patient
has them and no BASELINE is marked: each plasma sample is compared with the
last plasma sample collected up to the start of the treatment it is on
("Not evaluable (no pre-treatment sample)" if there is none); samples off
treatment are "Baseline" when a treatment's response is measured against
them, "No treatment ongoing" otherwise.

Resistance: plasma variants of serial patients matching
resistance_mutations.tsv - acquired-resistance mutations with their therapy
class and reference - are annotated. When the patient's treatments are
known, a match is kept only if the patient received a treatment of that
class before the sample (the others are listed in the report only);
without treatments the annotation just says which therapy it concerns.
"""

from __future__ import annotations

import ast
import json
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path

import pandas as pd
from loguru import logger

import clinical_timeline
import clonal_hematopoiesis as ch
import gene_panels
import specimen
from clinical_file import drop_sample_columns, set_sample_column
from config_loader import VARAN_DIR, config

DATE_COLUMN = "COLLECTION_DATE"
BASELINE_COLUMN = "BASELINE"
DATES_FILE = Path("service") / "collection_dates.tsv"
SUMMARY_FILE = Path("service") / "serial_ctdna.json"
TIMELINE_FILE = "data_timeline_specimen.txt"
TIMELINE_META = "meta_timeline_specimen.txt"
NAMESPACE = "SERIAL"
STATUS = "SERIAL.STATUS"
RESISTANCE = "SERIAL.RESISTANCE"

# Defaults, used when the [Serial] key is missing or blank.
DEFAULTS = {
    "ENABLED": "True",
    # ctMoniTR used the maximum VAF; "mean" is the other common choice.
    "METRIC": "max",
    # ctMoniTR TKI analysis (Clin Cancer Res 2025;31(11):2162): decrease
    # > 50%, increase > 20%, intermediate otherwise.
    "DECREASE_PCT": "50",
    "INCREASE_PCT": "20",
    # Leave CH confirmed/suspected variants out of the metrics.
    "EXCLUDE_CH": "True",
    # Germline-like VAF bands (heterozygous ~0.5, homozygous ~1.0), as
    # [SigMA] VAF_EXCLUDE_BANDS.
    "GERMLINE_VAF_BANDS": "[[0.45, 0.55], [0.90, 1.0]]",
    "RESISTANCE_LIST": "resistance_mutations.tsv",
}

SAMPLE_COLUMNS = {
    "DAYS_FROM_FIRST_SAMPLE": ("Days From First Sample",
                               "Days from the patient's first sample collection", "NUMBER"),
    "CTDNA_MAX_VAF_PCT": ("ctDNA Max VAF (%)",
                          "Highest VAF of the plasma sample's variants (CH and germline-like excluded)",
                          "NUMBER"),
    "CTDNA_MEAN_VAF_PCT": ("ctDNA Mean VAF (%)",
                           "Mean VAF of the plasma sample's variants (CH and germline-like excluded)",
                           "NUMBER"),
    "CTDNA_N_VARIANTS": ("ctDNA Variants",
                         "Variants counted in the ctDNA metrics", "NUMBER"),
    "CTDNA_CHANGE_PCT": ("ctDNA Change From Baseline (%)",
                         "Change of the ctDNA metric from the patient's baseline plasma sample",
                         "NUMBER"),
    "MOLECULAR_RESPONSE": ("Molecular Response",
                           "Baseline / Decrease / Intermediate / Increase / Not detected / "
                           "Not evaluable / No treatment ongoing (research use)", "STRING"),
    "CTDNA_N_NEW": ("ctDNA New Variants",
                    "Variants not detected in any earlier plasma sample of the patient", "NUMBER"),
    "CTDNA_N_LOST": ("ctDNA Lost Variants",
                     "Variants of the previous plasma sample no longer detected (gene covered)",
                     "NUMBER"),
    "RESISTANCE_VARIANTS": ("Resistance Variants",
                            "Known acquired-resistance mutations detected (see report)", "STRING"),
    "RESPONSE_BASELINE": ("Molecular Response Baseline",
                          "Plasma sample the molecular response is measured against", "STRING"),
    "TREATMENT_AT_SAMPLE": ("Treatment At Sample",
                            "Treatment the patient was on when the sample was collected", "STRING"),
    "DAYS_ON_TREATMENT": ("Days On Treatment",
                          "Days from the start of that treatment to the sample", "NUMBER"),
}

DISCLAIMER = (
    "Molecular response from serial ctDNA is a research tool: there is no "
    "consensus definition. Research use only.")


# ------------------------------------------------------------ settings

def setting(key: str) -> str:
    value = config.get("Serial", key, fallback="").strip()
    return value or DEFAULTS[key]


def enabled() -> bool:
    return setting("ENABLED").lower() in ("true", "t", "yes", "1")


# --------------------------------------------------------------- dates

def parse_date(value: object, where: str = "") -> date | None:
    """YYYY-MM-DD or DD/MM/YYYY; None for a blank value."""
    text = "" if value is None else str(value).strip()
    if text == "" or text.upper() in ("NA", "NAN"):
        return None
    for fmt in ("%Y-%m-%d", "%d/%m/%Y"):
        try:
            return datetime.strptime(text, fmt).date()
        except ValueError:
            pass
    msg = f"{where + ': ' if where else ''}'{text}' is not a date - use YYYY-MM-DD or DD/MM/YYYY."
    raise ValueError(msg)


def _is_yes(value: object) -> bool:
    return str(value).strip().lower() in ("yes", "y", "true", "t", "1", "x", "baseline")


def _column(df: pd.DataFrame, name: str) -> str | None:
    return next((c for c in df.columns if str(c).strip().upper() == name), None)


def collection_dates(sample_table: pd.DataFrame | None) -> pd.DataFrame:
    """SAMPLE_ID, PATIENT_ID, COLLECTION_DATE (ISO), BASELINE from a
    sample.tsv - empty when it has no COLLECTION_DATE column. Raises
    ValueError for a malformed date."""
    empty = pd.DataFrame(columns=["SAMPLE_ID", "PATIENT_ID", DATE_COLUMN, BASELINE_COLUMN])
    if sample_table is None:
        return empty
    date_col = _column(sample_table, DATE_COLUMN)
    if date_col is None:
        return empty
    base_col = _column(sample_table, BASELINE_COLUMN)
    rows = []
    for row in sample_table.fillna("").to_dict("records"):
        sample_id = str(row["SAMPLE_ID"]).strip()
        day = parse_date(row.get(date_col), f"{DATE_COLUMN} of sample {sample_id}")
        if day is None:
            continue
        rows.append({
            "SAMPLE_ID": sample_id,
            "PATIENT_ID": str(row.get("PATIENT_ID", "")).strip(),
            DATE_COLUMN: day.isoformat(),
            BASELINE_COLUMN: "yes" if base_col and _is_yes(row.get(base_col)) else "",
        })
    return pd.DataFrame(rows, columns=empty.columns)


def check_input(sample_table: pd.DataFrame, types: dict[str, str]) -> None:
    """Stop early on a malformed date, a BASELINE given to a non-plasma or
    undated sample, or two baselines for one patient."""
    dates = collection_dates(sample_table)
    base_col = _column(sample_table, BASELINE_COLUMN)
    if base_col is None:
        return
    dated = set(dates["SAMPLE_ID"])
    flagged = [str(r["SAMPLE_ID"]).strip() for r in sample_table.fillna("").to_dict("records")
               if _is_yes(r.get(base_col))]
    for sample_id in flagged:
        if types.get(sample_id) != specimen.PLASMA or sample_id not in dated:
            msg = (f"{BASELINE_COLUMN} marks sample {sample_id}, but a baseline must be a "
                   f"Plasma sample with a {DATE_COLUMN}.")
            raise ValueError(msg)
    per_patient = dates[dates[BASELINE_COLUMN] == "yes"].groupby("PATIENT_ID")["SAMPLE_ID"].apply(list)
    for patient, samples in per_patient.items():
        if len(samples) > 1:
            msg = f"Patient {patient} has more than one {BASELINE_COLUMN} sample: {', '.join(samples)}."
            raise ValueError(msg)


def write_collection_dates(folder: str | Path, sample_table: pd.DataFrame | None) -> None:
    """At creation: keep sample.tsv's collection dates in service/."""
    dates = collection_dates(sample_table)
    path = Path(folder) / DATES_FILE
    if dates.empty:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    dates.to_csv(path, sep="\t", index=False)


def carry_collection_dates(sources: list[str | Path], output: str | Path) -> None:
    """update/extract/remove: the output study gets its sources' dates
    (the later source wins for a sample in both); annotate_study() then
    keeps only its own samples."""
    frames = []
    for source in sources:
        path = Path(source) / DATES_FILE
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False))
    if not frames:
        return
    merged = pd.concat(frames).drop_duplicates("SAMPLE_ID", keep="last")
    out = Path(output) / DATES_FILE
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep="\t", index=False)


def read_collection_dates(folder: str | Path) -> pd.DataFrame | None:
    path = Path(folder) / DATES_FILE
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# ---------------------------------------------------------- resistance

@dataclass
class ResistanceRule:
    rule: ch.GeneRule
    label: str
    therapy: str
    drugs: str
    reference: str
    only_new: bool
    doi: str = ""
    note: str = ""

    @property
    def classes(self) -> list[str]:
        return clinical_timeline.split_classes(self.therapy)

    def text(self) -> str:
        return (f"{self.label}: resistance to {' or '.join(self.classes)} ({self.drugs})"
                f" - {self.reference}")


def resistance_list_path() -> Path:
    path = Path(setting("RESISTANCE_LIST"))
    return path if path.is_absolute() or path.exists() else VARAN_DIR / path


def load_resistance_list(path: Path | None = None) -> list[ResistanceRule]:
    path = path or resistance_list_path()
    table = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    known = {c.lower() for c in clinical_timeline.therapy_classes()}
    rules = []
    for row in table.to_dict("records"):
        unknown = [c for c in clinical_timeline.split_classes(row["therapy"]) if c.lower() not in known]
        if unknown:
            msg = (f"{path.name}, {row['label']}: therapy {', '.join(unknown)} isn't in "
                   f"{clinical_timeline.THERAPY_CLASSES.name}.")
            raise ValueError(msg)
        classes: set[str] = set()
        for name in filter(None, (c.strip().lower() for c in row["classes"].split(","))):
            classes |= ch.CLASSES[name]
        rule = ch.GeneRule(
            gene=row["gene"].strip(), tier=0, classes=classes,
            hotspots=[h.strip() for h in row["hotspots"].split(",") if h.strip()],
            exons={e.strip() for e in row["exons"].split(",") if e.strip()})
        rules.append(ResistanceRule(
            rule=rule, label=row["label"].strip() or row["gene"].strip(),
            therapy=row["therapy"].strip(), drugs=row["drugs"].strip(),
            reference=row["reference"].strip(), only_new=_is_yes(row["only_new"]),
            doi=row.get("doi", "").strip(), note=row.get("note", "").strip()))
    return rules


# ------------------------------------------------------------ analysis

def _bands() -> list[tuple[float, float]]:
    return [tuple(b) for b in ast.literal_eval(setting("GERMLINE_VAF_BANDS"))]


def _counted(maf: pd.DataFrame, vafs: pd.Series) -> pd.Series:
    """Rows that enter the ctDNA metrics: VAF known, not CH, not germline-like."""
    keep = vafs.notna()
    for low, high in _bands():
        keep &= ~((vafs >= low) & (vafs <= high))
    if setting("EXCLUDE_CH").lower() in ("true", "t", "yes", "1") and ch.EVIDENCE in maf.columns:
        keep &= ~maf[ch.EVIDENCE].fillna("").isin(ch.FLAGGED)
    return keep


def _category(change: float | None, baseline_value: float | None, value: float | None) -> str:
    if not baseline_value:
        return "Not evaluable"
    if value is None:
        return "Not detected"
    if change <= -float(setting("DECREASE_PCT")):
        return "Decrease"
    if change > float(setting("INCREASE_PCT")):
        return "Increase"
    return "Intermediate"


def annotate_study(folder: str | Path) -> dict | None:
    """Write the timeline, sample attributes and MAF columns described in
    the module docstring. Returns a summary for the report, or None when no
    sample has a collection date (everything serial is then removed)."""
    folder = Path(folder)
    _keys_cache.clear()
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    samples = specimen.read_study_samples(folder)
    dates = read_collection_dates(folder)
    present = set(samples["SAMPLE_ID"].astype(str).str.strip()) if samples is not None else set()
    if dates is not None:
        dates = dates[dates["SAMPLE_ID"].isin(present)]
        dates.to_csv(folder / DATES_FILE, sep="\t", index=False)

    maf_path = folder / ch.MAF_FILE
    maf = None
    if maf_path.exists():
        from filter_clinvar import read_maf
        maf = read_maf(maf_path)
        stale = [c for c in (STATUS, RESISTANCE) if c in maf.columns]
        if stale:
            maf = maf.drop(columns=stale)
            maf.to_csv(maf_path, sep="\t", index=False)

    types = specimen.study_types(folder)
    patient_of = {} if samples is None else dict(zip(
        samples["SAMPLE_ID"].astype(str).str.strip(), samples["PATIENT_ID"].astype(str).str.strip()))
    when = {} if dates is None else {
        r["SAMPLE_ID"]: date.fromisoformat(r[DATE_COLUMN]) for r in dates.to_dict("records")}
    # Treatments and status: their own timelines, and day 0 for all of them.
    timeline = clinical_timeline.annotate_study(folder, patient_of, when)

    if not enabled() or not when:
        drop_sample_columns(folder, list(SAMPLE_COLUMNS))
        (folder / TIMELINE_FILE).unlink(missing_ok=True)
        return None

    forced = {r["SAMPLE_ID"] for r in dates.to_dict("records") if r[BASELINE_COLUMN] == "yes"}

    by_patient: dict[str, list[str]] = {}
    for sample_id in sorted(when, key=lambda s: (when[s], s)):
        by_patient.setdefault(patient_of.get(sample_id, sample_id), []).append(sample_id)
    first_day = {p: when[s[0]] for p, s in by_patient.items()}
    days = {s: (when[s] - first_day[patient_of.get(s, s)]).days for s in when}

    _write_timeline(folder, by_patient, {
        s: timeline.day(patient_of.get(s, s), when[s]) for s in when}, types)

    # Per-sample ctDNA metrics (plasma samples with a date).
    metric_name = "mean" if setting("METRIC").lower() == "mean" else "max"
    values: dict[str, dict[str, float | None]] = {}
    detected: dict[str, set] = {}
    if maf is not None and not maf.empty:
        vafs = ch.vaf(maf)
        counted = _counted(maf, vafs)
        barcode = maf["Tumor_Sample_Barcode"].astype(str).str.strip()
        for sample_id in when:
            if types.get(sample_id) != specimen.PLASMA:
                continue
            rows = counted & (barcode == sample_id)
            v = vafs[rows]
            values[sample_id] = {
                "max": float(v.max()) if len(v) else None,
                "mean": float(v.mean()) if len(v) else None,
                "n": int(len(v)),
            }
            detected[sample_id] = {ch.maf_key(r) for r in maf[rows].to_dict("records")}
    else:
        for sample_id in when:
            if types.get(sample_id) == specimen.PLASMA:
                values[sample_id] = {"max": None, "mean": None, "n": 0}
                detected[sample_id] = set()

    resistance = load_resistance_list()
    attrs: dict[str, dict[str, object]] = {c: {} for c in SAMPLE_COLUMNS}
    status = pd.Series("", index=maf.index if maf is not None else [])
    resist = pd.Series("", index=maf.index if maf is not None else [])
    summary = {"metric": metric_name, "patients": {}, "resistance": [],
               "timeline_reference": timeline.reference,
               "treatments": bool(timeline.has_treatments())}
    coverage: dict[str, set[str] | None] = {}

    def covers(sample_id: str, gene: str) -> bool | None:
        if sample_id not in coverage:
            coverage[sample_id] = gene_panels.sample_genes(folder, sample_id, types.get(sample_id))
        genes = coverage[sample_id]
        return None if genes is None else gene in genes

    for sample_id, day in days.items():
        attrs["DAYS_FROM_FIRST_SAMPLE"][sample_id] = day
    for sample_id, v in values.items():
        attrs["CTDNA_MAX_VAF_PCT"][sample_id] = _pct(v["max"])
        attrs["CTDNA_MEAN_VAF_PCT"][sample_id] = _pct(v["mean"])
        attrs["CTDNA_N_VARIANTS"][sample_id] = v["n"]

    for patient, sample_ids in by_patient.items():
        plasma = [s for s in sample_ids if s in values]
        if len(plasma) < 2:
            continue
        treatments = timeline.treatments(patient)
        forced_base = next((s for s in plasma if s in forced), None)
        first = forced_base or plasma[0]
        rows_out = []
        for idx, sample_id in enumerate(plasma):
            value = values[sample_id][metric_name]
            current = clinical_timeline.in_force(treatments, when[sample_id])
            base, category = _reference(sample_id, plasma, when, forced_base, treatments, current)
            change = None
            if category is None:
                base_value = values[base][metric_name]
                change = ((value or 0.0) - base_value) / base_value * 100 if base_value else None
                category = _category(change, base_value, value)
            attrs["MOLECULAR_RESPONSE"][sample_id] = category
            attrs["CTDNA_CHANGE_PCT"][sample_id] = "NA" if change is None else round(change, 1)
            attrs["RESPONSE_BASELINE"][sample_id] = base or "NA"
            if treatments:
                attrs["TREATMENT_AT_SAMPLE"][sample_id] = current.label() if current else "None"
                attrs["DAYS_ON_TREATMENT"][sample_id] = (
                    (when[sample_id] - current.start).days if current else "NA")

            earlier = plasma[:idx]
            previous = earlier[-1] if earlier else None
            new = lost = 0
            if maf is not None:
                for i, r in maf[maf["Tumor_Sample_Barcode"].astype(str).str.strip() == sample_id].iterrows():
                    key, gene = ch.maf_key(r), str(r.get("Hugo_Symbol", "")).strip()
                    if sample_id == first or not earlier:
                        state = "Baseline" if sample_id == first else "First sample"
                    elif any(key in _all_keys(maf, e) for e in earlier):
                        state = "Persistent"
                    elif any(covers(e, gene) is not False for e in earlier):
                        state, new = "New", new + 1
                    else:
                        state = "New (not covered before)"
                    status[i] = state
                    hits = [res for res in resistance
                            if res.rule.gene == gene and res.rule.matches(r.to_dict())
                            and (not res.only_new or state.startswith("New"))]
                    kept = []
                    for res in hits:
                        given = clinical_timeline.received(treatments, res.classes, when[sample_id])
                        if not treatments:
                            received = "unknown (no treatments given)"
                        elif given:
                            received = (f"yes: {given.label()} from day "
                                        f"{timeline.day(patient, given.start)}")
                        else:
                            received = "no"
                        summary["resistance"].append({
                            "patient": patient, "sample": sample_id,
                            "day": timeline.day(patient, when[sample_id]),
                            "variant": f"{gene} {r.get('HGVSp_Short', '')}", "status": state,
                            "text": res.text(), "received": received,
                            "doi": res.doi, "note": res.note})
                        if received != "no":
                            kept.append(res.text() + (f" - received {received[5:]}"
                                                      if received.startswith("yes") else ""))
                    resist[i] = " | ".join(kept)
                if previous:
                    for key in _all_keys(maf, previous) - _all_keys(maf, sample_id):
                        gene = _gene_of(maf, previous, key)
                        if covers(sample_id, gene) is not False:
                            lost += 1
            attrs["CTDNA_N_NEW"][sample_id] = new
            attrs["CTDNA_N_LOST"][sample_id] = lost
            hits = sorted({h["variant"] for h in summary["resistance"]
                           if h["sample"] == sample_id and h["received"] != "no"})
            attrs["RESISTANCE_VARIANTS"][sample_id] = ", ".join(hits) or "None"
            rows_out.append({
                "sample": sample_id, "day": timeline.day(patient, when[sample_id]),
                "max": _pct(values[sample_id]["max"]), "mean": _pct(values[sample_id]["mean"]),
                "n": values[sample_id]["n"], "change": attrs["CTDNA_CHANGE_PCT"][sample_id],
                "category": category, "baseline": base or "NA",
                "treatment": current.label() if current else ("None" if treatments else "NA"),
                "new": new, "lost": lost})
        summary["patients"][patient] = rows_out

    for column, (short, long, datatype) in SAMPLE_COLUMNS.items():
        if attrs[column]:
            set_sample_column(folder, column, attrs[column], short, long, datatype)
        else:
            drop_sample_columns(folder, [column])

    if maf is not None and (status != "").any():
        maf[STATUS] = status
        maf[RESISTANCE] = resist
        maf.to_csv(maf_path, sep="\t", index=False)

    summary["settings"] = {k: setting(k) for k in DEFAULTS if k != "ENABLED"}
    (folder / SUMMARY_FILE).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info(
        f"Serial ctDNA: {len(summary['patients'])} patient(s) with 2+ dated plasma samples, "
        f"{len(summary['resistance'])} known resistance variant(s). {DISCLAIMER}")
    return summary


def _reference(sample_id: str, plasma: list[str], when: dict[str, date], forced_base: str | None,
               treatments: list, current) -> tuple[str | None, str | None]:
    """(baseline sample, fixed category) for a plasma sample; a None
    category means "compare with the baseline". See the module docstring."""
    if forced_base or not treatments:
        base = forced_base or plasma[0]
        if plasma.index(sample_id) < plasma.index(base):
            return base, "Before baseline"
        return base, "Baseline" if sample_id == base else None
    if current is None:
        bases = {_last_until(plasma, when, t.start) for t in treatments}
        return None, "Baseline" if sample_id in bases else "No treatment ongoing"
    base = _last_until(plasma, when, current.start)
    return base, None if base else "Not evaluable (no pre-treatment sample)"


def _last_until(plasma: list[str], when: dict[str, date], day: date) -> str | None:
    """The last plasma sample collected on or before `day`."""
    before = [s for s in plasma if when[s] <= day]
    return before[-1] if before else None


def _pct(value: float | None) -> object:
    return "NA" if value is None else round(value * 100, 2)


_keys_cache: dict[tuple[int, str], set] = {}


def _all_keys(maf: pd.DataFrame, sample_id: str) -> set:
    """Every variant of `sample_id` in the MAF (CH ones included: a
    variant's presence, not the metric)."""
    key = (id(maf), sample_id)
    if key not in _keys_cache:
        rows = maf[maf["Tumor_Sample_Barcode"].astype(str).str.strip() == sample_id]
        _keys_cache[key] = {ch.maf_key(r) for r in rows.to_dict("records")}
    return _keys_cache[key]


def _gene_of(maf: pd.DataFrame, sample_id: str, key: tuple) -> str:
    rows = maf[maf["Tumor_Sample_Barcode"].astype(str).str.strip() == sample_id]
    for r in rows.to_dict("records"):
        if ch.maf_key(r) == key:
            return str(r.get("Hugo_Symbol", "")).strip()
    return ""


def _write_timeline(folder: Path, by_patient: dict[str, list[str]], days: dict[str, int],
                    types: dict[str, str]) -> None:
    rows = []
    for patient, sample_ids in by_patient.items():
        for sample_id in sample_ids:
            rows.append({
                "PATIENT_ID": patient, "START_DATE": days[sample_id], "STOP_DATE": "",
                "EVENT_TYPE": "SPECIMEN", "SAMPLE_ID": sample_id,
                "SPECIMEN_TYPE": types.get(sample_id, "NA")})
    pd.DataFrame(rows).to_csv(folder / TIMELINE_FILE, sep="\t", index=False)


def read_summary(folder: str | Path) -> dict | None:
    """annotate_study()'s summary, for the report."""
    path = Path(folder) / SUMMARY_FILE
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def create_meta_timeline(project_id: str, folder: str | Path) -> None:
    if not (Path(folder) / TIMELINE_FILE).exists():
        return
    (Path(folder) / TIMELINE_META).write_text(
        f"cancer_study_identifier: {project_id}\n"
        "genetic_alteration_type: CLINICAL\n"
        "datatype: TIMELINE\n"
        f"data_filename: {TIMELINE_FILE}\n", encoding="utf-8")
    logger.info(f"{TIMELINE_META} created!")
