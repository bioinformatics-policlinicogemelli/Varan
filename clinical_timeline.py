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

"""Clinical timeline: treatments and clinical status of each patient.

Input, both optional:
- treatments.tsv, one row per event - found at conf.ini's [Timeline]
  TREATMENTS, else next to patient.tsv, else next to sample.tsv:
    PATIENT_ID, EVENT_TYPE (TREATMENT or STATUS), START_DATE, STOP_DATE
    (TREATMENT only; blank = ongoing), and
    for TREATMENT: AGENT (free text), THERAPY_CLASS (from
      therapy_classes.tsv, ";" between several), LINE (a number);
    for STATUS: STATUS (Progression / Response / Death) and RECIST (CR, PR,
      SD or PD, required for Response).
- patient.tsv's DIAGNOSIS_DATE column.
Dates: YYYY-MM-DD or DD/MM/YYYY. They are kept in the study's
service/clinical_events.tsv (like sample collection dates, see
serial_ctdna.py) and never reach cBioPortal, which gets days from day 0.

Day 0: the diagnosis when every patient with timeline events has a
DIAGNOSIS_DATE; otherwise, for every patient, the earliest of their sample
collection, treatment and status dates - one convention per study.

Output: data_timeline_treatment.txt and data_timeline_status.txt (with their
meta files), on the same day 0 as data_timeline_specimen.txt. serial_ctdna.py
uses the treatments for the molecular response baseline and to check that a
resistance mutation's therapy was received.

Stops the run (ValueError) on an unreadable date, a treatment ending before
it starts, a missing PATIENT_ID/EVENT_TYPE/START_DATE, or a value outside the
allowed ones (EVENT_TYPE, STATUS, RECIST, THERAPY_CLASS, LINE); anything
else that looks wrong is a warning.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

import pandas as pd
from loguru import logger

from clinical_file import PATIENT_FILE, drop_sample_columns
from config_loader import VARAN_DIR, config

TREATMENTS_FILE = "treatments.tsv"
DIAGNOSIS_COLUMN = "DIAGNOSIS_DATE"
EVENTS_FILE = Path("service") / "clinical_events.tsv"
SUMMARY_FILE = Path("service") / "clinical_timeline.json"
THERAPY_CLASSES = VARAN_DIR / "therapy_classes.tsv"
TREATMENT_TIMELINE = "data_timeline_treatment.txt"
STATUS_TIMELINE = "data_timeline_status.txt"
METAS = {TREATMENT_TIMELINE: "meta_timeline_treatment.txt",
         STATUS_TIMELINE: "meta_timeline_status.txt"}

TREATMENT = "TREATMENT"
STATUS = "STATUS"
DIAGNOSIS = "DIAGNOSIS"
STATUSES = ("Progression", "Response", "Death")
RECIST = ("CR", "PR", "SD", "PD")
COLUMNS = ["PATIENT_ID", "EVENT_TYPE", "START_DATE", "STOP_DATE", "AGENT",
           "THERAPY_CLASS", "LINE", "STATUS", "RECIST"]


# --------------------------------------------------------------- input

def treatments_path(sample_tsv: str | Path | None,
                    patient_tsv: str | Path | None) -> Path | None:
    """[Timeline] TREATMENTS, else treatments.tsv next to patient.tsv, else
    next to sample.tsv; None when there is none."""
    configured = config.get("Timeline", "TREATMENTS", fallback="").strip()
    if configured:
        path = Path(configured)
        if not path.exists():
            msg = f"[Timeline] TREATMENTS points to a missing file: {configured}"
            raise ValueError(msg)
        return path
    for neighbour in (patient_tsv, sample_tsv):
        if neighbour and str(neighbour).strip():
            path = Path(neighbour).parent / TREATMENTS_FILE
            if path.exists():
                return path
    return None


def therapy_classes() -> list[str]:
    table = pd.read_csv(THERAPY_CLASSES, sep="\t", comment="#", dtype=str,
                        keep_default_na=False)
    return [c.strip() for c in table["class"] if c.strip()]


def split_classes(value: str) -> list[str]:
    return [c.strip() for c in str(value).split(";") if c.strip()]


def _date(value: object, where: str) -> str:
    """ISO date, "" for blank; ValueError otherwise (same formats as
    COLLECTION_DATE)."""
    from serial_ctdna import parse_date
    day = parse_date(value, where)
    return day.isoformat() if day else ""


def _pick(row: dict, name: str) -> str:
    for key, value in row.items():
        if str(key).strip().upper() == name:
            return str(value).strip()
    return ""


def read_input(treatments: Path | None, patient_tsv: str | Path | None) -> pd.DataFrame:
    """Events (COLUMNS, dates as ISO) from treatments.tsv and patient.tsv's
    DIAGNOSIS_DATE, checked as the module docstring says."""
    rows: list[dict] = []
    if patient_tsv and str(patient_tsv).strip() and Path(patient_tsv).exists():
        patients = pd.read_csv(patient_tsv, sep="\t", dtype=str, keep_default_na=False)
        for row in patients.to_dict("records"):
            patient = _pick(row, "PATIENT_ID")
            day = _date(_pick(row, DIAGNOSIS_COLUMN), f"{DIAGNOSIS_COLUMN} of patient {patient}")
            if patient and day:
                rows.append({"PATIENT_ID": patient, "EVENT_TYPE": DIAGNOSIS, "START_DATE": day})
    if treatments is not None:
        table = pd.read_csv(treatments, sep="\t", dtype=str, keep_default_na=False, comment="#")
        known = {c.lower(): c for c in therapy_classes()}
        for n, row in enumerate(table.to_dict("records"), start=2):
            rows.append(_check_row(row, f"{Path(treatments).name} line {n}", known))
    return pd.DataFrame(rows, columns=COLUMNS).fillna("")


def _check_row(row: dict, where: str, known: dict[str, str]) -> dict:
    patient = _pick(row, "PATIENT_ID")
    kind = _pick(row, "EVENT_TYPE").upper()
    if not patient or not kind or not _pick(row, "START_DATE"):
        msg = f"{where}: PATIENT_ID, EVENT_TYPE and START_DATE are required."
        raise ValueError(msg)
    if kind not in (TREATMENT, STATUS):
        msg = f"{where}: EVENT_TYPE '{kind}' - use {TREATMENT} or {STATUS}."
        raise ValueError(msg)
    out = {c: "" for c in COLUMNS}
    out.update(PATIENT_ID=patient, EVENT_TYPE=kind,
               START_DATE=_date(_pick(row, "START_DATE"), f"{where}, START_DATE"),
               STOP_DATE=_date(_pick(row, "STOP_DATE"), f"{where}, STOP_DATE"))
    if kind == TREATMENT:
        if out["STOP_DATE"] and out["STOP_DATE"] < out["START_DATE"]:
            msg = f"{where}: STOP_DATE {out['STOP_DATE']} is before START_DATE {out['START_DATE']}."
            raise ValueError(msg)
        classes = []
        for name in split_classes(_pick(row, "THERAPY_CLASS")):
            if name.lower() not in known:
                msg = (f"{where}: THERAPY_CLASS '{name}' isn't in {THERAPY_CLASSES.name} - "
                       f"use one of: {', '.join(known.values())}.")
                raise ValueError(msg)
            classes.append(known[name.lower()])
        line = _pick(row, "LINE")
        if line and not line.isdigit():
            msg = f"{where}: LINE '{line}' must be a number (1, 2, ...)."
            raise ValueError(msg)
        out.update(AGENT=_pick(row, "AGENT"), THERAPY_CLASS="; ".join(classes), LINE=line)
        if not out["AGENT"] and not classes:
            msg = f"{where}: a TREATMENT needs an AGENT or a THERAPY_CLASS."
            raise ValueError(msg)
        if not classes:
            logger.warning(f"{where}: no THERAPY_CLASS - this treatment can't be linked to "
                           "resistance mutations.")
        return out
    status = _pick(row, "STATUS").capitalize()
    if status not in STATUSES:
        msg = f"{where}: STATUS '{_pick(row, 'STATUS')}' - use one of {', '.join(STATUSES)}."
        raise ValueError(msg)
    recist = _pick(row, "RECIST").upper()
    if status == "Response" and recist not in RECIST:
        msg = f"{where}: a Response needs RECIST = {', '.join(RECIST)} (got '{recist}')."
        raise ValueError(msg)
    if recist and status != "Response":
        logger.warning(f"{where}: RECIST is only read for STATUS Response - ignored.")
        recist = ""
    if out["STOP_DATE"]:
        logger.warning(f"{where}: STOP_DATE is only read for treatments - ignored.")
        out["STOP_DATE"] = ""
    out.update(STATUS=status, RECIST=recist)
    return out


def check_input(sample_tsv: str | Path | None, patient_tsv: str | Path | None,
                sample_table: pd.DataFrame, sample_dates: pd.DataFrame) -> pd.DataFrame:
    """Read and check the timeline input (ValueError on the errors listed
    in the module docstring) and warn about the rest. Returns the events."""
    events = read_input(treatments_path(sample_tsv, patient_tsv), patient_tsv)
    if events.empty:
        return events
    known = set(sample_table["PATIENT_ID"].astype(str).str.strip()) if "PATIENT_ID" in sample_table else set()
    unknown = sorted(set(events["PATIENT_ID"]) - known)
    if unknown:
        logger.warning(f"Timeline events of patients not in sample.tsv are left out: {', '.join(unknown)}")
    for patient, group in events.groupby("PATIENT_ID"):
        deaths = sorted(group.loc[group["STATUS"] == "Death", "START_DATE"])
        if len(deaths) > 1:
            logger.warning(f"Patient {patient} has {len(deaths)} Death events.")
        if deaths:
            later = group[(group["START_DATE"] > deaths[0]) | (group["STOP_DATE"] > deaths[0])]
            dated = sample_dates[(sample_dates["PATIENT_ID"] == patient)
                                 & (sample_dates["COLLECTION_DATE"] > deaths[0])]
            if len(later) or len(dated):
                logger.warning(f"Patient {patient} has events or samples after death ({deaths[0]}).")
        diagnosis = group.loc[group["EVENT_TYPE"] == DIAGNOSIS, "START_DATE"]
        if len(diagnosis) and (group["START_DATE"] < diagnosis.iloc[0]).any():
            logger.warning(f"Patient {patient} has events before the diagnosis date.")
    with_events = (set(events["PATIENT_ID"]) | set(sample_dates["PATIENT_ID"])) & known
    diagnosed = set(events.loc[events["EVENT_TYPE"] == DIAGNOSIS, "PATIENT_ID"]) & known
    if diagnosed and diagnosed != with_events:
        logger.warning(
            f"{DIAGNOSIS_COLUMN} is missing for patient(s) {', '.join(sorted(with_events - diagnosed))}: "
            "the timeline's day 0 is each patient's first event, for every patient.")
    logger.info(f"Clinical timeline: {len(events)} event(s) for {len(with_events)} patient(s).")
    return events


def write_events(folder: str | Path, sample_tsv: str | Path | None,
                 patient_tsv: str | Path | None) -> None:
    """At creation: keep the timeline input in service/ and leave
    DIAGNOSIS_DATE out of data_clinical_patient.txt."""
    drop_sample_columns(folder, [DIAGNOSIS_COLUMN], PATIENT_FILE)
    events = read_input(treatments_path(sample_tsv, patient_tsv), patient_tsv)
    path = Path(folder) / EVENTS_FILE
    path.unlink(missing_ok=True)
    if events.empty:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    events.to_csv(path, sep="\t", index=False)


def carry_events(sources: list[str | Path], output: str | Path) -> None:
    """update/extract/remove: the output study gets its sources' events. A
    patient's events come from the latest source that has any - a new
    batch's treatments.tsv lists all of that patient's events, so a
    correction replaces the old ones instead of adding to them."""
    merged: pd.DataFrame | None = None
    for source in sources:
        path = Path(source) / EVENTS_FILE
        if not path.exists():
            continue
        events = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if merged is not None:
            replaced = sorted(set(merged["PATIENT_ID"]) & set(events["PATIENT_ID"]))
            if replaced:
                logger.info(f"Timeline events of patient(s) {', '.join(replaced)} replaced by {source}'s.")
            merged = pd.concat([merged[~merged["PATIENT_ID"].isin(replaced)], events])
        else:
            merged = events
    if merged is None:
        return
    out = Path(output) / EVENTS_FILE
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep="\t", index=False)


# ------------------------------------------------------------ the study

@dataclass
class Treatment:
    start: date
    stop: date | None
    agent: str
    classes: list[str]
    line: str

    def label(self) -> str:
        name = self.agent or ", ".join(self.classes)
        return f"{name} (line {self.line})" if self.line else name


@dataclass
class Timeline:
    events: pd.DataFrame
    day_zero: dict[str, date] = field(default_factory=dict)
    reference: str = "first sample"

    def day(self, patient: str, when: date) -> int:
        return (when - self.day_zero[patient]).days

    def treatments(self, patient: str) -> list[Treatment]:
        rows = self.events[(self.events["PATIENT_ID"] == patient)
                           & (self.events["EVENT_TYPE"] == TREATMENT)]
        out = [Treatment(date.fromisoformat(r["START_DATE"]),
                         date.fromisoformat(r["STOP_DATE"]) if r["STOP_DATE"] else None,
                         r["AGENT"], split_classes(r["THERAPY_CLASS"]), r["LINE"])
               for r in rows.to_dict("records")]
        return sorted(out, key=lambda t: t.start)

    def has_treatments(self) -> bool:
        return (self.events["EVENT_TYPE"] == TREATMENT).any()


def in_force(treatments: list[Treatment], when: date) -> Treatment | None:
    """The treatment a sample collected on `when` is on: the latest one
    started before that day and not stopped before it. A sample taken on
    the start day counts as pre-treatment."""
    current = [t for t in treatments if t.start < when and (t.stop is None or when <= t.stop)]
    return current[-1] if current else None


def received(treatments: list[Treatment], classes: list[str], before: date) -> Treatment | None:
    """The first treatment of one of `classes` started before `before`."""
    wanted = {c.lower() for c in classes}
    return next((t for t in treatments
                 if t.start < before and wanted & {c.lower() for c in t.classes}), None)


def annotate_study(folder: str | Path, patient_of: dict[str, str],
                   sample_dates: dict[str, date]) -> Timeline:
    """Keep the study's own patients' events, set day 0 and write the
    treatment and status timelines (removed when there's nothing to show).
    `sample_dates` are the collection dates (serial_ctdna.py)."""
    folder = Path(folder)
    path = folder / EVENTS_FILE
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    events = pd.DataFrame(columns=COLUMNS)
    if path.exists():
        events = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        events = events[events["PATIENT_ID"].isin(set(patient_of.values()))]
        events.to_csv(path, sep="\t", index=False)

    firsts: dict[str, date] = {}
    for sample_id, when in sample_dates.items():
        patient = patient_of.get(sample_id, sample_id)
        firsts[patient] = min(when, firsts.get(patient, when))
    for r in events[events["EVENT_TYPE"] != DIAGNOSIS].to_dict("records"):
        when = date.fromisoformat(r["START_DATE"])
        firsts[r["PATIENT_ID"]] = min(when, firsts.get(r["PATIENT_ID"], when))
    diagnosis = {r["PATIENT_ID"]: date.fromisoformat(r["START_DATE"])
                 for r in events[events["EVENT_TYPE"] == DIAGNOSIS].to_dict("records")}
    timeline = Timeline(events)
    if diagnosis and set(firsts) <= set(diagnosis):
        timeline.day_zero, timeline.reference = diagnosis, "diagnosis"
    else:
        # A patient with nothing but a diagnosis date: that is their first event.
        timeline.day_zero = {**diagnosis, **firsts}
        timeline.reference = "first event" if len(events) else "first sample"

    _write(folder, events, timeline)
    if len(events):
        summary = {
            "reference": timeline.reference,
            "patients": int(events["PATIENT_ID"].nunique()),
            "treatments": int((events["EVENT_TYPE"] == TREATMENT).sum()),
            "status": int((events["EVENT_TYPE"] == STATUS).sum()),
            "diagnosis": len(diagnosis),
        }
        (folder / SUMMARY_FILE).write_text(json.dumps(summary, indent=2), encoding="utf-8")
        logger.info(f"Clinical timeline: day 0 = {timeline.reference}; {summary['treatments']} "
                    f"treatment(s), {summary['status']} status event(s).")
    return timeline


def _write(folder: Path, events: pd.DataFrame, timeline: Timeline) -> None:
    treatments, status = [], []
    for r in events.to_dict("records"):
        patient = r["PATIENT_ID"]
        start = timeline.day(patient, date.fromisoformat(r["START_DATE"]))
        if r["EVENT_TYPE"] == TREATMENT:
            stop = (timeline.day(patient, date.fromisoformat(r["STOP_DATE"]))
                    if r["STOP_DATE"] else "")
            treatments.append({
                "PATIENT_ID": patient, "START_DATE": start, "STOP_DATE": stop,
                "EVENT_TYPE": TREATMENT, "TREATMENT_TYPE": "Medical Therapy",
                "AGENT": r["AGENT"] or "NA", "THERAPY_CLASS": r["THERAPY_CLASS"] or "NA",
                "LINE": r["LINE"] or "NA"})
        else:
            status.append({
                "PATIENT_ID": patient, "START_DATE": start, "STOP_DATE": "",
                "EVENT_TYPE": STATUS,
                "STATUS": "Diagnosis" if r["EVENT_TYPE"] == DIAGNOSIS else r["STATUS"],
                "RECIST": r["RECIST"] or "NA"})
    for name, rows in ((TREATMENT_TIMELINE, treatments), (STATUS_TIMELINE, status)):
        if rows:
            pd.DataFrame(rows).sort_values(["PATIENT_ID", "START_DATE"]).to_csv(
                folder / name, sep="\t", index=False)
        else:
            (folder / name).unlink(missing_ok=True)


def read_summary(folder: str | Path) -> dict | None:
    path = Path(folder) / SUMMARY_FILE
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def create_meta_timelines(project_id: str, folder: str | Path) -> None:
    for data, meta in METAS.items():
        if not (Path(folder) / data).exists():
            continue
        (Path(folder) / meta).write_text(
            f"cancer_study_identifier: {project_id}\n"
            "genetic_alteration_type: CLINICAL\n"
            "datatype: TIMELINE\n"
            f"data_filename: {data}\n", encoding="utf-8")
        logger.info(f"{meta} created!")
