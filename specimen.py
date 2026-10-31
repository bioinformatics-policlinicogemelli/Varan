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

"""Specimen type of each sample: Tissue or Plasma.

The study-wide value is conf.ini's [Sample_Type] TYPE. An optional
SPECIMEN_TYPE column in sample.tsv overrides it for the samples that fill
it in (blank cells keep the study value), so one study can hold tissue and
plasma samples - e.g. after merging a tissue run and a plasma run with -u.
Every sample's type is written to data_clinical_sample.txt's SPECIMEN_TYPE
column, so it survives update/extract/remove.

Settings that depend on the specimen type can be overridden per type with a
"<section>.<type>" section in conf.ini, the same way [Vendor.guardant]
overrides the Guardant adapter's settings: [Filters.Plasma] t_VAF_min is
used for plasma samples instead of [Filters] t_VAF_min, and so on for
[MSI], [TMB], [SigMA] and [GenePanel]. See setting().

"Solid" and "Liquid", the values TYPE used to take, are still accepted
(Solid = Tissue, Liquid = Plasma) with a deprecation warning.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
from loguru import logger

from config_loader import config

TISSUE = "Tissue"
PLASMA = "Plasma"
SPECIMEN_TYPES = (TISSUE, PLASMA)
COLUMN = "SPECIMEN_TYPE"

_SYNONYMS = {"tissue": TISSUE, "plasma": PLASMA, "solid": TISSUE, "liquid": PLASMA}
_LEGACY = {"solid", "liquid"}

# Specimen type a [Vendor] PIPELINE value implies (see varan.py's
# reconcile_sample_type()). Vendor adapters declare theirs as a module-level
# SAMPLE_TYPE attribute instead.
ILLUMINA_PIPELINE_SAMPLE_TYPES = {
    "illumina_solid": TISSUE,
    "illumina_liquid": PLASMA,
}

_MISSING = object()
_warned_legacy = False


def normalize(value: object, where: str = "") -> str | None:
    """The specimen type `value` names (Tissue/Plasma), None for a blank
    value. Raises ValueError for anything else."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip().strip('"').strip("'")
    if text == "" or text.upper() in ("NA", "NAN"):
        return None
    try:
        return _SYNONYMS[text.lower()]
    except KeyError:
        msg = (f"{where + ': ' if where else ''}'{text}' is not a specimen type - "
               "use Tissue or Plasma.")
        raise ValueError(msg) from None


def study_specimen_type() -> str:
    """The study-wide specimen type, from conf.ini's [Sample_Type] TYPE."""
    global _warned_legacy
    raw = config.get("Sample_Type", "TYPE", fallback="")
    value = normalize(raw, "conf.ini [Sample_Type] TYPE")
    if value is None:
        msg = 'Please set [Sample_Type] TYPE in conf.ini to "Tissue" or "Plasma".'
        raise ValueError(msg)
    if raw.strip().strip('"').strip("'").lower() in _LEGACY and not _warned_legacy:
        _warned_legacy = True
        logger.warning(
            f"conf.ini [Sample_Type] TYPE = {raw.strip()} is deprecated: read as "
            f"'{value}'. Use Tissue or Plasma.")
    return value


def implied_by_pipeline(pipeline: str) -> str | None:
    """The specimen type a [Vendor] PIPELINE value implies, or None."""
    if pipeline in ILLUMINA_PIPELINE_SAMPLE_TYPES:
        return ILLUMINA_PIPELINE_SAMPLE_TYPES[pipeline]
    from vendor_adapters import ADAPTERS
    adapter = ADAPTERS.get(pipeline)
    declared = getattr(adapter, "SAMPLE_TYPE", None) if adapter else None
    return normalize(declared) if declared else None


def _column(df: pd.DataFrame, name: str) -> str | None:
    """The column of `df` named `name`, case-insensitively."""
    for col in df.columns:
        if str(col).strip().upper() == name.upper():
            return col
    return None


def specimen_types(sample_tsv: pd.DataFrame | None, log: bool = False) -> dict[str, str]:
    """{SAMPLE_ID: Tissue/Plasma} for the samples of a sample.tsv (or of
    data_clinical_sample.txt): the SPECIMEN_TYPE column where filled in,
    the study-wide type otherwise. With `log`, warns about the samples whose
    type overrides the study's. Raises ValueError for unknown values."""
    default = study_specimen_type()
    if sample_tsv is None or "SAMPLE_ID" not in sample_tsv.columns:
        return {}
    col = _column(sample_tsv, COLUMN)
    types: dict[str, str] = {}
    overridden: dict[str, list[str]] = {}
    for row in sample_tsv.to_dict("records"):
        sample_id = str(row["SAMPLE_ID"]).strip()
        given = normalize(row.get(col), f"{COLUMN} of sample {sample_id}") if col else None
        types[sample_id] = given or default
        if given and given != default:
            overridden.setdefault(given, []).append(sample_id)
    if log:
        for value, samples in overridden.items():
            logger.warning(
                f"{COLUMN} = {value} overrides [Sample_Type] TYPE = {default} for "
                f"{len(samples)} sample(s): {', '.join(samples)}.")
    return types


def read_sample_table(path: str | Path) -> pd.DataFrame | None:
    """sample.tsv as strings, or None when `path` isn't a file."""
    path = Path(path) if path else None
    if path is None or not path.is_file():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_study_samples(folder: str | Path) -> pd.DataFrame | None:
    """data_clinical_sample.txt's records (header rows skipped), or None."""
    path = Path(folder) / "data_clinical_sample.txt"
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t", dtype=str, comment="#", keep_default_na=False)
    except pd.errors.EmptyDataError:
        return None


def study_types(folder: str | Path) -> dict[str, str]:
    """{SAMPLE_ID: Tissue/Plasma} from a study's data_clinical_sample.txt.
    Samples without a (valid) SPECIMEN_TYPE are left out."""
    samples = read_study_samples(folder)
    if samples is None or "SAMPLE_ID" not in samples.columns:
        return {}
    col = _column(samples, COLUMN)
    if col is None:
        return {}
    types = {}
    for sample_id, value in zip(samples["SAMPLE_ID"], samples[col]):
        try:
            kind = normalize(value)
        except ValueError:
            kind = None
        if kind:
            types[str(sample_id).strip()] = kind
    return types


def check_pipeline(types: dict[str, str], pipeline: str) -> None:
    """Stop the run when a sample's specimen type contradicts the one its
    pipeline implies (e.g. a Guardant sample marked Tissue)."""
    implied = implied_by_pipeline(pipeline) if pipeline else None
    if implied is None:
        return
    wrong = sorted(s for s, t in types.items() if t != implied)
    if wrong:
        msg = (f"[Vendor] PIPELINE = {pipeline} is {implied}-only, but {COLUMN} marks "
               f"{len(wrong)} sample(s) as {types[wrong[0]]}: {', '.join(wrong)}. "
               f"Fix {COLUMN} in sample.tsv (or leave it blank).")
        raise ValueError(msg)


def setting(section: str, key: str, specimen_type: str | None = None,
            fallback: object = _MISSING) -> str:
    """conf.ini [section] key, taken from [section.<specimen_type>] when that
    section sets it."""
    if specimen_type:
        override = f"{section}.{specimen_type}"
        if config.has_option(override, key):
            return config.get(override, key)
    if fallback is _MISSING:
        return config.get(section, key)
    return config.get(section, key, fallback=fallback)


def overrides() -> dict[str, dict[str, str]]:
    """{"<section>.<type>": {key: value}} for every per-type section in
    conf.ini - shown in the report."""
    found = {}
    for section in config.sections():
        base, _, kind = section.rpartition(".")
        if base and kind in SPECIMEN_TYPES:
            found[section] = dict(config.items(section, raw=True))
    return found


def warn_missing_overrides(types: dict[str, str]) -> None:
    """When SPECIMEN_TYPE gives some samples a type other than the study's,
    warn if that type has no override section at all: its samples use the
    general settings, meant for the study's type."""
    present = set(types.values()) - {study_specimen_type()}
    defined = {s.rpartition(".")[2] for s in overrides()}
    for kind in sorted(present - defined):
        logger.warning(
            f"SPECIMEN_TYPE makes some samples {kind}, but conf.ini has no "
            f"[<section>.{kind}] section: they use the general "
            "[Filters]/[MSI]/[TMB]/[SigMA]/[GenePanel] settings.")


def backfill_specimen_type(folder: str | Path, sources: list[str | Path]) -> None:
    """After merging studies (update): give SPECIMEN_TYPE to the samples of
    a source study created before the column existed, from that study's
    provenance (its [Sample_Type] TYPE then: Solid -> Tissue, Liquid ->
    Plasma). Samples whose type can't be told stay "NA"."""
    from clinical_file import set_sample_column
    from write_report import read_provenance

    merged = read_study_samples(folder)
    if merged is None or "SAMPLE_ID" not in merged.columns:
        return
    types = study_types(folder)
    filled = []
    for source in sources:
        source_samples = read_study_samples(source)
        if (source_samples is None or "SAMPLE_ID" not in source_samples.columns
                or _column(source_samples, COLUMN)):
            continue
        try:
            kind = normalize((read_provenance(str(source)) or {}).get("sample_type"))
        except ValueError:
            kind = None
        if kind is None:
            logger.warning(
                f"{Path(source).name}: no {COLUMN} and no single sample type in its "
                f"provenance - its samples' {COLUMN} is left NA.")
            continue
        for sample_id in source_samples["SAMPLE_ID"].astype(str).str.strip():
            if sample_id not in types:
                types[sample_id] = kind
                filled.append(sample_id)
    if filled or _column(merged, COLUMN) is None:
        set_sample_column(
            folder, COLUMN, types, "Specimen Type",
            "Specimen type of the sequenced sample (Tissue or Plasma)")
    if filled:
        logger.info(f"{COLUMN} filled from the source studies' provenance for "
                    f"{len(filled)} sample(s).")


def describe(types: dict[str, str] | list[str]) -> str:
    """"Tissue", "Plasma" or "Mixed (3 Tissue, 2 Plasma)"."""
    values = list(types.values()) if isinstance(types, dict) else list(types)
    counts = Counter(values)
    if not counts:
        return ""
    if len(counts) == 1:
        return next(iter(counts))
    return "Mixed (" + ", ".join(f"{counts[t]} {t}" for t in SPECIMEN_TYPES if t in counts) + ")"
