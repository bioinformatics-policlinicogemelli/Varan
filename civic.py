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

"""CIViC (civicdb.org, CC0) clinical evidence for the study's alterations.

CIViC's TSV dumps (about 5 MB) are kept in [CIViC] DIR (default: the civic
folder in Varan's): downloaded on first use when missing (AUTO_DOWNLOAD),
or refreshed with `python civic.py --download`. release.json there records
when they were downloaded; the report shows it. No token is needed, and
matching is local: it takes seconds.

Matched: single-variant molecular profiles - protein changes (V600E),
codons (G12), exon deletions/insertions, amplifications, deletions and
fusions (both known partners, or the known one) - see study_alterations.py.
Gene-level profiles ("Mutation", "Loss-of-function", expression...) and
combinations ("A AND B") are left out. Evidence: accepted, not flagged,
"Supports" direction; levels from predictive evidence only (germline-origin
predictive items, e.g. BRCA1/2 and PARP inhibitors, are kept).

CIViC diseases aren't mapped to OncoTree: the levels are for ANY tumour
type, and the report and the MAF text name the disease of each item.

Output, separate from OncoKB's (no column shared):
- data_mutations_extended.txt, namespace CIVIC: CIVIC.SENSITIVITY and
  CIVIC.RESISTANCE (best predictive evidence level, A best ... E),
  CIVIC.EVIDENCE (level A-B predictive items: therapy, disease, level, EID),
  CIVIC.ASSERTIONS (AMP/ASCO/CAP tier assertions: tier, therapy, disease, AID).
- data_clinical_sample.txt: CIVIC_HIGHEST_LEVEL, CIVIC_RESISTANCE_LEVEL,
  CIVIC_N_ACTIONABLE / CIVIC_ACTIONABLE (sensitivity at an [Actionability]
  CIVIC_LEVELS level, default A and B), over mutations, fusions and CNAs.

References: Griffith et al., Nat Genet 2017, doi:10.1038/ng.3774; levels
from CIViC's evidence model (A validated, B clinical, C case study,
D preclinical, E inferential).
"""

from __future__ import annotations

import json
import sys
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd
from loguru import logger

import study_alterations as sa
from clinical_file import drop_sample_columns, set_sample_column
from config_loader import VARAN_DIR, config

URL = "https://civicdb.org/downloads/nightly/{name}"
FILES = {"evidence": "nightly-ClinicalEvidenceSummaries.tsv",
         "assertions": "nightly-AssertionSummaries.tsv",
         "variants": "nightly-VariantSummaries.tsv"}
RELEASE = "release.json"
SUMMARY_FILE = Path("service") / "civic.json"
NAMESPACE = "CIVIC"
MAF_COLUMNS = ["CIVIC.SENSITIVITY", "CIVIC.RESISTANCE", "CIVIC.EVIDENCE", "CIVIC.ASSERTIONS"]
LEVELS = ["A", "B", "C", "D", "E"]
COLUMNS = {
    "CIVIC_HIGHEST_LEVEL": ("CIViC Highest Level",
                           "Best CIViC predictive sensitivity evidence level of the sample's alterations "
                           "(any tumour type)", "STRING"),
    "CIVIC_RESISTANCE_LEVEL": ("CIViC Resistance Level",
                               "Best CIViC predictive resistance evidence level (any tumour type)", "STRING"),
    "CIVIC_N_ACTIONABLE": ("CIViC Actionable Alterations",
                           "Alterations with CIViC sensitivity evidence at [Actionability] CIVIC_LEVELS",
                           "NUMBER"),
    "CIVIC_ACTIONABLE": ("CIViC Actionable", "Those alterations, with their level", "STRING"),
}
MAX_ITEMS = 5


def _bool(key: str, default: str) -> bool:
    return config.get("CIViC", key, fallback=default).strip().lower() in ("true", "t", "yes", "1")


def data_dir() -> Path:
    value = config.get("CIViC", "DIR", fallback="").strip()
    return Path(value) if value else VARAN_DIR / "civic"


def actionable_levels() -> list[str]:
    value = config.get("Actionability", "CIVIC_LEVELS", fallback="").strip() or "A,B"
    return [v.strip().upper() for v in value.split(",") if v.strip()]


# ------------------------------------------------------------ the data

def download(directory: Path | None = None) -> Path:
    directory = Path(directory or data_dir())
    directory.mkdir(parents=True, exist_ok=True)
    for name in FILES.values():
        with urllib.request.urlopen(URL.format(name=name), timeout=120) as response:
            (directory / name).write_bytes(response.read())
    (directory / RELEASE).write_text(json.dumps({
        "downloaded": date.today().isoformat(), "source": URL.format(name="")}, indent=2),
        encoding="utf-8")
    logger.info(f"CIViC data downloaded to {directory}")
    return directory


def ensure_data() -> Path | None:
    """The CIViC folder, downloading it when missing and AUTO_DOWNLOAD is
    on; None (with a warning) when it can't be had."""
    directory = data_dir()
    if all((directory / n).exists() for n in FILES.values()):
        return directory
    if not _bool("AUTO_DOWNLOAD", "True"):
        logger.warning(f"CIViC data not found in {directory} and [CIViC] AUTO_DOWNLOAD is off: "
                       "run `python civic.py --download`. CIViC annotation skipped.")
        return None
    try:
        return download(directory)
    except Exception as err:  # noqa: BLE001 - no network: skip, don't fail the run
        logger.warning(f"CIViC data couldn't be downloaded ({err}): CIViC annotation skipped.")
        return None


def release(directory: Path) -> dict:
    path = directory / RELEASE
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


class Knowledge:
    """CIViC's single-variant profiles, as matchers, with their evidence
    and assertions."""

    def __init__(self, directory: Path) -> None:
        read = lambda key: pd.read_csv(directory / FILES[key], sep="\t", dtype=str,  # noqa: E731
                                       keep_default_na=False, quoting=3)
        variants, evidence, assertions = read("variants"), read("evidence"), read("assertions")
        self.mutations: dict[str, list[tuple[sa.Spec, str]]] = {}   # gene -> [(spec, profile)]
        self.cnas: dict[tuple[str, str], list[str]] = {}            # (gene, kind) -> [profile]
        self.fusions: list[tuple[str, str, str]] = []               # (5', 3', profile)
        for v in variants.to_dict("records"):
            profile = v.get("single_variant_molecular_profile_id", "")
            if not profile:
                continue
            if v["feature_type"] == "Fusion":
                five, three = v.get("5_prime_partner", ""), v.get("3_prime_partner", "")
                if five or three:
                    self.fusions.append((five, three, profile))
                continue
            if v["feature_type"] != "Gene":
                continue
            spec = sa.parse_one(v["variant"])
            if spec is None or (spec.kind == sa.MUTATION and not (spec.exact or spec.codon or spec.exon)):
                continue   # gene-level profiles: too broad for any tumour type
            gene = v.get("gene") or v.get("feature_name")
            if spec.kind == sa.MUTATION:
                self.mutations.setdefault(gene, []).append((spec, profile))
            elif spec.kind == sa.FUSION:
                self.fusions.append(("", gene, profile))
            else:
                self.cnas.setdefault((gene, spec.kind), []).append(profile)
        evidence = evidence[(evidence["evidence_direction"] == "Supports")
                            & (evidence["is_flagged"] != "true")]
        self.evidence = {k: g.to_dict("records") for k, g in evidence.groupby("molecular_profile_id")}
        assertions = assertions[(assertions["amp_category"] != "") & (assertions["is_flagged"] != "true")]
        self.assertions = {k: g.to_dict("records") for k, g in assertions.groupby("molecular_profile_id")}

    def for_mutation(self, row: dict) -> list[str]:
        return [p for spec, p in self.mutations.get(str(row.get("Hugo_Symbol", "")), [])
                if spec.matches_mutation(row)]

    def for_fusion(self, row: dict) -> list[str]:
        out = []
        for five, three, profile in self.fusions:
            known = [g for g in (five, three) if g]
            if all(sa.fusion_matches(g, row) for g in known):
                out.append(profile)
        return out

    def for_cna(self, call: dict) -> list[str]:
        kind = sa.AMPLIFICATION if call["value"] > 0 else sa.DELETION
        return self.cnas.get((call["gene"], kind), [])

    def describe(self, profiles: list[str]) -> dict:
        items = [e for p in profiles for e in self.evidence.get(p, [])]
        predictive = [e for e in items if e["evidence_type"] == "Predictive"]
        sens = [e["evidence_level"] for e in predictive if e["significance"] == "Sensitivity/Response"]
        res = [e["evidence_level"] for e in predictive if e["significance"] == "Resistance"]
        top = sorted((e for e in predictive if e["evidence_level"] in ("A", "B")),
                     key=lambda e: (e["evidence_level"], e["therapies"]))
        grouped: dict[tuple, list[str]] = {}
        for e in top:
            effect = "resistance" if e["significance"] == "Resistance" else "sensitivity"
            key = (e["therapies"] or "no therapy", effect, e["disease"], e["evidence_level"])
            grouped.setdefault(key, []).append(f"EID{e['evidence_id']}")
        texts = [f"{t} ({effect}) - {disease} [{level}, {', '.join(eids)}]"
                 for (t, effect, disease, level), eids in grouped.items()]
        asserted = list(dict.fromkeys(
            f"{a['amp_category']}: {a['therapies'] or a['assertion_type']} - {a['disease']} (AID{a['assertion_id']})"
            for p in profiles for a in self.assertions.get(p, [])))
        return {"sensitivity": _best(sens), "resistance": _best(res),
                "evidence": _cap(texts), "assertions": _cap(asserted)}


def _best(levels: list[str]) -> str:
    return next((lv for lv in LEVELS if lv in levels), "")


def _cap(items: list[str]) -> str:
    extra = len(items) - MAX_ITEMS
    return "; ".join(items[:MAX_ITEMS]) + (f"; +{extra} more" if extra > 0 else "")


# ------------------------------------------------------------- study

def annotate_study(folder: str | Path) -> dict | None:
    """Write the MAF columns and sample attributes of the module docstring;
    None (and those removed) when CIViC is off or its data is missing."""
    folder = Path(folder)
    (folder / SUMMARY_FILE).unlink(missing_ok=True)
    maf_path = folder / sa.MAF_FILE
    maf = sa.read_maf(folder)
    if maf is not None and any(c in maf.columns for c in MAF_COLUMNS):
        maf = maf.drop(columns=[c for c in MAF_COLUMNS if c in maf.columns])
        maf.to_csv(maf_path, sep="\t", index=False)
    directory = ensure_data() if _bool("ENABLED", "True") else None
    from actionability import _samples, cohort_summary
    samples = _samples(folder)
    if directory is None or samples is None:
        drop_sample_columns(folder, list(COLUMNS))
        return None
    knowledge = Knowledge(directory)

    found: dict[str, list[dict]] = {}
    if maf is not None and not maf.empty:
        described = [knowledge.describe(knowledge.for_mutation(r)) for r in maf.to_dict("records")]
        for column, key in zip(MAF_COLUMNS, ("sensitivity", "resistance", "evidence", "assertions")):
            maf[column] = [d[key] for d in described]
        maf.to_csv(maf_path, sep="\t", index=False)
        for r, d in zip(maf.to_dict("records"), described):
            if d["sensitivity"] or d["resistance"] or d["assertions"]:
                found.setdefault(str(r["Tumor_Sample_Barcode"]).strip(), []).append(
                    {"alteration": sa.mutation_label(r), **d})
    for r in sa.fusions(sa.read_sv(folder)).to_dict("records"):
        d = knowledge.describe(knowledge.for_fusion(r))
        if d["sensitivity"] or d["resistance"] or d["assertions"]:
            found.setdefault(str(r["Sample_Id"]).strip(), []).append({"alteration": sa.fusion_label(r), **d})
    for call in sa.read_cna(folder):
        d = knowledge.describe(knowledge.for_cna(call))
        if d["sensitivity"] or d["resistance"] or d["assertions"]:
            found.setdefault(call["sample"], []).append({"alteration": sa.cna_label(call), **d})

    wanted = actionable_levels()
    attrs: dict[str, dict[str, object]] = {c: {} for c in COLUMNS}
    per_sample = {}
    for sample in samples["SAMPLE_ID"]:
        hits = found.get(sample, [])
        best = _best([h["sensitivity"] for h in hits])
        actionable = sorted({f"{h['alteration']} ({h['sensitivity']})" for h in hits
                             if h["sensitivity"] in wanted})
        attrs["CIVIC_HIGHEST_LEVEL"][sample] = best or "None"
        attrs["CIVIC_RESISTANCE_LEVEL"][sample] = _best([h["resistance"] for h in hits]) or "None"
        attrs["CIVIC_N_ACTIONABLE"][sample] = len(actionable)
        attrs["CIVIC_ACTIONABLE"][sample] = ", ".join(actionable) or "None"
        per_sample[sample] = {"highest": best or "None", "actionable": actionable, "hits": hits}
    for column, (short, long, datatype) in COLUMNS.items():
        set_sample_column(folder, column, attrs[column], short, long, datatype)

    summary = cohort_summary(samples, per_sample, wanted, LEVELS)
    summary["release"] = release(directory)
    (folder / SUMMARY_FILE).parent.mkdir(parents=True, exist_ok=True)
    (folder / SUMMARY_FILE).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info(f"CIViC: {summary['n_actionable']}/{summary['n_samples']} sample(s) with level "
                f"{'/'.join(wanted)} evidence (any tumour type); data of "
                f"{summary['release'].get('downloaded', 'unknown date')}.")
    return summary


def read_summary(folder: str | Path) -> dict | None:
    path = Path(folder) / SUMMARY_FILE
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


if __name__ == "__main__":
    if "--download" in sys.argv:
        target = sys.argv[sys.argv.index("--download") + 1] if len(sys.argv) > sys.argv.index("--download") + 1 \
            else None
        download(Path(target) if target else VARAN_DIR / "civic")
    else:
        print("usage: python civic.py --download [folder]")
