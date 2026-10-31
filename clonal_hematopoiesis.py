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

"""Flag plasma variants that may come from clonal hematopoiesis (CH).

Most cell-free DNA comes from white blood cells (WBC), so a plasma sample
also carries the mutations of the patient's blood clones: over half of the
cfDNA mutations of patients with cancer had features of CH in Razavi et al.
2019 (Nat Med, doi:10.1038/s41591-019-0652-7). Only sequencing the matched
WBC can tell for sure; without it Varan can only flag. Research use only.

Plasma samples only (see specimen.py). Every variant of data_mutations_
extended.txt gets two columns, in cBioPortal's "CH" namespace so the portal
shows them in its mutation tables:

- CH.EVIDENCE: "CH confirmed (WBC)", "CH suspected", "CH compatible",
  "Tumor (tissue-supported)", or empty.
- CH.REASON: why, e.g. "DNMT3A (CH tier 1); absent in tissue S12
  (covered); VAF 0.4% vs sample max 8%".

Evidence, strongest first:
1. The variant is in the sample's matched WBC VCF (sample.tsv's optional
   wbc_path column) -> CH confirmed. With a WBC VCF, a variant missing from
   it isn't flagged.
2. The variant is in a Tissue sample of the same patient in the study ->
   Tumor (tissue-supported).
3. Gene and variant in ch_genes.tsv: tier 1 (typically CH) -> CH suspected;
   tier 2 (both CH and tumor) -> CH compatible, or CH suspected with one
   more piece of evidence: absent from a tissue sample of the patient that
   covers the gene, a low VAF relative to the sample's tumor variants
   ([CH] RELATIVE_VAF_MAX), or a stable VAF across the patient's plasma
   samples while their tumor variants change ([CH] SERIAL_FOLD).

The thresholds are heuristics, not validated against WBC sequencing.
Variants are never removed from the study; with [CH] EXCLUDE_FROM_DOWNSTREAM
the suspected/confirmed ones are left out of SigMA's input.

Each plasma sample also gets N_CH_SUSPECT and CH_SUSPECT_GENES in
data_clinical_sample.txt (CH confirmed + CH suspected).
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from loguru import logger

import specimen
from clinical_file import drop_sample_columns, set_sample_column
from config_loader import VARAN_DIR, config

MAF_FILE = "data_mutations_extended.txt"
NAMESPACE = "CH"
EVIDENCE = "CH.EVIDENCE"
REASON = "CH.REASON"

CONFIRMED = "CH confirmed (WBC)"
SUSPECTED = "CH suspected"
COMPATIBLE = "CH compatible"
TUMOR = "Tumor (tissue-supported)"
FLAGGED = (CONFIRMED, SUSPECTED)

WBC_COLUMN = "wbc_path"
N_COLUMN = "N_CH_SUSPECT"
GENES_COLUMN = "CH_SUSPECT_GENES"

DISCLAIMER = (
    "Without sequencing of matched white blood cells the origin of these "
    "variants cannot be confirmed. Research use only.")

TRUNCATING = {
    "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
    "Nonstop_Mutation", "Translation_Start_Site"}
MISSENSE = {"Missense_Mutation"}
CLASSES = {
    "truncating": TRUNCATING,
    "missense": MISSENSE,
    "any": TRUNCATING | MISSENSE | {"In_Frame_Del", "In_Frame_Ins"},
}


# ------------------------------------------------------------ settings

def enabled() -> bool:
    return _bool(config.get("CH", "ENABLED", fallback="True"))


def exclude_downstream() -> bool:
    return enabled() and _bool(config.get("CH", "EXCLUDE_FROM_DOWNSTREAM", fallback="False"))


def _bool(value: str) -> bool:
    return value.strip().lower() in ("true", "t", "yes", "1")


def _float(key: str, default: float) -> float:
    return float(config.get("CH", key, fallback=str(default)))


# ----------------------------------------------------------- gene list

@dataclass
class GeneRule:
    gene: str
    tier: int
    classes: set[str]
    hotspots: list[str]
    exons: set[str]

    def matches(self, row: dict) -> bool:
        protein = str(row.get("HGVSp_Short") or "")
        for spot in self.hotspots:
            if re.match(rf"^p\.{re.escape(spot)}(?![0-9])", protein):
                return True
        if not self.classes or row.get("Variant_Classification") not in self.classes:
            return False
        return not self.exons or _exon(row) in self.exons


def _exon(row: dict) -> str:
    for col in ("EXON", "Exon_Number"):
        value = str(row.get(col) or "")
        if value and value.lower() != "nan":
            return value.split("/")[0].strip()
    return ""


def gene_list_path() -> Path:
    path = Path(config.get("CH", "GENE_LIST", fallback="ch_genes.tsv").strip() or "ch_genes.tsv")
    return path if path.is_absolute() or path.exists() else VARAN_DIR / path


def load_gene_list(path: Path | None = None) -> dict[str, GeneRule]:
    table = pd.read_csv(
        path or gene_list_path(), sep="\t", comment="#", dtype=str, keep_default_na=False)
    rules = {}
    for row in table.to_dict("records"):
        classes: set[str] = set()
        for name in filter(None, (c.strip().lower() for c in row["classes"].split(","))):
            classes |= CLASSES[name]
        rules[row["gene"].strip()] = GeneRule(
            gene=row["gene"].strip(),
            tier=int(row["tier"]),
            classes=classes,
            hotspots=[h.strip() for h in row["hotspots"].split(",") if h.strip()],
            exons={e.strip() for e in row["exons"].split(",") if e.strip()},
        )
    return rules


# ------------------------------------------------------------ variants

def vaf(maf: pd.DataFrame) -> pd.Series:
    """Each variant's allele fraction: t_AF, t_VF, or t_alt_count/t_depth."""
    for col in ("t_AF", "t_VF"):
        if col in maf.columns:
            values = pd.to_numeric(maf[col], errors="coerce")
            if values.notna().any():
                return values
    if {"t_alt_count", "t_depth"} <= set(maf.columns):
        alt = pd.to_numeric(maf["t_alt_count"], errors="coerce")
        depth = pd.to_numeric(maf["t_depth"], errors="coerce")
        return alt / depth.where(depth > 0)
    return pd.Series(float("nan"), index=maf.index)


def _chrom(value: object) -> str:
    text = str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def maf_key(row: dict) -> tuple[str, str, str, str]:
    return (_chrom(row.get("Chromosome")), str(row.get("Start_Position")).strip(),
            str(row.get("Reference_Allele")).strip(), str(row.get("Tumor_Seq_Allele2")).strip())


def vcf_keys(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, str, str, str]:
    """A VCF record in MAF coordinates: indels lose their shared leading
    base(s), deletions start at the first deleted base, insertions stay on
    the base before them (as vcf2maf writes them)."""
    trimmed = 0
    while ref and alt and ref[0] == alt[0]:
        ref, alt, trimmed = ref[1:], alt[1:], trimmed + 1
    if not ref:
        start = pos + trimmed - 1
    else:
        start = pos + trimmed
    return (_chrom(chrom), str(start), ref or "-", alt or "-")


def read_wbc_variants(path: Path) -> set[tuple[str, str, str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    keys = set()
    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            for alt in fields[4].split(","):
                if alt and alt not in (".", "*"):
                    keys.add(vcf_keys(fields[0], int(fields[1]), fields[3], alt))
    return keys


# ------------------------------------------------------------ evidence

@dataclass
class Samples:
    types: dict[str, str]
    patients: dict[str, str]
    wbc: dict[str, Path]


def samples_from_table(table: pd.DataFrame | None, types: dict[str, str]) -> Samples:
    """Patients and WBC VCFs from a sample.tsv / data_clinical_sample.txt."""
    patients: dict[str, str] = {}
    wbc: dict[str, Path] = {}
    if table is not None and "SAMPLE_ID" in table.columns:
        wbc_col = next((c for c in table.columns if c.strip().lower() == WBC_COLUMN), None)
        for row in table.fillna("").to_dict("records"):
            sample_id = str(row["SAMPLE_ID"]).strip()
            patients[sample_id] = str(row.get("PATIENT_ID", "")).strip() or sample_id
            path = str(row.get(wbc_col, "")).strip() if wbc_col else ""
            if path:
                if Path(path).is_file():
                    wbc[sample_id] = Path(path)
                else:
                    logger.warning(f"{sample_id}: {WBC_COLUMN} {path} not found - ignored.")
    return Samples(types, patients, wbc)


def _covered_genes(folder: Path | None, sample_id: str, kind: str | None) -> set[str] | None:
    """Genes of the panel `sample_id` was sequenced with, None if unknown."""
    import gene_panels
    panel_id = None
    matrix = folder / gene_panels.MATRIX_FILE if folder else None
    if matrix and matrix.exists():
        table = gene_panels._read_matrix(matrix).set_index("SAMPLE_ID")
        if "mutations" in table.columns and sample_id in table.index:
            panel_id = table.at[sample_id, "mutations"]
    if not panel_id or panel_id == "NA":
        panel_id = gene_panels.panel_of(kind, "mutations")
    if not panel_id:
        return None
    genes = gene_panels.read_gene_list(panel_id)
    return set(genes) if genes is not None else None


def _pct(value: float) -> str:
    return f"{value * 100:.2g}%"


def classify(maf: pd.DataFrame, samples: Samples, folder: Path | None = None,
             rules: dict[str, GeneRule] | None = None) -> pd.DataFrame:
    """CH.EVIDENCE / CH.REASON for every row of `maf` (empty for non-plasma
    samples and unflagged variants)."""
    rules = rules if rules is not None else load_gene_list()
    relative_max = _float("RELATIVE_VAF_MAX", 0.25)
    serial_fold = _float("SERIAL_FOLD", 2.0)

    out = pd.DataFrame({EVIDENCE: "", REASON: ""}, index=maf.index)
    if maf.empty:
        return out
    rows = maf.to_dict("index")
    vafs = vaf(maf)
    sample_of = {i: str(r.get("Tumor_Sample_Barcode", "")).strip() for i, r in rows.items()}
    key_of = {i: maf_key(r) for i, r in rows.items()}
    rule_of = {}
    for i, r in rows.items():
        rule = rules.get(str(r.get("Hugo_Symbol", "")).strip())
        rule_of[i] = rule if rule and rule.matches(r) else None

    # Sample "tumor level": highest VAF among variants in genes outside the list.
    tumor_max: dict[str, float] = {}
    for i, r in rows.items():
        if str(r.get("Hugo_Symbol", "")).strip() in rules or pd.isna(vafs[i]):
            continue
        s = sample_of[i]
        tumor_max[s] = max(tumor_max.get(s, 0.0), float(vafs[i]))

    variants: dict[str, dict[tuple, float]] = {}
    for i in rows:
        variants.setdefault(sample_of[i], {})[key_of[i]] = vafs[i]

    by_patient: dict[str, list[str]] = {}
    for s in variants:
        by_patient.setdefault(samples.patients.get(s, s), []).append(s)

    wbc_cache: dict[str, set] = {}
    coverage_cache: dict[str, set[str] | None] = {}

    for i, r in rows.items():
        s = sample_of[i]
        if samples.types.get(s) != specimen.PLASMA:
            continue
        key, rule = key_of[i], rule_of[i]
        gene = str(r.get("Hugo_Symbol", "")).strip()

        if s in samples.wbc:
            if s not in wbc_cache:
                wbc_cache[s] = read_wbc_variants(samples.wbc[s])
            if key in wbc_cache[s]:
                out.loc[i] = [CONFIRMED, f"{gene}: present in matched WBC"]
            elif rule:
                out.loc[i] = ["", f"{gene}: not in matched WBC"]
            continue
        if rule is None:
            continue

        reasons = [f"{gene} (CH tier {rule.tier})"]
        evidence = 0
        patient = samples.patients.get(s, s)
        others = [o for o in by_patient.get(patient, []) if o != s]

        tissue = [o for o in others if samples.types.get(o) == specimen.TISSUE]
        in_tissue = [o for o in tissue if key in variants.get(o, {})]
        if in_tissue:
            out.loc[i] = [TUMOR, f"{gene}: also in tissue {', '.join(in_tissue)}"]
            continue
        for o in tissue:
            if o not in coverage_cache:
                coverage_cache[o] = _covered_genes(folder, o, specimen.TISSUE)
            covered = coverage_cache[o]
            if covered is None:
                reasons.append(f"absent in tissue {o} (coverage unknown)")
            elif gene in covered:
                reasons.append(f"absent in tissue {o} (covered)")
                evidence += 1
            else:
                reasons.append(f"not covered by tissue {o}'s panel")

        top = tumor_max.get(s)
        if top and not pd.isna(vafs[i]):
            reasons.append(f"VAF {_pct(float(vafs[i]))} vs sample max {_pct(top)}")
            if float(vafs[i]) / top < relative_max:
                evidence += 1

        plasma = [o for o in others if samples.types.get(o) == specimen.PLASMA
                  and key in variants.get(o, {})]
        series = [x for x in [s, *plasma] if tumor_max.get(x) and not pd.isna(variants[x][key])]
        if len(series) >= 2:
            v = [float(variants[x][key]) for x in series]
            t = [tumor_max[x] for x in series]
            if min(v) > 0 and max(v) / min(v) < serial_fold and max(t) / min(t) >= serial_fold:
                reasons.append(f"stable VAF across {len(series)} plasma samples while tumor VAF changed")
                evidence += 1

        label = SUSPECTED if rule.tier == 1 or evidence else COMPATIBLE
        out.loc[i] = [label, "; ".join(reasons)]
    return out


# --------------------------------------------------------------- study

def flag_study(folder: str | Path, sample_table: pd.DataFrame | None = None) -> None:
    """Write CH.EVIDENCE / CH.REASON into the study's data_mutations_
    extended.txt. Specimen types come from `sample_table` (sample.tsv, at
    creation) or from data_clinical_sample.txt (update). A variant already
    "CH confirmed (WBC)" keeps it when its WBC VCF isn't given again."""
    folder = Path(folder)
    maf_path = folder / MAF_FILE
    if not maf_path.exists():
        return
    from filter_clinvar import read_maf
    maf = read_maf(maf_path)

    if sample_table is not None:
        types = specimen.specimen_types(sample_table)
    else:
        types = specimen.study_types(folder)
        sample_table = specimen.read_study_samples(folder)
    previous = maf[[EVIDENCE, REASON]].copy() if {EVIDENCE, REASON} <= set(maf.columns) else None
    maf = maf.drop(columns=[c for c in (EVIDENCE, REASON) if c in maf.columns])

    if not enabled() or specimen.PLASMA not in types.values():
        if previous is not None:
            maf.to_csv(maf_path, sep="\t", index=False)
        return

    samples = samples_from_table(sample_table, types)
    flags = classify(maf, samples, folder)
    if previous is not None:
        keep = (previous[EVIDENCE] == CONFIRMED) & ~maf["Tumor_Sample_Barcode"].isin(samples.wbc)
        flags.loc[keep] = previous.loc[keep].to_numpy()
    maf[EVIDENCE] = flags[EVIDENCE]
    maf[REASON] = flags[REASON]
    maf.to_csv(maf_path, sep="\t", index=False)

    counts = flags[EVIDENCE].value_counts()
    summary = ", ".join(f"{counts[k]} {k}" for k in (CONFIRMED, SUSPECTED, COMPATIBLE, TUMOR)
                        if k in counts)
    logger.info(f"Clonal hematopoiesis flags on plasma variants: {summary or 'none'}. {DISCLAIMER}")


def write_sample_attributes(folder: str | Path) -> None:
    """N_CH_SUSPECT / CH_SUSPECT_GENES for every plasma sample ("NA" for
    the others), from data_mutations_extended.txt's CH.EVIDENCE. Dropped
    when the study has no flagged-variant columns or no plasma sample."""
    folder = Path(folder)
    maf_path = folder / MAF_FILE
    types = specimen.study_types(folder)
    maf = None
    if maf_path.exists():
        from filter_clinvar import read_maf
        maf = read_maf(maf_path)
    if maf is None or EVIDENCE not in maf.columns or specimen.PLASMA not in types.values():
        drop_sample_columns(folder, [N_COLUMN, GENES_COLUMN])
        return
    flagged = maf[maf[EVIDENCE].isin(FLAGGED)]
    counts, genes = {}, {}
    for sample_id, kind in types.items():
        if kind != specimen.PLASMA:
            continue
        rows = flagged[flagged["Tumor_Sample_Barcode"].astype(str).str.strip() == sample_id]
        counts[sample_id] = len(rows)
        genes[sample_id] = ",".join(sorted(set(rows["Hugo_Symbol"].astype(str)))) or "None"
    set_sample_column(
        folder, N_COLUMN, counts, "CH Suspected Variants",
        "Plasma variants flagged as suspected or confirmed clonal hematopoiesis",
        datatype="NUMBER")
    set_sample_column(
        folder, GENES_COLUMN, genes, "CH Suspected Genes",
        "Genes of the plasma variants flagged as suspected or confirmed clonal hematopoiesis")


def has_namespace(folder: str | Path) -> bool:
    """Whether data_mutations_extended.txt carries the CH namespace columns."""
    path = Path(folder) / MAF_FILE
    if not path.exists():
        return False
    with path.open(encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                return EVIDENCE in line.rstrip("\n").split("\t")
    return False


def drop_flagged(maf: pd.DataFrame, sample_id: str, kind: str | None,
                 wbc: Path | None = None) -> pd.DataFrame:
    """One sample's MAF without its suspected/confirmed CH variants - the
    evidence available for a single sample (WBC, gene, relative VAF)."""
    if kind != specimen.PLASMA or maf.empty:
        return maf
    tagged = maf.copy()
    tagged["Tumor_Sample_Barcode"] = sample_id
    samples = Samples({sample_id: kind}, {sample_id: sample_id},
                      {sample_id: wbc} if wbc else {})
    flags = classify(tagged, samples)
    keep = ~flags[EVIDENCE].isin(FLAGGED)
    if (~keep).any():
        logger.info(
            f"Sample {sample_id}: {int((~keep).sum())} suspected/confirmed CH variant(s) "
            "left out of SigMA's input ([CH] EXCLUDE_FROM_DOWNSTREAM).")
    return maf[keep]
