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

"""A study's alterations in one shape, for the modules that annotate them
(actionability.py, tumor_agnostic.py, civic.py, escat.py): mutations from
data_mutations_extended.txt, fusions from data_sv.txt, amplifications and
deep deletions (+2/-2) from data_cna.txt.

Also the alteration syntax those modules share (one item, or several
separated by ","):
  V600E          protein change (HGVSp_Short without "p.", any case)
  V600           any change at that codon        C77fs  frameshift at codon 77
  exon 19 deletion / exon 20 insertion (in frame) / exon 18 mutation / exon 7 frameshift
  mutation       any coding non-synonymous change  truncating  missense
  oncogenic      OncoKB Oncogenic / Likely Oncogenic (needs OncoKB annotation)
  fusion         a fusion involving the gene
  amplification  copy number +2     deletion (or loss)  copy number -2
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

import clonal_hematopoiesis as ch

MAF_FILE = "data_mutations_extended.txt"
SV_FILE = "data_sv.txt"
CNA_FILE = "data_cna.txt"

MUTATION = "mutation"
FUSION = "fusion"
AMPLIFICATION = "amplification"
DELETION = "deletion"

NON_SYNONYMOUS = ch.CLASSES["any"]
EXON_CLASSES = {
    # In frame, as EGFR exon 19 deletions and ERBB2/EGFR exon 20 insertions.
    "deletion": {"In_Frame_Del"},
    "insertion": {"In_Frame_Ins"},
    "frameshift": {"Frame_Shift_Del", "Frame_Shift_Ins"},
    "mutation": NON_SYNONYMOUS,
}
ONCOGENIC = {"Oncogenic", "Likely Oncogenic"}


# ------------------------------------------------------------ reading

def read_maf(folder: str | Path) -> pd.DataFrame | None:
    path = Path(folder) / MAF_FILE
    if not path.exists():
        return None
    from filter_clinvar import read_maf as _read
    return _read(path)


def read_sv(folder: str | Path) -> pd.DataFrame | None:
    path = Path(folder) / SV_FILE
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def fusions(sv: pd.DataFrame | None) -> pd.DataFrame:
    """data_sv.txt's fusion rows (splice variants left out)."""
    if sv is None or sv.empty:
        return pd.DataFrame(columns=["Sample_Id", "Site1_Hugo_Symbol", "Site2_Hugo_Symbol"])
    if "Class" in sv.columns:
        sv = sv[sv["Class"].str.upper() != "SPLICE"]
    return sv


def read_cna(folder: str | Path) -> list[dict]:
    """[{sample, gene, value}] for every +2/-2 call in data_cna.txt."""
    path = Path(folder) / CNA_FILE
    if not path.exists():
        return []
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    gene_col = "Hugo_Symbol" if "Hugo_Symbol" in table.columns else table.columns[0]
    calls = []
    for sample in (c for c in table.columns if c not in (gene_col, "Entrez_Gene_Id")):
        for gene, value in zip(table[gene_col], table[sample]):
            if value.strip() in ("2", "-2"):
                calls.append({"sample": sample, "gene": gene, "value": int(value)})
    return calls


def protein(row: dict) -> str:
    value = str(row.get("HGVSp_Short") or "").strip()
    return value[2:] if value.startswith("p.") else value


def mutation_label(row: dict) -> str:
    gene = str(row.get("Hugo_Symbol", "")).strip()
    return f"{gene} {protein(row) or row.get('Variant_Classification', '')}".strip()


def fusion_label(row: dict) -> str:
    a, b = row.get("Site1_Hugo_Symbol", ""), row.get("Site2_Hugo_Symbol", "")
    return f"{a}::{b} fusion" if a != b else f"{a} fusion"


def cna_label(call: dict) -> str:
    return f"{call['gene']} {'amplification' if call['value'] > 0 else 'deletion'}"


# ------------------------------------------------------------- syntax

@dataclass
class Spec:
    kind: str
    exact: str = ""          # protein change, upper case
    codon: str = ""          # e.g. "V600"
    frameshift: bool = False
    classes: set[str] = field(default_factory=set)
    exon: str = ""
    oncogenic: bool = False

    def matches_mutation(self, row: dict) -> bool:
        if self.kind != MUTATION:
            return False
        change = protein(row).upper()
        vc = str(row.get("Variant_Classification", ""))
        if self.exact:
            return change == self.exact
        if self.codon:
            if not re.match(rf"^{re.escape(self.codon)}(?![0-9])", change):
                return False
            return "FS" in change if self.frameshift else True
        if self.oncogenic:
            return str(row.get("ONCOGENIC", "")).strip() in ONCOGENIC
        if self.classes and vc not in self.classes:
            return False
        return not self.exon or ch._exon(row) == self.exon


def parse(text: str) -> list[Spec]:
    """The Specs of an alteration text (see the module docstring);
    ValueError for an unknown item."""
    specs = []
    for item in (i.strip() for i in str(text).split(",")):
        if item:
            spec = parse_one(item)
            if spec is None:
                msg = f"Alteration '{item}' isn't understood (see study_alterations.py)."
                raise ValueError(msg)
            specs.append(spec)
    return specs


def parse_one(item: str) -> Spec | None:
    text = re.sub(r"\s*\(.*\)$", "", item.strip())   # "Q195* (c.583C>T)" -> "Q195*"
    low = text.lower()
    if low in ("fusion", "fusions", "rearrangement"):
        return Spec(FUSION)
    if low in ("amplification", "amp"):
        return Spec(AMPLIFICATION)
    if low in ("deletion", "loss", "homozygous deletion", "deep deletion"):
        return Spec(DELETION)
    if low in ("mutation", "mutations", "any"):
        return Spec(MUTATION, classes=set(NON_SYNONYMOUS))
    if low in ("truncating", "truncating mutation"):
        return Spec(MUTATION, classes=set(ch.TRUNCATING))
    if low in ("missense", "missense mutation"):
        return Spec(MUTATION, classes=set(ch.MISSENSE))
    if low in ("oncogenic", "oncogenic mutation"):
        return Spec(MUTATION, oncogenic=True)
    m = re.fullmatch(r"exon\s*(\d+)\s+(deletion|insertion|frameshift|mutation)s?", low)
    if m:
        return Spec(MUTATION, classes=set(EXON_CLASSES[m.group(2)]), exon=m.group(1))
    up = text.upper()
    m = re.fullmatch(r"([A-Z])(\d+)(FS)?", up)
    if m:
        return Spec(MUTATION, codon=m.group(1) + m.group(2), frameshift=bool(m.group(3)))
    if re.fullmatch(r"[A-Z*]\d+[A-Z*_0-9]*(DEL|INS|DUP|DELINS|FS|\*)?[A-Z*0-9]*", up):
        return Spec(MUTATION, exact=up)
    return None


def fusion_matches(gene: str, row: dict, partner: str = "") -> bool:
    """A data_sv.txt row is a fusion of `gene` (and of `partner`, if given)."""
    genes = {row.get("Site1_Hugo_Symbol", ""), row.get("Site2_Hugo_Symbol", "")}
    return gene in genes and (not partner or partner in genes)
