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

"""DRAFT mapping from OncoTree code to SigMA's `tumor_type` parameter.

NOT YET VALIDATED - this is a first pass built from general OncoTree
knowledge, not cross-checked against a live query to oncotree.mskcc.org
(its API returned 403 to automated fetches while writing this). Every
entry needs a human check against https://oncotree.mskcc.org before this
is trusted in a clinical run - that is the explicit next step, not
something to skip.

SigMA's `tumor_type` parameter only accepts a fixed list of ~13 specific
values plus a generic "other" model (source: SigMA R/run.R and the
"Parameter choices" wiki page, both read 2026-07-19):
    eso, osteo, ovary, panc_ad, panc_en, prost, stomach, uterus, breast,
    bladder, bone_other, medullo, ewing, other

"other" is a real, usable pan-cancer fallback model (not a "skip" sentinel)
- any OncoTree code not covered by ONCOTREE_TO_SIGMA below still gets run
  through SigMA with tumor_type="other" rather than being excluded outright.
Whether that's the right call clinically (vs. only running the ~13
type-specific models and skipping everything else) is a decision to make
together, not baked into this file - see get_sigma_tumor_type()'s
`fallback_to_other` argument.

Also confirmed while researching this (SigMA "Parameter choices" wiki page):
`data="msk"` is documented as the intended baseline specifically "for
training models on larger gene panels exceeding 300 genes" - not just a
loose size analogy to TSO500 (523 genes), but SigMA's own stated use case.
"""

from __future__ import annotations

# DRAFT - verify every code against a live oncotree.mskcc.org lookup before
# trusting this in a clinical run. Parent/tissue-level codes are included
# alongside a few common child codes per lineage, not an exhaustive list.
ONCOTREE_TO_SIGMA: dict[str, str] = {
    # Breast
    "BREAST": "breast",
    "BRCA": "breast",
    "IDC": "breast",
    "ILC": "breast",
    "BRCANOS": "breast",
    # Ovary
    "OVARY": "ovary",
    "HGSOC": "ovary",
    "OCS": "ovary",
    "EOV": "ovary",
    # Prostate
    "PROSTATE": "prost",
    "PRAD": "prost",
    "PRSCC": "prost",
    # Bladder / urothelial
    "BLADDER": "bladder",
    "BLCA": "bladder",
    "UTUC": "bladder",
    # Stomach
    "STOMACH": "stomach",
    "STAD": "stomach",
    "TSTAD": "stomach",
    "DSTAD": "stomach",
    # Uterus
    "UTERUS": "uterus",
    "UCEC": "uterus",
    "UCS": "uterus",
    "UMEC": "uterus",
    # Esophagus
    "ESOPHAGUS": "eso",
    "ESCA": "eso",
    "ESCC": "eso",
    # Pancreas - adenocarcinoma vs. neuroendocrine are two different SigMA models
    "PAAD": "panc_ad",
    "PANET": "panc_en",
    # Bone - osteosarcoma has its own model, other bone sarcomas share "bone_other"
    "OS": "osteo",
    "CHS": "bone_other",
    "CHDM": "bone_other",
    "BONE": "bone_other",
    # Ewing sarcoma
    "ES": "ewing",
    # Medulloblastoma
    "MDB": "medullo",
}


def get_sigma_tumor_type(
    oncotree_code: str,
    fallback_to_other: bool = True) -> str | None:
    """Resolve an OncoTree code to a SigMA `tumor_type` value.

    Args:
        oncotree_code (str): The sample's ONCOTREE_CODE, any case.
        fallback_to_other (bool): If True (default), an unmapped code
            resolves to SigMA's generic "other" model rather than being
            excluded. See the module docstring - this is a real, usable
            model, not a null result, but whether to use it by default is
            a clinical judgment call, not something this function decides
            on its own; set False to only run the ~13 type-specific models
            and skip everything else.

    Returns:
        str | None: A valid SigMA `tumor_type` value, or None if the code
            is unmapped and fallback_to_other is False.

    """
    mapped = ONCOTREE_TO_SIGMA.get(oncotree_code.strip().upper())
    if mapped is not None:
        return mapped
    return "other" if fallback_to_other else None
