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

NOT YET FULLY VALIDATED - see SIGMA_INTEGRATION_FEASIBILITY.md at the repo
root for the complete verification writeup (round 2, 2026-08-02). Updated
after actually reading SigMA's R source (R/run.R, rdrr.io's rendered
run() man page) instead of relying on the round-1 draft's paraphrase.
oncotree.mskcc.org's API still returns 403 to automated fetches (tried
again, several endpoints/headers) - the child-level OncoTree codes below
(BRCANOS, HGSOC, PRSCC, UTUC, TSTAD, DSTAD, UMEC, ESCC, PANET, CHDM, MDB,
etc.) are still NOT independently confirmed against a live OncoTree query
and remain a human-verification TODO.

CORRECTED: SigMA's `tumor_type` accepts more values than round 1 claimed.
The official run() documentation (https://rdrr.io/github/parklab/SigMA/man/run.html)
lists these 17: bladder, bone_other, breast, crc, eso, gbm, lung, lymph,
medullo, osteo, ovary, panc_ad, panc_en, prost, stomach, thy, uterus.
Note "ewing" and "other" are each individually confirmed elsewhere in the
R source (see below) but do NOT appear in that rendered man-page list -
i.e. this is an inconsistency inside SigMA's own documentation, not
something introduced here.

IMPORTANT, verified straight from R/run.R's do_mva model-availability
check - this materially affects Varan, whose CNA/SNV pipeline is
panel-based (conf.ini [SigMA] DATA_PLATFORM = msk) rather than WES/WGS:

    if(do_mva & !custom & sum(tumor_type == names(gbm_models[[data]])) == 0){
      stop('No built-in MVA models for the tumor_type selected for
            targetted gene panels for "medullo" or "ewing" whole
            exome sequencing is available for others set do_mva to FALSE')
    }

Read plainly, this says: "medullo" and "ewing" only have trained MVA
(do_mva=TRUE) classifiers for whole-exome/WGS data, NOT for gene-panel
data ("msk", Varan's platform). Calling SigMA with tumor_type="medullo" or
"ewing" and data="msk" with do_mva=True (Varan's intended clinical use)
would very likely raise this stop() inside SigMA itself. MDB/ES are left
mapped below (they ARE valid tumor_type strings), but this is a SigMA-
internal limitation an end user has no way to know about - not a clinical
judgment call to expose as a toggle. get_sigma_call_params() therefore
forces do_mva=False for these two automatically and logs a warning,
rather than leaving the decision to whoever calls SigMA.

The panel-platform ("msk") model list is only explicitly named in that
same stop() message text as 10 values: eso, osteo, ovary, panc_ad,
panc_en, prost, stomach, uterus, breast, bladder - which is exactly
ONCOTREE_TO_SIGMA's non-bone/non-ewing/non-medullo coverage below. Common
oncology-panel tumor types with SigMA models for OTHER platforms - crc
(colorectal), gbm (glioblastoma), lung, lymph (lymphoma), thy (thyroid) -
are deliberately NOT added to ONCOTREE_TO_SIGMA: nothing found in the
public source confirms these have a panel/"msk" MVA model too (as opposed
to WES/WGS-only), so guessing them in would risk the same stop()-crash
class of bug as medullo/ewing. Confirming this needs a live R session
(`names(gbm_models[["msk"]])` after `library(SigMA)`), not another
doc-reading pass - flagged as an open question in
SIGMA_INTEGRATION_FEASIBILITY.md.

Also confirmed while researching this (SigMA "Parameter choices" wiki page,
and independently corroborated by run.R's `data` argument documentation on
rdrr.io): `data="msk"` is SigMA's own documented starting point "for
training models on larger gene panels exceeding 300 genes" - not just a
loose size analogy to TSO500 (523 genes), but SigMA's own stated use case.
"""

from __future__ import annotations

from loguru import logger

# DRAFT - verify every code against a live oncotree.mskcc.org lookup before
# trusting this in a clinical run. Parent/tissue-level codes are included
# alongside a few common child codes per lineage, not an exhaustive list.
#
# DO_MVA_UNSAFE_FOR_PANEL_DATA: SigMA tumor_type values whose MVA (do_mva=
# True) classifier is verified (R/run.R's gbm_models[[data]] check, see
# module docstring) to exist only for WES/WGS data, not for panel data
# (Varan's conf.ini DATA_PLATFORM = "msk"). Calling SigMA with one of these
# tumor_type values, data="msk", and do_mva=True is expected to raise
# SigMA's own stop() error. Any future caller must either set do_mva=False
# for samples mapped to one of these, route them to "other" instead, or
# skip SigMA for them outright - not decided here.
DO_MVA_UNSAFE_FOR_PANEL_DATA: frozenset[str] = frozenset({"medullo", "ewing"})

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
    # Ewing sarcoma - see DO_MVA_UNSAFE_FOR_PANEL_DATA: no confirmed panel
    # ("msk") MVA model, WES/WGS only per SigMA's own source.
    "ES": "ewing",
    # Medulloblastoma - see DO_MVA_UNSAFE_FOR_PANEL_DATA: same caveat as ES.
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


def get_sigma_call_params(
    oncotree_code: str,
    fallback_to_other: bool = True) -> tuple[str | None, bool]:
    """Resolve both the SigMA `tumor_type` and whether `do_mva` is safe.

    Varan always runs SigMA against panel data (conf.ini [SigMA]
    DATA_PLATFORM = "msk"). For tumor_type values in
    DO_MVA_UNSAFE_FOR_PANEL_DATA ("medullo", "ewing"), SigMA's own source
    has no trained do_mva=True classifier for panel data - only WES/WGS -
    and calling it that way is expected to raise SigMA's own stop() error
    (see module docstring). This isn't a choice to hand to whoever calls
    SigMA - they have no way to know about this SigMA-internal limitation
    - so do_mva is forced off here and a warning is logged, rather than
    exposing a toggle for it.

    Args:
        oncotree_code (str): The sample's ONCOTREE_CODE, any case.
        fallback_to_other (bool): See get_sigma_tumor_type().

    Returns:
        tuple[str | None, bool]: (tumor_type, do_mva) to pass straight
            into SigMA's run(). tumor_type is None only if
            fallback_to_other is False and the code is unmapped.

    """
    tumor_type = get_sigma_tumor_type(oncotree_code, fallback_to_other)
    if tumor_type in DO_MVA_UNSAFE_FOR_PANEL_DATA:
        logger.warning(
            f"SigMA tumor_type '{tumor_type}' (from OncoTree code "
            f"'{oncotree_code}') has no trained MVA classifier for panel "
            "data - forcing do_mva=False for this sample (its do_mva=True "
            "model is WES/WGS-only, not applicable to Varan's panel-based "
            "pipeline).")
        return tumor_type, False
    return tumor_type, True
