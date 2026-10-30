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

"""Mapping from OncoTree code to SigMA's `tumor_type` parameter.

VERIFIED 2026-10-30 (round 3): every code in ONCOTREE_TO_SIGMA checked
against the live OncoTree tree (oncotree.mskcc.org/api/tumorTypes - it
answered this time), and the MVA models per platform read from SigMA 1.6's
own gbm_models object in an R session (MVA_MODELS_BY_PLATFORM). Outcome:
MDB (no such OncoTree code) replaced by MBL and its subtypes, ESOPHAGUS
(no such code) dropped, colorectal codes added - SigMA does have a panel
("msk") MVA model for crc, unlike what round 2 could only assume - and
the do_mva decision now follows the configured DATA_PLATFORM. Everything
below this paragraph is the earlier rounds' history.

Round 2 notes - see SIGMA_INTEGRATION_FEASIBILITY.md at the repo root for
the complete verification writeup (round 2, 2026-08-02). Updated
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

# SigMA's trained MVA (do_mva=True) models, per `data` platform - read
# from SigMA 1.6's own gbm_models object (names(gbm_models[[data]]) in an
# R session, 2026-10-30). R/run.R stops with an error when do_mva=TRUE and
# tumor_type isn't among its platform's models (this includes the "other"
# fallback, bone_other, medullo and ewing on panel data), so do_mva is
# only turned on for these.
MVA_MODELS_BY_PLATFORM: dict[str, frozenset[str]] = {
    "msk": frozenset({
        "bladder", "breast", "crc", "eso", "osteo", "ovary", "panc_ad",
        "panc_en", "prost", "stomach", "uterus"}),
    "seqcap": frozenset({
        "bladder", "breast", "eso", "ewing", "medullo", "osteo", "ovary",
        "panc_ad", "panc_en", "prost", "stomach", "uterus"}),
    "seqcap_probe": frozenset({
        "bladder", "breast", "eso", "ewing", "medullo", "osteo", "ovary",
        "panc_ad", "panc_en", "prost", "stomach", "uterus"}),
    "tcga_mc3": frozenset({
        "adrenal", "biliary", "bladder", "breast", "cervix", "crc", "eso",
        "glioma", "head", "kidney", "lgg", "liver", "lung", "lymph", "meso",
        "ovary", "panc_ad", "prost", "skin", "stomach", "uterus"}),
    "wgs": frozenset({
        "biliary", "bladder", "bone_other", "breast", "cervix", "crc", "eso",
        "ewing", "kidney", "lms", "lung", "medullo", "osteo", "ovary",
        "panc_ad", "panc_en", "prost", "skin", "stomach", "uterus"}),
    "wgs_pancan": frozenset({"breast", "generic"}),
}

# Kept for callers that still import it: Varan's default platform's set.
PANEL_DATA_MVA_SAFE_TUMOR_TYPES = MVA_MODELS_BY_PLATFORM["msk"]


def is_mva_safe_for_panel_data(tumor_type: str, data_platform: str = "msk") -> bool:
    """Whether SigMA has a trained do_mva=True model for `tumor_type` on
    `data_platform` (conf.ini [SigMA] DATA_PLATFORM) - see
    MVA_MODELS_BY_PLATFORM. Use it for every resolved tumor_type, including
    the "other" fallback and a manual TUMOR_TYPE_OVERRIDE.

    Args:
        tumor_type (str): A SigMA `tumor_type` value.
        data_platform (str): SigMA's `data` argument (default "msk").

    Returns:
        bool: True if do_mva=True is safe for this tumor_type and platform.

    """
    return tumor_type in MVA_MODELS_BY_PLATFORM.get(data_platform, frozenset())


# OncoTree code -> SigMA tumor_type. Every code checked against the live
# OncoTree tree (oncotree.mskcc.org/api/tumorTypes, 2026-10-30); parent/
# tissue-level codes are included with their main child codes - not an
# exhaustive list of every OncoTree descendant.
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
    # Esophagus (OncoTree has no ESOPHAGUS code: its tissue node is STOMACH,
    # "Esophagus/Stomach", mapped to stomach above)
    "ESCA": "eso",
    "ESCC": "eso",
    # Colorectal - has a panel ("msk") MVA model
    "COADREAD": "crc",
    "COAD": "crc",
    "READ": "crc",
    "MACR": "crc",
    "SRCCR": "crc",
    # Pancreas - adenocarcinoma vs. neuroendocrine are two different SigMA models
    "PAAD": "panc_ad",
    "PANET": "panc_en",
    # Bone - osteosarcoma has its own model, other bone sarcomas share "bone_other"
    "OS": "osteo",
    "CHS": "bone_other",
    "CHDM": "bone_other",
    "BONE": "bone_other",
    # Ewing sarcoma - no panel ("msk") MVA model (seqcap/WGS only)
    "ES": "ewing",
    "ESST": "ewing",
    # Medulloblastoma - no panel ("msk") MVA model either. OncoTree's code is
    # MBL (the earlier MDB doesn't exist); its subtype codes map the same.
    "MBL": "medullo",
    "MBLNOS": "medullo",
    "MBLWNT": "medullo",
    "MBLSHH": "medullo",
    "MBLNWS": "medullo",
    "MBLG3": "medullo",
    "MBLG4": "medullo",
    "DMBL": "medullo",
    "DMBLNOS": "medullo",
    "DMBLSHH": "medullo",
    "AMBL": "medullo",
    "AMBLNOS": "medullo",
    "AMBLSHH": "medullo",
    "AMBLNWS": "medullo",
    "AMBLNWSG3": "medullo",
    "AMBLNWSG4": "medullo",
    "MBEN": "medullo",
    "MBENNOS": "medullo",
    "MBENSHH": "medullo",
    "MMBL": "medullo",
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
    fallback_to_other: bool = True,
    data_platform: str = "msk") -> tuple[str | None, bool]:
    """Resolve an OncoTree code to the (tumor_type, do_mva) pair walk.py
    needs to actually invoke SigMA for one sample.

    This is the single entry point walk.py's per-sample SigMA wiring
    should call - it wraps get_sigma_tumor_type() and additionally
    resolves do_mva via is_mva_safe_for_panel_data(), so a caller never
    has to remember to check tumor_type against the panel-safe whitelist
    itself (including the "other" fallback - see
    PANEL_DATA_MVA_SAFE_TUMOR_TYPES's docstring for why a plain
    medullo/ewing blacklist is not enough).

    Args:
        oncotree_code (str): The sample's ONCOTREE_CODE, any case.
        fallback_to_other (bool): Forwarded to get_sigma_tumor_type() -
            see there. Normally sourced from conf.ini's [SigMA]
            FALLBACK_TO_OTHER.

    Returns:
        tuple[str | None, bool]: (tumor_type, do_mva).
            tumor_type is None only when the code is unmapped and
            fallback_to_other is False - callers must treat that as
            "skip SigMA for this sample", not call SigMA with
            tumor_type=None.
            do_mva is always resolved for Varan's panel-data
            (DATA_PLATFORM = msk) use: True only for the confirmed
            panel-safe tumor_type values (see
            PANEL_DATA_MVA_SAFE_TUMOR_TYPES), False (with a logged
            warning) for everything else, including "other" and
            medullo/ewing. The raw Signature_3_ml/exp_sig3/rat_sig3
            features are still computed either way - only the
            GBM-combined Signature_3_mva score and pass_mva/
            pass_mva_strict calls are skipped when do_mva is False.

    """
    tumor_type = get_sigma_tumor_type(oncotree_code, fallback_to_other)
    if tumor_type is None:
        return None, False

    do_mva = is_mva_safe_for_panel_data(tumor_type, data_platform)
    if not do_mva:
        logger.warning(
            f"OncoTree code '{oncotree_code}' maps to SigMA tumor_type "
            f"'{tumor_type}', which has no trained MVA (do_mva) classifier "
            f"for DATA_PLATFORM={data_platform} - forcing do_mva=False for this sample. Only the "
            "raw Signature_3 likelihood/exposure features will be "
            "computed, not the Signature_3_mva score or pass_mva calls.")

    return tumor_type, do_mva
