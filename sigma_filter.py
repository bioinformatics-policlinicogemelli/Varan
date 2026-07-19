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

"""Build a MAF suitable for SigMA mutational-signature analysis.

Why this is not just filter_main() with a different flag combination:
SigMA's own make_matrix() (see the SigMA research notes in
project_future_implementations.md) filters to SNVs only, but does NOT
filter by VAF and does NOT distinguish germline from somatic - it trusts
the input is already clean. Varan's *clinical* MAF pipeline goes the other
direction: filter_main()'s 'o'/'i' flags narrow down to oncogenic/
high-impact variants only, which would leave far too few passenger
mutations for a trinucleotide-context spectrum to be meaningful. SigMA
needs a third tier: quality-filtered (PASS + population AF + VAF exclude
bands) but NOT clinically narrowed.

This reuses filter_clinvar.filter_vaf_exclude_bands() (the new 'g' filter
primitive) with SigMA's own conf.ini values ([SigMA] section - distinct
from [Filters]' pancancer VAF thresholds, per explicit request) rather than
going through filter_main()'s filters-string dispatch, since the exact
combination SigMA needs (PASS + population AF + VAF bands, explicitly
never 'o'/'i') isn't expressible as a single filters string without also
building a separate MAF_OncoKB-shaped output tree filter_main() assumes.

NOT YET WIRED into walk_folder()/Snakemake - this is the filtering
function only. Where it's called from (per-sample after vcf2maf, or once
on the concatenated unfiltered MAF) and where its output lands
(intermediate/sigma/, per the agreed design) is the next piece of this
branch's work.
"""

from __future__ import annotations

import ast

import pandas as pd
from loguru import logger

from config_loader import get_config
from filter_clinvar import check_bool, filter_vaf_exclude_bands

config = get_config()


def prepare_sigma_maf(maf_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """Filter one sample's unfiltered MAF down to what SigMA should see.

    Applies, in order:
    1. FILTER == "PASS" (if a FILTER column is present - vcf2maf output
       usually already is PASS-only, but don't assume it here).
    2. Population AF filter (germline removal) - reuses the same [Filters]
       AF/drop_NA_AF conf.ini values as filter_main()'s 'a' flag, since
       population-frequency-based germline filtering doesn't need
       SigMA-specific tuning the way VAF does.
    3. [SigMA] VAF_MIN (a simple lower bound, separate from [Filters]'
       t_VAF_min which is tuned for the general pancancer clinical MAF).
    4. [SigMA] VAF_EXCLUDE_BANDS via filter_vaf_exclude_bands() - the
       germline heterozygous/homozygous VAF-cluster scrub.

    Deliberately never applies filter_main()'s 'o' (oncogenic-only) or 'i'
    (impact-based) filters - those would drop the passenger mutations
    SigMA's trinucleotide-context spectrum actually needs.

    Args:
        maf_df (pd.DataFrame): One sample's unfiltered, annotated MAF
            (same file `filter_main()` reads from the `maf/` folder before
            any clinical filtering).
        sample_id (str): For logging only.

    Returns:
        pd.DataFrame: The SigMA-ready subset of maf_df.

    """
    n_start = len(maf_df)

    if "FILTER" in maf_df.columns:
        maf_df = maf_df[maf_df["FILTER"] == "PASS"]

    if "AF" in maf_df.columns:
        af = config.get("Filters", "AF")
        drop_na = check_bool(config.get("Filters", "drop_NA_AF"))
        af_numeric = pd.to_numeric(maf_df["AF"], errors="coerce")
        na_mask = af_numeric.isna()
        passes_af = eval(f"af_numeric {af}")  # noqa: S307 - same pattern as filter_main's 'a' flag
        keep_mask = passes_af.fillna(False) | (na_mask & (not drop_na))
        maf_df = maf_df[keep_mask]
    else:
        logger.warning(
            f"Sample {sample_id}: no AF column found - skipping the "
            "population-frequency germline filter for the SigMA MAF.")

    vaf_colname = ("t_AF"
        if "t_AF" in maf_df.columns and maf_df["t_AF"].notna().any()
        else "t_VF")
    if vaf_colname in maf_df.columns:
        vaf_min = float(config.get("SigMA", "VAF_MIN"))
        vaf_numeric = pd.to_numeric(maf_df[vaf_colname], errors="coerce")
        maf_df = maf_df[vaf_numeric.fillna(-1) >= vaf_min]

        exclude_bands = ast.literal_eval(config.get("SigMA", "VAF_EXCLUDE_BANDS"))
        maf_df = filter_vaf_exclude_bands(maf_df, vaf_colname, exclude_bands)
    else:
        logger.warning(
            f"Sample {sample_id}: neither t_AF nor t_VF column found - "
            "skipping the SigMA-specific VAF filters (min + exclude bands).")

    logger.info(
        f"Sample {sample_id}: SigMA MAF filtering kept {len(maf_df)}/{n_start} "
        "variants (PASS + population AF + VAF exclude bands, no "
        "oncogenic/impact restriction).")

    return maf_df
