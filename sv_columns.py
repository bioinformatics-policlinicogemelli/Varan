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

"""Column names of data_sv.txt (cBioPortal structural variant format).

The number of tumor reads supporting a fusion (TSO500 "Fusion Supporting
Reads") is written to Tumor_Variant_Count. Varan up to 2.1.3 wrote it to
Normal_Paired_End_Read_Count, which in the cBioPortal format means reads of
the *normal* sample; studies built with those versions are converted on
read (normalise_read_count_column), so old and new studies can be merged.
"""

from __future__ import annotations

import pandas as pd

READ_COUNT = "Tumor_Variant_Count"
LEGACY_READ_COUNT = "Normal_Paired_End_Read_Count"

# Columns written for fusions read from CombinedVariantOutput. Breakpoint
# chromosome/position are part of cBioPortal's uniqueness key for SV rows:
# without them, distinct breakpoints of the same gene pair (e.g. different
# EML4-ALK variants) look like duplicates.
FUSION_COLUMNS = [
    "Sample_Id", "SV_Status", "Class",
    "Site1_Hugo_Symbol", "Site1_Chromosome", "Site1_Position",
    "Site2_Hugo_Symbol", "Site2_Chromosome", "Site2_Position",
    READ_COUNT, "Event_Info", "RNA_Support"]


def normalise_read_count_column(df: pd.DataFrame) -> pd.DataFrame:
    """Rename the pre-2.1.4 read-count column to Tumor_Variant_Count.

    Only when the table has the legacy column and not the new one, so a
    table that already carries both is left untouched.

    Args:
        df (pd.DataFrame): A data_sv table (or a fusion input table).

    Returns:
        pd.DataFrame: The table with the read-count column renamed.

    """
    if LEGACY_READ_COUNT in df.columns and READ_COUNT not in df.columns:
        return df.rename(columns={LEGACY_READ_COUNT: READ_COUNT})
    return df


def normalise_header(columns: list[str]) -> list[str]:
    """Same as normalise_read_count_column, for a list of column names."""
    if LEGACY_READ_COUNT in columns and READ_COUNT not in columns:
        return [READ_COUNT if c == LEGACY_READ_COUNT else c for c in columns]
    return list(columns)
