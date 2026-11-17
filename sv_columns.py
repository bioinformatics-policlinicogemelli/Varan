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
the *normal* sample.

When studies are merged, a pre-2.1.4 data_sv.txt is converted only if it
is recognisably written by Varan (VARAN_LEGACY_HEADER): in any other table
Normal_Paired_End_Read_Count may really hold normal-sample reads, so it is
left untouched rather than relabelled.
"""

from __future__ import annotations

import pandas as pd
from loguru import logger

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

# The first columns of every data_sv.txt Varan wrote before 2.1.4 (same
# layout from 2024-03 to 2.1.3; OncoKB columns, if any, follow).
VARAN_LEGACY_HEADER = [
    "Sample_Id", "SV_Status", "Class",
    "Site1_Hugo_Symbol", "Site2_Hugo_Symbol",
    LEGACY_READ_COUNT, "Event_Info", "RNA_Support"]


def is_varan_legacy(columns: list[str]) -> bool:
    """Tell whether a data_sv header is the one written by Varan < 2.1.4."""
    return list(columns)[:len(VARAN_LEGACY_HEADER)] == VARAN_LEGACY_HEADER


def normalise_study_read_count(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Rename the read-count column of a pre-2.1.4 Varan data_sv.txt.

    Only for tables written by Varan (VARAN_LEGACY_HEADER), where the column
    always held tumor supporting reads. Any other table with
    Normal_Paired_End_Read_Count is left untouched, with a warning.

    Args:
        df (pd.DataFrame): A data_sv table read from an existing study.
        source (str): Where the table comes from, for the log.

    Returns:
        pd.DataFrame: The table, with the column renamed if applicable.

    """
    if LEGACY_READ_COUNT not in df.columns or READ_COUNT in df.columns:
        return df
    if is_varan_legacy(df.columns):
        logger.info(f"{source}: data_sv.txt written by Varan < 2.1.4, read "
                    f"count column {LEGACY_READ_COUNT} renamed to {READ_COUNT}")
        return df.rename(columns={LEGACY_READ_COUNT: READ_COUNT})
    logger.warning(f"{source}: data_sv.txt has {LEGACY_READ_COUNT} but was not "
                   f"written by Varan: the column is kept as it is (it may hold "
                   f"normal-sample reads), not converted to {READ_COUNT}")
    return df


def normalise_input_read_count(columns: list[str]) -> list[str]:
    """Rename Normal_Paired_End_Read_Count in a fusion INPUT table header.

    Varan has always used that column of its fusion input tables as the
    tumor supporting-read count (it applied the minimum read count to it),
    so input tables keep being accepted with the old name.

    Args:
        columns (list[str]): Header of the fusion input table.

    Returns:
        list[str]: The header with the read-count column renamed if needed.

    """
    if LEGACY_READ_COUNT in columns and READ_COUNT not in columns:
        logger.info(f"Fusion input table: {LEGACY_READ_COUNT} read as "
                    f"{READ_COUNT} (tumor supporting reads)")
        return [READ_COUNT if c == LEGACY_READ_COUNT else c for c in columns]
    return list(columns)
