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
Reads") is written to Tumor_Variant_Count, and the minimum read count is
applied to that column only.

Verified on DRAGEN TSO 500 v2.5 output: "Fusion Supporting Reads" in
CombinedVariantOutput.tsv equals Alt Split Dedup + Alt Pair Dedup of the
same fusion in the sample's *_AllFusions.csv (e.g. EML4-ALK 487 + 279 = 766),
i.e. all deduplicated tumor RNA reads supporting the fusion, split and
paired together. That is cBioPortal's Tumor_Variant_Count ("number of reads
of the tumor tissue that have the variant"); Tumor_Split_Read_Count or
Tumor_Paired_End_Read_Count alone would each be only part of it.

Varan up to 2.1.3 wrote that number to Normal_Paired_End_Read_Count, which
in the cBioPortal format means reads of the *normal* sample. Existing
tables are never renamed or rewritten: when an old study is merged with a
new one, both columns are kept, each with its own values.
"""

from __future__ import annotations

READ_COUNT = "Tumor_Variant_Count"
LEGACY_READ_COUNT = "Normal_Paired_End_Read_Count"

# Read-count columns: rows that differ only in these are the same event.
READ_COUNT_COLUMNS = [READ_COUNT, LEGACY_READ_COUNT]

# Columns written for fusions read from CombinedVariantOutput. Breakpoint
# chromosome/position are part of cBioPortal's uniqueness key for SV rows:
# without them, distinct breakpoints of the same gene pair (e.g. different
# EML4-ALK variants) look like duplicates.
FUSION_COLUMNS = [
    "Sample_Id", "SV_Status", "Class",
    "Site1_Hugo_Symbol", "Site1_Chromosome", "Site1_Position",
    "Site2_Hugo_Symbol", "Site2_Chromosome", "Site2_Position",
    READ_COUNT, "Event_Info", "RNA_Support"]
