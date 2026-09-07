"""Vendor-agnostic building blocks for converting a third-party panel-
sequencing vendor's raw run output into the sample.tsv / patient.tsv /
{run_id}_fusions.tsv shape Varan's own `varan.py -i` flag consumes.

This module holds ONLY the pieces that are not specific to any one vendor's
file formats: generic S3 listing helpers, the ONCOTREE dict.csv lookup
(a Varan-side concept, not a vendor one), the shared `SampleRow` result
shape and sample.tsv writer, and the fusion-table writer for Varan's
`FUSIONS/*.tsv` ingestion path (`fill_fusion_from_temp()` in walk.py) -
that path and its column shape belong to Varan, not to any single vendor,
so every vendor adapter should reuse the same writer here rather than
re-implementing it.

Vendor-specific parsing (VCF dialects, metadata XML/JSON shapes, MSI/CNV
report formats, etc.) belongs in that vendor's own module
(e.g. `vendor_adapters/guardant.py`), not here. See
`vendor_adapters/__init__.py` for the adapter registry and the interface
each vendor module is expected to implement, and
`MULTIVENDOR_INTEGRATION_NOTES.md` for the full design rationale.
"""

import csv
import os
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from loguru import logger

# Varan's own data_sv.txt column shape (Sample_Id/SV_Status/Site1_Hugo_Symbol/
# Site2_Hugo_Symbol are the ones fill_fusion_from_temp() in walk.py actually
# requires; the rest are carried through as-is). This is Varan's ingestion
# contract, not a vendor's - shared by every adapter.
FUSION_TABLE_HEADER = ("Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\t"
                       "Site2_Hugo_Symbol\tNormal_Paired_End_Read_Count\t"
                       "Event_Info\tRNA_Support\n")

# Varan's Templates/sample.tsv column order (see Templates/sample.tsv /
# Templates/sample_guardant_example.tsv). Any additional column a vendor
# wants to pass through (e.g. a future "golden" BRCA Reversion or AutoQC
# field) is appended after these, per-row, and left blank for rows/vendors
# that don't populate it - write_clinical_sample() in Update_functions.py
# already carries unknown extra sample.tsv columns straight through into
# data_clinical_sample.txt with no extra code required on the Varan side.
SAMPLE_TSV_COLUMNS = [
    "SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path",
    "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR",
]


@dataclass
class SampleRow:
    """One sample's worth of data destined for a row of sample.tsv.

    `run_id` is bookkeeping only (used to name the output report / fusion
    table file) - it is deliberately NOT one of Varan's own sample.tsv
    columns and is never written out.

    `extra` holds any additional passthrough columns a vendor wants to
    populate (e.g. future BRCA_REVERSION), keyed by column name. Columns
    present in `extra` for at least one row in a batch are added to the
    written sample.tsv header; rows that don't set a given key get a blank
    value for it, matching Varan's own "blank means not-applicable, not
    zero" convention for MSI/TMB.
    """

    sample_id: str
    patient_id: str
    run_id: str
    oncotree_code: str
    snv_path: str = ""
    cnv_path: str = ""
    comb_path: str = ""
    msi: str = ""
    tmb: str = ""
    msi_thr: str = ""
    tmb_thr: str = ""
    extra: Dict[str, str] = field(default_factory=dict)


def run_cmd(cmd: str) -> Optional[str]:
    """Run a shell command, returning stripped stdout on success or None.

    NOTE: kept as `shell=True` string-concatenation on purpose here since
    it only ever runs fixed `aws s3 ls`/`aws s3 cp` invocations built from
    already-validated folder/file names in practice - this is a separate
    (and much lower-risk) code path from the vcf-query shell=True bug fixed
    elsewhere in Varan's own vcf2maf_constructor, which took untrusted
    sample-derived strings straight into a shell string. Still, callers
    should avoid passing untrusted input through `cmd`.
    """
    result = subprocess.run(
        cmd, shell=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.strip() if result.returncode == 0 else None


def list_s3_files(s3_folder: str) -> List[str]:
    """List file names (not full paths) directly under an S3 folder."""
    output = run_cmd(f"aws s3 ls {s3_folder.rstrip('/')}/")
    if not output:
        return []
    return [line.split()[-1] for line in output.splitlines()]


def load_oncotree_dict(dict_path: str) -> Tuple[Dict[str, str], set]:
    """Load a vendor-diagnosis -> ONCOTREE_CODE lookup from a dict.csv with
    `name`/`code` columns.

    This is a Varan-side concept (ONCOTREE codes), not vendor-specific, so
    every vendor adapter that needs a free-text-diagnosis -> ONCOTREE_CODE
    mapping should reuse this rather than re-implementing its own CSV
    lookup. The mapping file's actual content/path is vendor- and
    deployment-specific and is passed in by the caller.
    """
    name_to_code: Dict[str, str] = {}
    valid_codes = set()
    if dict_path and os.path.exists(dict_path):
        # "utf-8-sig" transparently strips a leading BOM if present (e.g.
        # a dict.csv saved from Excel as "CSV UTF-8") and behaves exactly
        # like plain "utf-8" otherwise. Without this, a BOM sticks to the
        # first header name (turning "code" into "﻿code"), silently
        # emptying the whole mapping with no error - every diagnosis would
        # then come back NOT_IN_DICT even though the file "looks" right.
        with open(dict_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "name" in row and "code" in row:
                    name = row["name"].strip().lower()
                    code = row["code"].strip().upper()
                    name_to_code[name] = code
                    valid_codes.add(code)
        if not name_to_code:
            logger.warning(
                f"{dict_path} exists but no usable name/code rows were "
                "read from it - check that its header row is exactly "
                "'name,code' (or 'code,name') and that it isn't empty.")
    return name_to_code, valid_codes


def append_fusions_to_table(fusion_table_path: str, sample_id: str,
                             fusions: Sequence[Dict[str, str]]) -> None:
    """Append this sample's confirmed fusions to the shared run-level
    fusions.tsv file, in Varan's own data_sv.txt column shape.

    `fusions` is a list of dicts with `gene_a`/`gene_b`/`supp` keys
    (supp = supporting read/molecule count, vendor-defined units).

    Fusions are written through this plain `FUSIONS/*.tsv` path (read by
    `fill_fusion_from_temp()` in walk.py) rather than through a synthesized
    CombinedVariantOutput.tsv / comb_path file. That distinction matters:
    `_walk_setup()` in walk.py detects "does this run have CombinedOutput"
    by folder existence, not per-sample - so as soon as ANY sample in a
    batch has a non-empty comb_path, the WHOLE run's MSI/TMB handling
    silently switches from the sample.tsv-driven VALUE+THR precedence path
    (`fill_from_file`) to the native-CombinedOutput path
    (`fill_from_combined`), defeating any vendor's own MSI/TMB columns for
    every sample in the batch, not just the one needing fusions. Originally
    found and fixed for Guardant; the same hazard applies to any future
    vendor, which is why this lives here instead of in guardant.py.
    """
    is_new = not os.path.exists(fusion_table_path)
    with open(fusion_table_path, "a") as f:
        if is_new:
            f.write(FUSION_TABLE_HEADER)
        for fus in fusions:
            f.write(
                f"{sample_id}\tSOMATIC\tFUSION\t{fus['gene_a']}\t"
                f"{fus['gene_b']}\t{fus['supp']}\t"
                f"{fus['gene_a']}-{fus['gene_b']} Fusion\tYes\n")


def get_incremental_report_path(base_dir: str, run_id: str) -> str:
    """Pick `{run_id}_VARAN.tsv`, or `{run_id}_1_VARAN.tsv`,
    `{run_id}_2_VARAN.tsv`, ... - whichever doesn't already exist - so
    re-running a batch never silently overwrites a previous report."""
    filename = f"{run_id}_VARAN.tsv"
    path = os.path.join(base_dir, filename)
    if not os.path.exists(path):
        return path
    counter = 1
    while True:
        filename = f"{run_id}_{counter}_VARAN.tsv"
        path = os.path.join(base_dir, filename)
        if not os.path.exists(path):
            return path
        counter += 1


def write_sample_tsv(rows: Sequence[SampleRow], path: str) -> None:
    """Write `rows` out in Varan's Templates/sample.tsv column order,
    appending any vendor-populated `extra` columns (union across all rows,
    stable-sorted) after the standard ten, blank where a given row doesn't
    set that key.
    """
    extra_cols: List[str] = []
    seen = set()
    for row in rows:
        for key in row.extra:
            if key not in seen:
                seen.add(key)
                extra_cols.append(key)

    header = SAMPLE_TSV_COLUMNS + extra_cols
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            base = [
                row.sample_id, row.patient_id, row.oncotree_code,
                row.snv_path, row.cnv_path, row.comb_path,
                row.msi, row.tmb, row.msi_thr, row.tmb_thr,
            ]
            extras = [row.extra.get(col, "") for col in extra_cols]
            f.write("\t".join(map(str, base + extras)) + "\n")
