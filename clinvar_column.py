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

"""CLIN_SIG from the ClinVar file set in conf.ini (CLINV).

vcf2maf passes CLINV to VEP as a custom annotation, but it copies a custom
annotation into the MAF only when it is listed in --retain-ann. Without it,
the MAF's CLIN_SIG column comes from the ClinVar data bundled in the VEP
cache (for cache 116, ClinVar 2023-06) and the CLINV file has no effect.

Varan therefore asks vcf2maf to keep the custom column (e.g. ClinVar_CLNSIG)
and then moves its values into CLIN_SIG, dropping the extra column: the MAF
keeps the same columns as before, and CLIN_SIG holds the clinical
significance from the ClinVar release set in conf.ini (ClinVar's own
spelling, e.g. "Benign/Likely_benign"). If CLINV is not set, CLIN_SIG still
comes from the VEP cache.
"""

from __future__ import annotations

from pathlib import Path

from loguru import logger


def custom_column(clinv: str) -> str | None:
    """Return the MAF column VEP fills from the CLINV custom annotation.

    CLINV uses VEP's --custom syntax:
    file,short_name,format,type,force_report_coordinates,field[,field...]
    and VEP names the annotation <short_name>_<field>, e.g. ClinVar_CLNSIG.

    Args:
        clinv (str): The CLINV value from conf.ini.

    Returns:
        str or None: The column of the first field, or None if CLINV is empty
        or does not name a field.

    """
    parts = [p.strip() for p in clinv.split(",")]
    if len(parts) < 6 or not parts[1] or not parts[5]:  # noqa: PLR2004
        return None
    return f"{parts[1]}_{parts[5]}"


def use_clinvar_file(maf_path: str | Path, column: str) -> None:
    """Replace CLIN_SIG with the values of the ClinVar custom column.

    Lines starting with "#" (e.g. "#version 2.4") are kept unchanged, and so
    is every other column; only CLIN_SIG gets new values and the custom
    column is removed.

    Args:
        maf_path (str or Path): The MAF written by vcf2maf.
        column (str): The custom column, as returned by custom_column().

    """
    maf_path = Path(maf_path)
    if not maf_path.is_file():
        return
    lines = maf_path.read_text(encoding="utf-8").split("\n")
    out = []
    header = None
    for line in lines:
        if header is None and (line.startswith("#") or not line):
            out.append(line)
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            if "CLIN_SIG" not in header or column not in header:
                logger.warning(f"{maf_path.name}: no CLIN_SIG or {column} "
                "column, CLIN_SIG left as annotated by the VEP cache.")
                return
            i_sig = header.index("CLIN_SIG")
            i_custom = header.index(column)
        elif line:
            fields += [""] * (len(header) - len(fields))
            fields[i_sig] = fields[i_custom]
        if line:
            del fields[i_custom]
        out.append("\t".join(fields))
    maf_path.write_text("\n".join(out), encoding="utf-8")
