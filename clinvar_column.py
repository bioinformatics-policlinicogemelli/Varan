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

"""ClinVar columns from the ClinVar file set in conf.ini (CLINV).

vcf2maf passes CLINV to VEP as a custom annotation, but it copies a custom
annotation into the MAF only when it is listed in --retain-ann. Without it,
the MAF only had VEP's CLIN_SIG column, filled from the ClinVar data bundled
in the VEP cache (for cache 116, ClinVar 2023-06), and the CLINV file had no
effect.

Since 2024 ClinVar keeps three separate classifications: CLNSIG (germline
classification), ONC (oncogenicity) and SCI (somatic clinical impact); a
somatic hotspot may have ONC and SCI but no CLNSIG. With CLINV set, Varan
annotates all three and keeps them in the MAF as <name>_CLNSIG, <name>_ONC
and <name>_SCI (ClinVar_CLNSIG, ClinVar_ONC, ClinVar_SCI with the usual
CLINV), and removes the older CLIN_SIG column, so that each MAF has a single,
up-to-date ClinVar source. Without CLINV the MAF keeps CLIN_SIG from the VEP
cache, as before.
"""

from __future__ import annotations

from pathlib import Path

from loguru import logger

FIELDS = ("CLNSIG", "ONC", "SCI")
CACHE_COLUMN = "CLIN_SIG"


def custom_annotation(clinv: str) -> tuple[str, list[str]]:
    """Return CLINV with the ClinVar fields Varan needs, and their columns.

    CLINV uses VEP's --custom syntax:
    file,short_name,format,type,force_report_coordinates,field[,field...]
    VEP names each annotation <short_name>_<field>. CLNSIG, ONC and SCI are
    added to the field list if missing, so older configurations that list
    only CLNSIG get all three.

    Args:
        clinv (str): The CLINV value from conf.ini.

    Returns:
        tuple: The CLINV string to pass to VEP, and the MAF columns to keep
        (empty if CLINV does not have the five leading parts).

    """
    parts = [p.strip() for p in clinv.split(",")]
    if len(parts) < 5 or not parts[1]:  # noqa: PLR2004
        return clinv, []
    fields = parts[5:] + [f for f in FIELDS if f not in parts[5:]]
    return (",".join(parts[:5] + fields),
            [f"{parts[1]}_{f}" for f in fields])


def clnsig_column(clinv: str) -> str | None:
    """Return the column with ClinVar's germline classification (CLNSIG).

    Args:
        clinv (str): The CLINV value from conf.ini.

    Returns:
        str or None: e.g. "ClinVar_CLNSIG", or None if CLINV is not usable.

    """
    _, columns = custom_annotation(clinv)
    return next((c for c in columns if c.endswith("_CLNSIG")), None)


def drop_cache_column(maf_path: str | Path, clnsig: str) -> None:
    """Remove the VEP-cache CLIN_SIG column from a MAF that has <clnsig>.

    Lines starting with "#" (e.g. "#version 2.4") and every other column are
    kept unchanged.

    Args:
        maf_path (str or Path): The MAF written by vcf2maf.
        clnsig (str): The ClinVar CLNSIG column, as from clnsig_column().

    """
    maf_path = Path(maf_path)
    if not maf_path.is_file():
        return
    lines = maf_path.read_text(encoding="utf-8").split("\n")
    out = []
    drop = None
    for line in lines:
        if drop is None and (line.startswith("#") or not line):
            out.append(line)
            continue
        fields = line.split("\t")
        if drop is None:
            if clnsig not in fields:
                logger.warning(f"{maf_path.name}: no {clnsig} column, "
                f"{CACHE_COLUMN} left as annotated by the VEP cache.")
                return
            if CACHE_COLUMN not in fields:
                return
            drop = fields.index(CACHE_COLUMN)
        if line and drop < len(fields):
            del fields[drop]
        out.append("\t".join(fields))
    maf_path.write_text("\n".join(out), encoding="utf-8")
