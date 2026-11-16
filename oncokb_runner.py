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

"""Run the bundled oncokb-annotator scripts safely and record the OncoKB
data release they annotated against.

Two problems this module exists to close:

1. Silent partial annotation. When an OncoKB API call fails for anything
   other than a 401, AnnotatorCore's pull_*_info() functions print
   "Error on annotating the url ...", write that alteration out with no
   annotation at all, and the script still exits 0 - so a network hiccup or
   an OncoKB outage would produce a study with unannotated variants and no
   error. run_annotator() treats that marker as a failure.

2. Traceability. OncoKB's knowledge base changes between releases, so the
   same variant can get a different oncogenicity/level on different days.
   The data release used is fetched once per process (at the first
   annotation) and written to <output_folder>/service/oncokb_info.json, next
   to provenance.json, so the report can show it. A file rather than a
   module-level value because the Snakemake path (walk_stage.py) runs the
   annotation stages and the report in separate processes.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import requests
from loguru import logger

# Same default base URL AnnotatorCore uses (Varan never passes -u).
ONCOKB_INFO_URL = "https://www.oncokb.org/api/v1/info"
INFO_TIMEOUT_S = 60

# Printed (to stdout) by AnnotatorCore.pull_*_info() for every alteration
# whose API call failed and was left unannotated.
ANNOTATION_ERROR_MARKER = "Error on annotating the url"

INFO_FILE = Path("service") / "oncokb_info.json"

_recorded_for: set[str] = set()


def fetch_data_version() -> str | None:
    """Return the current OncoKB data release, e.g. "v6.2 (2026-10-30)".

    Returns:
        str | None: The release and its date, or None if OncoKB's /info
            endpoint could not be reached or parsed.

    """
    try:
        response = requests.get(ONCOKB_INFO_URL, timeout=INFO_TIMEOUT_S)
        response.raise_for_status()
        data_version = response.json()["dataVersion"]
        return f"{data_version['version']} ({data_version['date']})"
    except (requests.RequestException, KeyError, TypeError, ValueError) as err:
        logger.warning(f"Could not read the OncoKB data version from "
                       f"{ONCOKB_INFO_URL}: {err}")
        return None


def record_data_version(output_folder: str | Path) -> None:
    """Write the OncoKB data release to <output_folder>/service, once per process.

    Args:
        output_folder (str | Path): The study's output folder for this run.

    Returns:
        None

    """
    key = str(Path(output_folder).resolve())
    if key in _recorded_for:
        return
    _recorded_for.add(key)

    data_version = fetch_data_version()
    info_path = Path(output_folder) / INFO_FILE
    info_path.parent.mkdir(parents=True, exist_ok=True)
    # Per-process temp name: Snakemake can run the CNV/SNV/fusion stages
    # (each of which may annotate) at the same time.
    tmp_path = info_path.with_suffix(f".json.{os.getpid()}.tmp")
    with tmp_path.open("w") as f:
        json.dump({"oncokb_data_version": data_version}, f, indent=2)
    tmp_path.replace(info_path)
    logger.info(f"OncoKB data version: {data_version or 'unavailable'}")


def read_data_version(output_folder: str | Path) -> str | None:
    """Read the OncoKB data release recorded for this output folder, if any.

    Args:
        output_folder (str | Path): The study's output folder.

    Returns:
        str | None: The recorded release, "unavailable" if annotation ran but
            the release could not be fetched, or None if no annotation ran.

    """
    info_path = Path(output_folder) / INFO_FILE
    if not info_path.exists():
        return None
    try:
        with info_path.open() as f:
            return json.load(f).get("oncokb_data_version") or "unavailable"
    except (json.JSONDecodeError, OSError):
        return "unavailable"


def run_annotator(cmd: list[str], description: str,
                  output_folder: str | Path) -> None:
    """Run one oncokb-annotator command and fail loudly on any annotation error.

    Args:
        cmd (list[str]): The full annotator command line.
        description (str): What is being annotated, for log/error messages
            (e.g. "sample X (MAF)").
        output_folder (str | Path): The study's output folder, where the
            OncoKB data release is recorded.

    Returns:
        None

    Raises:
        RuntimeError: If the annotator exits non-zero, or if it reports that
            any alteration could not be annotated.

    """
    record_data_version(output_folder)

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    for line in (result.stdout + result.stderr).splitlines():
        if line.strip():
            logger.debug(f"[oncokb-annotator] {line}")

    if result.returncode != 0:
        tail = "\n".join(result.stderr.strip().splitlines()[-10:])
        msg = (f"OncoKB annotation failed for {description} "
               f"(exit code {result.returncode}).\n{tail}")
        raise RuntimeError(msg)

    failed = [line for line in result.stdout.splitlines()
              if ANNOTATION_ERROR_MARKER in line]
    if failed:
        msg = (f"OncoKB could not annotate {len(failed)} alteration(s) for "
               f"{description}: the annotator left them without annotation. "
               f"First failed request: {failed[0]}")
        raise RuntimeError(msg)
