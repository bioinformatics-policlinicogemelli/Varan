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

"""Locate the programs of the VEP installation that conf.ini points to.

vcf2maf runs VEP with the same perl that runs vcf2maf itself, and looks up
samtools and tabix on PATH. Calling plain "perl" therefore makes VEP depend
on whichever conda environment happens to be active, so VEP and Varan had
to live in one environment. When VEP_PATH is the bin/ folder of a conda
environment that holds VEP (it then also holds perl, samtools and tabix),
Varan uses those programs instead, and VEP can sit in its own environment.
Anything not found there falls back to PATH, as before.
"""

from __future__ import annotations

import os
from pathlib import Path


def vep_tool(vep_path: str | Path, name: str) -> str | None:
    """Return the executable `name` inside VEP_PATH, or None if absent.

    Args:
        vep_path (str | Path): conf.ini [Paths] VEP_PATH.
        name (str): Program name, e.g. "perl", "samtools", "tabix".

    Returns:
        str | None: Full path of the program, or None if VEP_PATH has no
            executable with that name.

    """
    candidate = Path(vep_path) / name
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def vep_environment(vep_path: str | Path) -> dict[str, str] | None:
    """Return an environment for running VEP with VEP_PATH first on PATH.

    VEP also looks up some programs on PATH by itself (for instance it runs
    `which tabix` when it loads its variant cache, printing "which: no
    tabix in ..." if it is missing). Putting VEP's own bin/ folder first
    makes it find the programs of its own environment, as if that
    environment were active.

    Args:
        vep_path (str | Path): conf.ini [Paths] VEP_PATH.

    Returns:
        dict[str, str] | None: A copy of the current environment with
            VEP_PATH prepended to PATH, or None (= inherit the current
            environment unchanged) when VEP_PATH is not a conda bin/ folder
            with its own perl.

    """
    if not vep_tool(vep_path, "perl"):
        return None
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(vep_path), env.get("PATH", "")])
    return env


def vep_perl(vep_path: str | Path) -> str:
    """Return the perl to run VEP and vcf2maf with.

    Args:
        vep_path (str | Path): conf.ini [Paths] VEP_PATH.

    Returns:
        str: VEP_PATH/perl if it exists, otherwise "perl" from PATH.

    """
    return vep_tool(vep_path, "perl") or "perl"
