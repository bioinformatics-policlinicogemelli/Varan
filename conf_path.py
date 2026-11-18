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

"""Paths of Varan's own files and of the conf.ini file it reads.

Varan's support files (cancer_list.txt, styles.css, the OncoKB annotator,
the cBioPortal validator...) are looked up in the folder where Varan is
installed, not in the current directory, so Varan can be run from any
folder (e.g. from a module that only puts it on the PATH).

By default Varan reads the conf.ini in that same folder. Setting the
VARAN_CONF environment variable points it to a configuration kept
elsewhere, so the code can be installed read-only and shared, while the
configuration stays outside it. An environment variable rather than a
command-line option because every module reads the configuration as soon
as it is imported, before varan.py parses its arguments.
"""

from __future__ import annotations

import os
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent


def code_path(*parts: str) -> str:
    """Return the path of a file of the Varan installation.

    Args:
        *parts (str): Path components relative to the Varan folder,
            e.g. ("oncokb-annotator", "MafAnnotator.py").

    Returns:
        str: The absolute path.

    """
    return str(CODE_DIR.joinpath(*parts))


def conf_path() -> str:
    """Return the conf.ini path to read.

    Returns:
        str: The VARAN_CONF environment variable if set, else the conf.ini
            in the Varan folder.

    Raises:
        FileNotFoundError: If VARAN_CONF is set but the file does not exist
            (configparser would otherwise silently read nothing).

    """
    path = os.environ.get("VARAN_CONF", "").strip()
    if not path:
        return code_path("conf.ini")
    if not Path(path).is_file():
        msg = f"VARAN_CONF points to a file that does not exist: {path}"
        raise FileNotFoundError(msg)
    return path
