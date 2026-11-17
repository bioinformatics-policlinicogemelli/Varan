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

"""Path of the conf.ini file Varan reads.

By default Varan reads ./conf.ini, i.e. the one in the code folder it is
run from. Setting the VARAN_CONF environment variable points it to a
configuration kept elsewhere, so the code can be installed read-only and
shared, while the configuration (paths and the OncoKB token) stays outside
it. An environment variable rather than a command-line option because
every module reads the configuration as soon as it is imported, before
varan.py parses its arguments.
"""

from __future__ import annotations

import os
from pathlib import Path


def conf_path() -> str:
    """Return the conf.ini path to read.

    Returns:
        str: The VARAN_CONF environment variable if set, else "conf.ini".

    Raises:
        FileNotFoundError: If VARAN_CONF is set but the file does not exist
            (configparser would otherwise silently read nothing).

    """
    path = os.environ.get("VARAN_CONF", "").strip()
    if not path:
        return "conf.ini"
    if not Path(path).is_file():
        msg = f"VARAN_CONF points to a file that does not exist: {path}"
        raise FileNotFoundError(msg)
    return path
