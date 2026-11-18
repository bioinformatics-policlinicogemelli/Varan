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

"""OncoKB API token.

Varan reads the token only from the ONCOKB_API_TOKEN environment variable.
It is never stored in conf.ini (a file that gets copied, shared and
committed) and never passed on a command line, where any user on the same
node could read it with `ps`: oncokb_run.py hands it to the annotator inside
its own process. Varan never writes the token to logs or reports.
"""

from __future__ import annotations

import os

TOKEN_ENV = "ONCOKB_API_TOKEN"


def get_token() -> str:
    """Return the OncoKB API token, or "" if it is not set.

    Returns:
        str: The value of the ONCOKB_API_TOKEN environment variable, stripped.

    """
    return os.environ.get(TOKEN_ENV, "").strip()
