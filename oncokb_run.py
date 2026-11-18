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

"""Run an oncokb-annotator script without the token on the command line.

Usage: python3 oncokb_run.py oncokb-annotator/<Script>.py [options]

The annotator scripts take the token only as "-b <token>". Passing it that
way would show it to every user of the node in `ps`. This wrapper reads it
from the ONCOKB_API_TOKEN environment variable and adds "-b <token>" to the
script's arguments inside this process, then runs the script unchanged, as
if it had been launched directly. Only the other arguments are visible.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

from oncokb_token import TOKEN_ENV, get_token


def main() -> None:
    """Run the annotator script given as first argument with the token."""
    if len(sys.argv) < 2:  # noqa: PLR2004
        sys.exit("usage: oncokb_run.py <annotator script> [options]")
    token = get_token()
    if not token:
        sys.exit(f"{TOKEN_ENV} is not set: cannot run the OncoKB annotation.")
    script = Path(sys.argv[1]).resolve()
    # Exactly as when the script is run directly: its folder, not Varan's,
    # is sys.path[0] (the annotators import AnnotatorCore from there) and
    # sys.argv[0] is the script itself.
    sys.path[0] = str(script.parent)
    sys.argv = [str(script), *sys.argv[2:], "-b", token]
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
