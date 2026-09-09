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

"""Shared logging helper used by varan.py's own startup and by
Update_script.py/Delete_script.py/ExtractSamples_script.py, each of which
resolves its own versioned output folder independently (there's no single
place upstream, like walk_folder() for a normal create run, where all
three converge) and wants a copy of the log written there too.
"""

from __future__ import annotations

from pathlib import Path

from loguru import logger


def attach_file_log_sink(directory: Path, logfile: str) -> Path | None:
    """Attach a loguru file sink under `directory`, creating it if needed.

    Never raises. A shared cluster working directory's `Logs/` folder can
    easily end up owned by whichever user ran Varan there first, with
    permissions that then block every other user's own runs from writing
    a log at all - a logging setup problem should never be what crashes an
    otherwise-fine run, so any failure here (permission denied, read-only
    filesystem, whatever) is swallowed and logging just falls back to
    stderr only (already attached separately, unaffected by this).

    Returns the path actually used, or None if it couldn't be created.
    """
    try:
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / logfile
        logger.add(
            path,
            format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",
            mode="w")
        return path
    except OSError as err:
        logger.warning(f"Could not set up a log file under {directory}: {err}")
        return None
