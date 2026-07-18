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

"""Shared conf.ini loader, so `--config <path>` on the CLI actually works.

Every Varan module used to build its own `ConfigParser` and call
`config.read("conf.ini")` at import time (module-level code), reading a
literal, cwd-relative path before varan.py had even parsed its arguments.
That made a `--config` CLI flag pointless: by the time the flag was read,
every already-imported module had already read (or failed to find)
"./conf.ini" on its own.

The fix is two-part:
1. This module holds the *one* ConfigParser instance for the whole process,
   built lazily on first use via `get_config()`.
2. `varan.py` calls `set_config_path()` immediately after parsing `--config`,
   and only *then* imports the other Varan modules - so their first call to
   `get_config()` (still at their own module level, that part didn't need to
   change) reads the path the user actually asked for.

Any module that today does `config = ConfigParser(); config.read("conf.ini")`
should instead do `from config_loader import get_config` and use
`get_config()` wherever it used to use `config` - either assigning it to a
module-level name right after the import (safe *only* as long as varan.py
keeps importing that module after `set_config_path()`, as it now does), or
calling `get_config()` directly inside each function for extra safety.
"""

from __future__ import annotations

from configparser import ConfigParser

_config_path = "conf.ini"
_config: ConfigParser | None = None


def set_config_path(path: str) -> None:
    """Set the conf.ini path this run should use.

    Must be called before the first `get_config()` call actually reads the
    file - in practice, before any other Varan module is imported, since
    most of them read the config at import time.

    Args:
        path (str): Path to the conf.ini file to use for this run.

    Returns:
        None

    """
    global _config_path, _config
    _config_path = path
    _config = None  # force a re-read on the next get_config() call


def get_config() -> ConfigParser:
    """Return the process-wide ConfigParser, reading it lazily on first use.

    Returns:
        ConfigParser: The parsed configuration. If the file at the configured
            path doesn't exist, this is an empty ConfigParser (configparser's
            own behavior - `.read()` never raises for a missing file), and a
            warning is logged so a missing/misspelled --config path doesn't
            fail silently.

    """
    global _config
    if _config is None:
        _config = ConfigParser()
        found = _config.read(_config_path)
        if not found:
            from loguru import logger
            logger.warning(
                f"conf.ini not found at '{_config_path}' - continuing with "
                "no configuration loaded. Use --config <path> to point at "
                "the right file if this wasn't expected.")
    return _config
