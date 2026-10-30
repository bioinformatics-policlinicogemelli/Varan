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

This module holds the *one* ConfigParser instance for the whole process,
built lazily on first use via `get_config()`, from the path set by
`set_config_path()` (varan.py, walk_stage.py and pipeline_stage.py call it
right after parsing `-C/--config`).

Modules use either `get_config()` inside their functions, or the module-
level `config` proxy below (`from config_loader import config`), which
resolves `get_config()` on every attribute access. Neither reads conf.ini
at import time, so importing a Varan module before `set_config_path()`
(or before varan.py's [Sample_Type] reconciliation) is harmless: nothing
is read until a function actually runs.
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


class _LazyConfig:
    """Stand-in for get_config()'s ConfigParser, looked up on every use
    (`config.get(...)`, `config.has_section(...)`, ...), so a module can keep
    a module-level `config` name without reading conf.ini at import time."""

    def __getattr__(self, name: str):
        return getattr(get_config(), name)

    # Special methods are looked up on the class, never via __getattr__:
    # delegate the ones ConfigParser users rely on ("x" in config,
    # config["Section"], iterating sections) explicitly.
    def __contains__(self, key: object) -> bool:
        return key in get_config()

    def __getitem__(self, key: str):
        return get_config()[key]

    def __iter__(self):
        return iter(get_config())

    def __len__(self) -> int:
        return len(get_config())


config = _LazyConfig()
