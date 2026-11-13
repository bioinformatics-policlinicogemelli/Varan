"""Settings kept between sessions, and the files shipped with the app."""

from __future__ import annotations

import json
import os
import sys
from dataclasses import asdict, dataclass, fields
from pathlib import Path

DEFAULT_IMAGE = "ghcr.io/bioinformatics-policlinicogemelli/varan:latest"


def resource(name: str) -> Path:
    """A file shipped with the app: next to this module, or in the
    PyInstaller bundle (cancer_list.txt comes from Varan's own folder)."""
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    for folder in (base, base / "varan_gui", base.parent):
        if (folder / name).exists():
            return folder / name
    return base / name


def app_dir() -> Path:
    """Where settings and run folders live (created on first use)."""
    if sys.platform == "win32":
        root = Path(os.environ.get("APPDATA") or Path.home() / "AppData" / "Roaming")
        folder = root / "Varan"
    elif sys.platform == "darwin":
        folder = Path.home() / "Library" / "Application Support" / "Varan"
    else:
        root = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config")
        folder = root / "varan"
    folder.mkdir(parents=True, exist_ok=True)
    return folder


@dataclass
class Settings:
    engine: str = "docker"
    image: str = DEFAULT_IMAGE
    vep_cache: str = ""
    vep_version: str = ""
    ref_fasta: str = ""
    oncokb_token: str = ""
    conf_ini: str = ""
    last_dir: str = ""

    @classmethod
    def load(cls) -> Settings:
        path = app_dir() / "settings.json"
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return cls()
        known = {f.name for f in fields(cls)}
        return cls(**{k: str(v) for k, v in data.items() if k in known})

    def save(self) -> None:
        path = app_dir() / "settings.json"
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")


def vep_cache_root(folder: str | Path) -> Path | None:
    """The folder holding homo_sapiens/, given it, homo_sapiens itself or
    the folder one level above (e.g. vep_cache/homo_sapiens_vep_111_GRCh37)."""
    folder = Path(folder)
    if not folder.is_dir():
        return None
    if folder.name == "homo_sapiens":
        return folder.parent
    if (folder / "homo_sapiens").is_dir():
        return folder
    inner = [c for c in folder.iterdir() if c.is_dir() and (c / "homo_sapiens").is_dir()]
    return inner[0] if len(inner) == 1 else None


def vep_cache_versions(folder: str | Path) -> list[str]:
    """GRCh37 cache versions in a VEP cache folder (homo_sapiens/111_GRCh37
    gives "111"), newest first."""
    root = vep_cache_root(folder)
    if root is None:
        return []
    versions = set()
    for child in (root / "homo_sapiens").iterdir():
        head = child.name.split("_")[0]
        if child.is_dir() and head.isdigit() and "GRCh37" in child.name:
            versions.add(head)
    return sorted(versions, key=int, reverse=True)


def cancer_types() -> list[tuple[str, str]]:
    """(id, name) pairs from Varan's cancer_list.txt, for the -c menu."""
    try:
        lines = resource("cancer_list.txt").read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    out = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) >= 2 and cells[0].strip():
            out.append((cells[0].strip(), cells[1].strip()))
    return out
