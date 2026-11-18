"""From the user's files to what the container sees.

The container only sees the folders bind-mounted into it, under Linux
paths. A run gets its own folder (the "staging" folder, under the app's
folder) holding conf.ini and copies of the input tables whose file paths
are rewritten to container paths; every folder those paths, the reference
data and the output point into is mounted.
"""

from __future__ import annotations

import configparser
import csv
import io
import os
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath

RUN_MOUNT = "/varan/run"
DATA_ROOT = "/varan/data"

# Varan's own tools inside the image (see the Dockerfile).
IMAGE_PATHS = {
    "VCF2MAF": "/mskcc-vcf2maf/vcf2maf.pl",
    "VEP_PATH": "/opt/vep/src/ensembl-vep",
}
IMAGE_CLINVAR = "/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG"


def _norm(path: Path) -> str:
    return os.path.normcase(str(path))


def _inside(path: Path, folder: Path) -> bool:
    try:
        path.relative_to(folder)
    except ValueError:
        return os.path.normcase(str(path)).startswith(os.path.normcase(str(folder)) + os.sep)
    return True


@dataclass
class Mounts:
    """Host folders mounted into the container, each under /varan/data/<n>."""

    folders: list[Path] = field(default_factory=list)

    def add(self, folder: Path) -> None:
        folder = Path(os.path.abspath(folder))
        if any(_inside(folder, f) for f in self.folders):
            return
        self.folders = [f for f in self.folders if not _inside(f, folder)] + [folder]

    def container_path(self, host: str | Path) -> str:
        host = Path(os.path.abspath(host))
        for n, folder in enumerate(self.folders):
            if _inside(host, folder):
                rel = os.path.relpath(host, folder)
                parts = [] if rel == "." else Path(rel).parts
                return str(PurePosixPath(DATA_ROOT, str(n), *parts))
        msg = f"{host} is not inside a mounted folder"
        raise ValueError(msg)

    def args(self, staging: Path) -> list[str]:
        out = _mount_arg(staging, RUN_MOUNT)
        for n, folder in enumerate(self.folders):
            out += _mount_arg(folder, f"{DATA_ROOT}/{n}")
        return out


def _mount_arg(source: Path, target: str) -> list[str]:
    # --mount instead of -v: a Windows drive letter's ":" can't be confused
    # with -v's separator. Its value is CSV, so quote a source with "," in it.
    src = f"source={source}"
    if "," in src or '"' in src:
        src = '"' + src.replace('"', '""') + '"'
    return ["--mount", f"type=bind,{src},target={target}"]


def as_host_path(value: str, base: Path) -> Path | None:
    """The host file or folder a table cell / conf value points at, if it
    is one: absolute, or relative to the table's own folder."""
    text = value.strip().strip('"').strip("'")
    if not text or "\n" in text or len(text) > 1000:
        return None
    path = Path(text)
    if not path.is_absolute():
        # A relative path needs a folder or an extension, so that ids and
        # codes ("S1", "LUAD", "27") are never taken for files.
        if text in (".", "..") or not any(c in text for c in ("/", "\\", ".")):
            return None
        path = base / text
    try:
        if path.exists():
            return Path(os.path.abspath(path))
    except OSError:
        return None
    return None


def path_cells(table: Path) -> list[Path]:
    """Every existing file or folder a TSV table points at."""
    rows = _read_tsv(table)
    found = []
    for row in rows[1:]:
        for cell in row:
            host = as_host_path(cell, table.parent)
            if host is not None:
                found.append(host)
    return found


# The columns Varan reads files from: sample.tsv's, and the Guardant
# adapter's selection.tsv aliases (vendor_adapters/guardant.py,
# _SELECTION_PATH_COLUMNS). Other columns, even named *_path (e.g. a
# RUN_PATH kept for bookkeeping), are never opened, so never checked.
FILE_COLUMNS = {
    "snv_path", "cnv_path", "comb_path",
    "xml_path", "metadata_path", "finalmetadata_path", "vcf_path", "snv_vcf_path",
    "msi_path", "msi_call_path", "cnv_call_path", "fusion_path", "fus_path", "fusion_call_path",
}


def missing_path_cells(table: Path) -> list[str]:
    """Cells of the file columns (FILE_COLUMNS) that point at nothing."""
    rows = _read_tsv(table)
    if not rows:
        return []
    header = [h.strip().lower() for h in rows[0]]
    cols = [i for i, h in enumerate(header) if h in FILE_COLUMNS]
    out = []
    for row in rows[1:]:
        for i in cols:
            if i < len(row) and row[i].strip() and as_host_path(row[i], table.parent) is None:
                out.append(row[i].strip())
    return out


def _read_tsv(table: Path) -> list[list[str]]:
    # utf-8-sig: tables saved by Excel start with a byte order mark.
    text = table.read_text(encoding="utf-8-sig", errors="replace")
    return [row for row in csv.reader(io.StringIO(text), delimiter="\t") if any(c.strip() for c in row)]


def rewrite_table(table: Path, mounts: Mounts, dest: Path) -> None:
    """Copy `table` to `dest` with its file paths as container paths, Unix
    line endings and no byte order mark."""
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    for n, row in enumerate(_read_tsv(table)):
        if n > 0:
            row = [_translate(cell, table.parent, mounts) for cell in row]
        writer.writerow(row)
    dest.write_text(out.getvalue(), encoding="utf-8", newline="")


def _translate(cell: str, base: Path, mounts: Mounts) -> str:
    host = as_host_path(cell, base)
    if host is None:
        return cell
    try:
        return mounts.container_path(host)
    except ValueError:
        return cell


def read_conf(path: Path | None) -> configparser.RawConfigParser:
    conf = configparser.RawConfigParser(interpolation=None)
    conf.optionxform = str
    if path is not None:
        conf.read(path, encoding="utf-8")
    return conf


def merge_conf(base: configparser.RawConfigParser,
               user: configparser.RawConfigParser) -> configparser.RawConfigParser:
    """`base` with every key the user's conf.ini sets."""
    for section in user.sections():
        if not base.has_section(section):
            base.add_section(section)
        for key, value in user.items(section):
            base.set(section, key, value)
    return base


def conf_host_paths(conf: configparser.RawConfigParser, base: Path) -> list[tuple[str, str, Path]]:
    """(section, key, host path) of conf values naming an existing file or
    folder, outside [Paths] (handled on its own)."""
    out = []
    for section in conf.sections():
        if section == "Paths":
            continue
        for key, value in conf.items(section):
            host = as_host_path(value, base) if os.path.isabs(value.strip()) else None
            if host is not None:
                out.append((section, key, host))
    return out


def write_conf(conf: configparser.RawConfigParser, dest: Path) -> None:
    with dest.open("w", encoding="utf-8", newline="\n") as handle:
        conf.write(handle)
