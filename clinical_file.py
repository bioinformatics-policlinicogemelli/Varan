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

"""Add, overwrite or drop one attribute (column) of an existing
data_clinical_sample.txt, keeping its four "#" header rows in step - the
same way walk_clinical.py layers the BRCA exon and HRD columns onto it."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

SAMPLE_FILE = "data_clinical_sample.txt"


def _read(path: Path) -> tuple[list[list[str]], pd.DataFrame]:
    with path.open(encoding="utf-8") as f:
        header = [next(f).rstrip("\n").split("\t") for _ in range(4)]
    body = pd.read_csv(path, sep="\t", header=4, dtype=str, keep_default_na=False)
    return header, body


def _write(path: Path, header: list[list[str]], body: pd.DataFrame) -> None:
    with path.open("w", encoding="utf-8") as f:
        for row in header:
            f.write("\t".join(row) + "\n")
        # lineterminator="\n": f is a text-mode handle that already turns
        # "\n" into the platform's line ending (see walk_clinical.py).
        body.to_csv(f, sep="\t", index=False, lineterminator="\n")


def set_sample_column(
    folder: str | Path, column: str, values: dict[str, object],
    short: str, long: str, datatype: str = "STRING", default: str = "NA") -> bool:
    """Write `column` of data_clinical_sample.txt from {SAMPLE_ID: value}
    (`default` for samples not in it). A new column gets the given header
    rows; an existing one keeps its own. False if there's no such file."""
    path = Path(folder) / SAMPLE_FILE
    if not path.exists():
        return False
    header, body = _read(path)
    new_values = [
        str(values.get(str(s).strip(), default)) for s in body["SAMPLE_ID"]]
    if column not in body.columns:
        for row, value in zip(header, (short, long, datatype, "1")):
            row.append(value)
    body[column] = new_values
    _write(path, header, body)
    return True


def drop_sample_columns(folder: str | Path, columns: list[str]) -> None:
    """Remove `columns` (those present) from data_clinical_sample.txt."""
    path = Path(folder) / SAMPLE_FILE
    if not path.exists():
        return
    header, body = _read(path)
    present = [c for c in columns if c in body.columns]
    if not present:
        return
    keep = [i for i, c in enumerate(body.columns) if c not in present]
    header = [[row[i] for i in keep] for row in header]
    _write(path, header, body.drop(columns=present))
