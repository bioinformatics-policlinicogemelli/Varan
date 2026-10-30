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

"""Threshold conditions from conf.ini (e.g. THRESHOLD_FUSION=>=15, AF=<0.0004,
TMB's {'Low':'<=5', ...}): an operator followed by a number.

Replaces building a Python expression out of the condition (and sometimes
out of the input value itself) and eval()-ing it: a malformed threshold now
fails with a message naming it, and nothing read from conf.ini or an input
file can be executed as code.
"""

from __future__ import annotations

import operator
import re
from collections.abc import Callable

import pandas as pd

_OPERATORS: dict[str, Callable] = {
    ">=": operator.ge, "<=": operator.le, "==": operator.eq,
    "!=": operator.ne, ">": operator.gt, "<": operator.lt,
}
_CONDITION = re.compile(
    r"^\s*(>=|<=|==|!=|>|<)\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*$")


def parse_threshold(condition: str) -> tuple[Callable, float]:
    """Split a condition like ">=15" into (operator function, 15.0)."""
    match = _CONDITION.match(str(condition))
    if not match:
        msg = (f"Invalid threshold '{condition}' in conf.ini: expected an "
               "operator (>=, <=, >, <, ==, !=) followed by a number, e.g. '>=15'.")
        raise ValueError(msg)
    return _OPERATORS[match.group(1)], float(match.group(2))


def meets_threshold(value: object, condition: str) -> bool:
    """Whether `value` (a number, or text holding one) satisfies `condition`.
    Raises ValueError if `value` isn't a number."""
    compare, threshold = parse_threshold(condition)
    return bool(compare(float(value), threshold))


def meets_threshold_series(values: pd.Series, condition: str) -> pd.Series:
    """meets_threshold() over a numeric Series (NaN never satisfies it)."""
    compare, threshold = parse_threshold(condition)
    return compare(values, threshold)
