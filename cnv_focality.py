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

"""Classify CNV segments as focal or broad (arm-level), hg19/GRCh37 only.

Follows the GISTIC2 convention: a segment is "broad" when it covers at
least `threshold` (default 0.5) of the chromosome arm it overlaps;
otherwise it's "focal". A segment spanning the centromere - so already
covering more than a single arm - is broad by definition regardless of
`threshold`.

This is informational only (see classify_focality's return value) and is
not meant to feed any filtering/threshold logic - callers should treat it
as an extra column, not a gate.
"""

from __future__ import annotations

# UCSC hg19 (GRCh37) chromosome sizes (chrom.sizes), autosomes + X/Y only -
# this module doesn't classify segments on other contigs (e.g. chrM).
CHROM_LENGTH_HG19: dict[str, int] = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
}

# UCSC hg19 centromere gap boundaries (gap.txt, type="centromere"), used to
# split each chromosome into its p/q arms.
CENTROMERE_HG19: dict[str, tuple[int, int]] = {
    "1": (121535434, 124535434), "2": (92326171, 95326171),
    "3": (90504854, 93504854), "4": (49660117, 52660117),
    "5": (46405641, 49405641), "6": (58830166, 61830166),
    "7": (58054331, 61054331), "8": (43838887, 46838887),
    "9": (47367679, 50367679), "10": (39254935, 42254935),
    "11": (51644205, 54644205), "12": (34856694, 37856694),
    "13": (16000000, 19000000), "14": (16000000, 19000000),
    "15": (17000000, 19000000), "16": (35335801, 38335801),
    "17": (22263006, 25263006), "18": (15460898, 18460898),
    "19": (24681782, 27681782), "20": (26369569, 29369569),
    "21": (11288129, 14288129), "22": (13000000, 16000000),
    "X": (58632012, 61632012), "Y": (10104553, 13104553),
}


def classify_focality(
    chrom: str, start: int, end: int, threshold: float = 0.5,
    ) -> str | None:
    """Classify a CNV segment as "focal" or "broad" (arm-level).

    Parameters
    ----------
    chrom : str
        Chromosome name, with or without a "chr" prefix (e.g. "7", "chr7",
        "X").
    start : int
        Segment start position (1-based, as in the CNV table).
    end : int
        Segment end position.
    threshold : float
        Minimum fraction of the overlapping arm's length a segment must
        cover to be called "broad" - 0.5 (50%) is the GISTIC2 convention.

    Returns
    -------
    str | None
        "focal" or "broad", or None if `chrom` isn't a known hg19
        chromosome (e.g. an alt contig) and can't be classified.

    """
    chrom = str(chrom).strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]

    if chrom not in CHROM_LENGTH_HG19:
        return None

    chrom_length = CHROM_LENGTH_HG19[chrom]
    cent_start, cent_end = CENTROMERE_HG19[chrom]
    start, end = int(start), int(end)

    if start < cent_start and end > cent_end:
        # Spans the centromere - already more than one arm.
        return "broad"
    if end <= cent_start:
        arm_length = cent_start
    elif start >= cent_end:
        arm_length = chrom_length - cent_end
    else:
        # Overlaps the centromere gap without clearing it on both sides.
        return "broad"

    seg_length = end - start
    fraction = seg_length / arm_length if arm_length > 0 else 1.0
    return "broad" if fraction >= threshold else "focal"
