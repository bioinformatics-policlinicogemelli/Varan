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


"""Module for generating bar plots summarizing data across versions.

This script reads clinical sample and mutation data from different versioned
folders and produces visual summaries:
- Bar plot comparing the number of samples and patients across the latest 5 versions.
- Horizontal bar plots showing counts of SNVs, CNVs, and SVs per version.

Expected folder structure:
Each versioned folder (e.g., 'project_v1', 'project_v2') should contain:
- data_clinical_sample.txt
- data_mutations_extended.txt
- data_cna.txt
- data_sv.txt

Plots are saved in an `img` subdirectory under the specified output folder.

"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle


def load_clinical_data(outputfolderpath: Path, folder: str) -> tuple:
    """Load clinical data and return counts for unique samples and patients."""
    clin_sam_path = outputfolderpath / folder / "data_clinical_sample.txt"
    if clin_sam_path.exists():
        clin_sam_df = pd.read_csv(clin_sam_path, sep="\t", header=4)
        unique_sam = len(set(clin_sam_df.iloc[:, 0]))
        unique_pat = len(set(clin_sam_df.iloc[:, 1]))
    else:
        unique_sam = unique_pat = 0
    return unique_sam, unique_pat

# Fixed labels/colors for the 4 non-neutral discrete CNA values, following
# cBioPortal's own oncoprint convention (blue = loss, red = gain) so the
# stacked bar reads consistently with the rest of the report.
CNA_CATEGORIES = {
    -2: ("Deep Deletion", "#2166AC"),
    -1: ("Shallow Deletion", "#92C5DE"),
    1: ("Gain", "#F4A582"),
    2: ("Amplification", "#B2182B"),
}

# Common MAF Variant_Classification values get a fixed color so the legend
# stays stable across reports; anything else (rarer classifications) falls
# back to a cycling qualitative palette assigned on first appearance.
SNV_CATEGORY_COLORS = {
    "Missense_Mutation": "#4DAF4A",
    "Nonsense_Mutation": "#000000",
    "Frame_Shift_Ins": "#984EA3",
    "Frame_Shift_Del": "#A65628",
    "In_Frame_Ins": "#FF7F00",
    "In_Frame_Del": "#FFD92F",
    "Splice_Site": "#377EB8",
    "Translation_Start_Site": "#F781BF",
    "Nonstop_Mutation": "#999999",
}
_SNV_FALLBACK_PALETTE = plt.get_cmap("tab20").colors


def load_genomic_data(outputfolderpath: Path, folder: str) -> tuple:
    """Load genomic data (SNV, CNV, SV) and return per-category breakdowns.

    Returns
    -------
    tuple
        (snv_breakdown, cnv_breakdown, sv_number): snv_breakdown maps
        Variant_Classification -> row count; cnv_breakdown maps each
        non-zero discrete CNA value (-2, -1, 1, 2) -> call count across all
        genes/samples; sv_number is a plain row count (SV isn't broken down
        further here).

    """
    snv_breakdown: dict = {}
    cnv_breakdown: dict = {}
    sv_number = 0

    # Load SNV data
    data_snv_path = outputfolderpath / folder / "data_mutations_extended.txt"
    if data_snv_path.exists():
        data_snv_df = pd.read_csv(data_snv_path, sep="\t", header=0, low_memory=False)
        if "Variant_Classification" in data_snv_df.columns:
            snv_breakdown = data_snv_df["Variant_Classification"].value_counts().to_dict()
        elif len(data_snv_df):
            snv_breakdown = {"Unknown": len(data_snv_df)}

    # Load CNV data
    data_cnv_path = outputfolderpath / folder / "data_cna.txt"
    if data_cnv_path.exists():
        data_cnv_df = pd.read_csv(data_cnv_path, sep="\t", header=0, index_col=0)
        calls = data_cnv_df.to_numpy().flatten()
        for value in (-2, -1, 1, 2):
            count = int((calls == value).sum())
            if count:
                cnv_breakdown[value] = count

    # Load SV data
    data_sv_path = outputfolderpath / folder / "data_sv.txt"
    if data_sv_path.exists():
        data_sv_df = pd.read_csv(data_sv_path, sep="\t", header=0)
        sv_number = len(data_sv_df)

    return snv_breakdown, cnv_breakdown, sv_number


def _snv_color(category: str, assigned: dict) -> str:
    """Return a stable color for an SNV Variant_Classification category."""
    if category in SNV_CATEGORY_COLORS:
        return SNV_CATEGORY_COLORS[category]
    if category not in assigned:
        assigned[category] = _SNV_FALLBACK_PALETTE[
            len(assigned) % len(_SNV_FALLBACK_PALETTE)]
    return assigned[category]


def _draw_stacked_barh(
    ax: Axes, studies: list, breakdowns: list[dict],
    category_order: list, color_of: dict) -> None:
    """Draw one horizontal stacked bar per study, segmented by category.

    Bar positions/order match the plain (non-stacked) bars used elsewhere in
    this report - only the coloring changes, splitting each bar into its
    per-category segments instead of one solid color.

    Args:
        ax (Axes): Subplot to draw into.
        studies (list): Version folder names, one bar each, in draw order.
        breakdowns (list[dict]): Per-study {category: count} dicts, same
            order as `studies`.
        category_order (list): Categories to stack, in stacking order.
        color_of (dict): category -> matplotlib color.

    Returns:
        None

    """
    y_pos = np.arange(len(studies))
    left = np.zeros(len(studies))

    for category in category_order:
        widths = np.array([b.get(category, 0) for b in breakdowns], dtype=float)
        if not widths.any():
            continue
        ax.barh(y_pos, widths, left=left, color=color_of[category],
                label=str(category), height=0.8)
        left += widths

    ax.set_yticks(y_pos)
    ax.set_yticklabels(studies)

def create_general_plot(limited_dic: dict, output_folder: str) -> None:
    """Generate and save the general plot for samples and patients."""
    values = list(limited_dic.values())
    keys = list(limited_dic.keys())
    samples = [count[0] for count in values]
    patients = [count[1] for count in values]

    x = np.arange(len(keys))
    width = 0.3

    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - width/2, samples, width, label="Samples", color="#008080")
    rects2 = ax.bar(x + width/2, patients, width, label="Patients", color="#FF6F61")

    ax.set_ylabel("Count")
    ax.set_xticks(x)
    ax.set_xticklabels(keys, rotation=45, ha="right")
    ax.legend()

    autolabel_ver(ax, rects1)
    autolabel_ver(ax, rects2)

    fig.tight_layout()
    plt.show()
    plt.savefig(Path(output_folder) / "img" / "general.png")

def create_barplots(output_folder: str) -> int:
    """Generate bar plots of clinical and genomic data from multiple dataset versions.

    Parameters
    ----------
    output_folder : str
        Path to the current output folder. The function will look for
        other versions in the same parent directory and save plots
        in an 'img' subdirectory under this folder.

    Returns
    -------
    int: Number of versioned folders processed.

    """
    (Path(output_folder) / "img").mkdir(parents=True, exist_ok=True)
    outputfolderpath = Path(output_folder).parent
    output_folder_base = re.sub(r"_v\d+$", "", Path(output_folder).name)
    old_versions = {\
        file.name for file in Path(outputfolderpath).resolve().iterdir()\
            if re.match(rf"^{re.escape(output_folder_base)}_v[0-9]+$", file.name)}

    total_dic = {}

    for folder in old_versions:
        unique_sam, unique_pat = load_clinical_data(outputfolderpath, folder)
        total_dic[folder] = [unique_sam, unique_pat]

    n = 5
    sorted_total = dict(sorted(total_dic.items(), key=lambda item: \
        int(re.search(r"_v(\d+)$", item[0]).group(1))))
    limited_dic = dict(list(sorted_total.items())[-n:])

    create_general_plot(limited_dic, output_folder)

    total_genes = {}

    for folder in old_versions:
        snv_breakdown, cnv_breakdown, sv_number = load_genomic_data(
            outputfolderpath, folder)
        total_genes[folder] = [snv_breakdown, cnv_breakdown, sv_number]

    sorted_total = dict(sorted(total_genes.items(), \
        key=lambda item: int(re.search(r"_v(\d+)$", item[0]).group(1))))
    limited_dic = dict(list(sorted_total.items())[-n:])

    studies = list(limited_dic.keys())
    values = list(limited_dic.values())

    snv_breakdowns = [count[0] for count in values]
    cnv_breakdowns = [{CNA_CATEGORIES[k][0]: v for k, v in count[1].items()}
                       for count in values]
    sv = [count[2] for count in values]

    fig, axes = plt.subplots(
        1, 3,
        figsize=(10, 5),
        sharey=True,
        gridspec_kw={"wspace": 0.1},
        constrained_layout=True,
    )

    # SNV: stacked by Variant_Classification, one color per category.
    snv_categories = sorted({cat for b in snv_breakdowns for cat in b})
    snv_colors: dict = {}
    snv_color_of = {cat: _snv_color(cat, snv_colors) for cat in snv_categories}
    _draw_stacked_barh(axes[0], studies, snv_breakdowns, snv_categories, snv_color_of)
    axes[0].set_title("SNV", fontsize=10)
    axes[0].invert_yaxis()
    if snv_categories:
        axes[0].legend(fontsize=6, loc="lower right")

    # CNV: stacked by discrete CNA value, fixed cBioPortal-style colors.
    cnv_categories = [label for _, (label, _) in CNA_CATEGORIES.items()]
    cnv_color_of = {label: color for _, (label, color) in CNA_CATEGORIES.items()}
    _draw_stacked_barh(axes[1], studies, cnv_breakdowns, cnv_categories, cnv_color_of)
    axes[1].set_title("CNV", fontsize=10)
    if any(cnv_breakdowns):
        axes[1].legend(fontsize=6, loc="lower right")

    # SV: unchanged, single count per version (not broken down further).
    bars = axes[2].barh(studies, sv, color="green")
    axes[2].set_title("SV", fontsize=10)

    for i, breakdowns in enumerate([snv_breakdowns, cnv_breakdowns, None]):
        axes[i].set_xlabel("Counts", fontsize=9)
        if breakdowns is not None:
            totals = [sum(b.values()) for b in breakdowns]
        else:
            totals = sv
        max_value = max(totals) if totals else 0
        axes[i].set_xlim([0, max_value * 1.15] if max_value > 0 else [0, 10])
        y_pos = np.arange(len(studies))
        for y, total in zip(y_pos, totals):
            axes[i].text(total + max_value * 0.02, y, f"{total}",
                         va="center", ha="left", fontsize=8)

    plt.savefig(Path(output_folder) / "img" / "genes.png")

    return len(old_versions)

def autolabel_ver(ax: Axes, rects: Rectangle) -> None:
    """Attach a text label above each bar in a bar chart.

    Parameters
    ----------
    ax : Axes
        The axes object to which the bar chart is attached. It is used to
        position and place the text labels.
    rects : Rectangle
        A collection of Rectangle objects (bars) in the bar chart. Each
        Rectangle represents a bar in the chart, and the function will
        label each bar with its height value.

    """
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f"{height}",
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center", va="bottom")
