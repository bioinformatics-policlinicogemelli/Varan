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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from pathlib import Path


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
    output_folder_base = re.sub(r"_v\d+$", "", output_folder)
    old_versions = {\
        file.name for file in Path(outputfolderpath).resolve().iterdir()\
            if re.match(rf"^{re.escape(output_folder_base)}_v[0-9]+$", file.name)}

    total_dic = {}

    for folder in old_versions:
        clin_sam_path = Path(outputfolderpath).resolve() \
            / folder / "data_clinical_sample.txt"
        if clin_sam_path.exists():
            clin_sam_df = pd.read_csv(clin_sam_path, sep="\t", header = 4)
            unique_sam = len(set(clin_sam_df.iloc[:, 0]))
            unique_pat = len(set(clin_sam_df.iloc[:, 1]))
        else:
            unique_sam = 0
            unique_pat = 0

        total_dic[folder] = [unique_sam, unique_pat]

    n = 5
    sorted_total = dict(sorted(total_dic.items(), key=lambda item: \
        int(re.search(r"_v(\d+)$", item[0]).group(1))))
    limited_dic = dict(list(sorted_total.items())[-n:])

    values = list(limited_dic.values())
    keys = [*limited_dic]
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


    total_genes = {}

    for folder in old_versions:
        data_snv_path = Path(outputfolderpath).resolve() \
            / folder / "data_mutations_extended.txt"
        data_cnv_path = Path(outputfolderpath).resolve() \
            / folder / "data_cna.txt"
        data_sv_path = Path(outputfolderpath).resolve() \
            / folder / "data_sv.txt"

        if data_snv_path.exists():
            data_snv_df = pd.read_csv(data_snv_path, \
                sep="\t", header=0, low_memory=False)
            snv_number = len(data_snv_df)
        else:
            snv_number = 0

        if data_cnv_path.exists():
            data_cnv_df = pd.read_csv(data_cnv_path, \
                sep="\t", header=0, index_col=0)
            cnv_number = (data_cnv_df != 0).sum().sum()
        else:
            cnv_number = 0

        if data_sv_path.exists():
            data_sv_df = pd.read_csv(data_sv_path, \
                sep="\t", header=0)
            sv_number = len(data_sv_df)
        else:
            sv_number = 0

        total_genes[folder] = [snv_number, cnv_number, sv_number]


    sorted_total = dict(sorted(total_genes.items(), \
        key=lambda item: int(re.search(r"_v(\d+)$", item[0]).group(1))))
    limited_dic = dict(list(sorted_total.items())[-n:])

    studies = list(limited_dic.keys())
    values = list(limited_dic.values())


    snv = [count[0] for count in values]
    cnv = [count[1] for count in values]
    sv = [count[2] for count in values]

    categories = {"SNV": snv, "CNV": cnv, "SV": sv}
    colors = {"SNV": "blue", "CNV": "orange", "SV": "green"}

    fig, axes = plt.subplots(
        1, len(categories),
        figsize=(10, 5),
        sharey=True,
        gridspec_kw={"wspace": 0.1},
        constrained_layout=True,
    )

    for i, (cat, values) in enumerate(categories.items()):
        bars = axes[i].barh(studies, values, color=colors[cat])
        axes[i].set_title(cat, fontsize=10)
        axes[i].set_xlabel("Counts", fontsize=9)
        if i == 0:
            axes[i].invert_yaxis()

        max_value = max(values) if values else 0
        if max_value > 0:
            axes[i].set_xlim([0, max_value * 1.15])
        else:
            axes[i].set_xlim([0, 10])


        for bar in bars:
            width = bar.get_width()
            axes[i].text(width + max_value * 0.02, bar.get_y() + bar.get_height()/2,
                         f"{width}", va="center", ha="left", fontsize=8)

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
