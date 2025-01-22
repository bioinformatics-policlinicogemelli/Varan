#####################################
# NAME: create_graphs.py
# Date: 08/01/2025
version = "1.0"
# ===================================

import os 
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def create_barplots(output_folder):
    outputfolderpath = os.path.dirname(output_folder)
    output_folder_base = re.sub(r'_v\d+$', '', output_folder)
    old_versions = [file for file in os.listdir(os.path.realpath(outputfolderpath)) if re.split(fr'{output_folder_base}_v[0-9]+$', os.path.basename(output_folder_base))[0] + "_v" in file]
    
    old_versions.append(output_folder)
    old_versions = set(old_versions)
    
    total_dic = {}

    for folder in old_versions:
        clin_sam_path = os.path.join(os.path.realpath(outputfolderpath), folder, "data_clinical_sample.txt")
        if os.path.exists(clin_sam_path):
            clin_sam_df = pd.read_csv(clin_sam_path, sep="\t", header = 4)
            unique_sam = len(set(clin_sam_df.iloc[:, 0]))
            unique_pat = len(set(clin_sam_df.iloc[:, 1]))
        else:
            unique_sam = 0
            unique_pat = 0

        total_dic[folder] = [unique_sam, unique_pat]

    n = 5
    sorted_total = dict(sorted(total_dic.items(), key=lambda item: int(re.search(r'_v(\d+)$', item[0]).group(1))))
    limited_dic = dict(list(sorted_total.items())[:n])

    values = list(limited_dic.values())
    keys = [*limited_dic]
    samples = [count[0] for count in values]
    patients = [count[1] for count in values]

    x = np.arange(len(keys))
    width = 0.3

    fig, ax = plt.subplots(figsize=(10, 6))
    
    rects1 = ax.bar(x - width/2, samples, width, label='Samples', color='#008080')
    rects2 = ax.bar(x + width/2, patients, width, label='Patients', color='#FF6F61')


    ax.set_ylabel('Count')
    ax.set_xticks(x)
    ax.set_xticklabels(keys, rotation=45, ha='right') 
    ax.legend()

    autolabel_ver(ax, rects1)
    autolabel_ver(ax, rects2)

    fig.tight_layout()

    plt.show()
    plt.savefig(os.path.join(output_folder, "img", 'general.png'))


    total_genes = {}

    for folder in old_versions:
        data_snv_path = os.path.join(os.path.realpath(outputfolderpath), folder, "data_mutations_extended.txt")
        data_cnv_path = os.path.join(os.path.realpath(outputfolderpath), folder, "data_cna.txt")
        data_sv_path = os.path.join(os.path.realpath(outputfolderpath), folder, "data_sv.txt")

        if os.path.exists(data_snv_path):
            data_snv_df = pd.read_csv(data_snv_path, sep="\t", header=0)
            snv_number = len(data_snv_df)
        else:
            snv_number = 0

        if os.path.exists(data_cnv_path):
            data_cnv_df = pd.read_csv(data_cnv_path, sep="\t", header=0)
            cnv_number = len(data_cnv_df)
        else:
            cnv_number = 0

        if os.path.exists(data_sv_path):
            data_sv_df = pd.read_csv(data_sv_path, sep="\t", header=0)
            sv_number = len(data_sv_df)
        else:
            sv_number = 0

        total_genes[folder] = [snv_number, cnv_number, sv_number]


    sorted_total = dict(sorted(total_genes.items(), key=lambda item: int(re.search(r'_v(\d+)$', item[0]).group(1))))
    limited_dic = dict(list(sorted_total.items())[:n])

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
        gridspec_kw={'wspace': 0.1}, 
        constrained_layout=True
    )

    for i, (cat, values) in enumerate(categories.items()):
        axes[i].barh(studies, values, color=colors[cat])
        axes[i].set_title(cat, fontsize=10)
        axes[i].set_xlabel("Counts", fontsize=9)
        if i == 0:  
            axes[i].invert_yaxis()

    plt.savefig(os.path.join(output_folder, "img", 'genes.png'))






    











    # # Organizza i dati per categorie (SNV, CNV, SV)
    # categories = ['SNV', 'CNV', 'SV']
    # category_values = [snv, cnv, sv]

    # # Imposta il colore per ciascun studio (basato su un colormap)
    # colors = plt.cm.tab10(np.linspace(0, 1, len(studies)))

    # # Crea il grafico
    # fig, ax = plt.subplots(figsize=(10, 6))

    # # Larghezza delle barre
    # barWidth = 0.25

    # # Posizioni delle barre lungo l'asse delle x per ciascuna categoria
    # r = np.arange(len(studies))

    # # Creazione delle barre per ciascuna categoria
    # for i, (category, values) in enumerate(zip(categories, category_values)):
    #     ax.barh(r + i * barWidth, values, color=colors[i], edgecolor='grey', height=barWidth, label=category)

    # # Imposta le etichette dell'asse y (le categorie)
    # ax.set_yticks(r + barWidth)
    # ax.set_yticklabels(studies)

    # # Etichette e titolo
    # ax.set_xlabel('Total Count')
    # ax.set_title('Comparison of SNV, CNV, and SV Across Studies')

    # # Aggiungi la legenda
    # ax.legend()

    # # Aggiungi le etichette sui valori delle barre
    # for i, category in enumerate(categories):
    #     for j in range(len(studies)):
    #         ax.text(category_values[i][j] + 10, r[j] + i * barWidth, str(category_values[i][j]), va='center')

    # # Mostra il grafico
    # plt.tight_layout()    
    # plt.savefig(os.path.join(output_folder, "img", 'genes.png'))

def autolabel_ver(ax, rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom')