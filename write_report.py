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

import ast
import os
import re
import shutil
import sys
from configparser import ConfigParser
from datetime import datetime
from operator import index

import pandas as pd

import versioning
import walk
from filter_clinvar import check_bool

config = ConfigParser()
config.read("conf.ini")

annotations_list = ast.literal_eval(config.get("Annotations", "ANNOTATIONS"))

def extract_sample_list(filecase):
    with open(filecase, "r") as meta:
        for line in meta:
            if line.startswith("case_list_ids:"):
               sample_part = line.split(":")[1]
               samples = sample_part.split("\t")
               sample_list = [sample.strip() for sample in samples]
    return sample_list


def ghost_sample(output_folder):
    new_samples = set()

    for file in ["data_mutations_extended.txt", "data_cna.txt", "data_sv.txt"]:
        new_samples = get_samples(file, new_samples, output_folder)        

    all_path = os.path.join(output_folder, "data_clinical_sample.txt")
    if os.path.exists(all_path):
        all_df = pd.read_csv(all_path, sep="\t", header=4)
        all_samples = set(all_df["SAMPLE_ID"])
    else:
        all_samples = set()

    ghosts = all_samples - new_samples

    return list(ghosts)


def get_samples(file, sample_list, output_folder):
    path = os.path.join(output_folder, file)
    if os.path.exists(path):
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if file == "data_mutations_extended.txt":
            samples = set(df["Tumor_Sample_Barcode"])
        elif file == "data_cna.txt":
            samples = set(df.columns[1:])
        elif file == "data_sv.txt":
            samples = set(df["Sample_Id"])
        sample_list = sample_list.union(samples)

    return sample_list


###########################
#           Main          #
###########################

def write_report_main(output_folder, cancer, oncoKB, filters, number_for_graph):
    os.system("cp " + "styles.css" + " " + os.path.join(output_folder, "img", "styles.css"))
    now = datetime.now()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")

    img_path = os.path.join("img", "logo_VARAN.png")
    general_graph_path = os.path.join("img", "general.png")
    genes_graph_path = os.path.join("img", "genes.png")

    my_filters = write_filters_report(output_folder)
    cancer = cancer.capitalize()

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > 4:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(output_folder)

    clin_sam_path = os.path.join(output_folder, "data_clinical_sample.txt")
    if os.path.exists(clin_sam_path):
        clin_sam = pd.read_csv(clin_sam_path, sep="\t", header=4)

        new_samples = set(clin_sam["SAMPLE_ID"])
        new_patients = set(clin_sam["PATIENT_ID"])
        new_smpl_nr = len(new_samples)
        new_pt_nr = len(new_patients)

    if versioning.old_version_exists:
        name = re.search(r"^(.+_v)[0-9]+$", os.path.basename(output_folder)).group(1)
        actual_version = int(re.search(r"^.+_v([0-9]+)$", os.path.basename(output_folder)).group(1))

        old_name = name + str(actual_version - 1)
        old_file = os.path.join(os.path.dirname(output_folder), old_name, "data_clinical_sample.txt")

        if actual_version != 1:
            if os.path.exists(old_file):
                old_clin_sam = pd.read_csv(old_file, sep="\t", header=4)

                old_samples = set(old_clin_sam["SAMPLE_ID"])
                old_patients = set(old_clin_sam["PATIENT_ID"])

                only_new_sam = new_samples - old_samples
                only_new_pat = new_patients - old_patients

                only_old_sam = old_samples - new_samples
                only_old_pat = old_patients - new_patients

            else:
                only_new_sam = new_samples
                only_new_pat = new_patients

                only_old_sam = []
                only_old_pat = []

    html_content = f"""
    <!DOCTYPE html>
    <html lang="it">
    <head>
        <link rel="stylesheet" type="text/css" href="{os.path.join("img", "styles.css")}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan Report</title>
    </head>
    <body>
        <header>
            <img src="{img_path}" alt="Logo Varan">
            <h1>VARAN</h1>
        </header>

        <h2>Report generate on {date}</h2>

        <div class="container">
            <section class="general-info">
                <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>STUDY NAME:</strong> {os.path.basename(os.path.normpath(output_folder))}</p>
                    <p><strong>CANCER TYPE:</strong> {cancer}</p>
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>TOTAL SAMPLE(S):</strong> {new_smpl_nr}</p>
                    <p><strong>TOTAL PATIENT(S):</strong> {new_pt_nr}</p>
                </div>"""


    if ghosts:
        html_content += f"""
        <div class="content">
            <p><strong><span>&#9888;</span></strong> The following samples are not present in cnv, snv and fusions after filtering: {ghosts}
        </div>"""

    html_content += """</section>"""

    if versioning.old_version_exists:
        if actual_version != 1:
            html_content += f"""
            <section class="comparison">
                <div class="section-title">Comparison with Previous Version (v_{actual_version - 1})</div>
                <div class="content">
                    <p><strong>ADDED PATIENT(S):</strong> {len(only_new_pat)}</p>
                    <p><strong>ADDED SAMPLE(S):</strong> {len(only_new_sam)}</p>
                    <p><strong>REMOVED PATIENT(S):</strong> {len(only_old_pat)}</p>
                    <p><strong>REMOVED SAMPLE(S):</strong> {len(only_old_sam)}</p>
                </div>
            </section>
            """

    html_content += f"""
        <section class="filters">
            <div class="section-title">Filters & Configuration</div>"""

    if any(letter in filters for letter in "d"):
        html_content += f"""
            <div class="subtitle">VCF Filters</div>"""

    if "d" in filters:
        html_content += f"""
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(letter in filters for letter in "va"):
        html_content += f"""
            <div class="subtitle">VAF Filters</div>"""
    if "v" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>T_VAF_MIN:</strong> {extract_key_value(my_filters, "T_VAF_MIN")}</p>
            <p><strong>T_VAF_MAX:</strong> {extract_key_value(my_filters, "T_VAF_MAX")}</p>
        </div>"""
        if "n" in filters:
            html_content += f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL:</strong> {extract_key_value(my_filters, "T_VAF_MIN_NOVEL")}</p>
            </div>"""

    if "a" in filters:
        drop_NA_AF = config.get("Filters", "drop_NA_AF")
        drop_NA_AF = check_bool(drop_NA_AF)
        if drop_NA_AF:
            excl_or_incl = "exclude"
        else:
            excl_or_incl = "include"
        html_content += f"""
                <div class="content">
                <p><strong>AF:</strong> {extract_key_value(my_filters, "AF")} & {excl_or_incl} NA</p>
            </div>"""

    if any(letter in filters for letter in "oqycbi"):
        html_content += f"""
            <div class="subtitle">MAF Filters</div>"""

    if oncoKB and "o" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>ONCOKB:</strong> include {", ".join([item.strip() for item in extract_key_value(my_filters, "ONCOKB_FILTER").strip("[]").replace('"', "").split(",")])}</p>
            </div>"""

    if "q" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>CONSEQUENCES:</strong> include {", ".join([item.strip() for item in extract_key_value(my_filters, "CONSEQUENCES").strip("[]").replace('"', "").split(",")])}</p>
            </div>"""

    if "y" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>POLYPHEN:</strong> include {", ".join([item.strip() for item in extract_key_value(my_filters, "POLYPHEN").strip("[]").replace('"', "").split(",")])}</p>
            </div>"""

    if "c" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>CLIN_SIG:</strong> exclude {", ".join([item.strip() for item in extract_key_value(my_filters, "CLIN_SIG").strip("[]").replace('"', "").split(",")])}</p>
            </div>"""

    if "i" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>IMPACT:</strong> exclude {", ".join([item.strip() for item in extract_key_value(my_filters, "IMPACT").strip("[]").replace('"', "").split(",")])}</p>
            </div>"""

    if "s" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>SIFT:</strong> include {", ".join([item.strip() for item in extract_key_value(my_filters, "SIFT").strip("[]").replace('"', "").split(",")])}</p>
            </div>"""

    if "p" in filters:
        html_content += f"""
                <div class="content">
                 <p><strong>FILTER:</strong> = PASS
            </div>"""

    TMB_section = re.sub(r"[{}']", "", extract_section(my_filters, "TMB"))
    TMB_section = re.sub(r"([,:])(?=\S)", r"\1 ", TMB_section)

    html_content += f"""
            <div class="subtitle">Copy Number Alterations (CNA)</div>
            <div class="content">
                {extract_section(my_filters, "Cna")}
            </div>

            <div class="subtitle">Tumor Mutational Burden (TMB)</div>
            <div class="content">
                {TMB_section}
            </div>

            <div class="subtitle">Microsatellite Instability (MSI)</div>
            <div class="content">
                {extract_section(my_filters, "MSI")}
            </div>

            <div class="subtitle">Fusions</div>
            <div class="content">
                {extract_section(my_filters, "FUSION")}
            </div>
        </section>"""

    if os.path.exists(os.path.join(output_folder, general_graph_path)) or os.path.exists(os.path.join(output_folder, genes_graph_path)):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">Graphical Overview {graph_expression}</div>"""

    if os.path.exists(os.path.join(output_folder, general_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""


    if os.path.exists(os.path.join(output_folder, genes_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if os.path.exists(os.path.join(output_folder, general_graph_path)) or os.path.exists(os.path.join(output_folder, genes_graph_path)):
        html_content += f"""</section>"""

    if annotations_list != []:
        html_content += f"""
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += f"""
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += f"""
    </body>
    </html>
    """

    with open(os.path.join(output_folder, "report_VARAN.html"), "w") as f:
        f.write(html_content)


def write_filters_report(output_folder):
    sections_to_include = {
            "Filters": ["BENIGN", "CLIN_SIG", "CONSEQUENCES", "ONCOKB_FILTER", 
                        "t_VAF_min", "t_VAF_min_novel", "t_VAF_max", 
                        "AF", "POLYPHEN", "IMPACT", "SIFT"],
            "Cna": ["HEADER_CNV", "PLOIDY", "CNVkit"],
            "TMB": ["THRESHOLD_TMB"],
            "MSI": ["THRESHOLD_SITES", "THRESHOLD_MSI"],
            "FUSION": ["THRESHOLD_FUSION"]
        }

    conf_content = []

    for section, keys in sections_to_include.items():
        if section in config:
            conf_content.append(f"[{section}]")  
            for key in keys:
                if key in config[section]:  
                    value = config[section][key]
                    conf_content.append(f"{key.upper()} = {value}")

    return conf_content


def extract_section(content, section_name):
    content = "\n".join(content)
    match = re.search(rf"\[{section_name}\](.*?)(?=\n\[[^\]]+\]|\Z)", content, re.DOTALL | re.IGNORECASE)

    if not match:
        return "N/A"
    section_content = match.group(1).strip()
    formatted_content = re.sub(r"^\s*(\S+)\s*=\s*", r"<strong>\1</strong> = ", section_content, flags=re.MULTILINE)

    return formatted_content.replace("\n", "<br>")


def extract_key_value(filters, key_name):
    pattern = rf"{key_name}\s*=\s*(.+)"
    for line in filters:
        match = re.search(pattern, line)
        if match:
            return match.group(1).strip().strip('"').strip("'")
    return "N/A"


###########################
#          Update         #
###########################

def write_report_update(original_study, updating_with, new_study, number_for_graph):
    old_img_path = os.path.join(original_study, "img", "logo_VARAN.png")
    general_graph_path = os.path.join("img", "general.png")
    genes_graph_path = os.path.join("img", "genes.png")

    if os.path.exists(old_img_path):
        img_output_dir = os.path.join(new_study, "img")
        os.makedirs(img_output_dir, exist_ok=True)
        new_img_path = os.path.join(img_output_dir, "logo_VARAN.png")
        shutil.copy(old_img_path, new_img_path)

    os.system("cp " + "styles.css" + " " + os.path.join(new_study, "img", "styles.css"))

    now = datetime.now()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > 4:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(new_study)

    output_file = os.path.join(new_study, "report_VARAN.html")

    case_list1 = os.path.join(original_study, "case_lists")
    case_list2 = os.path.join(updating_with, "case_lists")
    cna_1 = os.path.join(case_list1, "cases_cna.txt")
    cna_2 = os.path.join(case_list2, "cases_cna.txt")
    sequenced_1 = os.path.join(case_list1, "cases_sequenced.txt")
    sequenced_2 = os.path.join(case_list2, "cases_sequenced.txt")
    sv_1 = os.path.join(case_list1, "cases_sv.txt")
    sv_2 = os.path.join(case_list2, "cases_sv.txt")

    clin_sam_old_df = pd.read_csv(os.path.join(original_study, "data_clinical_sample.txt"), sep="\t")
    clin_sam_new_df = pd.read_csv(os.path.join(updating_with, "data_clinical_sample.txt"), sep="\t")

    updated_clin_pat = set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
    updated_clin_sample = set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])
    added_clin_pat = set(clin_sam_new_df.iloc[4:, 1]).difference(clin_sam_old_df.iloc[4:, 1])
    added_clin_sample = set(clin_sam_new_df.iloc[4:, 0]).difference(clin_sam_old_df.iloc[4:, 0])

    updated_samples_cna, added_samples_cna, total_patients, total_samples = compare_sample_file_update(cna_1, cna_2, new_study)
    updated_samples_sequenced, added_samples_sequenced, _, _ = compare_sample_file_update(sequenced_1, sequenced_2, new_study)
    updated_samples_sv, added_samples_sv, _, _ = compare_sample_file_update(sv_1, sv_2, new_study)

    old_report = os.path.join(original_study, "report_VARAN.html")
    updating_report = os.path.join(updating_with, "report_VARAN.html")

    if os.path.exists(old_report) and os.path.exists(updating_report):
        filters1 = extract_filters_from_html(old_report)
        filters2 = extract_filters_from_html(updating_report)
        cancer_type1 = extract_cancer_type_from_html(old_report)
        cancer_type2 = extract_cancer_type_from_html(updating_report)
    else:
        filters1 = dict()
        filters2 = dict()
        cancer_type1 = None
        cancer_type2 = None

    order = ["T_VAF_MIN", "T_VAF_MIN_NOVEL", "T_VAF_MAX", "AF", "ONCOKB", "IMPACT", "CLIN_SIG", "CONSEQUENCES", "POLYPHEN", "SIFT", \
    "HEADER_CNV", "PLOIDY", "CNVKIT", "THRESHOLD_TMB", "THRESHOLD_SITES", "THRESHOLD_MSI", "THRESHOLD_FUSION"]

    changed_filters = []

    for filter_name in order:
        if filter_name not in (filters1.keys() | filters2.keys()):
            continue 
        value1 = filters1.get(filter_name, "Not Present")
        value2 = filters2.get(filter_name, "Not Present")
        if value1 != value2:
            changed_filters.append(filter_name)

    common_filters = list(set(order) - set(changed_filters))
    common_filters = list(set(common_filters) & set(filters1.keys()))

    html_content = f"""
    <!DOCTYPE html>
    <html lang="it">
    <head>
        <link rel="stylesheet" type="text/css" href="{os.path.join("img", "styles.css")}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan - Update</title>
    </head>
    <body>
        <header>
            <img src="{os.path.join("img", "logo_VARAN.png")}" alt="Logo Varan">
            <h1>VARAN - Update</h1>
        </header>

        <h2>Report generate on {date}</h2>
        <div class="container">
            <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>ORIGINAL STUDY:</strong> {os.path.basename(os.path.normpath(original_study))}</p>
                    <p><strong>UPDATING WITH:</strong> {os.path.basename(os.path.normpath(updating_with))}</p>
                    <p><strong>NEW STUDY:</strong> {os.path.basename(os.path.normpath(new_study))}</p>"""

    if (cancer_type1 == cancer_type2) and (cancer_type1 != None):
        html_content += f"""<p><strong>CANCER TYPE:</strong> {cancer_type1}</p>"""
    else:
        html_content += f"""<p><strong>CANCER TYPE:</strong> Mixed</p>"""

    html_content += f"""
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>Total Patients:</strong> {total_patients}</p>
                    <p><strong>Total Samples:</strong> {total_samples}</p>
                </div>
    """

    if ghosts:
        html_content += f"""
        <div class="content">
            <p><strong><span>&#9888;</span></strong> The following samples are not present in cnv, snv and fusions after filtering: {ghosts}
        </div>"""

    if updated_clin_sample:
        html_content += f"""
            <div class="content">
            <p><strong>UPDATED:</strong></p>
            <p>&emsp;<strong>{len(updated_clin_pat)} Patients:</strong> {", ".join(updated_clin_pat)}</p>
            <p>&emsp;<strong>{len(updated_clin_sample)}Samples:</strong> {", ".join(updated_clin_sample)}</p>
            </div>
            """

    if added_clin_sample:
        html_content += f"""
            <div class="content">
            <p><strong>ADDED:</strong></p>
            <p>&emsp;<strong>{len(added_clin_pat)} Patients:</strong> {", ".join(added_clin_pat)}</p>
            <p>&emsp;<strong>{len(added_clin_pat)} Samples:</strong> {", ".join(added_clin_sample)}</p>
            </div>
            """

    if updated_samples_cna or updated_samples_sequenced or updated_samples_sv or added_samples_cna or added_samples_sequenced or added_samples_sv:
        html_content += f"""
        <div class="section-title">Detailed Overview</div>
            <div class="content">"""

    if updated_samples_cna:
        html_content += f"""
        <p><strong>Cases_CNA:</strong></p>
        <p>&emsp;<strong>Updated:</strong> {len(updated_samples_cna)} samples ({", ".join(updated_samples_cna)})</p>
        <p>&emsp;<strong>Added:</strong> {len(added_samples_cna)} samples ({", ".join(added_samples_cna)})</p>"""

    if updated_samples_sequenced:
        html_content += f"""
        <p><strong>Cases_sequenced:</strong></p>
        <p>&emsp;<strong>Updated:</strong> {len(updated_samples_sequenced)} samples ({", ".join(updated_samples_sequenced)})</p>
        <p>&emsp;<strong>Added:</strong> {len(added_samples_sequenced)} samples ({", ".join(added_samples_sequenced)})</p>"""

    if updated_samples_sv:
        html_content += f"""
        <p><strong>Cases_sv:</strong></p>
        <p>&emsp;<strong>Updated:</strong> {len(updated_samples_sv)} samples ({", ".join(updated_samples_sv)})</p>
        <p>&emsp;<strong>Added:</strong> {len(added_samples_sv)} samples ({", ".join(added_samples_sv)})</p>"""

    if updated_samples_cna or updated_samples_sequenced or updated_samples_sv or added_samples_cna or added_samples_sequenced or added_samples_sv:
        html_content += f"""
        </div>
        """

    if common_filters != []:
        html_content += f"""
            <section class="filters">
                <div class="section-title">Filters & Configuration</div>"""

    if any(filter in common_filters for filter in ["ALT"]):
        html_content += f"""
            <div class="subtitle">VCF Filters</div>"""

    if "ALT" in common_filters:
        html_content += f"""
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(filter in common_filters for filter in ["T_VAF_MIN", "T_VAF_MAX", "AF"]):
        html_content += f"""
            <div class="subtitle">VAF Filters</div>"""

    if any(filter in common_filters for filter in ["T_VAF_MIN", "T_VAF_MAX"]):
        html_content += f"""
            <div class="content">
            <p><strong>T_VAF_MIN</strong>: {filters1["T_VAF_MIN"]}</p>
            <p><strong>T_VAF_MAX</strong>: {filters1["T_VAF_MAX"]}</p>
        </div>"""

        if "T_VAF_MIN_NOVEL" in common_filters:
            html_content += f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL</strong>: {filters1["T_VAF_MIN_NOVEL"]}</p>
            </div>"""

    if "AF" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>AF</strong>: {filters1["AF"]}</p>
            </div>"""

    if any(filter in common_filters for filter in ["ONCOKB", "CONSEQUENCES", "POLYPHEN", "CLIN_SIG", "IMPACT", "SIFT", "FILTER"]):
        html_content += f"""
            <div class="subtitle">MAF Filters</div>"""

    if "ONCOKB" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>ONCOKB</strong>: {filters1["ONCOKB"]}</p>
            </div>"""

    if "CONSEQUENCES" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>CONSEQUENCES</strong>: {filters1["CONSEQUENCES"]}</p>
            </div>"""

    if "POLYPHEN" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>POLYPHEN</strong>: {filters1["POLYPHEN"]}</p>
            </div>"""

    if "CLIN_SIG" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>CLIN_SIG</strong>: {filters1["CLIN_SIG"]}</p>
            </div>"""

    if "IMPACT" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>IMPACT</strong>: {filters1["IMPACT"]}</p>
            </div>"""

    if "SIFT" in common_filters:
        html_content += f"""
                <div class="content">
                <p><strong>SIFT</strong>: {filters1["SIFT"]}</p>
            </div>"""

    if "FILTER" in common_filters:
        html_content += f"""
                <div class="content">
                 <p><strong>FILTER</strong>: {filters1["FILTER"]}</p>
            </div>"""

    keys_to_check = {"HEADER_CNV", "PLOIDY", "CNVKIT", "THRESHOLD_TMB", "THRESHOLD_SITES", "THRESHOLD_MSI", "THRESHOLD_FUSION"}
    if common_filters != {}:
        if all(key in filters1 for key in keys_to_check):
            html_content += f"""
                    <div class="subtitle">Copy Number Alterations (CNA)</div>
                    <div class="content">
                        <p><strong>HEADER_CNV</strong> = {ast.literal_eval(filters1["HEADER_CNV"])}</p>
                        <p><strong>PLOIDY</strong> = {filters1["PLOIDY"]}</p>
                        <p><strong>CNVKIT</strong> = {filters1["CNVKIT"]}</p>
                    </div>

                    <div class="subtitle">Tumor Mutational Burden (TMB)</div>
                    <div class="content">
                        <p><strong>THRESHOLD_TMB</strong>: {filters1["THRESHOLD_TMB"]}</p>
                    </div>

                    <div class="subtitle">Microsatellite Instability (MSI)</div>
                    <div class="content">
                        <p><strong>THRESHOLD_SITES</strong>: {filters1["THRESHOLD_SITES"]}</p>
                        <p><strong>THRESHOLD_MSI</strong>: {filters1["THRESHOLD_MSI"]}</p>
                    </div>

                    <div class="subtitle">Fusions</div>
                    <div class="content">
                        <p><strong>THRESHOLD_FUSION</strong>: {filters1["THRESHOLD_FUSION"]}</p>
                    </div>
                </section>
            </section>"""

    if any(filters1.get(f) != filters2.get(f) for f in filters1.keys() | filters2.keys()):
        html_content += f"""
            <div class="section-title">Differences in Filters & Configurations</div>
            <div class="content">
            <br />
            <table>
                <thead>
                    <tr>
                        <th>Filter</th>
                        <th>{os.path.basename(os.path.normpath(original_study))}</th>
                        <th>{os.path.basename(os.path.normpath(updating_with))}</th>
                    </tr>
                </thead>
                <tbody>
        """

        for filter_name in order:
            if filter_name in changed_filters:
                value1 = filters1.get(filter_name, "Not Present")
                value2 = filters2.get(filter_name, "Not Present")
                if value1 != value2:
                    html_content += f"""
                        <tr>
                            <td>{filter_name}</td>
                            <td>{value1}</td>
                            <td>{value2}</td>
                        </tr>"""
        html_content += """
                </tbody>
            </table>
            </div>
        """

    if os.path.exists(os.path.join(new_study, general_graph_path)) or os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">Graphical Overview {graph_expression}</div>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""

    if os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)) or os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
            </section>"""

    if annotations_list != []:
        html_content += f"""
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += f"""
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += f"""
    </body>
    </html>
    """

    with open(output_file, "w") as f:
        f.write(html_content)


def compare_sample_file_update(file1, file2, outputfolder):
    clin_sam_output_df = pd.read_csv(os.path.join(outputfolder, "data_clinical_sample.txt"), sep="\t")

    if os.path.exists(file1) and os.path.exists(file2):
        samples_file1 = extract_sample_list(file1)  
        samples_file2 = extract_sample_list(file2)

        updated_samples = [sample for sample in samples_file1 if sample in samples_file2 and sample != ""]
        new_samples = [sample for sample in samples_file2 if not sample in samples_file1 and sample != ""]

    else:
        updated_samples = None
        new_samples = None

    return updated_samples, new_samples, len(set(clin_sam_output_df.iloc[4:, 1])), len(set(clin_sam_output_df.iloc[4:, 0]))


###########################
#          Extract        #
###########################

def write_report_extract(original_study, new_study, number_for_graph):
    old_img_path = os.path.join(original_study, "img", "logo_VARAN.png")
    general_graph_path = os.path.join("img", "general.png")
    genes_graph_path = os.path.join("img", "genes.png")

    if os.path.exists(old_img_path):
        img_output_dir = os.path.join(new_study, "img")
        os.makedirs(img_output_dir, exist_ok=True)
        new_img_path = os.path.join(img_output_dir, "logo_VARAN.png")
        shutil.copy(old_img_path, new_img_path)   

    os.system("cp " + "styles.css" + " " + os.path.join(new_study, "img", "styles.css"))

    now = datetime.now()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > 4:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(new_study)

    old_report = os.path.join(original_study, "report_VARAN.html")
    if os.path.exists(old_report):
        filters = extract_filters_from_html(old_report)
        cancer_type = extract_cancer_type_from_html(old_report)
    else:
        filters = dict()
        cancer_type = None

    case_list1 = os.path.join(original_study, "case_lists")
    case_list2 = os.path.join(new_study, "case_lists")
    cna_1 = os.path.join(case_list1, "cases_cna.txt")
    cna_2 = os.path.join(case_list2, "cases_cna.txt")
    sequenced_1 = os.path.join(case_list1, "cases_sequenced.txt")
    sequenced_2 = os.path.join(case_list2, "cases_sequenced.txt")
    sv_1 = os.path.join(case_list1, "cases_sv.txt")
    sv_2 = os.path.join(case_list2, "cases_sv.txt")

    extracted_samples_cna, total_patients, total_samples = compare_sample_file_extract(cna_1, cna_2, original_study, new_study)
    extracted_samples_sequenced, _, _ = compare_sample_file_extract(sequenced_1, sequenced_2, original_study, new_study)
    extracted_samples_sv, _, _ = compare_sample_file_extract(sv_1, sv_2, original_study, new_study)

    html_content = f"""
    <!DOCTYPE html>
    <html lang="it">
    <head>
        <link rel="stylesheet" type="text/css" href="{os.path.join("img", "styles.css")}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan - Extract</title>
    </head>
    <body>
        <header>
            <img src="{os.path.join("img", "logo_VARAN.png")}" alt="Logo Varan">
            <h1>VARAN - Extract</h1>
        </header>

        <h2>Generate on {date}</h2>

        <div class="container">
            <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>ORIGINAL STUDY:</strong> {os.path.basename(os.path.normpath(original_study))}</p>
                    <p><strong>NEW STUDY:</strong> {os.path.basename(os.path.normpath(new_study))}</p>"""

    if cancer_type:
        html_content += f"""<p><strong>CANCER TYPE:</strong> {cancer_type}</p>"""

    html_content += f""" 
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>Total Patients:</strong> {total_patients}</p>
                    <p><strong>Total Samples:</strong> {total_samples}</p>
                </div>"""

    if ghosts:
        html_content += f"""
        <div class="content">
            <p><strong><span>&#9888;</span></strong> The following samples are not present in cnv, snv and fusions after filtering: {ghosts}
        </div>"""

    html_content += """
            <div class="section-title">Detailed Overview</div>
                <div class="content">
                    <p><strong>Cases_CNA:</strong></p>
    """

    if extracted_samples_cna:
        html_content += f"""
            <p>&emsp;<strong>Extracted:</strong> {len(extracted_samples_cna)} samples ({", ".join(extracted_samples_cna)})</p>
        """
    else:
        html_content += f"""
            <p>&emsp;None of the extracted samples was in the original study ({os.path.basename(os.path.normpath(original_study))})</p>
        """

    html_content += f"""
        <p><strong>Cases_sequenced:</strong></p>
        """

    if extracted_samples_sequenced:
        html_content += f"""
            <p>&emsp;<strong>Extracted:</strong> {len(extracted_samples_sequenced)} samples ({", ".join(extracted_samples_sequenced)})</p>
        """
    else:
        html_content += f"""
            <p>&emsp;None of the extracted samples was in the original study ({os.path.basename(os.path.normpath(original_study))})</p>
        """

    html_content += f"""
        <p><strong>Cases_sv:</strong></p>
        """

    if extracted_samples_sv:
        html_content += f"""
            <p>&emsp;<strong>Extracted:</strong> {len(extracted_samples_sv)} samples ({", ".join(extracted_samples_sv)})</p>
        """
    else:
        html_content += f"""
            <p>&emsp;None of the extracted samples was in the original study ({os.path.basename(os.path.normpath(original_study))})</p>
        """

    if filters != {}:
        html_content += f"""
            <section class="filters">
                <div class="section-title">Filters & Configuration</div>"""

    if any(filter in filters.keys() for filter in ["ALT"]):
        html_content += f"""
            <div class="subtitle">VCF Filters</div>"""

    if "ALT" in filters.keys():
        html_content += f"""
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(filter in filters.keys() for filter in ["T_VAF_MIN", "T_VAF_MAX", "AF"]):
        html_content += f"""
            <div class="subtitle">VAF Filters</div>"""

    if any(filter in filters.keys() for filter in ["T_VAF_MIN", "T_VAF_MAX"]):
        html_content += f"""
            <div class="content">
            <p><strong>T_VAF_MIN</strong>: {filters["T_VAF_MIN"]}</p>
            <p><strong>T_VAF_MAX</strong>: {filters["T_VAF_MAX"]}</p>
        </div>"""
        if "T_VAF_MIN_NOVEL" in filters.keys():
            html_content += f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL</strong>: {filters["T_VAF_MIN_NOVEL"]}</p>
            </div>"""

    if "AF" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>AF</strong>: {filters["AF"]}</p>
            </div>"""

    if any(filter in filters.keys() for filter in ["ONCOKB", "CONSEQUENCES", "POLYPHEN", "CLIN_SIG", "IMPACT", "SIFT", "FILTER"]):
        html_content += f"""
            <div class="subtitle">MAF Filters</div>"""

    if "ONCOKB" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>ONCOKB</strong>: {filters["ONCOKB"]}</p>
            </div>"""

    if "CONSEQUENCES" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>CONSEQUENCES</strong>: {filters["CONSEQUENCES"]}</p>
            </div>"""

    if "POLYPHEN" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>POLYPHEN</strong>: {filters["POLYPHEN"]}</p>
            </div>"""

    if "CLIN_SIG" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>CLIN_SIG</strong>: {filters["CLIN_SIG"]}</p>
            </div>"""

    if "IMPACT" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>IMPACT</strong>: {filters["IMPACT"]}</p>
            </div>"""

    if "SIFT" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>SIFT</strong>: {filters["SIFT"]}</p>
            </div>"""

    if "FILTER" in filters.keys():
        html_content += f"""
                <div class="content">
                 <p><strong>FILTER</strong>: {filters["FILTER"]}</p>
            </div>"""

    if filters != {}:
        html_content += f"""
                <div class="subtitle">Copy Number Alterations (CNA)</div>
                <div class="content">
                    <p><strong>HEADER_CNV</strong> = {ast.literal_eval(filters["HEADER_CNV"])}</p>
                    <p><strong>PLOIDY</strong> = {filters["PLOIDY"]}</p>
                    <p><strong>CNVKIT</strong> = {filters["CNVKIT"]}</p>
                </div>

                <div class="subtitle">Tumor Mutational Burden (TMB)</div>
                <div class="content">
                    <p><strong>THRESHOLD_TMB</strong>: {filters["THRESHOLD_TMB"]}</p>
                </div>

                <div class="subtitle">Microsatellite Instability (MSI)</div>
                <div class="content">
                    <p><strong>THRESHOLD_SITES</strong>: {filters["THRESHOLD_SITES"]}</p>
                    <p><strong>THRESHOLD_MSI</strong>: {filters["THRESHOLD_MSI"]}</p>
                </div>

                <div class="subtitle">Fusions</div>
                <div class="content">
                    <p><strong>THRESHOLD_FUSION</strong>: {filters["THRESHOLD_FUSION"]}</p>
                </div>
            </section>"""

        if filters != {}:
            html_content += f"""
            </section>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)) or os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">Graphical Overview {graph_expression}</div>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""


    if os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)) or os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""</section>"""

    if annotations_list != []:
        html_content += f"""
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += f"""
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += f"""
    </body>
    </html>
    """

    with open(os.path.join(new_study, "report_VARAN.html"), "w") as f:
        f.write(html_content)


def compare_sample_file_extract(file1, file2, input_folder, outputfolder):
    clin_sam_old_df = pd.read_csv(os.path.join(input_folder, "data_clinical_sample.txt"), sep="\t")
    clin_sam_new_df = pd.read_csv(os.path.join(outputfolder, "data_clinical_sample.txt"), sep="\t")

    clin_pat = set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
    clin_sample = set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])

    if os.path.exists(file1) and os.path.exists(file2):
        samples_file1 = extract_sample_list(file1)  
        samples_file2 = extract_sample_list(file2)

        extracted_samples = [sample for sample in samples_file2 if sample in samples_file1 and sample != ""]
    else:
        extracted_samples = None

    return extracted_samples, len(set(clin_sam_new_df.iloc[4:, 1])), len(set(clin_sam_new_df.iloc[4:, 0]))


###########################
#          Remove         #
###########################

def write_report_remove(original_study, new_study, number_for_graph):
    old_img_path = os.path.join(original_study, "img", "logo_VARAN.png")
    general_graph_path = os.path.join("img", "general.png")
    genes_graph_path = os.path.join("img", "genes.png")

    if os.path.exists(old_img_path):
        img_output_dir = os.path.join(new_study, "img")
        os.makedirs(img_output_dir, exist_ok=True)
        new_img_path = os.path.join(img_output_dir, "logo_VARAN.png")
        shutil.copy(old_img_path, new_img_path)

    os.system("cp " + "styles.css" + " " + os.path.join(new_study, "img", "styles.css"))

    now = datetime.now()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > 4:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(new_study)

    old_report = os.path.join(original_study, "report_VARAN.html")
    if os.path.exists(old_report):
        filters = extract_filters_from_html(old_report)
        cancer_type = extract_cancer_type_from_html(old_report)
    else:
        filters = dict()
        cancer_type = None

    case_list1 = os.path.join(original_study, "case_lists")
    case_list2 = os.path.join(new_study, "case_lists")
    cna_1 = os.path.join(case_list1, "cases_cna.txt")
    cna_2 = os.path.join(case_list2, "cases_cna.txt")
    sequenced_1 = os.path.join(case_list1, "cases_sequenced.txt")
    sequenced_2 = os.path.join(case_list2, "cases_sequenced.txt")
    sv_1 = os.path.join(case_list1, "cases_sv.txt")
    sv_2 = os.path.join(case_list2, "cases_sv.txt")

    removed_samples_cna, left_samples_cna, total_patients, total_samples = compare_sample_file_remove(cna_1, cna_2, original_study, new_study)
    removed_samples_sequenced, left_samples_sequenced, _, _ = compare_sample_file_remove(sequenced_1, sequenced_2, original_study, new_study)
    removed_samples_sv, left_samples_sv, _, _ = compare_sample_file_remove(sv_1, sv_2, original_study, new_study)

    html_content = f"""
    <!DOCTYPE html>
    <html lang="it">
    <head>
        <link rel="stylesheet" type="text/css" href="{os.path.join("img", "styles.css")}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan - Remove</title>
    </head>
    <body>
        <header>
            <img src="{os.path.join("img", "logo_VARAN.png")}" alt="Logo Varan">
            <h1>VARAN - Remove</h1>
        </header>

        <h2>Generate on {date}</h2>

        <div class="container">
            <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>ORIGINAL STUDY:</strong> {os.path.basename(os.path.normpath(original_study))}</p>
                    <p><strong>NEW STUDY:</strong> {os.path.basename(os.path.normpath(new_study))}</p>"""

    if cancer_type:
        html_content += f"""<p><strong>CANCER TYPE:</strong> {cancer_type}</p>"""

    html_content += f"""
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>Total Patients:</strong> {total_patients}</p>
                    <p><strong>Total Samples:</strong> {total_samples}</p>
                </div>"""

    if ghosts:
        html_content += f"""
        <div class="content">
            <p><strong><span>&#9888;</span></strong> The following samples are not present in cnv, snv and fusions after filtering: {ghosts}
        </div>"""

    html_content += """
            <div class="section-title">Detailed Overview</div>
                <div class="content">
                    <p><strong>Cases_CNA:</strong></p>"""

    if not removed_samples_cna:
        html_content += f"""
            <p>&emsp;{os.path.basename(os.path.normpath(original_study))}'s cases_cna was empty.</p>"""
    else:
        html_content += f"""
            <p>&emsp;<strong>Removed:</strong> {len(removed_samples_cna)} samples ({", ".join(removed_samples_cna)})</p>
            <p>&emsp;There are now {left_samples_cna} samples in cases_cna.</p>"""

    html_content += f"""
        <p><strong>Cases_sequenced:</strong></p>"""

    if not removed_samples_sequenced:
        html_content += f"""
            <p>&emsp;{os.path.basename(os.path.normpath(original_study))}'s cases_sequenced was empty.</p>"""
    else:
        html_content += f"""
            <p>&emsp;<strong>Removed:</strong> {len(removed_samples_sequenced)} samples ({", ".join(removed_samples_sequenced)})</p>
            <p>&emsp;There are now {left_samples_sequenced} samples in cases_sequenced.</p>"""

    html_content += f"""
        <p><strong>Cases_sv:</strong></p>"""

    if not removed_samples_sv:
        html_content += f"""
            <p>&emsp;{os.path.basename(os.path.normpath(original_study))}'s cases_sv was empty.</p>"""
    else:
        html_content += f"""
            <p>&emsp;<strong>Removed:</strong> {len(removed_samples_sv)} samples ({", ".join(removed_samples_sv)})</p>
            <p>&emsp;There are now {left_samples_sv} samples in cases_sv.</p>"""

    if filters != {}:
        html_content += f"""
            <section class="filters">
                <div class="section-title">Filters & Configuration</div>"""

    if any(filter in filters.keys() for filter in ["ALT"]):
        html_content += f"""
            <div class="subtitle">VCF Filters</div>"""

    if "ALT" in filters.keys():
        html_content += f"""
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(filter in filters.keys() for filter in ["T_VAF_MIN", "T_VAF_MAX", "AF"]):
        html_content += f"""
            <div class="subtitle">VAF Filters</div>"""

    if any(filter in filters.keys() for filter in ["T_VAF_MIN", "T_VAF_MAX"]):
        html_content += f"""
            <div class="content">
            <p><strong>T_VAF_MIN</strong>: {filters["T_VAF_MIN"]}</p>
            <p><strong>T_VAF_MAX</strong>: {filters["T_VAF_MAX"]}</p>
        </div>"""
        if "T_VAF_MIN_NOVEL" in filters.keys():
            html_content += f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL</strong>: {filters["T_VAF_MIN_NOVEL"]}</p>
            </div>"""

    if "AF" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>AF</strong>: {filters["AF"]}</p>
            </div>"""

    if any(filter in filters.keys() for filter in ["ONCOKB", "CONSEQUENCES", "POLYPHEN", "CLIN_SIG", "IMPACT", "SIFT", "FILTER"]):
        html_content += f"""
            <div class="subtitle">MAF Filters</div>"""

    if "ONCOKB" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>ONCOKB</strong>: {filters["ONCOKB"]}</p>
            </div>"""

    if "CONSEQUENCES" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>CONSEQUENCES</strong>: {filters["CONSEQUENCES"]}</p>
            </div>"""

    if "POLYPHEN" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>POLYPHEN</strong>: {filters["POLYPHEN"]}</p>
            </div>"""

    if "CLIN_SIG" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>CLIN_SIG</strong>: {filters["CLIN_SIG"]}</p>
            </div>"""

    if "IMPACT" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>IMPACT</strong>: {filters["IMPACT"]}</p>
            </div>"""

    if "SIFT" in filters.keys():
        html_content += f"""
                <div class="content">
                <p><strong>SIFT</strong>: {filters["SIFT"]}</p>
            </div>"""

    if "FILTER" in filters.keys():
        html_content += f"""
                <div class="content">
                 <p><strong>FILTER</strong>: {filters["FILTER"]}</p>
            </div>"""

    if filters != {}:
        html_content += f"""
                <div class="subtitle">Copy Number Alterations (CNA)</div>
                <div class="content">
                    <p><strong>HEADER_CNV</strong> = {ast.literal_eval(filters["HEADER_CNV"])}</p>
                    <p><strong>PLOIDY</strong> = {filters["PLOIDY"]}</p>
                    <p><strong>CNVKIT</strong> = {filters["CNVKIT"]}</p>
                </div>

                <div class="subtitle">Tumor Mutational Burden (TMB)</div>
                <div class="content">
                    <p><strong>THRESHOLD_TMB</strong>: {filters["THRESHOLD_TMB"]}</p>
                </div>

                <div class="subtitle">Microsatellite Instability (MSI)</div>
                <div class="content">
                    <p><strong>THRESHOLD_SITES</strong>: {filters["THRESHOLD_SITES"]}</p>
                    <p><strong>THRESHOLD_MSI</strong>: {filters["THRESHOLD_MSI"]}</p>
                </div>

                <div class="subtitle">Fusions</div>
                <div class="content">
                    <p><strong>THRESHOLD_FUSION</strong>: {filters["THRESHOLD_FUSION"]}</p>
                </div>
            </section>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)) or os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">Graphical Overview {graph_expression}</div>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""

    if os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if os.path.exists(os.path.join(new_study, general_graph_path)) or os.path.exists(os.path.join(new_study, genes_graph_path)):
        html_content += f"""</section>"""

    if annotations_list != []:
        html_content += f"""
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += f"""
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += f"""
    </body>
    </html>
    """

    with open(os.path.join(new_study, "report_VARAN.html"), "w") as f:
        f.write(html_content)


def compare_sample_file_remove(file1, file2, input_folder, outputfolder):
    clin_sam_old_df = pd.read_csv(os.path.join(input_folder, "data_clinical_sample.txt"), sep="\t")
    clin_sam_new_df = pd.read_csv(os.path.join(outputfolder, "data_clinical_sample.txt"), sep="\t")

    if os.path.exists(file1) and os.path.exists(file2):
        samples_file1 = extract_sample_list(file1)  
        samples_file2 = extract_sample_list(file2)

        removed_samples = [sample for sample in samples_file1 if not sample in samples_file2 and sample != ""]
        left_samples = len(set(samples_file1)-set(removed_samples))

    elif os.path.exists(file1) and not os.path.exists(file2):
        removed_samples = extract_sample_list(file1)
        left_samples = 0

    else:
        removed_samples = None
        left_samples = None

    return removed_samples, left_samples, len(set(clin_sam_new_df.iloc[4:, 1])), len(set(clin_sam_new_df.iloc[4:, 0]))


def extract_filters_from_html(report):
    filters = {}
    with open(report, "r", encoding="utf-8") as file:
        html_content = file.read()

    fs_match = re.search(r"<section class="filters">.*?</section>", html_content, re.DOTALL)
    if fs_match:
        section_text = fs_match.group()
        p_items = re.findall(r"<p>\s*<strong>([^<]+?)[:]*</strong>\s*[:=]?\s*(.*?)\s*</p>", section_text, re.DOTALL)
        for key, value in p_items:
            filters[key.strip().rstrip(":")] = value.strip()
        other_items = re.findall(r"<strong>([^<]+?)</strong>\s*[:=]?\s*(.*?)(?=<br>|</div>|</section>)", section_text, re.DOTALL)
        for key, value in other_items:
            k = key.strip().rstrip(":")
            if k not in filters:
                filters[k] = value.strip()

    vcf_match = re.search(r"<div class="subtitle">VCF Filters</div>\s*<div class="content">(.*?)</div>", html_content, re.DOTALL)
    if vcf_match:
        content = vcf_match.group(1)
        alt_match = re.search(r'Keep only rows where\s*<strong>\s*ALT\s*</strong>\s*is different from\s*"(.*?)"', content, re.DOTALL)
        if alt_match:
            alt_val = alt_match.group(1).strip()
            filters["ALT"] = '!="{}"'.format(alt_val)

    return filters


def extract_cancer_type_from_html(report):
    cancer_type = None

    with open(report, "r", encoding="utf-8") as file:
        html_content = file.read()

    cancer_type_match = re.search(r"<p><strong>\s*CANCER TYPE:\s*</strong>\s*(.*?)\s*</p>", html_content, re.DOTALL | re.IGNORECASE)

    if cancer_type_match:
        cancer_type = cancer_type_match.group(1).strip()

    return cancer_type
