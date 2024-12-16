#####################################
# NAME: write_report.py
# Date: 12/12/2024
version = "1.0"
# ===================================

from operator import index
import pandas as pd
import os
import re
import sys
from datetime import datetime
from configparser import ConfigParser
import versioning


config = ConfigParser()
config.read('conf.ini')

def write_infos_report(output_folder):
    clin_sam_path = os.path.join(output_folder, "data_clinical_sample.txt")
    if os.path.exists(clin_sam_path):
        clin_sam = pd.read_csv(clin_sam_path, sep="\t", header=4)

        new_samples = set(clin_sam["SAMPLE_ID"])
        new_patients = set(clin_sam["PATIENT_ID"])
        new_smpl_nr = len(new_samples)
        new_pt_nr = len(new_patients)

    if versioning.old_version_exists:

        name = re.search(r'^(.+_v)[0-9]+$', os.path.basename(output_folder)).group(1)
        actual_version = int(re.search(r'^.+_v([0-9]+)$', os.path.basename(output_folder)).group(1))

        if actual_version != 1:
            old_name = name + str(actual_version - 1)
            old_file = os.path.join(old_name, "data_clinical_sample.txt")

            if os.path.exists(old_file):
                old_clin_sam = pd.read_csv(old_file, sep="\t", header=4)

                old_samples = set(old_clin_sam["SAMPLE_ID"])
                old_patients = set(old_clin_sam["PATIENT_ID"])
                old_sample_nr = len(old_samples)
                old_patient_nr = len(old_patients)

                only_new_sam = new_samples - old_samples
                only_new_pat = new_patients - old_patients

                only_old_sam = old_samples - new_samples
                only_old_pat = old_patients - new_patients
            

    with open(os.path.join(output_folder, "report_VARAN.txt"), "w") as file:
        now = datetime.now()
        date = now.strftime("%d/%m/%Y, %H:%M:%S")
        file.write(f"Varan run - {date}\n")

        file.write(f"\nSTUDY NAME: {output_folder}")
        file.write(f"\nTOTAL SAMPLE(S): {len(new_samples)}")
        file.write(f"\nTOTAL PATIENT(S): {len(new_patients)}")

        if actual_version != 1 and os.path.exists(old_file):
            if len(only_new_sam) > 0 or len(only_old_sam) > 0:
                file.write(f"\n\nComparing this version (v_{actual_version}) with the previous existing one (v_{actual_version - 1}).\n")

            if len(only_new_sam) > 0:
                only_new_samples = ", ".join(str(e) for e in only_new_sam)
                file.write(f"\nADDED SAMPLE(S): {len(only_new_sam)} ({only_new_samples})")
                only_new_patients = ", ".join(str(e) for e in only_new_pat)
                file.write(f"\nADDED PATIENT(S): {len(only_new_pat)} ({only_new_patients})")

            if len(only_old_sam) > 0:
                only_old_samples = ", ".join(str(e) for e in only_old_sam)
                file.write(f"\nREMOVED SAMPLE(S): {len(only_old_sam)} ({only_old_samples})")
                only_old_patients = ", ".join(str(e) for e in only_old_pat)
                file.write(f"\nREMOVED PATIENT(S): {len(only_old_pat)} ({only_old_patients})")



def write_filters_report(output_folder):
    sections_to_include = {
            "Filters": ["BENIGN", "CLIN_SIG", "CONSEQUENCES", "ONCOKB_FILTER", 
                        "t_VAF_min", "t_VAF_min_novel", "t_VAF_max", 
                        "AF", "POLYPHEN", "IMPACT", "SIFT"],
            "Cna": ["HEADER_CNV", "PLOIDY"],
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

    conf_content = "\n".join(conf_content)

    report = os.path.join(output_folder, "report_VARAN.txt")
    with open(report, "a") as file:
        file.write(f"\n\nThe following configuration and filters have been used:\n")
        file.write(conf_content)



def convert_txt_to_html(input_file, output_file):
    with open(input_file, 'r') as f:
        content = f.read()

    run_date_match = re.search(r"Varan run - (\d{2}/\d{2}/\d{4}, \d{2}:\d{2}:\d{2})", content)
    study_name = re.search(r"STUDY NAME:\s*(\S+)", content)
    total_samples_match = re.search(r"TOTAL SAMPLE\(S\): (\d+)", content)
    total_patients_match = re.search(r"TOTAL PATIENT\(S\): (\d+)", content)

    run_date = run_date_match.group(1) if run_date_match else "N/A"

    current_version_match = re.search(r"this version \(v_(\d+)\)", content)
    previous_version_match = re.search(r"previous existing one \(v_(\d+)\)", content)

    current_version = current_version_match.group(1) if current_version_match else "N/A"
    previous_version = previous_version_match.group(1) if previous_version_match else "N/A"

    added_samples = extract_changes(content, "ADDED SAMPLE\(S\)")
    removed_samples = extract_changes(content, "REMOVED SAMPLE\(S\)")
    added_patients = extract_changes(content, "ADDED PATIENT\(S\)")
    removed_patients = extract_changes(content, "REMOVED PATIENT\(S\)")

    comparison_section = ""
    if added_samples or removed_samples or added_patients or removed_patients:
        comparison_section = f"""
        <section class=\"comparison\">
            <div class=\"section-title\">Comparison with Previous Version (v_{previous_version})</div>
            <div class=\"content\">
                {added_samples}
                {removed_samples}
                {added_patients}
                {removed_patients}
            </div>
        </section>
        """

    html_content = f"""
    <!DOCTYPE html>
    <html lang=\"it\">
    <head>
        <meta charset=\"UTF-8\">
        <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
        <title>Varan Run Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                background-color: #e6f0ff;
                color: #333;
                margin: 0;
                padding: 0;
                text-align: center;
            }}
            header {{
                background-color: #003366;
                color: white;
                padding: 20px;
                display: flex;
                align-items: center;
                justify-content: flex-start;
            }}
            header img {{
                max-width: 80px;
                height: auto;
                margin-right: 20px;
            }}
            h1 {{
                font-size: 3em;
                margin: 0;
                text-transform: uppercase;
            }}
            h2 {{
                font-size: 1.4em;
                margin: 5px 0 0;
                color: #b0c4de;
                font-weight: normal;
            }}
            .container {{
                max-width: 1000px;
                margin: 30px auto;
                padding: 20px;
                background-color: rgba(255, 255, 255, 0.85);
                border-radius: 8px;
                box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
                position: relative;
            }}
            .section-title {{
                font-size: 1.6em;
                color: white;
                background-color: #003366;
                padding: 10px;
                margin-bottom: 20px;
                border-radius: 5px;
            }}
            .content {{
                text-align: left;
                font-size: 1.1em;
                line-height: 1.6;
                margin-bottom: 20px;
            }}
            .filters, .cna, .tmb, .msi, .fusion {{
                margin-bottom: 40px;
            }}
        </style>
    </head>
    <body>
        <header>
            <img src=\"logo_VARAN.png\" alt=\"Logo Varan\">
            <h1>VARAN</h1>
        </header>

        <h2>Run on {run_date}</h2>

        <div class=\"container\">
            <section class=\"general-info\">
                <div class=\"section-title\">General Information</div>
                <div class=\"content\">
                    <strong>STUDY NAME:</strong> {study_name.group(1) if study_name else "N/A"}<br>
                    <strong>TOTAL SAMPLE(S):</strong> {total_samples_match.group(1) if total_samples_match else "N/A"}<br>
                    <strong>TOTAL PATIENT(S):</strong> {total_patients_match.group(1) if total_patients_match else "N/A"}<br>
                </div>
            </section>

            {comparison_section}

            <section class=\"filters\">
                <div class=\"section-title\">Filters Applied</div>
                <div class=\"content\">
                    {extract_section(content, "Filters")}
                </div>
            </section>

            <section class=\"cna\">
                <div class=\"section-title\">Copy Number Alterations (CNA)</div>
                <div class=\"content\">
                    {extract_section(content, "Cna")}
                </div>
            </section>

            <section class=\"tmb\">
                <div class=\"section-title\">Tumor Mutational Burden (TMB)</div>
                <div class=\"content\">
                    {extract_section(content, "TMB")}
                </div>
            </section>

            <section class=\"msi\">
                <div class=\"section-title\">Microsatellite Instability (MSI)</div>
                <div class=\"content\">
                    {extract_section(content, "MSI")}
                </div>
            </section>

            <section class=\"fusion\">
                <div class=\"section-title\">Fusions</div>
                <div class=\"content\">
                    {extract_section(content, "FUSION")}
                </div>
            </section>
        </div>
    </body>
    </html>
    """

    with open(output_file, 'w') as f:
        f.write(html_content)

    os.remove(input_file)


def extract_section(content, section_name):
    section_match = re.search(rf"\[{section_name}\](.*?)(\n\[|$)", content, re.DOTALL)
    section_content = section_match.group(1).strip() if section_match else "No data available."
    return section_content.replace("\n", "<br>")


def extract_changes(content, change_type):
    change_match = re.search(rf"{change_type}: (.+?)\n", content)
    if change_match:
        cleaned_type = change_type.replace(r"\(S\)", "(S)")
        return f"<strong>{cleaned_type}:</strong> {change_match.group(1)}<br>"
    return ""


