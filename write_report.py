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

"""Module for generating the VARAN HTML report.

This script collects and organizes information about the current study,
as well as previous studies, and presents the data in a structured and
readable HTML format. The generated report is intended to facilitate
data review and interpretation.

"""
from __future__ import annotations

import ast
import json
import os
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
import subprocess

import pandas as pd

import versioning
from config_loader import get_config
from filter_clinvar import check_bool
from versioning import get_git_version

config = get_config()

annotations_list = ast.literal_eval(config.get("Annotations", "ANNOTATIONS"))
vep_cache_version = int(config.get("Paths", "CACHE"))
vep_path = config.get("Paths", "VEP_PATH")
sample_type = config.get("Sample_Type", "TYPE")
ref_genome = config.get("Paths", "REF_FASTA")


def get_python_version() -> str:
    """Return the current Python version."""
    return sys.version.split()[0]


def get_clinvar_update_date(config: configparser.ConfigParser) -> str | None:
    """Extract ClinVar update date from configuration file.

    The function retrieves the CLINV entry from the configuration file,
    extracts the date in YYYYMMDD format from the filename (e.g. 
    clinvar_20260201.vcf.gz), and converts it to ISO format (YYYY-MM-DD).

    Parameters
    ----------
    config : ConfigParser
        Configuration object containing ClinVar path definition.

    Returns
    -------
    update_date : str or None
        ClinVar update date formatted as YYYY-MM-DD if successfully
        extracted, otherwise None.
    """
    if not config.has_option("Paths", "CLINV"):
        return None

    raw_value = config.get("Paths", "CLINV")

    clinvar_path = raw_value.split(",")[0].strip()

    filename = os.path.basename(clinvar_path)

    match = re.search(r"clinvar_(\d{8})", filename)
    if not match:
        return None

    date_str = match.group(1)

    try:
        formatted_date = datetime.strptime(date_str, "%Y%m%d").strftime("%Y-%m-%d")
        return formatted_date
    except ValueError:
        return None


def get_vep_version(vep_executable: str = "vep") -> str | None:
    """Retrieve the VEP version by parsing the help output.

    Some VEP installations do not accept the '--version' flag. In this case,
    we extract the version from the first line of 'vep --help' output.

    Parameters
    ----------
    vep_executable : str, optional
        Path or name of the VEP executable. Default is "vep".

    Returns
    -------
    version : str or None
        The version string of VEP if successfully retrieved, otherwise None.
    """
    try:
        result = subprocess.run(
            [vep_executable, "--help"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        # Prendi la prima riga non vuota e cerca un numero di versione
        for line in result.stdout.splitlines():
            if "ensembl-vep" in line.lower() or "vep" in line.lower():
                # estrae la prima cosa che assomiglia a x.y.z
                import re
                match = re.search(r"\d+\.\d+(\.\d+)?", line)
                if match:
                    return match.group(0)
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error running VEP: {e.stderr}")
        return None
    except FileNotFoundError:
        print(f"Executable '{vep_executable}' not found")
        return None


def write_provenance(output_folder: str, **fields: str | None) -> None:
    """Write service/provenance.json: software versions and run metadata.

    Lives in a "service" folder alongside "img", not in the study root, so it
    doesn't show up as one more file the user has to make sense of - it's
    machine-readable bookkeeping, not part of the clinical output. extract and
    remove read this back instead of recomputing/re-claiming tool versions for
    tools they may not actually invoke (e.g. VEP isn't run during an extract).

    Args:
        output_folder (str): The study's output folder for this run.
        **fields: Arbitrary provenance fields (varan_version, python_version,
            clinvar_update, vep_version, analysis_start_time, ...).

    Returns:
        None

    """
    service_dir = Path(output_folder) / "service"
    service_dir.mkdir(parents=True, exist_ok=True)
    with (service_dir / "provenance.json").open("w") as f:
        json.dump(fields, f, indent=2, default=str)


def read_provenance(study_folder: str) -> dict | None:
    """Read service/provenance.json from a study folder, if present.

    Args:
        study_folder (str): Path to a study folder that may hold provenance
            written by an earlier `write_provenance` call.

    Returns:
        dict | None: The provenance fields, or None if the file is missing or
            unreadable (e.g. an older study created before this existed).

    """
    path = Path(study_folder) / "service" / "provenance.json"
    if not path.exists():
        return None
    try:
        with path.open() as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def render_provenance_html(label: str, provenance: dict | None) -> str:
    """Render a compact "Software Environment" block from a provenance dict.

    Used by extract/remove/update to show the versions of the study they
    derive from, instead of recomputing/re-claiming tool usage they didn't
    actually invoke themselves (e.g. VEP isn't run during an extract).

    Args:
        label (str): Heading for this block, e.g. "Source study" or
            "Original study" / "Incoming data" (update shows both).
        provenance (dict | None): Fields as written by `write_provenance`, or
            None if the source study predates provenance tracking.

    Returns:
        str: The HTML block, or "" if there is nothing to show.

    """
    if not provenance:
        return ""
    field_labels = {
        "varan_version": "Varan Version",
        "python_version": "Python Version",
        "clinvar_update": "ClinVar Last Update",
        "vep_version": "VEP Version",
        "analysis_start_time": "Analysis Started",
        "report_generated": "Source Report Generated",
    }
    rows = [
        f"<p><strong>{field_label}:</strong> {provenance[key]}</p>"
        for key, field_label in field_labels.items()
        if provenance.get(key)
    ]
    if not rows:
        return ""
    return (
        f'<div class="content"><p><strong>{label}</strong></p>'
        + "".join(rows) + "</div>")


def extract_sample_list(filecase: str) -> list[str]:
    """Extract a list of sample IDs from a metadata file.

    The function opens a metadata file, searches for a line starting with
    "case_list_ids:", and extracts the sample IDs listed in that line,
    assuming they are tab-separated.

    Parameters
    ----------
    filecase : str
        Path to the metadata file from which to extract sample IDs.

    Returns
    -------
    sample_list : list of str
        List of sample IDs extracted from the "case_list_ids:" line.

    """
    sample_list = []
    file_path = Path(filecase)

    with file_path.open() as meta:
        for line in meta:
            if line.startswith("case_list_ids:"):
                sample_part = line.split(":")[1]
                samples = sample_part.split("\t")
                sample_list = [sample.strip() for sample in samples]

    return sample_list


def ghost_sample(output_folder: str) -> list[str]:
    """Identify samples missing from mutation, CNA, or SV files.

    This function reads sample IDs from three molecular data files
    (`data_mutations_extended.txt`, `data_cna.txt`, `data_sv.txt`) within the given
    output folder, and compares them to the sample IDs in `data_clinical_sample.txt`.
    Samples that are present in the clinical file but not in any of the three molecular
    files are considered "ghost" samples.

    Parameters
    ----------
    output_folder : str
        Path to the folder containing the input files.

    Returns
    -------
    list of str
        A list of ghost sample IDs.

    """
    new_samples = set()
    output_path = Path(output_folder)

    for file in ["data_mutations_extended.txt", "data_cna.txt", "data_sv.txt"]:
        new_samples = get_samples(file, new_samples, output_folder)

    all_path = output_path / "data_clinical_sample.txt"
    if all_path.exists():
        all_df = pd.read_csv(all_path, sep="\t", header=4)
        all_samples = set(all_df["SAMPLE_ID"])
    else:
        all_samples = set()

    ghosts = all_samples - new_samples

    return list(ghosts)


def build_warnings_section(
    output_folder: str,
    ghosts: list[str],
    validation_status: str | None = None,
) -> str:
    """Build a single HTML "Warnings" section collecting everything worth flagging.

    Centralizes what used to be scattered across an inline mention (ghost samples)
    and log files nobody reads unless something already looks wrong
    (noParsed_snv.log, noParsed_cnv.log, sampleID_dup.log). The cBioPortal
    validation result is shown here too, as a simple status - full validator
    detail lives in its own dedicated report_validate.html, linked from here.

    Args:
        output_folder (str): The study's output folder for this run.
        ghosts (list[str]): Sample IDs present in clinical data but with no
            detected SNV/CNA/SV alterations (see `ghost_sample`).
        validation_status (str | None): "PASS" / "REVIEW" / "FAIL", or None if
            this report type doesn't run cBioPortal validation (update/extract/
            remove don't - only walk does).

    Returns:
        str: The HTML for the section, or "" if there is nothing to warn about.

    """
    output_path = Path(output_folder)
    items = []

    if validation_status is not None:
        link = ""
        if (output_path / "report_validate.html").exists():
            link = ' — <a href="report_validate.html">validation report</a>'
        items.append(
            f"<p><strong>&#9888; cBioPortal validation:</strong> "
            f"{validation_status}{link}</p>")

    if ghosts:
        items.append(
            f"<p><strong>&#9888; Ghost samples</strong> (present in clinical data, "
            f"no SNV/CNA/SV alterations found after filtering): {ghosts}</p>")

    for log_name, label in [
        ("noParsed_snv.log", "SNV VCF files that failed to parse"),
        ("noParsed_cnv.log", "CNV VCF files that failed to parse"),
        ("sampleID_dup.log", "Duplicate sample IDs encountered"),
    ]:
        log_path = output_path / log_name
        if log_path.exists():
            content = log_path.read_text().strip()
            if content:
                entries = "; ".join(line for line in content.splitlines() if line)
                items.append(
                    f"<p><strong>&#9888; {label}</strong> "
                    f"({log_name}): {entries}</p>")

    if not items:
        return ""

    return (
        '<section class="warnings">'
        '<div class="section-title">Warnings</div>'
        '<div class="content">' + "".join(items) + "</div>"
        "</section>")


def get_samples(file: str, sample_list: set[str], output_folder: str) -> set[str]:
    """Extract sample id from genomic data file and update the existing sample set.

    Depending on the file type (`data_mutations_extended.txt`, `data_cna.txt`,
    or `data_sv.txt`),
    the function reads the appropriate column(s) to retrieve sample IDs and adds them
    to the provided sample list.

    Parameters
    ----------
    file : str
        Name of the file to process. Supported values are:
        'data_mutations_extended.txt', 'data_cna.txt', and 'data_sv.txt'.
    sample_list : set of str
        Existing set of sample identifiers to be updated.
    output_folder : str
        Path to the directory containing the input files.

    Returns
    -------
    set of str
        Updated set of sample identifiers including those extracted from the specified
        file.

    """
    path = Path(output_folder) / file
    if path.exists():
        dataframe = pd.read_csv(path, sep="\t", low_memory=False)
        if file == "data_mutations_extended.txt":
            samples = set(dataframe["Tumor_Sample_Barcode"])
        elif file == "data_cna.txt":
            samples = set(dataframe.columns[1:])
        elif file == "data_sv.txt":
            samples = set(dataframe["Sample_Id"])
        else:
            samples = set()
        sample_list = sample_list.union(samples)

    return sample_list


###########################
#           Main          #
###########################

def write_report_main(
    output_folder: str,
    cancer: str,
    filters: str,
    number_for_graph: int,
    oncokb: bool = False,
    start_time: str = "",
    validation_status: str | None = None,
) -> None:
    """Generate an HTML report summarizing the results of a VARAN analysis.

    The report includes visual graphs, filter configurations, and sample annotations,
    depending on availability and input parameters. It also identifies "ghost samples"
    (samples present in the clinical file but not in SNV, CNV, or fusion files).

    Parameters
    ----------
    output_folder : str
        Path to the directory containing analysis outputs and where the report will
        be saved.

    cancer : str
        Type of cancer being analyzed. Used for labeling purposes in the report.

    filters : str
        A string of filter codes used during variant processing (e.g., 'dvaoqc...')
        that determine which filtering criteria are applied and displayed in the report.

    number_for_graph : int
        Number of versions to consider when generating graphical representations and
        captions.

    oncokb : bool
        Flag indicating whether OncoKB annotations are included in the analysis.

    Returns
    -------
    None
        The function saves the HTML report to `output_folder` as 'report_VARAN.html'.

    """
    shutil.copy("styles.css", Path(output_folder) / "img" / "styles.css")
    now = datetime.now().astimezone()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")
    start_time_html = (
        f"<h3>Analysis started on {start_time}</h3>" if start_time else "")

    output_folder = Path(output_folder)

    img_path = Path("img") / "logo_VARAN.png"
    general_graph_path = Path("img") / "general.png"
    genes_graph_path = Path("img") / "genes.png"

    my_filters = write_filters_report()
    cancer = cancer.capitalize()

    graph_expression = get_graph_expression(number_for_graph)

    ghosts = ghost_sample(output_folder)
    new_samples, new_patients, new_smpl_nr, new_pt_nr = parse_clinical_sample(
        output_folder)

    varan_version = get_git_version()
    python_version = get_python_version()
    clinvar_update = get_clinvar_update_date(config)
    vep_version = get_vep_version(Path(vep_path, "vep"))

    write_provenance(
        output_folder,
        varan_version=varan_version,
        python_version=python_version,
        clinvar_update=clinvar_update,
        vep_version=vep_version,
        analysis_start_time=start_time,
        report_generated=datetime.now().astimezone().strftime("%d/%m/%Y, %H:%M:%S"))

    if versioning.old_version_exists:
        name = re.search(r"^(.+_v)[0-9]+$", output_folder.name).group(1)
        actual_version = int(re.search(r"^.+_v([0-9]+)$",\
        Path(output_folder).name).group(1))

        old_name = name + str(actual_version - 1)
        old_file = Path(output_folder).parent / old_name / "data_clinical_sample.txt"

        if actual_version != 1:
            if Path(old_file).exists():
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
        <link rel="stylesheet" type="text/css" href="{Path('img') / 'styles.css'}">
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
        {start_time_html}

        <div class="container">
            <section class="general-info">
                <div class="section-title">General Information</div>
                <div class="subtitle">Execution Information</div>
                <div class="content">
                    <p><strong>Command Line:</strong> {" ".join(sys.argv)}</p>
                <div class="subtitle">Software Environment</div>"""

    if isinstance(varan_version, str) and varan_version.strip():
        html_content += f"""
            <p><strong>Varan Version:</strong> {varan_version}</p>"""

    if isinstance(python_version, str) and python_version.strip():
        html_content += f"""
            <p><strong>Python Version:</strong> {python_version}</p>"""

    if isinstance(clinvar_update, str) and clinvar_update:
        html_content += f"""
            <p><strong>ClinVar Last Update:</strong> {clinvar_update}</p>"""
    
    if isinstance(vep_version, str) and vep_version.strip():
        html_content += f"""
            <p><strong>VEP Version:</strong> {vep_version}</p>"""

    if isinstance(vep_cache_version, int):
        html_content += f"""
            <p><strong>VEP Cache:</strong> {vep_cache_version}</p>"""

    html_content += f"""
                    <div class="subtitle">Study Metadata</div>
                    <p><strong>Study Name:</strong> {Path(output_folder).name}</p>
                    <p><strong>Sample Type:</strong> {sample_type}</p>
                    <p><strong>Cancer Type:</strong> {cancer}</p>"""

    if isinstance(ref_genome, str) and ref_genome.strip():
        html_content += f"""
            <p><strong>Reference Genome:</strong> {ref_genome}</p>"""

    html_content += f"""
                    <p><strong>Total Sample(s):</strong> {new_smpl_nr}</p>
                    <p><strong>Total Patient(s):</strong> {new_pt_nr}</p>
                </div>"""

    html_content += """</section>"""

    html_content += build_warnings_section(output_folder, ghosts, validation_status)

    html_content += """
        <section class="filters">
            <div class="section-title">Filters & Configuration</div>"""

    if any(letter in filters for letter in "d"):
        html_content += """
            <div class="subtitle">VCF Filters</div>"""

    if "d" in filters:
        html_content += """
                <div class="content">
                <p>Keep only VCF rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only VCF rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(letter in filters for letter in "va"):
        html_content += """
            <div class="subtitle">VAF Filters</div>"""
    if "v" in filters:
        html_content += f"""
            <div class="content">
            <p>
                <strong>T_VAF_MIN:</strong> {extract_key_value(my_filters, "T_VAF_MIN")}
            </p>
            <p>
                <strong>T_VAF_MAX:</strong> {extract_key_value(my_filters, "T_VAF_MAX")}
            </p>
        </div>"""

        if "n" in filters:
            html_content += (f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL:</strong>
                {extract_key_value(my_filters, "T_VAF_MIN_NOVEL")}
                </p>
            </div>""")

    if "a" in filters:
        drop_na_af = config.get("Filters", "drop_NA_AF")
        drop_na_af = check_bool(drop_na_af)
        excl_or_incl = "exclude" if drop_na_af else "include"

        html_content += f"""
                <div class="content">
                <p><strong>AF:</strong> {extract_key_value(my_filters, "AF")}
                & {excl_or_incl} NA</p>
            </div>"""

    if any(letter in filters for letter in "oqpycbi"):
        html_content += """
            <div class="subtitle">MAF Filters</div>"""

    if oncokb and "o" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>ONCOKB:</strong> include {', '.join([
                item.strip()
                for item in extract_key_value(my_filters, "ONCOKB_FILTER")
                .strip('[]')
                .replace('"', '')
                .split(',')
            ])}</p>
        </div>"""

    if "q" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>CONSEQUENCES:</strong> include {', '.join([
                item.strip()
                for item in extract_key_value(my_filters, "CONSEQUENCES")
                .strip('[]')
                .replace('"', '')
                .split(',')
            ])}</p>
        </div>"""

    if "y" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>POLYPHEN:</strong> include {', '.join([
                item.strip()
                for item in extract_key_value(my_filters, "POLYPHEN")
                .strip('[]')
                .replace('"', '')
                .split(',')
            ])}</p>
        </div>"""

    if "c" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>CLIN_SIG:</strong> exclude {', '.join([
                item.strip()
                for item in extract_key_value(my_filters, "CLIN_SIG")
                .strip('[]')
                .replace('"', '')
                .split(',')
            ])}</p>
        </div>"""

    if "i" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>IMPACT:</strong> exclude {', '.join([
                item.strip()
                for item in extract_key_value(my_filters, "IMPACT")
                .strip('[]')
                .replace('"', '')
                .split(',')
            ])}</p>
        </div>"""

    if "s" in filters:
        html_content += f"""
            <div class="content">
            <p><strong>SIFT:</strong> include {', '.join([
                item.strip()
                for item in extract_key_value(my_filters, "SIFT")
                .strip('[]')
                .replace('"', '')
                .split(',')
            ])}</p>
        </div>"""

    if "p" in filters:
        html_content += """
                <div class="content">
                 <p><strong>FILTER:</strong> = PASS
            </div>"""

    tmb_section = re.sub(r"[{}']", "", extract_section(my_filters, "TMB"))
    tmb_section = re.sub(r"([,:])(?=\S)", r"\1 ", tmb_section)

    html_content += f"""
            <div class="subtitle">Copy Number Alterations (CNA)</div>
            <div class="content">
                {extract_section(my_filters, "Cna")}
            </div>

            <div class="subtitle">Tumor Mutational Burden (TMB)</div>
            <div class="content">
                {tmb_section}
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

    if versioning.old_version_exists and actual_version != 1:
            html_content += f"""
            <section class="comparison">
                <div class="section-title">
                    Comparison with Previous Version (_v{actual_version - 1})
                </div>
                <div class="content">
                    <p><strong>ADDED PATIENT(S):</strong> {len(only_new_pat)}</p>
                    <p><strong>ADDED SAMPLE(S):</strong> {len(only_new_sam)}</p>
                    <p><strong>REMOVED PATIENT(S):</strong> {len(only_old_pat)}</p>
                    <p><strong>REMOVED SAMPLE(S):</strong> {len(only_old_sam)}</p>
                </div>
            </section>
            """

    if (Path(output_folder) / general_graph_path).exists()\
    or (Path(output_folder) / genes_graph_path).exists():
        html_content += f"""
            <section class="graphs">
                <div class="section-title">
                    Graphical Overview {graph_expression}
                </div>"""

    if (Path(output_folder) / general_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""

    if (Path(output_folder) / genes_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if (Path(output_folder) / general_graph_path).exists()\
    or (Path(output_folder) / genes_graph_path).exists():
        html_content += """</section>"""

    if annotations_list != []:
        html_content += """
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += """
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += """
    </body>
    </html>
    """

    output_path = Path(output_folder) / "report_VARAN.html"
    with output_path.open("w") as f:
        f.write(html_content)


def parse_clinical_sample(output_folder: str) -> tuple[set[str], set[str], int, int]:
    """Parse clinical sample data from a predefined file.

    This function attempts to read a tab-delimited clinical sample file named
    'data_clinical_sample.txt' located in the given output folder. If the file
    exists, it extracts the SAMPLE_ID and PATIENT_ID columns to determine
    the set of unique sample and patient identifiers.

    Parameters
    ----------
    output_folder : str
        The path to the folder containing the clinical sample file.

    Returns
    -------
    tuple of (set of str, set of str, int, int)
        A tuple containing:
        - Set of unique sample IDs
        - Set of unique patient IDs
        - Count of unique samples
        - Count of unique patients

    """
    path = Path(output_folder) / "data_clinical_sample.txt"
    if not path.exists():
        return set(), set(), 0, 0

    clin_sam = pd.read_csv(path, sep="\t", header=4)
    samples = set(clin_sam["SAMPLE_ID"])
    patients = set(clin_sam["PATIENT_ID"])
    return samples, patients, len(samples), len(patients)


def get_graph_expression(number_for_graph: int) -> str:
    """Indicate how many graph versions to display.

    This function returns an empty string for a single graph version,
    a custom message if the number exceeds a predefined limit (4),
    or a formatted string indicating how many recent versions are included.

    Parameters
    ----------
    number_for_graph : int
        The number of versions to consider for graph display.

    Returns
    -------
    str
        A string describing the graph version selection.

    """
    limit = 4
    if number_for_graph == 1:
        return ""
    if number_for_graph > limit:
        return "(last 5 versions)"
    return f"(last {number_for_graph} versions)"


def write_filters_report() -> list[str]:
    """Extract selected filtering and parameters from the configuration.

    This function reads predefined keys from specific sections such as Filters, CNA,
    TMB, MSI, and FUSION in the global `config` object and compiles them into a list
    of lines in `[SECTION]` and `KEY = VALUE` format. The output is used for
    reporting purposes, such as populating filter summaries in HTML reports.

    Returns
    -------
    list of str
        A list of strings representing selected configuration settings in INI-style
        format.

    """
    sections_to_include = {
            "Filters": ["BENIGN", "CLIN_SIG", "CONSEQUENCES", "ONCOKB_FILTER",
                        "t_VAF_min", "t_VAF_min_novel", "t_VAF_max",
                        "AF", "POLYPHEN", "IMPACT", "SIFT"],
            "Cna": ["PLOIDY", "CNVKIT_algorithm"],
            "TMB": ["THRESHOLD_TMB"],
            "MSI": ["THRESHOLD_SITES", "THRESHOLD_MSI"],
            "FUSION": ["THRESHOLD_FUSION"],
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


def extract_section(content: list[str], section_name: str) -> str:
    """Extract and formats a specific configuration section from a list of strings.

    Given a list of strings representing a configuration file (e.g., INI-style content),
    this function locates the specified section by its name, extracts its key-value
    pairs, and returns them as an HTML-formatted string. Each key is bolded using
    <strong> tags, and entries are separated by <br> line breaks.

    Parameters
    ----------
    content : list of str
        The configuration content, typically generated by `write_filters_report`,
        as a list of lines (e.g., ["[Section]", "KEY = value", ...]).
    section_name : str
        The name of the section to extract (case-insensitive).

    Returns
    -------
    str
        The HTML-formatted string of the section's contents. If the section is not
        found, returns "N/A".

    """
    content = "\n".join(content)
    match = re.search(
        rf"\[{section_name}\](.*?)(?=\n\[[^\]]+\]|\Z)",
        content,
        re.DOTALL | re.IGNORECASE,
    )

    if not match:
        return "N/A"
    section_content = match.group(1).strip()
    formatted_content = re.sub(
    r"^\s*(\S+)\s*=\s*",
    r"<strong>\1</strong> = ",
    section_content,
    flags=re.MULTILINE,
    )

    return formatted_content.replace("\n", "<br>")


def extract_key_value(filters: list[str], key_name: str) -> str | None:
    """Extract the value associated with a given key from a list of filter lines.

    Searches through a list of strings (each representing a line of a configuration or
    filter) for a line containing the specified key in the format "key = value".
    Returns the value as a stripped string without surrounding quotes. If the key is
    not found, returns "NA".

    Parameters
    ----------
    filters : list of str
        List of strings representing configuration lines, e.g. ["KEY = value", ...].
    key_name : str
        The key whose value needs to be extracted (case-sensitive).

    Returns
    -------
    str
        The extracted value as a string without quotes, or "N/A" if the key is not
        found.

    """
    pattern = rf"{key_name}\s*=\s*(.+)"
    for line in filters:
        match = re.search(pattern, line)
        if match:
            return match.group(1).strip().strip('"').strip("'")
    return "N/A"


###########################
#          Update         #
###########################

def write_report_update(original_study: Path, updating_with: Path,
new_study: Path, number_for_graph: int, start_time: str = "") -> None:
    """Create an HTML update report comparing an original study with a new update.

    Parameters
    ----------
    original_study : Path
        Path to the folder of the original study.
    updating_with : Path
        Path to the folder containing the update.
    new_study : Path
        Path to the folder where the updated study has been saved.
    number_for_graph : int
        Number of samples used in generating the bar plot.

    Returns
    -------
    None
        This function writes the HTML report directly to disk.

    """
    old_img_path = original_study / "img" / "logo_VARAN.png"
    general_graph_path = Path("img") / "general.png"
    genes_graph_path = Path("img") / "genes.png"

    if old_img_path.exists():
        img_output_dir = Path(new_study) / "img"
        img_output_dir.mkdir(parents=True, exist_ok=True)
        new_img_path = img_output_dir / "logo_VARAN.png"
        shutil.copy(old_img_path, new_img_path)

    shutil.copy("styles.css", Path(new_study) / "img" / "styles.css")

    now = datetime.now().astimezone()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")
    start_time_html = (
        f"<h3>Analysis started on {start_time}</h3>" if start_time else "")

    limit = 4

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > limit:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(new_study)

    new_samples, new_patients, new_smpl_nr, new_pt_nr = parse_clinical_sample(
    new_study)
    if versioning.old_version_exists:
        name = re.search(r"^(.+_v)[0-9]+$", new_study.name).group(1)
        actual_version = int(re.search(r"^.+_v([0-9]+)$",\
        Path(new_study).name).group(1))

        old_name = name + str(actual_version - 1)
        old_file = Path(new_study).parent / old_name / "data_clinical_sample.txt"

        if actual_version != 1:
            if Path(old_file).exists():
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

    output_file = Path(new_study) / "report_VARAN.html"

    case_list1 = Path(original_study) / "case_lists"
    case_list2 = Path(updating_with) / "case_lists"
    cna_1 = Path(case_list1) / "cases_cna.txt"
    cna_2 = Path(case_list2) / "cases_cna.txt"
    sequenced_1 = Path(case_list1) / "cases_sequenced.txt"
    sequenced_2 = Path(case_list2) / "cases_sequenced.txt"
    sv_1 = Path(case_list1) / "cases_sv.txt"
    sv_2 = Path(case_list2) / "cases_sv.txt"

    clin_sam_old_df = pd.read_csv(Path(original_study) / "data_clinical_sample.txt",\
    sep="\t")
    clin_sam_new_df = pd.read_csv(Path(updating_with) / "data_clinical_sample.txt",\
    sep="\t")

    updated_clin_pat = (
        set(clin_sam_old_df.iloc[4:, 1]) & set(clin_sam_new_df.iloc[4:, 1])
    )

    updated_clin_sample = (
        set(clin_sam_old_df.iloc[4:, 0]) & set(clin_sam_new_df.iloc[4:, 0])
    )

    added_clin_pat = (
        set(clin_sam_new_df.iloc[4:, 1]).difference(clin_sam_old_df.iloc[4:, 1])
    )

    added_clin_sample = (
        set(clin_sam_new_df.iloc[4:, 0]).difference(clin_sam_old_df.iloc[4:, 0])
    )

    updated_samples_cna, added_samples_cna, total_patients, total_samples = (
        compare_sample_file_update(cna_1, cna_2, new_study)
    )
    updated_samples_sequenced, added_samples_sequenced, _, _ =(
        compare_sample_file_update(sequenced_1, sequenced_2, new_study)
    )

    updated_samples_sv, added_samples_sv, _, _ = (
        compare_sample_file_update(sv_1, sv_2, new_study)
    )

    old_report = Path(original_study) / "report_VARAN.html"
    updating_report = Path(updating_with) / "report_VARAN.html"

    if old_report.exists() and updating_report.exists():
        filters1 = extract_filters_from_html(old_report)
        filters2 = extract_filters_from_html(updating_report)
        cancer_type1 = extract_cancer_type_from_html(old_report)
        cancer_type2 = extract_cancer_type_from_html(updating_report)
        sample_type1 = extract_sample_type_from_html(old_report)
        sample_type2 = extract_sample_type_from_html(updating_report)

    else:
        filters1 = {}
        filters2 = {}
        cancer_type1 = None
        cancer_type2 = None
        sample_type1 = None
        sample_type2 = None

    order = ["T_VAF_MIN", "T_VAF_MIN_NOVEL", "T_VAF_MAX", "AF", "ONCOKB", "IMPACT",\
    "CLIN_SIG", "CONSEQUENCES", "POLYPHEN", "SIFT", "PLOIDY", "CNVKIT_algorithm",\
    "THRESHOLD_TMB", "THRESHOLD_SITES", "THRESHOLD_MSI", "THRESHOLD_FUSION"]

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
        <link rel="stylesheet" type="text/css" href="{Path('img') / 'styles.css'}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan - Update</title>
    </head>
    <body>
        <header>
            <img src="{Path('img') / 'logo_VARAN.png'}" alt="Logo Varan">
            <h1>VARAN - Update</h1>
        </header>

        <h2>Report generate on {date}</h2>
        {start_time_html}
        <div class="container">
            <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>ORIGINAL STUDY:</strong> {Path(original_study).name}</p>
                    <p><strong>UPDATING WITH:</strong> {Path(updating_with).name}</p>
                    <p><strong>NEW STUDY:</strong> {Path(new_study).name}</p>"""

    if (sample_type1 == sample_type2) and (sample_type1 is not None):
        html_content += f"""<p><strong>SAMPLE TYPE:</strong> {sample_type1}</p>"""
    else:
        html_content += """<p><strong>SAMPLE TYPE:</strong> Mixed</p>"""

    if (cancer_type1 == cancer_type2) and (cancer_type1 is not None):
        html_content += f"""<p><strong>CANCER TYPE:</strong> {cancer_type1}</p>"""
    else:
        html_content += """<p><strong>CANCER TYPE:</strong> Mixed</p>"""

    html_content += f"""
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>Total Patients:</strong> {total_patients}</p>
                    <p><strong>Total Samples:</strong> {total_samples}</p>
                </div>
    """

    original_provenance = read_provenance(original_study)
    incoming_provenance = read_provenance(updating_with)
    write_provenance(
        new_study,
        original_study=original_provenance,
        incoming_data=incoming_provenance)
    html_content += render_provenance_html("Original study", original_provenance)
    html_content += render_provenance_html("Incoming data", incoming_provenance)

    html_content += build_warnings_section(new_study, ghosts)

    if updated_clin_sample:
        html_content += (
            f"""
            <div class="content">
                <p><strong>UPDATED:</strong></p>
                <p>&emsp;<strong>{len(updated_clin_pat)} Patients:</strong>
                {", ".join(updated_clin_pat)}</p>
                <p>&emsp;<strong>{len(updated_clin_sample)} Samples:</strong>
                {", ".join(updated_clin_sample)}</p>
            </div>
            """
        )

    if added_clin_sample:
        html_content += (
            f"""
            <div class="content">
                <p><strong>ADDED:</strong></p>
                <p>&emsp;<strong>{len(added_clin_pat)} Patients:</strong>
                {", ".join(added_clin_pat)}</p>
                <p>&emsp;<strong>{len(added_clin_sample)} Samples:</strong>
                {", ".join(added_clin_sample)}</p>
            </div>
            """
        )

    if (updated_samples_cna or updated_samples_sequenced or updated_samples_sv
    or added_samples_cna or added_samples_sequenced or added_samples_sv):
        html_content += """
        <div class="section-title">Detailed Overview</div>
            <div class="content">"""

    if updated_samples_cna or added_samples_cna:
        html_content += f"""
    <p><strong>Cases_CNA:</strong></p>"""

    if updated_samples_cna:
        html_content += f"""
    <p>&emsp;<strong>Updated:</strong> {len(updated_samples_cna)} samples (
        {', '.join(updated_samples_cna)})</p>"""

    if added_samples_cna:
        html_content += f"""
    <p>&emsp;<strong>Added:</strong> {len(added_samples_cna)} samples (
        {', '.join(added_samples_cna)})</p>"""

    if updated_samples_sequenced or added_samples_sequenced:
        html_content += f"""
    <p><strong>Cases_sequenced:</strong></p>"""

    if updated_samples_sequenced:
        html_content += f"""
    <p>&emsp;<strong>Updated:</strong> {len(updated_samples_sequenced)} samples (
        {', '.join(updated_samples_sequenced)})</p>"""

    if added_samples_sequenced:
        html_content += f"""
    <p>&emsp;<strong>Added:</strong> {len(added_samples_sequenced)} samples (
        {', '.join(added_samples_sequenced)})</p>"""

    if updated_samples_sv or added_samples_sv:
        html_content += f"""
    <p><strong>Cases_sv:</strong></p>"""

    if updated_samples_sv:
        html_content += f"""
    <p>&emsp;<strong>Updated:</strong> {len(updated_samples_sv)} samples (
        {', '.join(updated_samples_sv)})</p>"""

    if added_samples_sv:
        html_content += f"""
    <p>&emsp;<strong>Added:</strong> {len(added_samples_sv)} samples (
        {', '.join(added_samples_sv)})</p>"""

    if (updated_samples_cna or updated_samples_sequenced or updated_samples_sv
    or added_samples_cna or added_samples_sequenced or added_samples_sv):
        html_content += """
        </div>
        """

    if common_filters != []:
        html_content += """
            <section class="filters">
                <div class="section-title">Filters & Configuration</div>"""

    if any(filt in common_filters for filt in ["ALT"]):
        html_content += """
            <div class="subtitle">VCF Filters</div>"""

    if "ALT" in common_filters:
        html_content += """
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(filt in common_filters for filt in ["T_VAF_MIN", "T_VAF_MAX", "AF"]):
        html_content += """
            <div class="subtitle">VAF Filters</div>"""

    if any(filt in common_filters for filt in ["T_VAF_MIN", "T_VAF_MAX"]):
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

    if any(filt in common_filters for filt in ["ONCOKB", "CONSEQUENCES", "POLYPHEN",
    "CLIN_SIG", "IMPACT", "SIFT"]):
        html_content += """
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

    keys_to_check = {"PLOIDY", "CNVKIT_algorithm", "THRESHOLD_TMB",
    "THRESHOLD_SITES", "THRESHOLD_MSI", "THRESHOLD_FUSION"}

    if common_filters != {} and all(key in filters1 for key in keys_to_check):
        html_content += f"""
                <div class="subtitle">Copy Number Alterations (CNA)</div>
                <div class="content">
                    <p><strong>PLOIDY</strong> = {filters1["PLOIDY"]}</p>
                    <p><strong>CNVKIT Algorithm</strong> = {filters1["CNVKIT_algorithm"]}</p>
                </div>

                <div class="subtitle">Tumor Mutational Burden (TMB)</div>
                <div class="content">
                    <p><strong>THRESHOLD_TMB</strong>: {filters1["THRESHOLD_TMB"]}</p>
                </div>

                <div class="subtitle">Microsatellite Instability (MSI)</div>
                <div class="content">
                    <p><strong>
                        THRESHOLD_SITES</strong>: {filters1["THRESHOLD_SITES"]}
                    </p>
                    <p><strong>THRESHOLD_MSI</strong>: {filters1["THRESHOLD_MSI"]}</p>
                </div>

                <div class="subtitle">Fusions</div>
                <div class="content">
                    <p><strong>
                        THRESHOLD_FUSION</strong>: {filters1["THRESHOLD_FUSION"]}
                    </p>
                </div>
            </section>
        </section>"""

    all_keys = filters1.keys() | filters2.keys()
    if any(filters1.get(f) != filters2.get(f) for f in all_keys):
        html_content += f"""
            <div class="section-title">Differences in Filters & Configurations</div>
            <div class="content">
            <br />
            <table>
                <thead>
                    <tr>
                        <th>Filter</th>
                        <th>{Path(original_study).name}</th>
                        <th>{Path(updating_with).name}</th>
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

    if versioning.old_version_exists and actual_version != 1:
            html_content += f"""
            <section class="comparison">
                <div class="section-title">
                    Comparison with Previous Version (_v{actual_version - 1})
                </div>
                <div class="content">
                    <p><strong>ADDED PATIENT(S):</strong> {len(only_new_pat)}</p>
                    <p><strong>ADDED SAMPLE(S):</strong> {len(only_new_sam)}</p>
                    <p><strong>REMOVED PATIENT(S):</strong> {len(only_old_pat)}</p>
                    <p><strong>REMOVED SAMPLE(S):</strong> {len(only_old_sam)}</p>
                </div>
            </section>
            """

    new_study_path = Path(new_study)
    if ((new_study_path / general_graph_path).exists()
    or (new_study_path / genes_graph_path).exists()):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">
                    Graphical Overview {graph_expression}
                </div>"""

    if (Path(new_study) / general_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""

    if (Path(new_study) / genes_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if (
        (Path(new_study) / general_graph_path).exists()
        or (Path(new_study) / genes_graph_path).exists()
    ):
        html_content += """</section>"""

    if annotations_list != []:
        html_content += """
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += """
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += """
    </body>
    </html>
    """

    output_file = Path(output_file)
    with output_file.open("w") as f:
        f.write(html_content)


def compare_sample_file_update(
    file1: str | Path,
    file2: str | Path,
    outputfolder: str | Path,
) -> tuple[list[str] | None, list[str] | None, int, int]:
    """Compare two sample files and identifies updated and new samples.

    Args:
        file1 (Union[str, Path]): Path to the first sample file.
        file2 (Union[str, Path]): Path to the second sample file.
        outputfolder (Union[str, Path]): Path to the folder containing the clinical
        sample data file.

    Returns:
        Tuple:
            - updated_samples (Optional[list[str]]): A list of samples that are present
            in both files.
            - new_samples (Optional[list[str]]): A list of samples that are present only
            in the second file.
            - total_patients (int): The number of unique patients in the clinical sample
            data.
            - total_samples (int): The number of unique samples in the clinical sample
            data.

    """
    clin_sam_output_df = pd.read_csv(
        Path(outputfolder) / "data_clinical_sample.txt", sep="\t")

    if Path(file1).exists() and Path(file2).exists():
        samples_file1 = extract_sample_list(file1)
        samples_file2 = extract_sample_list(file2)

        updated_samples = [
            sample for sample in samples_file1
            if sample in samples_file2 and sample != ""
        ]

        new_samples = [
            sample for sample in samples_file2
            if sample not in samples_file1 and sample != ""
        ]

    else:
        updated_samples = None
        new_samples = None

    return (
        updated_samples,
        new_samples,
        len(set(clin_sam_output_df.iloc[4:, 1])),
        len(set(clin_sam_output_df.iloc[4:, 0])),
    )


###########################
#          Extract        #
###########################

def write_report_extract(original_study: str, new_study: str,
    number_for_graph: int, start_time: str = "") -> None:
    """Generate an HTML report comparing original and new studies.

    Copies logo and styles, extracts sample data and filters,
    summarizes sample changes, and includes graphical and filter info.

    Args:
        original_study (str): Path to the original study folder.
        new_study (str): Path to the new study folder.
        number_for_graph (int): Number of recent versions for graph titles.

    Writes:
        report_VARAN.html in the new_study folder.

    """
    old_img_path = Path(original_study) / "img" / "logo_VARAN.png"
    general_graph_path = Path("img") / "general.png"
    genes_graph_path = Path("img") / "genes.png"

    if old_img_path.exists():
        img_output_dir = Path(new_study) / "img"
        img_output_dir.mkdir(parents=True, exist_ok=True)
        new_img_path = img_output_dir / "logo_VARAN.png"
        shutil.copy(old_img_path, new_img_path)

    shutil.copy("styles.css", Path(new_study) / "img" / "styles.css")

    now = datetime.now().astimezone()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")
    start_time_html = (
        f"<h3>Analysis started on {start_time}</h3>" if start_time else "")

    limit = 4

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > limit:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(new_study)

    new_samples, new_patients, new_smpl_nr, new_pt_nr = parse_clinical_sample(
    new_study)
    if versioning.old_version_exists:
        name = re.search(r"^(.+_v)[0-9]+$", new_study.name).group(1)
        actual_version = int(re.search(r"^.+_v([0-9]+)$",\
        Path(new_study).name).group(1))

        old_name = name + str(actual_version - 1)
        old_file = Path(new_study).parent / old_name / "data_clinical_sample.txt"

        if actual_version != 1:
            if Path(old_file).exists():
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

    old_report = Path(original_study) / "report_VARAN.html"
    if Path(old_report).exists():
        filters = extract_filters_from_html(old_report)
        cancer_type = extract_cancer_type_from_html(old_report)
        sample_type = extract_sample_type_from_html(old_report)
    else:
        filters = {}
        cancer_type = None
        sample_type = None

    case_list1 = Path(original_study) / "case_lists"
    case_list2 = Path(new_study) / "case_lists"
    cna_1 = case_list1 / "cases_cna.txt"
    cna_2 = case_list2 / "cases_cna.txt"
    sequenced_1 = case_list1 / "cases_sequenced.txt"
    sequenced_2 = case_list2 / "cases_sequenced.txt"
    sv_1 = case_list1 / "cases_sv.txt"
    sv_2 = case_list2 / "cases_sv.txt"

    extracted_samples_cna, total_patients, total_samples = compare_sample_file_extract(
        cna_1, cna_2, new_study)
    extracted_samples_sequenced, _, _ = compare_sample_file_extract(
        sequenced_1, sequenced_2, new_study)
    extracted_samples_sv, _, _ = compare_sample_file_extract(
        sv_1, sv_2, new_study)

    html_content = f"""
    <!DOCTYPE html>
    <html lang="it">
    <head>
        <link rel="stylesheet" type="text/css" href="{Path('img') / 'styles.css'}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan - Extract</title>
    </head>
    <body>
        <header>
            <img src="{Path('img') / 'logo_VARAN.png'}" alt="Logo Varan">
            <h1>VARAN - Extract</h1>
        </header>

        <h2>Generate on {date}</h2>
        {start_time_html}

        <div class="container">
            <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>ORIGINAL STUDY:</strong> {Path(original_study).name}</p>
                    <p><strong>NEW STUDY:</strong> {Path(new_study).name}</p>"""

    if sample_type:
        html_content += f"""<p><strong>SAMPLE TYPE:</strong> {sample_type}</p>"""

    if cancer_type:
        html_content += f"""<p><strong>CANCER TYPE:</strong> {cancer_type}</p>"""

    html_content += f"""
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>Total Patients:</strong> {total_patients}</p>
                    <p><strong>Total Samples:</strong> {total_samples}</p>
                </div>"""

    source_provenance = read_provenance(original_study)
    if source_provenance:
        write_provenance(new_study, **source_provenance)
    html_content += render_provenance_html("Source study", source_provenance)

    html_content += build_warnings_section(new_study, ghosts)

    html_content += """
            <div class="section-title">Detailed Overview</div>
                <div class="content">
                    <p><strong>Cases_CNA:</strong></p>
    """

    if extracted_samples_cna:
        html_content += f"""
            <p>&emsp;<strong>Extracted:</strong> {len(extracted_samples_cna)}
            samples ({', '.join(extracted_samples_cna)})</p>
        """
    else:
        html_content += f"""
            <p>&emsp;None of the extracted samples was in the original
            study ({Path(original_study).name})</p>
        """

    html_content += """
        <p><strong>Cases_sequenced:</strong></p>
        """

    sample_count = len(extracted_samples_sequenced)
    sample_list = ", ".join(extracted_samples_sequenced)
    if extracted_samples_sequenced:
        html_content += (
            f"""
            <p>&emsp;<strong>Extracted:</strong> {sample_count} samples
            ({sample_list})</p>
            """
        )

    else:
        html_content += f"""
            <p>&emsp;None of the extracted samples was in the original study
            ({Path(original_study).name})</p>
        """

    html_content += """
        <p><strong>Cases_sv:</strong></p>
        """

    if extracted_samples_sv:
        html_content += f"""
            <p>&emsp;<strong>Extracted:</strong> {len(extracted_samples_sv)} samples
            ({', '.join(extracted_samples_sv)})</p>
        """
    else:
        html_content += f"""
            <p>&emsp;None of the extracted samples was in the original study
            ({Path(original_study).name})</p>
        """

    if filters != {}:
        html_content += """
            <section class="filters">
                <div class="section-title">Filters & Configuration</div>"""

    if any(filt in filters for filt in ["ALT"]):
        html_content += """
            <div class="subtitle">VCF Filters</div>"""

    if "ALT" in filters:
        html_content += """
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(filt in filters for filt in ["T_VAF_MIN", "T_VAF_MAX", "AF"]):
        html_content += """
            <div class="subtitle">VAF Filters</div>"""

    if any(filt in filters for filt in ["T_VAF_MIN", "T_VAF_MAX"]):
        html_content += f"""
            <div class="content">
            <p><strong>T_VAF_MIN</strong>: {filters["T_VAF_MIN"]}</p>
            <p><strong>T_VAF_MAX</strong>: {filters["T_VAF_MAX"]}</p>
        </div>"""
        if "T_VAF_MIN_NOVEL" in filters:
            html_content += f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL</strong>: {filters["T_VAF_MIN_NOVEL"]}</p>
            </div>"""

    if "AF" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>AF</strong>: {filters["AF"]}</p>
            </div>"""

    if any(
        filt in filters
        for filt in [
            "ONCOKB", "CONSEQUENCES", "POLYPHEN",
            "CLIN_SIG", "IMPACT", "SIFT", "FILTER",
        ]
    ):
        html_content += """
            <div class="subtitle">MAF Filters</div>"""

    if "ONCOKB" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>ONCOKB</strong>: {filters["ONCOKB"]}</p>
            </div>"""

    if "CONSEQUENCES" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>CONSEQUENCES</strong>: {filters["CONSEQUENCES"]}</p>
            </div>"""

    if "POLYPHEN" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>POLYPHEN</strong>: {filters["POLYPHEN"]}</p>
            </div>"""

    if "CLIN_SIG" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>CLIN_SIG</strong>: {filters["CLIN_SIG"]}</p>
            </div>"""

    if "IMPACT" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>IMPACT</strong>: {filters["IMPACT"]}</p>
            </div>"""

    if "SIFT" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>SIFT</strong>: {filters["SIFT"]}</p>
            </div>"""

    keys_to_check = {"PLOIDY", "CNVKIT_algorithm", "THRESHOLD_TMB",
    "THRESHOLD_SITES", "THRESHOLD_MSI", "THRESHOLD_FUSION"}

    if filters != {} and all(key in filters for key in keys_to_check):
        html_content += f"""
                <div class="subtitle">Copy Number Alterations (CNA)</div>
                <div class="content">
                    <p><strong>PLOIDY</strong> = {filters["PLOIDY"]}</p>
                    <p><strong>CNVKIT Algorithm</strong> = {filters["CNVKIT_algorithm"]}</p>
                </div>

                <div class="subtitle">Tumor Mutational Burden (TMB)</div>
                <div class="content">
                    <p><strong>THRESHOLD_TMB</strong>: {filters["THRESHOLD_TMB"]}</p>
                </div>

                <div class="subtitle">Microsatellite Instability (MSI)</div>
                <div class="content">
                    <p><strong>
                        THRESHOLD_SITES</strong>: {filters["THRESHOLD_SITES"]}
                    </p>
                    <p><strong>THRESHOLD_MSI</strong>: {filters["THRESHOLD_MSI"]}</p>
                </div>

                <div class="subtitle">Fusions</div>
                <div class="content">
                    <p><strong>
                        THRESHOLD_FUSION</strong>: {filters["THRESHOLD_FUSION"]}
                    </p>
                </div>
            </section>"""

        if filters != {}:
            html_content += """
            </section>"""

    if versioning.old_version_exists and actual_version != 1:
            html_content += f"""
            <section class="comparison">
                <div class="section-title">
                    Comparison with Previous Version (_v{actual_version - 1})
                </div>
                <div class="content">
                    <p><strong>ADDED PATIENT(S):</strong> {len(only_new_pat)}</p>
                    <p><strong>ADDED SAMPLE(S):</strong> {len(only_new_sam)}</p>
                    <p><strong>REMOVED PATIENT(S):</strong> {len(only_old_pat)}</p>
                    <p><strong>REMOVED SAMPLE(S):</strong> {len(only_old_sam)}</p>
                </div>
            </section>
            """

    if ((Path(new_study) / general_graph_path).exists()
        or (Path(new_study) / genes_graph_path).exists()):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">
                    Graphical Overview {graph_expression}
                </div>"""

    if (Path(new_study) / general_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""


    if (Path(new_study) / genes_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if (
        (Path(new_study) / general_graph_path).exists()
        or (Path(new_study) / genes_graph_path).exists()
    ):
        html_content += """</section>"""

    if annotations_list != []:
        html_content += """
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += """
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += """
    </body>
    </html>
    """

    report_path = Path(new_study) / "report_VARAN.html"
    with report_path.open("w") as f:
        f.write(html_content)


def compare_sample_file_extract(
    file1: str | Path,
    file2: str | Path,
    outputfolder: str | Path,
) -> tuple[list[str] | None, int, int]:
    """Compare two sample files and associated clinical data between two studies.

    Args:
        file1 (str | Path): Path to the original sample file.
        file2 (str | Path): Path to the updated sample file.
        outputfolder (str | Path): Folder containing the updated clinical sample data.

    Returns:
        tuple[list[str] | None, int, int]: A tuple containing:
            - The list of overlapping sample names between the files
            (or None if files do not exist).
            - The count of unique patients in the updated clinical file.
            - The count of unique samples in the updated clinical file.

    """
    outputfolder = Path(outputfolder)
    clin_sam_new_df = pd.read_csv(outputfolder / "data_clinical_sample.txt", sep="\t")

    if Path(file1).exists() and Path(file2).exists():
        samples_file1 = extract_sample_list(file1)
        samples_file2 = extract_sample_list(file2)

        extracted_samples = [
            sample
            for sample in samples_file2
            if sample in samples_file1 and sample != ""
        ]
    else:
        extracted_samples = None

    num_patients = len(set(clin_sam_new_df.iloc[4:, 1]))
    num_samples = len(set(clin_sam_new_df.iloc[4:, 0]))

    return extracted_samples, num_patients, num_samples


###########################
#          Remove         #
###########################

def write_report_remove(
    original_study: str | Path,
    new_study: str | Path,
    number_for_graph: int,
    start_time: str = "",
) -> None:
    """Generate an HTML report summarizing sample removals and study updates.

    This function copies necessary images and stylesheets, extracts information
    from the original study, compares sample lists between the original and new studies,
    and creates a detailed HTML report that includes general information, removed
    samples, filters applied, and graphical overviews if available.

    Args:
        original_study (str): Path to the directory of the original study.
        new_study (str): Path to the directory of the updated/new study.
        number_for_graph (int): Number indicating how many recent versions to display in
        graphs.

    Returns:
        None

    """
    old_img_path = original_study / "img" / "logo_VARAN.png"
    general_graph_path = Path("img") / "general.png"
    genes_graph_path = Path("img") / "genes.png"

    if old_img_path.exists():
        img_output_dir = Path(new_study) / "img"
        img_output_dir.mkdir(parents=True, exist_ok=True)
        new_img_path = img_output_dir / "logo_VARAN.png"
        shutil.copy(old_img_path, new_img_path)

    shutil.copy("styles.css", Path(new_study) / "img" / "styles.css")

    now = datetime.now().astimezone()
    date = now.strftime("%d/%m/%Y, %H:%M:%S")
    start_time_html = (
        f"<h3>Analysis started on {start_time}</h3>" if start_time else "")

    limit = 4

    if number_for_graph == 1:
        graph_expression = ""
    elif number_for_graph > limit:
        graph_expression = "(last 5 versions)"
    else:
        graph_expression = f"(last {number_for_graph} versions)"

    ghosts = ghost_sample(new_study)

    new_samples, new_patients, new_smpl_nr, new_pt_nr = parse_clinical_sample(
    new_study)
    if versioning.old_version_exists:
        name = re.search(r"^(.+_v)[0-9]+$", new_study.name).group(1)
        actual_version = int(re.search(r"^.+_v([0-9]+)$",\
        Path(new_study).name).group(1))

        old_name = name + str(actual_version - 1)
        old_file = Path(new_study).parent / old_name / "data_clinical_sample.txt"

        if actual_version != 1:
            if Path(old_file).exists():
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

    old_report = Path(original_study) / "report_VARAN.html"
    if Path(old_report).exists():
        filters = extract_filters_from_html(old_report)
        cancer_type = extract_cancer_type_from_html(old_report)
        sample_type = extract_sample_type_from_html(old_report)
    else:
        filters = {}
        cancer_type = None
        sample_type = None

    case_list1 = Path(original_study) / "case_lists"
    case_list2 = Path(new_study) / "case_lists"
    cna_1 = case_list1 / "cases_cna.txt"
    cna_2 = case_list2 / "cases_cna.txt"
    sequenced_1 = case_list1 / "cases_sequenced.txt"
    sequenced_2 = case_list2 / "cases_sequenced.txt"
    sv_1 = case_list1 / "cases_sv.txt"
    sv_2 = case_list2 / "cases_sv.txt"

    (
        old_samples_cna,
        removed_samples_cna,
        left_samples_cna,
        total_patients,
        total_samples,
    ) = compare_sample_file_remove(cna_1, cna_2, new_study)

    old_samples_sequenced, removed_samples_sequenced, left_samples_sequenced, _, _ = (
        compare_sample_file_remove(sequenced_1, sequenced_2, new_study))

    old_samples_sv, removed_samples_sv, left_samples_sv, _, _ = (
        compare_sample_file_remove(sv_1, sv_2, new_study))


    html_content = f"""
    <!DOCTYPE html>
    <html lang="it">
    <head>
        <link rel="stylesheet" type="text/css" href="{Path('img') / 'styles.css'}">
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Varan - Remove</title>
    </head>
    <body>
        <header>
            <img src="{Path('img') / 'logo_VARAN.png'!s}" alt="Logo Varan">
            <h1>VARAN - Remove</h1>
        </header>

        <h2>Generate on {date}</h2>
        {start_time_html}

        <div class="container">
            <div class="section-title">General Information</div>
                <div class="content">
                    <p><strong>COMMAND LINE:</strong> {" ".join(sys.argv)}</p>
                    <p><strong>ORIGINAL STUDY:</strong> {Path(original_study).name}</p>
                    <p><strong>NEW STUDY:</strong> {Path(new_study).name}</p>"""

    if sample_type:
        html_content += f"""<p><strong>SAMPLE TYPE:</strong> {sample_type}</p>"""

    if cancer_type:
        html_content += f"""<p><strong>CANCER TYPE:</strong> {cancer_type}</p>"""

    html_content += f"""
                    <hr width="100%" size="2" color="#003366" noshade>
                    <p><strong>Total Patients:</strong> {total_patients}</p>
                    <p><strong>Total Samples:</strong> {total_samples}</p>
                </div>"""

    source_provenance = read_provenance(original_study)
    if source_provenance:
        write_provenance(new_study, **source_provenance)
    html_content += render_provenance_html("Source study", source_provenance)

    html_content += build_warnings_section(new_study, ghosts)

    html_content += """
            <div class="section-title">Detailed Overview</div>
                <div class="content">
                    <p><strong>Cases_CNA:</strong></p>"""

    if not old_samples_cna:
        html_content += f"""
            <p>&emsp;{Path(original_study).name}'s cases_cna was empty.</p>"""
    else:
        html_content += f"""
            <p>&emsp;<strong>Removed:</strong> {len(removed_samples_cna)} samples
            ({", ".join(removed_samples_cna)})</p>
            <p>&emsp;There is/are now {left_samples_cna} samples in cases_cna.</p>"""

    html_content += """
        <p><strong>Cases_sequenced:</strong></p>"""

    if not old_samples_sequenced:
        html_content += f"""
            <p>&emsp;{Path(original_study).name}'s cases_sequenced was empty.</p>"""
    else:
        html_content += f"""
        <p>&emsp;<strong>Removed:</strong> {len(removed_samples_sequenced)} samples
        ({', '.join(removed_samples_sequenced)})</p>
        <p>&emsp;There is/are now {left_samples_sequenced}
        samples in cases_sequenced.</p>
    """

    html_content += """
        <p><strong>Cases_sv:</strong></p>"""

    if not old_samples_sv:
        html_content += f"""
            <p>&emsp;{Path(original_study).name}'s cases_sv was empty.</p>"""
    else:
        html_content += (
            f"""
            <p>&emsp;<strong>Removed:</strong> {len(removed_samples_sv)} samples
            ({", ".join(removed_samples_sv)})</p>
            <p>&emsp;There is/are now {left_samples_sv} samples in cases_sv.</p>""")

    if filters != {}:
        html_content += """
            <section class="filters">
                <div class="section-title">Filters & Configuration</div>"""

    if any(filt in filters for filt in ["ALT"]):
        html_content += """
            <div class="subtitle">VCF Filters</div>"""

    if "ALT" in filters:
        html_content += """
                <div class="content">
                <p>Keep only rows where <strong>ALT</strong> is different from "."</p>
                <p>Keep only rows where <strong>FILTER</strong> = "PASS"</p>
            </div>"""

    if any(filt in filters for filt in ["T_VAF_MIN", "T_VAF_MAX", "AF"]):
        html_content += """
            <div class="subtitle">VAF Filters</div>"""

    if any(filt in filters for filt in ["T_VAF_MIN", "T_VAF_MAX"]):
        html_content += f"""
            <div class="content">
            <p><strong>T_VAF_MIN</strong>: {filters["T_VAF_MIN"]}</p>
            <p><strong>T_VAF_MAX</strong>: {filters["T_VAF_MAX"]}</p>
        </div>"""
        if "T_VAF_MIN_NOVEL" in filters:
            html_content += f"""
                <div class="content">
                <p><strong>T_VAF_MIN_NOVEL</strong>: {filters["T_VAF_MIN_NOVEL"]}</p>
            </div>"""

    if "AF" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>AF</strong>: {filters["AF"]}</p>
            </div>"""

    if any(
        filt in filters
        for filt in [
            "ONCOKB", "CONSEQUENCES", "POLYPHEN",
            "CLIN_SIG", "IMPACT", "SIFT",
        ]
    ):
        html_content += """
            <div class="subtitle">MAF Filters</div>"""

    if "ONCOKB" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>ONCOKB</strong>: {filters["ONCOKB"]}</p>
            </div>"""

    if "CONSEQUENCES" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>CONSEQUENCES</strong>: {filters["CONSEQUENCES"]}</p>
            </div>"""

    if "POLYPHEN" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>POLYPHEN</strong>: {filters["POLYPHEN"]}</p>
            </div>"""

    if "CLIN_SIG" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>CLIN_SIG</strong>: {filters["CLIN_SIG"]}</p>
            </div>"""

    if "IMPACT" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>IMPACT</strong>: {filters["IMPACT"]}</p>
            </div>"""

    if "SIFT" in filters:
        html_content += f"""
                <div class="content">
                <p><strong>SIFT</strong>: {filters["SIFT"]}</p>
            </div>"""

    keys_to_check = {"PLOIDY", "CNVKIT_algorithm", "THRESHOLD_TMB",
    "THRESHOLD_SITES", "THRESHOLD_MSI", "THRESHOLD_FUSION"}

    if filters != {} and all(key in filters for key in keys_to_check):
        html_content += f"""
                <div class="subtitle">Copy Number Alterations (CNA)</div>
                <div class="content">
                    <p><strong>PLOIDY</strong> = {filters["PLOIDY"]}</p>
                    <p><strong>CNVKIT Algorithm</strong> = {filters["CNVKIT_algorithm"]}</p>
                </div>

                <div class="subtitle">Tumor Mutational Burden (TMB)</div>
                <div class="content">
                    <p><strong>THRESHOLD_TMB</strong>: {filters["THRESHOLD_TMB"]}</p>
                </div>

                <div class="subtitle">Microsatellite Instability (MSI)</div>
                <div class="content">
                    <p><strong>
                        THRESHOLD_SITES</strong>: {filters["THRESHOLD_SITES"]}
                    </p>
                    <p><strong>THRESHOLD_MSI</strong>: {filters["THRESHOLD_MSI"]}</p>
                </div>

                <div class="subtitle">Fusions</div>
                <div class="content">
                    <p><strong>
                        THRESHOLD_FUSION</strong>: {filters["THRESHOLD_FUSION"]}
                    </p>
                </div>
            </section>"""

    if versioning.old_version_exists and actual_version != 1:
            html_content += f"""
            <section class="comparison">
                <div class="section-title">
                    Comparison with Previous Version (_v{actual_version - 1})
                </div>
                <div class="content">
                    <p><strong>ADDED PATIENT(S):</strong> {len(only_new_pat)}</p>
                    <p><strong>ADDED SAMPLE(S):</strong> {len(only_new_sam)}</p>
                    <p><strong>REMOVED PATIENT(S):</strong> {len(only_old_pat)}</p>
                    <p><strong>REMOVED SAMPLE(S):</strong> {len(only_old_sam)}</p>
                </div>
            </section>
            """

    if ((Path(new_study) / general_graph_path).exists()
        or (Path(new_study) / genes_graph_path).exists()):
        html_content += f"""
            <section class="graphs">
                <div class="section-title">
                    Graphical Overview {graph_expression}
                </div>"""

    if (Path(new_study) / general_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{general_graph_path}" alt="Samples and Patients barchart">
            </div>"""

    if (Path(new_study) / genes_graph_path).exists():
        html_content += f"""
                <div class="content">
                <img src="{genes_graph_path}" alt="SNV, CNV and Fusions barchart">
            </div>"""

    if ((Path(new_study) / general_graph_path).exists()
        or (Path(new_study) / genes_graph_path).exists()):
        html_content += """</section>"""

    if annotations_list != []:
        html_content += """
            <section class="annotations">
                <div class="section-title">Annotations</div>
                    <div class="content">
                        <ul>"""

        for annotation in annotations_list:
            html_content += f"""
                         <li><p>{annotation}</p></li>"""

    if annotations_list != []:
        html_content += """
                    </ul>
                </div>
            </div>
        </section>"""

    html_content += """
    </body>
    </html>
    """

    report_path = Path(new_study) / "report_VARAN.html"
    with report_path.open("w") as f:
        f.write(html_content)


def compare_sample_file_remove(
    file1: str | Path,
    file2: str | Path,
    outputfolder: str | Path,
) -> tuple[list[str] | None, list[str] | None, list[str] | int | None, int, int]:
    """Identify removed samples in two studies and summarize the updated clinical data.

    Args:
        file1 (Union[str, Path]): Path to the original (old) sample file.
        file2 (Union[str, Path]): Path to the new sample file.
        outputfolder (Union[str, Path]): Folder containing the new clinical sample data
        file.

    Returns:
        Tuple:
            - samples_file1 (list[str] | None)
            - removed_samples (list[str] | None): Samples that were in the old file
            but not in the new file,
              or None if both files are missing.
            - left_samples (Optional[int]): Number of samples left in the new file, or
            None if both files are missing.
            - total_patients (int): Number of unique patients in the new clinical sample
            data.
            - total_samples (int): Number of unique samples in the new clinical sample
            data.

    """
    outputfolder = Path(outputfolder)
    clin_sam_new_df = pd.read_csv(outputfolder / "data_clinical_sample.txt", sep="\t")

    if Path(file1).exists() and Path(file2).exists():
        samples_file1 = extract_sample_list(file1)
        samples_file2 = extract_sample_list(file2)

        removed_samples = [
            sample for sample in samples_file1
            if sample not in samples_file2 and sample != ""
        ]
        left_samples = len(set(samples_file1) - set(removed_samples))

    elif Path(file1).exists() and not Path(file2).exists():
        samples_file1 = extract_sample_list(file1)
        removed_samples = extract_sample_list(file1)
        left_samples = 0

    else:
        samples_file1 = None
        removed_samples = None
        left_samples = None

    return (
        samples_file1,
        removed_samples,
        left_samples,
        len(set(clin_sam_new_df.iloc[4:, 1])),
        len(set(clin_sam_new_df.iloc[4:, 0])),
    )


def extract_filters_from_html(report: str) -> dict | None:
    """Extract filtering criteria from a given HTML report file.

    This function parses the specified HTML report to extract filtering parameters
    listed under the filters section and VCF filters subsection. It returns these
    filters as a dictionary where keys are filter names and values are their criteria.

    Args:
        report (str): Path to the HTML report file.

    Returns:
        dict: A dictionary containing filter names as keys and their corresponding
              filter criteria as values.

    """
    filters = {}
    report_path = Path(report)
    with report_path.open(encoding="utf-8") as file:
        html_content = file.read()

    fs_match = re.search(
        r'<section class="filters">.*?</section>', html_content, re.DOTALL)
    if fs_match:
        section_text = fs_match.group()
        p_items = re.findall(
            r"<p>\s*<strong>([^<]+?)[:]*</strong>\s*[:=]?\s*(.*?)\s*</p>",
            section_text, re.DOTALL)
        for key, value in p_items:
            filters[key.strip().rstrip(":")] = value.strip()
        other_items = re.findall(
            r"<strong>([^<]+?)</strong>\s*[:=]?\s*(.*?)(?=<br>|</div>|</section>)",
            section_text, re.DOTALL)
        for key, value in other_items:
            k = key.strip().rstrip(":")
            if k not in filters:
                filters[k] = value.strip()

    vcf_match = re.search(
        r'<div class="subtitle">VCF Filters</div>\s*<div class="content">(.*?)</div>',
        html_content, re.DOTALL)
    if vcf_match:
        content = vcf_match.group(1)
        alt_match = re.search(
            (r'Keep only rows where\s*<strong>\s*ALT\s*</strong>\s*'
            r'is different from\s*"(.*?)"'),
            content, re.DOTALL)
        if alt_match:
            alt_val = alt_match.group(1).strip()
            filters["ALT"] = f'!="{alt_val}"'

    return filters


def extract_cancer_type_from_html(report: str | Path) -> str | None:
    """Extract the cancer type information from an HTML report file.

    This function searches the HTML report for the section indicating the cancer type
    and returns its value if found.

    Args:
        report (str): Path to the HTML report file.

    Returns:
        str or None: The extracted cancer type as a string, or None if not found.

    """
    cancer_type = None

    report = Path(report)
    with report.open(encoding="utf-8") as file:
        html_content = file.read()

    cancer_type_match = re.search(
    r"<p><strong>\s*CANCER TYPE:\s*</strong>\s*(.*?)\s*</p>",
    html_content, re.DOTALL | re.IGNORECASE)

    if cancer_type_match:
        cancer_type = cancer_type_match.group(1).strip()

    return cancer_type


def extract_sample_type_from_html(report: str | Path) -> str | None:
    """Extract the sample type information from an HTML report file.

    This function searches the HTML report for the section indicating the sample type
    and returns its value if found.

    Args:
        report (str | Path): Path to the HTML report file.

    Returns:
        str or None: The extracted sample type as a string, or None if not found.
    """
    sample_type = None
    report = Path(report)

    try:
        with report.open(encoding="utf-8") as file:
            html_content = file.read()

        sample_type_match = re.search(
            r"<p><strong>\s*SAMPLE TYPE:\s*</strong>\s*(.*?)\s*</p>",
            html_content, 
            re.DOTALL | re.IGNORECASE
        )

        if sample_type_match:
            sample_type = sample_type_match.group(1).strip()

    except Exception as e:
        return None

    return sample_type
