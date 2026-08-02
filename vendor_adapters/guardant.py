"""Guardant360 (FPG360) vendor adapter: converts a per-sample VCF +
`.cnv_call.hdr.tsv` + `.msi_call.hdr.tsv` + `_finalmetadata.xml` (fetched
from an S3 run folder) into the sample.tsv / {run_id}_fusions.tsv shape
Varan's own `varan.py -i` flag consumes.

This is the first (reference) vendor adapter under `vendor_adapters/`,
moved here from the formerly-standalone `create_Varan_input.py` script.
Everything vendor-specific (Guardant's VCF dialect, its XML metadata
shape, its MSI/CNV report columns, its S3 layout conventions incl. the
"_2" backup-folder fallback) lives in this module. Generic, reusable
pieces (S3 listing, the ONCOTREE dict.csv lookup, the sample.tsv writer,
the Varan FUSIONS/*.tsv writer) live in `vendor_adapters/common.py` - see
that module's docstring for why the split is drawn there.

Every documented, verified bug fix from the original review is preserved
here unchanged (see MULTIVENDOR_INTEGRATION_NOTES.md for the full list):
  1. SNV VCF only carries FILTER=PASS rows (germline contamination fix).
  2. The dead `cn_value == 2.0` check is gone (never fired; the real
     filtering already happens via the `call` column in process_vcf()).
  3/4. SVTYPE=BND (fusion breakend) rows are explicitly skipped before the
     structural/SNV split, so they neither leak into the SNV VCF nor
     desync the sequential gene<->copy-number index pairing.
  5. MSI carries the real numeric msi_score; MSI_THR is derived from
     msi_status via substring match, left blank (not guessed) for any
     unrecognized status string.
  6. Fusions are written via `common.append_fusions_to_table()` into a
     plain FUSIONS/*.tsv file, never through a synthesized
     CombinedVariantOutput.tsv / comb_path (comb_path is always blank).

Still-open questions carried forward from the original notes (NOT resolved
here - see MULTIVENDOR_INTEGRATION_NOTES.md):
  - `fill_fusion_from_temp()`'s hardcoded `min_read_count=15` in walk.py is
    independent of both conf.ini's THRESHOLD_FUSION and Guardant's own
    `call=1` confidence flag - a real, Guardant-confirmed fusion with <15
    supporting molecules is silently dropped. Needs a human call on
    whether 15 is the right cutoff for Guardant's molecule-count scale.
  - ONCOTREE_CODE mapping via dict.csv: unchanged from the original
    script's logic, not reviewed further.
  - "Golden" extra fields (BRCA Reversion, AutoQC metrics): the mechanism
    exists (see common.SampleRow.extra / write_sample_tsv), but the real
    board-summary/AutoQC file to verify column names against was never
    available.
"""

import csv
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from vendor_adapters.common import (
    SampleRow,
    append_fusions_to_table,
    get_incremental_report_path,
    list_s3_files,
    load_oncotree_dict,
    run_cmd,
    write_sample_tsv,
)

NAME = "guardant"

# --- DEFAULT PATHS (deployment-specific; override via run()'s kwargs for
# a different environment/test setup rather than editing these) ---
DICT_PATH = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/Preprocessing_Files/dict.csv"
REPORT_BASE_DIR = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/Varan_Input"
VCF_BASE_DIR = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/processed_VCF"
TEMP_LOCAL_DIR = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/tmp"

# --- VCF HEADER TEMPLATES ---
COMMON_HEADER = [
    "##fileformat=VCFv4.2\n", "##reference=file://hashtable/reference.bin\n",
    "##contig=<ID=chr1,length=249250621>\n", "##contig=<ID=chr2,length=243199373>\n",
    "##contig=<ID=chr3,length=198022430>\n", "##contig=<ID=chr4,length=191154276>\n",
    "##contig=<ID=chr5,length=180915260>\n", "##contig=<ID=chr6,length=171115067>\n",
    "##contig=<ID=chr7,length=159138663>\n", "##contig=<ID=chr8,length=146364022>\n",
    "##contig=<ID=chr9,length=141213431>\n", "##contig=<ID=chr10,length=135534747>\n",
    "##contig=<ID=chr11,length=135006516>\n", "##contig=<ID=chr12,length=133851895>\n",
    "##contig=<ID=chr13,length=115169878>\n", "##contig=<ID=chr14,length=107349540>\n",
    "##contig=<ID=chr15,length=102531392>\n", "##contig=<ID=chr16,length=90354753>\n",
    "##contig=<ID=chr17,length=81195210>\n", "##contig=<ID=chr18,length=78077248>\n",
    "##contig=<ID=chr19,length=59128983>\n", "##contig=<ID=chr20,length=63025520>\n",
    "##contig=<ID=chr21,length=48129895>\n", "##contig=<ID=chr22,length=51304566>\n",
    "##contig=<ID=chrX,length=155270560>\n", "##contig=<ID=chrY,length=59373566>\n",
    "##contig=<ID=chrM,length=16569>\n",
    "##ALT=<ID=CNV,Description='Copy number variant region'>\n",
    "##ALT=<ID=DEL,Description='Deletion relative to the reference'>\n",
    "##ALT=<ID=DUP,Description='Region of elevated copy number relative to the reference'>\n",
]
SNV_SPECIFIC = [
    "##INFO=<ID=DP,Number=1,Type=Integer,Description='Approximate read depth'>\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description='Allelic depths'>\n",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description='Allele fractions'>\n",
    "##FILTER=<ID=PASS,Description='Pass'>\n",
]
CNV_SPECIFIC = [
    "##INFO=<ID=REFLEN,Number=1,Type=Integer,Description='REF length'>\n",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description='Type of structural variant'>\n",
    "##INFO=<ID=END,Number=1,Type=Integer,Description='End position'>\n",
    "##INFO=<ID=SEGID,Number=1,Type=String,Description='Segment ID'>\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n",
    "##FORMAT=<ID=CN,Number=1,Type=Float,Description='Copy Number'>\n",
    "##FORMAT=<ID=SM,Number=1,Type=Float,Description='Fold Change (mapped to SM for compatibility)'>\n",
]


def get_xml_data(xml_path: str, onco_info: Tuple[Dict[str, str], set]) -> Tuple[str, str]:
    """Parse Guardant's `_finalmetadata.xml`/`_metadata.xml` for
    SubjectId/AccessionId and Diagnosis, mapping Diagnosis to an ONCOTREE
    code via the shared dict.csv lookup. Unchanged from the original
    script's logic - not reviewed further this pass (see module docstring
    open questions)."""
    name_to_code, valid_codes = onco_info
    p_id, o_code = "N/A", "UNKNOWN"
    try:
        root = ET.parse(xml_path).getroot()
        for param in root.iter():
            if any(c.tag.split("}")[-1] == "Name" and c.text == "accessionID" for c in param):
                for c in param:
                    if c.tag.split("}")[-1] == "Value":
                        p_id = c.text
                        break
        for diag in root.iter():
            if diag.tag.split("}")[-1] == "Diagnosis" and diag.text:
                raw_diag = diag.text.strip()
                diag_lower = raw_diag.lower()
                diag_upper = raw_diag.upper()
                if diag_lower in name_to_code:
                    o_code = name_to_code[diag_lower]
                elif diag_upper in valid_codes:
                    o_code = diag_upper
                else:
                    o_code = f"NOT_IN_DICT ({raw_diag})"
                break
    except Exception:
        pass
    return p_id, o_code


def get_msi_data(msi_path: Optional[str]) -> Dict[str, str]:
    """Read Guardant's .msi_call.hdr.tsv - has a REAL numeric score, not
    just a category.

    Returns the numeric msi_score as-is, plus msi_status normalized to
    Varan's own Stable/Unstable vocabulary (substring match against
    MSI-H/MSS/MSI-L/STABLE/UNSTABLE, case-insensitive; left blank rather
    than guessed for anything unrecognized - see module docstring).
    """
    data = {"run_id": "UNKNOWN_RUN", "msi_score": "", "msi_thr": ""}
    if msi_path and os.path.exists(msi_path):
        try:
            with open(msi_path, "r") as f:
                r = csv.DictReader(f, delimiter="\t")
                for row in r:
                    row = {k.strip(): v.strip() for k, v in row.items() if k}
                    data["run_id"] = row.get("runid", "UNKNOWN_RUN")
                    data["msi_score"] = row.get("msi_score", "")
                    status = row.get("msi_status", "").upper()
                    if "MSI-H" in status or "UNSTABLE" in status:
                        data["msi_thr"] = "Unstable"
                    elif "MSS" in status or "MSI-L" in status or "STABLE" in status:
                        data["msi_thr"] = "Stable"
                    # else: leave msi_thr empty rather than guess - an
                    # unrecognized status string should not silently
                    # produce a wrong clinical call.
                    break
        except Exception:
            pass
    return data


def load_fusions(fusion_path: Optional[str]) -> List[Dict[str, str]]:
    """Read Guardant's .fusion_call.hdr.tsv, keeping only `call == "1"`
    (Guardant-confirmed) rows."""
    fusions: List[Dict[str, str]] = []
    if fusion_path and os.path.exists(fusion_path):
        try:
            with open(fusion_path, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if row.get("call") == "1":
                        fusions.append({
                            "gene_a": row["gene_a"],
                            "gene_b": row["gene_b"],
                            "supp": row.get("fusion_molecule_count_ab", "0"),
                        })
        except Exception:
            pass
    return fusions


def load_cnv_tsv_ordered(tsv_path: Optional[str]) -> List[Dict]:
    """Read Guardant's .cnv_call.hdr.tsv, one row per gene, in file order
    (order matters: process_vcf() pairs these back to VCF structural rows
    by sequential index).

    The original `if cn_value == 2.0: continue` filter is intentionally
    gone: copy_number is a continuous value (e.g. 2.07, 1.84, 3.17) that's
    essentially never exactly 2.0, so it never excluded anything - the
    real filtering already happens in process_vcf() via the `call` column
    (0 = no significant call, 1/2 = deletion/amplification).
    """
    cnv_list: List[Dict] = []
    if not tsv_path or not os.path.exists(tsv_path):
        return cnv_list
    with open(tsv_path, "r") as f:
        lines = f.readlines()
        start = next(
            (i for i, l in enumerate(lines)
             if "gene" in l.lower() and "copy_number" in l.lower()), -1)
        if start == -1:
            return cnv_list
        f.seek(0)
        [next(f) for _ in range(start)]
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                cnv_list.append({
                    "gene": row["gene"].strip(),
                    "cn": float(row["copy_number"]),
                    "call": row["call"].strip(),
                })
            except Exception:
                continue
    return cnv_list


def process_vcf(vcf_in: str, cnv_tsv_in: Optional[str], snv_out: str,
                 cnv_out: str, sample_id: str) -> None:
    """Split a Guardant VCF into a PASS-only SNV VCF and a per-gene CNV
    VCF, pairing structural rows to `.cnv_call.hdr.tsv` rows by sequential
    index.

    Preserves (unchanged from the reviewed/corrected version):
      - FILTER=PASS-only SNV rows (germline-contamination fix).
      - Explicit SVTYPE=BND skip, checked BEFORE the structural/SNV split,
        so breakend/fusion-junction rows neither leak into the SNV VCF as
        pseudo point-mutations nor silently consume a cnv_index slot meant
        for a real CNV row.
      - is_struct narrowed to SVTYPE=CNV/DUP/DEL (not the original overly
        broad "SVTYPE=" in info_str, which also matched BND).
    """
    ordered_cnv_data = load_cnv_tsv_ordered(cnv_tsv_in)
    cnv_index = 0
    col_header = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n"

    with open(vcf_in, "r") as f_in, open(snv_out, "w") as f_snv, open(cnv_out, "w") as f_cnv:
        f_snv.writelines(COMMON_HEADER + SNV_SPECIFIC + [col_header])
        f_cnv.writelines(COMMON_HEADER + CNV_SPECIFIC + [col_header])

        for line in f_in:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue

            chrom = cols[0]
            if not chrom.startswith("chr"):
                chrom = "chrM" if chrom == "MT" else f"chr{chrom}"
            cols[0] = chrom

            while len(cols) < 10:
                cols.append(".")
            alt, info_str = cols[4], cols[7]

            is_struct = any(x in alt for x in ["<CNV>", "<DUP>", "<DEL>"]) or any(
                svtype in info_str for svtype in ("SVTYPE=CNV", "SVTYPE=DUP", "SVTYPE=DEL"))

            # BND (breakend) rows are fusion junction evidence, already
            # captured properly via load_fusions()/.fusion_call.hdr.tsv -
            # neither a CNV call nor a point mutation. Must be checked
            # BEFORE the is_struct branch below.
            if "SVTYPE=BND" in info_str:
                continue

            if is_struct:
                if cnv_index < len(ordered_cnv_data):
                    d = ordered_cnv_data[cnv_index]
                    fc = round(d["cn"] / 2.0, 4)

                    cols[5] = "1"

                    call_val = str(d["call"]).strip()
                    if call_val in ["1", "2"]:
                        cols[4] = "<DUP>"
                        cols[6] = "PASS"
                    else:
                        cols[4] = "."
                        cols[6] = "FAIL"

                    original_end = next(
                        (x.split("=")[1] for x in info_str.split(";") if x.startswith("END=")),
                        str(int(cols[1]) + 1))
                    cols[7] = f"SVTYPE=CNV;END={original_end};SEGID={d['gene']}"
                    cols[8], cols[9] = "GT:CN:SM", f"0/1:{d['cn']}:{fc}"
                    cnv_index += 1
                else:
                    cols[4] = "."
                    cols[5] = "1"
                    cols[6] = "FAIL"
                f_cnv.write("\t".join(cols) + "\n")
            else:
                if cols[6] != "PASS":
                    continue
                try:
                    format_keys = cols[8].split(":")
                    sample_vals = cols[9].split(":")
                    if "AD" in format_keys:
                        idx_ad = format_keys.index("AD")
                        ad_parts = sample_vals[idx_ad].split(",")
                        ref_count, alt_count = int(ad_parts[0]), int(ad_parts[1])
                        depth = ref_count + alt_count
                        af = round(alt_count / depth, 4) if depth > 0 else 0
                        cols[8], cols[9] = "GT:AD:AF:DP", f"{sample_vals[0]}:{ref_count},{alt_count}:{af}:{depth}"
                except Exception:
                    pass
                f_snv.write("\t".join(cols) + "\n")


def process_single_sample(
    sid: str, s3_folder: str, run_id_default: str,
    onco_dict: Tuple[Dict[str, str], set], fusion_table_path: str,
    xml_folder: Optional[str] = None, *,
    vcf_base_dir: str = VCF_BASE_DIR, temp_local_dir: str = TEMP_LOCAL_DIR,
) -> Optional[SampleRow]:
    """Fetch one Guardant sample's files from S3, convert them, and return
    a `SampleRow` (or None if the required VCF/metadata can't be found)."""
    if xml_folder is None:
        xml_folder = s3_folder

    s3_files = list_s3_files(s3_folder)
    xml_files = list_s3_files(xml_folder) if xml_folder != s3_folder else s3_files

    clean_sid = sid.rstrip("_")

    def find_file_in_list(suffix, files_list):
        for f in files_list if files_list else []:
            if f.endswith(suffix) and clean_sid in f and not f.startswith("AIO"):
                return f
        return None

    xml_name = find_file_in_list("_finalmetadata.xml", xml_files)

    if not xml_name and xml_folder == s3_folder:
        backup_folder = s3_folder.rstrip("/") + "_2"
        print(f"[{sid}] Metadata not found in {s3_folder}. Trying backup path: {backup_folder}")
        backup_files = list_s3_files(backup_folder)
        xml_name = find_file_in_list("_finalmetadata.xml", backup_files)
        if xml_name:
            xml_folder = backup_folder
            print(f"[{sid}] Metadata FOUND in backup folder: {xml_name}")
        else:
            print(f"[{sid}] Metadata not found in any path.")
            return None

    # NOTE: no "comb" (comb_path/CombinedVariantOutput) target - fusions
    # are routed through common.append_fusions_to_table() instead, see
    # module docstring.
    file_targets = {
        "xml": (xml_name, xml_folder),
        "vcf": (find_file_in_list(".vcf", s3_files), s3_folder),
        "msi": (find_file_in_list(".msi_call.hdr.tsv", s3_files), s3_folder),
        "cnv": (find_file_in_list(".cnv_call.hdr.tsv", s3_files), s3_folder),
        "fus": (find_file_in_list(".fusion_call.hdr.tsv", s3_files), s3_folder),
    }

    local_files = {}
    for key, target in file_targets.items():
        if target:
            fname, folder_path = target
            if fname:
                local_p = os.path.join(temp_local_dir, f"{sid}_{fname}")
                if run_cmd(f"aws s3 cp {folder_path.rstrip('/')}/{fname} {local_p}"):
                    local_files[key] = local_p

    if "vcf" not in local_files:
        print(f"VCF missing for {sid}")
        return None

    p_id, o_code = get_xml_data(local_files["xml"], onco_dict)
    m_info = get_msi_data(local_files.get("msi"))
    run_id_final = m_info["run_id"] if m_info["run_id"] != "UNKNOWN_RUN" else run_id_default

    v_out = os.path.join(vcf_base_dir, run_id_final, sid)
    os.makedirs(v_out, exist_ok=True)

    snv_f, cnv_f = os.path.join(v_out, f"{sid}.snv.vcf"), os.path.join(v_out, f"{sid}.cnv.vcf")

    process_vcf(local_files["vcf"], local_files.get("cnv"), snv_f, cnv_f, sid)

    fusions_list = load_fusions(local_files.get("fus"))
    append_fusions_to_table(fusion_table_path, sid, fusions_list)

    for p in local_files.values():
        if os.path.exists(p):
            os.remove(p)

    # comb_path is always blank and TMB/TMB_THR are always blank: Guardant360
    # CDx has no TMB field anywhere in either metadata XML variant checked.
    return SampleRow(
        sample_id=sid, patient_id=p_id, run_id=run_id_final,
        oncotree_code=o_code, snv_path=snv_f, cnv_path=cnv_f, comb_path="",
        msi=m_info["msi_score"], tmb="", msi_thr=m_info["msi_thr"], tmb_thr="",
    )


def run(
    folder: Optional[str] = None, selection: Optional[str] = None, *,
    dict_path: str = DICT_PATH, report_base_dir: str = REPORT_BASE_DIR,
    vcf_base_dir: str = VCF_BASE_DIR, temp_local_dir: str = TEMP_LOCAL_DIR,
) -> Optional[Dict[str, object]]:
    """Convert one Guardant run (`folder`, an S3 run folder) or a batch of
    samples pulled from arbitrary run folders (`selection`, a TSV with
    `sample_id`/`s3_path_run` columns) into a sample.tsv + fusions.tsv
    pair.

    Plain function, explicit inputs/outputs (no argparse/sys.argv/module
    globals involved beyond the path defaults above) so this can be called
    directly - from a test, a notebook, or eventually a Snakemake rule -
    without going through the CLI wrapper in `create_Varan_input.py`.

    Returns a dict with `report_path`, `fusion_table_path`, and
    `report_data` (the list of `SampleRow`s written), or None if no sample
    produced any data.
    """
    onco_dict = load_oncotree_dict(dict_path)
    os.makedirs(temp_local_dir, exist_ok=True)
    os.makedirs(report_base_dir, exist_ok=True)

    report_data: List[SampleRow] = []
    main_run_id = "MANUAL_BATCH"
    fusion_table_path = None

    if folder:
        s3_folder = folder.rstrip("/")
        print(f"Analyzing folder: {s3_folder}")
        main_run_id = os.path.basename(s3_folder).split(".")[0]
        s3_files = list_s3_files(s3_folder)
        xml_files = [f for f in s3_files if f.endswith("_finalmetadata.xml") and not f.startswith("AIO")]

        xml_folder = s3_folder
        if not xml_files:
            backup_folder = s3_folder + "_2"
            print(f"No metadata found in {s3_folder}. Searching backup folder: {backup_folder}...")
            backup_files = list_s3_files(backup_folder)
            xml_files = [f for f in backup_files if f.endswith("_finalmetadata.xml") and not f.startswith("AIO")]
            if xml_files:
                xml_folder = backup_folder
                print(f"Found {len(xml_files)} metadata file(s) in the backup folder.")

        fusion_table_path = os.path.join(report_base_dir, f"{main_run_id}_fusions.tsv")
        for xml_f in xml_files:
            sid = xml_f.replace("_finalmetadata.xml", "").split("_")[-1]
            print(f"Processing {sid}...")
            res = process_single_sample(
                sid, s3_folder, main_run_id, onco_dict, fusion_table_path, xml_folder,
                vcf_base_dir=vcf_base_dir, temp_local_dir=temp_local_dir)
            if res:
                report_data.append(res)

    elif selection:
        print(f"Analyzing selection: {selection}")
        main_run_id = Path(selection).stem
        fusion_table_path = os.path.join(report_base_dir, f"{main_run_id}_fusions.tsv")
        with open(selection, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sid, s3_folder = row["sample_id"], row["s3_path_run"].rstrip("/")
                run_id_tmp = os.path.basename(s3_folder).split(".")[0]
                res = process_single_sample(
                    sid, s3_folder, run_id_tmp, onco_dict, fusion_table_path,
                    vcf_base_dir=vcf_base_dir, temp_local_dir=temp_local_dir)
                if res:
                    report_data.append(res)

    if not report_data:
        print("\nNo data collected.")
        return None

    r_path = get_incremental_report_path(report_base_dir, main_run_id)
    write_sample_tsv(report_data, r_path)
    print(f"\nReport generated: {r_path}")
    print(f"Fusions (if any) in: {fusion_table_path}")
    print("Run Varan with: python varan.py -i "
          f"{r_path} <patient.tsv or \"\"> {fusion_table_path} -o ... -c ...")
    return {
        "report_path": r_path,
        "fusion_table_path": fusion_table_path,
        "report_data": report_data,
    }
