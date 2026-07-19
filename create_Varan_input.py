import os
import sys
import csv
import argparse
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

# --- CONFIGURAZIONE PERCORSI ---
DICT_PATH = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/Preprocessing_Files/dict.csv"
REPORT_BASE_DIR = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/Varan_Input"
VCF_BASE_DIR = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/processed_VCF"
TEMP_LOCAL_DIR = "/data/data_storage/novaseq_results/research/CbioPortal/FPG360/tmp"

# --- HEADER VCF ---
COMMON_HEADER = ["##fileformat=VCFv4.2\n", "##reference=file://hashtable/reference.bin\n", "##contig=<ID=chr1,length=249250621>\n", "##contig=<ID=chr2,length=243199373>\n", "##contig=<ID=chr3,length=198022430>\n", "##contig=<ID=chr4,length=191154276>\n", "##contig=<ID=chr5,length=180915260>\n", "##contig=<ID=chr6,length=171115067>\n", "##contig=<ID=chr7,length=159138663>\n", "##contig=<ID=chr8,length=146364022>\n", "##contig=<ID=chr9,length=141213431>\n", "##contig=<ID=chr10,length=135534747>\n", "##contig=<ID=chr11,length=135006516>\n", "##contig=<ID=chr12,length=133851895>\n", "##contig=<ID=chr13,length=115169878>\n", "##contig=<ID=chr14,length=107349540>\n", "##contig=<ID=chr15,length=102531392>\n", "##contig=<ID=chr16,length=90354753>\n", "##contig=<ID=chr17,length=81195210>\n", "##contig=<ID=chr18,length=78077248>\n", "##contig=<ID=chr19,length=59128983>\n", "##contig=<ID=chr20,length=63025520>\n", "##contig=<ID=chr21,length=48129895>\n", "##contig=<ID=chr22,length=51304566>\n", "##contig=<ID=chrX,length=155270560>\n", "##contig=<ID=chrY,length=59373566>\n", "##contig=<ID=chrM,length=16569>\n", "##ALT=<ID=CNV,Description='Copy number variant region'>\n", "##ALT=<ID=DEL,Description='Deletion relative to the reference'>\n", "##ALT=<ID=DUP,Description='Region of elevated copy number relative to the reference'>\n"]
SNV_SPECIFIC = ["##INFO=<ID=DP,Number=1,Type=Integer,Description='Approximate read depth'>\n", "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n", "##FORMAT=<ID=AD,Number=R,Type=Integer,Description='Allelic depths'>\n", "##FORMAT=<ID=AF,Number=A,Type=Float,Description='Allele fractions'>\n", "##FILTER=<ID=PASS,Description='Pass'>\n"]
CNV_SPECIFIC = ["##INFO=<ID=REFLEN,Number=1,Type=Integer,Description='REF length'>\n", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description='Type of structural variant'>\n", "##INFO=<ID=END,Number=1,Type=Integer,Description='End position'>\n", "##INFO=<ID=SEGID,Number=1,Type=String,Description='Segment ID'>\n", "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n", "##FORMAT=<ID=CN,Number=1,Type=Float,Description='Copy Number'>\n", "##FORMAT=<ID=SM,Number=1,Type=Float,Description='Fold Change (mapped to SM for compatibility)'>\n"]

# Varan's own data_sv.txt column shape (Sample_Id/SV_Status/Site1_Hugo_Symbol/
# Site2_Hugo_Symbol are the ones fill_fusion_from_temp() in walk.py actually
# requires; the rest are carried through as-is). Fusions now go through
# Varan's FUSIONS/*.tsv ingestion path, NOT a synthesized CombinedVariantOutput
# comb_path file - see the note on MSI/comb_path below for why.
FUSION_TABLE_HEADER = ("Sample_Id\tSV_Status\tClass\tSite1_Hugo_Symbol\t"
                       "Site2_Hugo_Symbol\tNormal_Paired_End_Read_Count\t"
                       "Event_Info\tRNA_Support\n")

def run_cmd(cmd):
    result = subprocess.run(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.strip() if result.returncode == 0 else None

def list_s3_files(s3_folder):
    output = run_cmd(f"aws s3 ls {s3_folder.rstrip('/')}/")
    if not output: return []
    return [line.split()[-1] for line in output.splitlines()]

def load_oncotree_dict():
    name_to_code = {}
    valid_codes = set()
    if os.path.exists(DICT_PATH):
        with open(DICT_PATH, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'name' in row and 'code' in row:
                    name = row['name'].strip().lower()
                    code = row['code'].strip().upper()
                    name_to_code[name] = code
                    valid_codes.add(code)
    return name_to_code, valid_codes

def get_xml_data(xml_path, onco_info):
    name_to_code, valid_codes = onco_info
    p_id, o_code = "N/A", "UNKNOWN"
    try:
        root = ET.parse(xml_path).getroot()
        for param in root.iter():
            if any(c.tag.split('}')[-1] == "Name" and c.text == "accessionID" for c in param):
                for c in param:
                    if c.tag.split('}')[-1] == "Value":
                        p_id = c.text
                        break
        for diag in root.iter():
            if diag.tag.split('}')[-1] == "Diagnosis" and diag.text:
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
    except:
        pass
    return p_id, o_code

def get_msi_data(msi_path):
    """Read Guardant's .msi_call.hdr.tsv - has a REAL numeric score, not just
    a category.

    FIX (2026-07-19): the old version only kept a binary SUM_JSD=0/999999
    encoding for Varan's now-removed workaround, throwing away the real
    msi_score. Varan no longer needs that encoding - it accepts a real MSI
    VALUE and/or a pre-computed MSI_THR directly (see sample.tsv's MSI/
    MSI_THR columns). This returns both: the numeric msi_score as MSI, and
    msi_status normalized to Varan's own Stable/Unstable vocabulary as
    MSI_THR - both get written straight into sample.tsv, no synthesized
    CombinedVariantOutput file needed for MSI anymore.

    NOTE ON "MSS/MSI-L" matching: the real example file used to verify this
    (20250630_30ng_22.msi_call.hdr.tsv) only had an "MSI-H" example - the
    exact string for a stable call was never seen. This checks for "MSI-H"
    or "MSS" or "STABLE" as substrings (case-insensitive) rather than one
    hardcoded compound string, but still flag the first real MSS/MSI-L
    example you get to confirm the exact text.
    """
    data = {"run_id": "UNKNOWN_RUN", "msi_score": "", "msi_thr": ""}
    if msi_path and os.path.exists(msi_path):
        try:
            with open(msi_path, 'r') as f:
                r = csv.DictReader(f, delimiter='\t')
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
        except: pass
    return data

def load_fusions(fusion_path):
    fusions = []
    if fusion_path and os.path.exists(fusion_path):
        try:
            with open(fusion_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    if row.get("call") == "1":
                        fusions.append({
                            "gene_a": row["gene_a"],
                            "gene_b": row["gene_b"],
                            "supp": row.get("fusion_molecule_count_ab", "0"),
                        })
        except: pass
    return fusions

def append_fusions_to_table(fusion_table_path, sample_id, fusions):
    """Append this sample's real (call=1) fusions to the shared FUSIONS/
    Fusions.tsv file, in Varan's own data_sv.txt column shape.

    FIX (2026-07-19): replaces create_combined_file()'s synthesized
    CombinedVariantOutput.tsv with a [MSI]/[Fusions] section. That approach
    had a real side effect: as soon as ANY sample's comb_path is non-empty,
    Varan's walk_folder routes the WHOLE run's MSI/TMB through the
    CombinedOutput-driven path instead of honoring sample.tsv's own MSI/
    MSI_THR columns (see _walk_setup/_walk_write_clinical_tables in
    walk.py) - which would have silently defeated the new MSI VALUE+THR
    precedence logic for every sample in the batch, not just broken
    Fusions. Writing fusions into a plain FUSIONS/*.tsv file instead (passed
    as the third path to `varan.py -i sample.tsv patient.tsv fusions.tsv`)
    uses a separate ingestion path (fill_fusion_from_temp) that doesn't
    touch CombinedOutput at all.
    """
    is_new = not os.path.exists(fusion_table_path)
    with open(fusion_table_path, 'a') as f:
        if is_new:
            f.write(FUSION_TABLE_HEADER)
        for fus in fusions:
            f.write(
                f"{sample_id}\tSOMATIC\tFUSION\t{fus['gene_a']}\t"
                f"{fus['gene_b']}\t{fus['supp']}\t"
                f"{fus['gene_a']}-{fus['gene_b']} Fusion\tYes\n")

def load_cnv_tsv_ordered(tsv_path):
    """Read Guardant's .cnv_call.hdr.tsv, one row per gene.

    FIX (2026-07-19): the old `if cn_value == 2.0: continue` filter never
    actually excluded anything - copy_number is a continuous value (e.g.
    2.07, 1.84, 3.17 in the real example file), essentially never exactly
    2.0, so every gene's row survived this filter regardless of its `call`
    value. Removed rather than "fixed", since process_vcf() below already
    correctly gates on the real `call` column (0 = no significant call,
    1/2 = deletion/amplification) when deciding PASS vs FAIL for the CNV
    VCF it writes - that check was already doing the real filtering, this
    dead line was just misleading.
    """
    cnv_list = []
    if not tsv_path or not os.path.exists(tsv_path): return cnv_list
    with open(tsv_path, 'r') as f:
        lines = f.readlines()
        start = next((i for i, l in enumerate(lines) if "gene" in l.lower() and "copy_number" in l.lower()), -1)
        if start == -1: return cnv_list
        f.seek(0); [next(f) for _ in range(start)]
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                cnv_list.append({
                    'gene': row['gene'].strip(),
                    'cn': float(row['copy_number']),
                    'call': row['call'].strip()})
            except: continue
    return cnv_list

def process_vcf(vcf_in, cnv_tsv_in, snv_out, cnv_out, sample_id):
    ordered_cnv_data = load_cnv_tsv_ordered(cnv_tsv_in)
    cnv_index = 0
    col_header = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n"

    with open(vcf_in, 'r') as f_in, open(snv_out, 'w') as f_snv, open(cnv_out, 'w') as f_cnv:
        f_snv.writelines(COMMON_HEADER + SNV_SPECIFIC + [col_header])
        f_cnv.writelines(COMMON_HEADER + CNV_SPECIFIC + [col_header])

        for line in f_in:
            line = line.strip()
            if not line or line.startswith("#"): continue
            cols = line.split("\t")
            if len(cols) < 5: continue

            chrom = cols[0]
            if not chrom.startswith("chr"):
                chrom = "chrM" if chrom == "MT" else f"chr{chrom}"
            cols[0] = chrom

            while len(cols) < 10: cols.append(".")
            alt, info_str = cols[4], cols[7]
            # FIX (2026-07-19): was `"SVTYPE=" in info_str`, which also
            # matches "SVTYPE=BND" (fusion breakend candidate rows) - those
            # got misclassified as CNV-structural and consumed a slot in
            # ordered_cnv_data via cnv_index, which happened to be harmless
            # in the one real file checked (all BND rows came after every
            # real CNV row) but is not guaranteed in general: if a future
            # file interleaves BND rows among CNV rows, gene<->copy-number
            # pairing (done purely by sequential index, see below) would
            # silently misattribute copy numbers to the wrong genes.
            is_struct = any(x in alt for x in ["<CNV>", "<DUP>", "<DEL>"]) or any(
                svtype in info_str for svtype in ("SVTYPE=CNV", "SVTYPE=DUP", "SVTYPE=DEL"))
            # BND (breakend) rows are fusion junction evidence, already
            # captured properly via load_fusions()/.fusion_call.hdr.tsv -
            # neither a CNV call nor a point mutation. Checked against the
            # real example VCF: 4 BND rows exist, all FILTER=PASS, and
            # without this explicit skip they fell through into the SNV
            # branch below (is_struct=False for them) and got written into
            # the SNV VCF as if they were point mutations - REF=N, ALT in
            # breakend notation (e.g. "]6:117647236]N"), which vcf2maf
            # would either choke on or misinterpret.
            if "SVTYPE=BND" in info_str:
                continue

            if is_struct:
                if cnv_index < len(ordered_cnv_data):
                    d = ordered_cnv_data[cnv_index]
                    fc = round(d['cn'] / 2.0, 4)

                    # Forziamo QUAL (cols[5]) a "1"
                    cols[5] = "1"

                    # Logica specifica richiesta: ALT = <DUP> e FILTER = PASS solo per call 1 o 2.
                    # Altrimenti ALT = . e FILTER = FAIL.
                    call_val = str(d['call']).strip()
                    if call_val in ["1", "2"]:
                        cols[4] = "<DUP>"
                        cols[6] = "PASS"
                    else:
                        cols[4] = "."
                        cols[6] = "FAIL"

                    original_end = next((x.split("=")[1] for x in info_str.split(";") if x.startswith("END=")), str(int(cols[1]) + 1))
                    cols[7] = f"SVTYPE=CNV;END={original_end};SEGID={d['gene']}"
                    cols[8], cols[9] = "GT:CN:SM", f"0/1:{d['cn']}:{fc}"
                    cnv_index += 1
                else:
                    # Se i dati del TSV finiscono ma ci sono righe strutturali nel VCF
                    cols[4] = "."
                    cols[5] = "1"
                    cols[6] = "FAIL"
                f_cnv.write("\t".join(cols) + "\n")
            else:
                # FIX (2026-07-19): this branch used to write EVERY non-
                # structural row to the SNV VCF regardless of cols[6]
                # (FILTER). Checked against the real example VCF: 225 of
                # 343 non-structural rows are FILTER=FAIL with GT=1/1
                # (homozygous) and a dbSNP rsID - classic germline common-
                # SNP signatures, vs. 118 real FILTER=PASS somatic calls
                # (GT=0/1, low VAF, COSMIC IDs). Writing the FAIL rows
                # through meant germline variants nearly 2x outnumbering
                # real somatic calls were flowing into Varan's SNV
                # pipeline unfiltered - a real contamination risk for both
                # the clinical MAF and any downstream signature analysis.
                if cols[6] != "PASS":
                    continue
                try:
                    format_keys = cols[8].split(':')
                    sample_vals = cols[9].split(':')
                    if "AD" in format_keys:
                        idx_ad = format_keys.index("AD")
                        ad_parts = sample_vals[idx_ad].split(',')
                        ref_count, alt_count = int(ad_parts[0]), int(ad_parts[1])
                        depth = ref_count + alt_count
                        af = round(alt_count / depth, 4) if depth > 0 else 0
                        cols[8], cols[9] = "GT:AD:AF:DP", f"{sample_vals[0]}:{ref_count},{alt_count}:{af}:{depth}"
                except: pass
                f_snv.write("\t".join(cols) + "\n")

def get_incremental_report_path(base_dir, run_id):
    filename = f"{run_id}_VARAN.tsv"
    path = os.path.join(base_dir, filename)
    if not os.path.exists(path): return path
    counter = 1
    while True:
        filename = f"{run_id}_{counter}_VARAN.tsv"
        path = os.path.join(base_dir, filename)
        if not os.path.exists(path): return path
        counter += 1

def process_single_sample(sid, s3_folder, run_id_default, onco_dict, fusion_table_path, xml_folder=None):
    if xml_folder is None:
        xml_folder = s3_folder

    s3_files = list_s3_files(s3_folder)
    xml_files = list_s3_files(xml_folder) if xml_folder != s3_folder else s3_files

    clean_sid = sid.rstrip('_')

    def find_file_in_list(suffix, files_list):
        for f in files_list if files_list else []:
            if f.endswith(suffix) and clean_sid in f and not f.startswith("AIO"):
                return f
        return None

    # --- RICERCA DEI METADATI XML CON FALLBACK SPECULATIVO ---
    xml_name = find_file_in_list("_finalmetadata.xml", xml_files)

    if not xml_name and xml_folder == s3_folder:
        backup_folder = s3_folder.rstrip("/") + "_2"
        print(f"[{sid}] Metadata non trovato in {s3_folder}. Tentativo nel percorso di backup: {backup_folder}")
        backup_files = list_s3_files(backup_folder)
        xml_name = find_file_in_list("_finalmetadata.xml", backup_files)
        if xml_name:
            xml_folder = backup_folder
            print(f"[{sid}] Metadata TROVATO nella cartella di backup: {xml_name}")
        else:
            print(f"[{sid}] Metadata non trovato in nessun percorso.")
            return None

    # Configurazione target per il download dei differenti file
    # NOTE: "comb" (comb_path/CombinedVariantOutput) target removed - see
    # append_fusions_to_table()'s docstring for why fusions no longer go
    # through a synthesized CombinedOutput file.
    file_targets = {
        "xml": (xml_name, xml_folder),
        "vcf": (find_file_in_list(".vcf", s3_files), s3_folder),
        "msi": (find_file_in_list(".msi_call.hdr.tsv", s3_files), s3_folder),
        "cnv": (find_file_in_list(".cnv_call.hdr.tsv", s3_files), s3_folder),
        "fus": (find_file_in_list(".fusion_call.hdr.tsv", s3_files), s3_folder)
    }

    local_files = {}
    for key, target in file_targets.items():
        if target:
            fname, folder_path = target
            if fname:
                local_p = os.path.join(TEMP_LOCAL_DIR, f"{sid}_{fname}")
                if run_cmd(f"aws s3 cp {folder_path.rstrip('/')}/{fname} {local_p}"):
                    local_files[key] = local_p

    if "vcf" not in local_files:
        print(f"VCF mancante per {sid}"); return None

    p_id, o_code = get_xml_data(local_files["xml"], onco_dict)
    m_info = get_msi_data(local_files.get("msi"))
    run_id_final = m_info["run_id"] if m_info["run_id"] != "UNKNOWN_RUN" else run_id_default

    v_out = os.path.join(VCF_BASE_DIR, run_id_final, sid)
    os.makedirs(v_out, exist_ok=True)

    snv_f, cnv_f = os.path.join(v_out, f"{sid}.snv.vcf"), os.path.join(v_out, f"{sid}.cnv.vcf")

    # --- CHIAMATA VCF ATTIVATA ---
    process_vcf(local_files["vcf"], local_files.get("cnv"), snv_f, cnv_f, sid)

    fusions_list = load_fusions(local_files.get("fus"))
    append_fusions_to_table(fusion_table_path, sid, fusions_list)

    for p in local_files.values():
        if os.path.exists(p): os.remove(p)

    # No more comb_path (blank string) and no more SUM_JSD-encoded MSI:
    # real msi_score as MSI, msi_status mapped to Varan's Stable/Unstable
    # vocabulary as MSI_THR - Varan's own VALUE+THR precedence keeps this
    # THR as-is without recomputing from conf.ini (see walk.py
    # _resolve_biomarker_threshold). TMB/TMB_THR left blank - Guardant
    # doesn't report TMB (checked: no TMB field anywhere in
    # _finalmetadata.xml or _metadata.xml for the real example sample).
    return [sid, p_id, run_id_final, o_code, snv_f, cnv_f, "",
            m_info["msi_score"], "", m_info["msi_thr"], ""]

def main():
    parser = argparse.ArgumentParser(description="FPG360 Varan Input Generator")
    parser.add_argument("-f", "--folder", help="Path S3 della run")
    parser.add_argument("-s", "--selection", help="File TSV con lista campioni")

    if len(sys.argv) == 1:
        parser.print_help(); sys.exit(1)

    args = parser.parse_args()
    onco_dict = load_oncotree_dict()
    os.makedirs(TEMP_LOCAL_DIR, exist_ok=True)
    os.makedirs(REPORT_BASE_DIR, exist_ok=True)

    report_data = []
    main_run_id = "MANUAL_BATCH"

    if args.folder:
        s3_folder = args.folder.rstrip("/")
        print(f"Analizzando folder: {s3_folder}")
        main_run_id = os.path.basename(s3_folder).split(".")[0]
        s3_files = list_s3_files(s3_folder)
        xml_files = [f for f in s3_files if f.endswith("_finalmetadata.xml") and not f.startswith("AIO")]

        # --- SCOPERTA CAMPIONI CON FALLBACK NELLA CARTELLA DI BACKUP _2 ---
        xml_folder = s3_folder
        if not xml_files:
            backup_folder = s3_folder + "_2"
            print(f"Nessun metadata trovato in {s3_folder}. Ricerca nella cartella di backup: {backup_folder}...")
            backup_files = list_s3_files(backup_folder)
            xml_files = [f for f in backup_files if f.endswith("_finalmetadata.xml") and not f.startswith("AIO")]
            if xml_files:
                xml_folder = backup_folder
                print(f"Trovati {len(xml_files)} file di metadata nella cartella di backup.")

        fusion_table_path = os.path.join(REPORT_BASE_DIR, f"{main_run_id}_fusions.tsv")
        for xml_f in xml_files:
            sid = xml_f.replace("_finalmetadata.xml", "").split("_")[-1]
            print(f"Elaborazione {sid}...")
            res = process_single_sample(sid, s3_folder, main_run_id, onco_dict, fusion_table_path, xml_folder)
            if res: report_data.append(res)

    elif args.selection:
        print(f"Analizzando selezione: {args.selection}")
        main_run_id = Path(args.selection).stem
        fusion_table_path = os.path.join(REPORT_BASE_DIR, f"{main_run_id}_fusions.tsv")
        with open(args.selection, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sid, s3_folder = row['sample_id'], row['s3_path_run'].rstrip("/")
                run_id_tmp = os.path.basename(s3_folder).split(".")[0]
                res = process_single_sample(sid, s3_folder, run_id_tmp, onco_dict, fusion_table_path)
                if res: report_data.append(res)

    if report_data:
        r_path = get_incremental_report_path(REPORT_BASE_DIR, main_run_id)
        with open(r_path, 'w') as f:
            # Matches Varan's real Templates/sample.tsv column order exactly
            # (RUN_ID dropped - it was never a Varan column, only used
            # internally above for path bookkeeping; MSI_SCORE dropped as a
            # separate column since it now lives in MSI itself).
            f.write("\t".join(["SAMPLE_ID", "PATIENT_ID", "ONCOTREE_CODE", "snv_path", "cnv_path", "comb_path", "MSI", "TMB", "MSI_THR", "TMB_THR"]) + "\n")
            for row in report_data:
                sid, p_id, run_id_final, o_code, snv_f, cnv_f, comb_f, msi, tmb, msi_thr, tmb_thr = row
                f.write("\t".join(map(str, [sid, p_id, o_code, snv_f, cnv_f, comb_f, msi, tmb, msi_thr, tmb_thr])) + "\n")
        print(f"\nReport generato: {r_path}")
        print(f"Fusioni (se presenti) in: {fusion_table_path}")
        print("Lancia Varan con: python varan.py -i "
              f"{r_path} <patient.tsv o \"\"> {fusion_table_path} -o ... -c ...")
    else:
        print("\nNessun dato raccolto.")

if __name__ == "__main__":
    main()
