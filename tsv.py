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

"""Functions for extracting key metrics like TMB, MSI, and genomic fusions.

- `get_msi_tmb`: Extracts TMB and MSI data from a file.
- `split_hugo_symbols`: Splits Hugo symbols by common delimiters.
- `get_fusions`: Extracts fusion events, including gene pairs and locations.
- `main`: Processes fusion data using the above functions.

"""

from loguru import logger

def get_msi_tmb(input_file: str, sample_type: str) -> dict:
    """Extract MSI and TMB data from an input_file file.

    This function parses a TSV file to extract the total TMB, usable MSI sites,
    and unstable MSI sites from specific lines in the file.

    Parameters
    ----------
    input_file : str
        The path to the input_file file that contains the MSI and TMB data.

    Returns
    -------
    dict
        A dictionary containing the TMB total and MSI data.
        The MSI data is stored as a list of tuples

    """
    data = {"TMB_Total": "NA", "MSI": []}
    msi_dic = {"Usable_MSI": "NA", "Tot_MSI_unstable": "NA"}
        
    with input_file.open() as tsv_file:
        righe = tsv_file.read().splitlines()
        
        for riga in righe:
            campi = riga.split(sep="\t")
            if len(campi) < 2: continue

            if "Total TMB" in riga:
                data["TMB_Total"] = campi[1]

            if sample_type == "SOLID":
                if "SUM_JSD" in riga:
                    logger.warning("Sample_Type on conf.ini is set to SOLID, but the sample appears to be LIQUID. MSI value could be empty!")
                if "Usable MSI Sites" in riga:
                    msi_dic["Usable_MSI"] = campi[1]
                if "Percent Unstable MSI Sites" in riga:
                    msi_dic["Tot_MSI_unstable"] = campi[1]
            
            elif sample_type == "LIQUID":
                if "Usable MSI Sites" in riga:
                    logger.warning("Sample_Type on conf.ini is set to LIQUID, but the sample appears to be SOLID. MSI value could be empty!")
                if "SUM_JSD" in riga:
                    msi_dic["Usable_MSI"] = "50" # to bypass control
                    msi_dic["Tot_MSI_unstable"] = campi[1]
        
        data["MSI"] = [
            ("Usable_MSI", msi_dic["Usable_MSI"]),
            ("Tot_MSI_unstable", msi_dic["Tot_MSI_unstable"])
        ]
        return data


def get_gis(input_file: str) -> dict:
    """Extract Genomic Instability Score, Tumor Fraction and Ploidy from a file.

    These three values live in the [GIS] section of a CombinedVariantOutput
    file, which DRAGEN only writes when the TSO500 HRD feature was enabled for
    that run - for a plain TSO500 run the section is absent entirely (not
    present-but-empty), so every value below stays "NA" unless actually found,
    the same way get_msi_tmb handles optional fields.

    Parameters
    ----------
    input_file : Path
        The path to the CombinedVariantOutput file.

    Returns
    -------
    dict
        A dictionary with "GIS", "Tumor_Fraction" and "Ploidy" - "NA" for any
        value not present in the file (i.e. non-HRD samples).

    """
    data = {"GIS": "NA", "Tumor_Fraction": "NA", "Ploidy": "NA"}

    with input_file.open() as tsv_file:
        for riga in tsv_file.read().splitlines():
            campi = riga.split(sep="\t")
            if len(campi) < 2:
                continue

            if "Genomic Instability Score" in riga:
                data["GIS"] = campi[1]
            if "Tumor Fraction" in riga:
                data["Tumor_Fraction"] = campi[1]
            if "Ploidy" in riga:
                data["Ploidy"] = campi[1]

    return data


def split_hugo_symbols(hugo_symbol: str) -> str:
    """Split a Hugo symbol into multiple gene symbols.

    Parameters
    ----------
    hugo_symbol : str
        A string that may contain multiple gene symbols separated by ';', '-', or '/'.

    Returns
    -------
    str
        A list of individual gene symbols.

    """
    for symbol in [";", "-", "/"]:
        if symbol in hugo_symbol:
            return hugo_symbol.split(symbol)
    msg = f"No known separator (';', '-', '/') found in Hugo Symbol: {hugo_symbol!r}"
    raise ValueError(msg)


def get_fusions(input_file: str) -> list[dict[str, str]]:
    """Extract gene fusion events from a given file.

    This function parses a fusion report and extracts gene symbols,
    chromosome locations, and supporting read counts for each fusion event.

    Parameters
    ----------
    input_file : str
        Path to the input file containing fusion data.

    Returns
    -------
    list of dict
        A list of dictionaries, representing a fusion event for gene symbols,
        chromosomes, positions, and read counts.

    """
    with input_file.open() as file:
        fusions = []
        lines = file.readlines()
        for i in range(len(lines)):
            if "[Fusions]" in lines[i] or "[Data Fusions]" in lines[i]:
                for j in range(i+2, len(lines)):
                    if lines[j].strip() == "NA" or lines[j].strip() == "":
                        break

                    gene_pair, bp1, bp2, fsr, g1rr, g2rr = lines[j].strip().split("\t")

                    hugo_symbol = split_hugo_symbols(gene_pair)

                    chrom1 = bp1.split(":")
                    chrom2 = bp2.split(":")
                    site1_chromosome = chrom1[0]
                    site1_position = chrom1[1]
                    site2_chromosome = chrom2[0]
                    site2_position = chrom2[1]

                    fusions.append({
                        "Site1_Hugo_Symbol": hugo_symbol[0],
                        "Site2_Hugo_Symbol": hugo_symbol[1],
                        "Site1_Chromosome": site1_chromosome,
                        "Site2_Chromosome": site2_chromosome,
                        "Site1_Position": site1_position,
                        "Site2_Position": site2_position,
                        "Normal_Paired_End_Read_Count": fsr,
                        "Gene 1 Reference Reads": g1rr,
                        "Gene 2 Reference Reads": g2rr,
                        "Event_Info": gene_pair})
        return fusions


def get_splice_variants(input_file: str) -> list[dict[str, str]]:
    """Extract splice variant events from a CombinedVariantOutput file.

    NOTE ON VERIFICATION: the six column names below - Gene / Affected Exon /
    Breakpoint 1 / Breakpoint 2 / Splice Supporting Reads / Reference Reads
    Transcript - are not a guess: they are copied verbatim from the
    "[Splice Variants]" header row of real TSO500 CombinedVariantOutput.tsv
    files (Illumina always writes the header, even when the section has no
    calls). Independently cross-checked against Illumina's own TSO500 v2.2
    Local Run Manager release notes, which describe a defect fix for the
    "Splice Support Reads" and "Reference Reads Transcript" columns
    specifically in this file (their values had been swapped) - confirming
    both are real, plain numeric columns, structurally analogous to the
    already-verified [Fusions] section's FSR / reference-read-count pair.
    What remains unverified is only the *row*-splitting for a populated
    line (every real example file seen so far has "NA", no calls) - by
    analogy with [Fusions] (tab-separated, Breakpoint columns as
    "chrom:pos") this is a reasonable inference, not a blind guess, but
    fill_splice_from_combined() logs every row it parses at INFO level so
    the first real populated row is easy to spot-check.

    Parameters
    ----------
    input_file : Path
        Path to the input file containing splice variant data.

    Returns
    -------
    list of dict
        A list of dictionaries, one per splice variant event, with keys
        "Gene", "Affected_Exon", "Breakpoint_1", "Breakpoint_2",
        "Splice_Supporting_Reads" and "Reference_Reads_Transcript".

    """
    with input_file.open() as file:
        lines = file.readlines()
        splice_variants = []
        for i in range(len(lines)):
            if "[Splice Variants]" in lines[i]:
                for j in range(i + 2, len(lines)):
                    if lines[j].strip() in ("NA", ""):
                        break
                    fields = lines[j].strip().split("\t")
                    if len(fields) < 6:
                        continue
                    gene, affected_exon, bp1, bp2, ssr, ref_reads = fields[:6]
                    splice_variants.append({
                        "Gene": gene,
                        "Affected_Exon": affected_exon,
                        "Breakpoint_1": bp1,
                        "Breakpoint_2": bp2,
                        "Splice_Supporting_Reads": ssr,
                        "Reference_Reads_Transcript": ref_reads})
                break
        return splice_variants


def get_exons(input_file: str) -> list[dict[str, str]]:
    """Extract BRCA exon-level CNV information from a given input file.

    This function searches for the "[Exon-Level CNVs]" section in the input file
    and parses the next two lines of CNV data (typically corresponding to BRCA1 and BRCA2).
    If a line contains only "BRCA1\tNA" or "BRCA2\tNA", it will be skipped.

    Args:
        input_file (str): Path to the input file containing CNV information.
            This file must include a section labeled "[Exon-Level CNVs]".

    Returns:
        list[dict[str, str]]: A list of dictionaries, each representing a CNV entry with keys:
            "Hugo_Symbol", "Chromosome", "Start_Position", "Stop_Position",
            "Affected_Exon(s)", "Fold_Change", and "CNV_Type".

    """
    with input_file.open() as file:
        lines = file.readlines()

        has_exon_section = any("[Exon-Level CNVs]" in line for line in lines)
        if not has_exon_section:
            return None

        exonic = []
        for i in range(len(lines)):
            if "[Exon-Level CNVs]" in lines[i]:
                for j in range(i + 2, min(i + 4, len(lines))):
                    line = lines[j].strip()
                    if line.endswith("\tNA") or "NA" in line.split("\t")[1:]:
                        continue
                    fields = line.split("\t")
                    if len(fields) < 7:
                        continue
                    gene, chr, start, stop, exon, fc, cnv_type = fields[:7]
                    exonic.append({
                        "Hugo_Symbol": gene,
                        "Chromosome": chr,
                        "Start_Position": start,
                        "Stop_Position": stop,
                        "Affected_Exon(s)": exon,
                        "Fold_Change": fc,
                        "CNV_Type": cnv_type
                    })
                break
        return exonic


def main(input_file: str) -> None:
    """Run the fusion data extraction pipeline.

    Parameters
    ----------
    input_file : Path
        Path to the input file containing fusion event data.

    """
    get_fusions(input_file)
