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


def get_msi_tmb(input_file: str) -> dict:
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
    data = {}
    with input_file.open() as tsv_file:
        righe = tsv_file.read().splitlines()
        msi_dic = {}
        for riga in righe:
            if("Total TMB" in riga):
                campi = riga.split(sep="\t")
                tmb_total = campi[1]
                data["TMB_Total"] = tmb_total
            if("Usable MSI Sites" in riga):
                campi = riga.split(sep="\t")
                usable_msi = campi[1]
                msi_dic["Usable_MSI"] = usable_msi
            if("Percent Unstable MSI Sites" in riga):
                campi = riga.split(sep="\t")
                tot_msi_unstable = campi[1]
                msi_dic["Tot_MSI_unstable"] = tot_msi_unstable
            if("SUM_JSD" in riga):
                msi_dic["Usable_MSI"] = 50
                campi = riga.split(sep="\t")
                msi_dic["Tot_MSI_unstable"] = campi[1]
                usable_msi = campi[1]
        data["MSI"] = list(msi_dic.items())
        return(data)


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
            split_symbols = hugo_symbol.split(symbol)
    return split_symbols


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


def main(input_file: str) -> None:
    """Run the fusion data extraction pipeline.

    Parameters
    ----------
    input_file : Path
        Path to the input file containing fusion event data.

    """
    get_fusions(input_file)
