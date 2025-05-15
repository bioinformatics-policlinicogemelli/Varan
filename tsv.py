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



def get_msi_tmb(INPUT):
    data = {}
    with open(INPUT) as tsvFile:
        righe = tsvFile.read().splitlines()
        MSI_dic = {}
        MSI = []
        for riga in righe:
            if("Total TMB" in riga):
                campi = riga.split(sep="\t")
                TMB_Total = campi[1]
                data["TMB_Total"] = TMB_Total
            if("Usable MSI Sites" in riga):
                campi = riga.split(sep="\t")
                Usable_MSI = campi[1]
                MSI_dic["Usable_MSI"] = Usable_MSI
            if("Percent Unstable MSI Sites" in riga):
                campi = riga.split(sep="\t")
                Tot_MSI_unstable = campi[1]
                MSI_dic["Tot_MSI_unstable"] = Tot_MSI_unstable
            if("SUM_JSD" in riga):
                MSI_dic["Usable_MSI"] = 50
                campi = riga.split(sep="\t")
                MSI_dic["Tot_MSI_unstable"] = campi[1]
                Usable_MSI = campi[1]
        data["MSI"] = [(key, value) for key, value in MSI_dic.items()]
        return(data)


def split_hugo_symbols(hugo_symbol):
    for simbolo in [";", "-", "/"]:
        if simbolo in hugo_symbol:
            split_symbols = hugo_symbol.split(simbolo)
            return split_symbols
    return split_symbols


def get_fusions(INPUT):
    with open(INPUT) as file:
        fusioni = []
        lines = file.readlines()
        for i in range(len(lines)):
            if ("[Fusions]" or "[Data Fusions]") in lines[i]:
                for j in range(i+2, len(lines)):
                    if lines[j].strip() == "NA" or lines[j].strip() == "":
                        break

                    gene_pair, bp1, bp2, fsr, g1rr, g2rr = lines[j].strip().split("\t")

                    Hugo_Symbol = split_hugo_symbols(gene_pair)

                    chrom1 = bp1.split(":")
                    chrom2 = bp2.split(":")
                    Site1_Chromosome = chrom1[0]
                    Site1_Position = chrom1[1]
                    Site2_Chromosome = chrom2[0]
                    Site2_Position = chrom2[1]

                    fusioni.append({"Site1_Hugo_Symbol": Hugo_Symbol[0], "Site2_Hugo_Symbol": Hugo_Symbol[1], "Site1_Chromosome": Site1_Chromosome,
                                    "Site2_Chromosome": Site2_Chromosome, "Site1_Position": Site1_Position, "Site2_Position": Site2_Position,
                                    "Normal_Paired_End_Read_Count": fsr, "Gene 1 Reference Reads": g1rr,
                                    "Gene 2 Reference Reads": g2rr, "Event_Info": gene_pair})
        return fusioni


def main(INPUT):
    get_fusions(INPUT)
