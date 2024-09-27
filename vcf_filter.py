#funzione per filtrare le colonne ALT e FILTER di un file vcf e creare un nuovo
#vcf filtrato

import argparse
import pandas as pd

def parsing_vcf(file_input,file_output):
    
    #lista contenente i dati da #CHROM in poi (no meta-information)
    with open (file_input) as f:
        correct_data=[line.strip().split("\t") for line in f if not line.startswith("##")]

    #creazione dataframe contenente i dati
    data=pd.DataFrame(correct_data[1:],columns=correct_data[0])
    
    #filtro per dati in cui ALT!="."
    data=data[data["ALT"]!="."]

    #filtro per dati in cui FILTER == "PASS"
    data=data[data["FILTER"]=="PASS"]
    #per salvare il file così filtrato
    data.to_csv(file_output,sep="\t", mode="a", index=False)


def write_header_lines(input_vcf, output_vcf):
    with open(input_vcf, 'r') as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            if line.startswith("##"):
                f_out.write(line)



def main(INPUT, OUTPUT):
    write_header_lines(INPUT,OUTPUT)
    parsing_vcf(INPUT,OUTPUT)

if __name__ == '__main__':


# parse arguments
    parser = argparse.ArgumentParser(description="Filter for a vcf",
                                            epilog="Version: 1.0\n\
                                            Author: Bioinformatics Facility GSTeP'\n\
                                            email: luciano.giaco@policlinicogemelli.it")

    # arguments
    parser.add_argument('-i', '--input', help="<input.vcf>\
                                            VCF file to filter",
                                            required=True)
    parser.add_argument('-o', '--output', help="<output-file.tab>\
                                            file path of the Table output",
                                            required=True)

    args = parser.parse_args()
    INPUT = args.input
    OUTPUT = args.output

    main(INPUT, OUTPUT)

