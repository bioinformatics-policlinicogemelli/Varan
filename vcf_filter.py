import argparse
import pandas as pd

def parsing_vcf(file_input,file_output):
    with open (file_input) as f:
        correct_data=[line.strip().split("\t") for line in f if not line.startswith("##")]

    data=pd.DataFrame(correct_data[1:],columns=correct_data[0])

    data=data[data["ALT"]!="."]

    data=data[data["FILTER"]=="PASS"]

    data.to_csv(file_output,sep="\t", mode="a", index=False)


def write_header_lines(input_vcf, output_vcf):
    with open(input_vcf, 'r') as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            if line.startswith("##"):
                f_out.write(line)


def main(input, output):
    write_header_lines(input, output)
    parsing_vcf(input, output)