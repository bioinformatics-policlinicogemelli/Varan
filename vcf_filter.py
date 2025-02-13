import argparse
import pandas as pd

def parsing_vcf(file_input,file_output, filters):

    with open (file_input) as f:
        correct_data=[line.strip().split("\t") for line in f if not line.startswith("##")]

    data=pd.DataFrame(correct_data[1:],columns=correct_data[0])
    
    if "d" in filters:
        data=data[data["ALT"]!="."]
        data=data[data["FILTER"]=="PASS"]

    if "f" in filters:
        data['FIR'] = data['INFO'].str.extract("FractionInformativeReads=([0-9]{1}\.[0-9]{3})")
        data['DP'] = data['INFO'].str.extract("^DP=([0-9]+);")
        data["FIR"] = data["FIR"].apply(lambda x: float(x))
        data["DP"] = data["DP"].apply(lambda x: float(x))
        data["AASR"] = data['DP'] * data['FIR']

        data = data[data["AASR"] > 10]
        data.drop(["FIR", "DP", "AASR"], inplace=True, axis=1)

    if "z" in filters:
        format_columns = data['FORMAT'].iloc[0].split(":")
        sample_col = data.columns[data.columns.get_loc("FORMAT") + 1]

        # gt_index = format_columns.index('GT')
        # data['GT_extracted'] = data[sample_col].str.split(":").str[gt_index]
        # data = data[data['GT_extracted'] != '0/0']

        ad_index = format_columns.index("AD")
        ad_values = data[sample_col].str.split(":").str[ad_index]
        alt_allele_reads = ad_values.str.split(",").str[1]

        data["ALT_Reads"] = alt_allele_reads.apply(lambda x: float(x))
        data = data[data["ALT_Reads"] > 10]
        data.drop("ALT_Reads", inplace=True, axis=1)
    
    data.to_csv(file_output, sep="\t", mode="a", index=False)


def write_header_lines(input_vcf, output_vcf):
    with open(input_vcf, 'r') as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            if line.startswith("##"):
                f_out.write(line)



def main(input, output, filters):
    write_header_lines(input, output)
    parsing_vcf(input, output, filters)















# if __name__ == '__main__':


# # parse arguments
#     parser = argparse.ArgumentParser(description="Filter for a vcf",
#                                             epilog="Version: 1.0\n\
#                                             Author: Bioinformatics Facility GSTeP'\n\
#                                             email: luciano.giaco@policlinicogemelli.it")

#     # arguments
#     parser.add_argument('-i', '--input', help="<input.vcf>\
#                                             VCF file to filter",
#                                             required=True)
#     parser.add_argument('-o', '--output', help="<output-file.tab>\
#                                             file path of the Table output",
#                                             required=True)

#     args = parser.parse_args()
#     INPUT = args.input
#     OUTPUT = args.output

#     main(INPUT, OUTPUT)

