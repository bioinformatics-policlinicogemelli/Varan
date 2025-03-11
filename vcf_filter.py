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