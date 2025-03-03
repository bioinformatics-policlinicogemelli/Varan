# VARAN

[![DOI](https://zenodo.org/badge/788270006.svg)](https://zenodo.org/doi/10.5281/zenodo.12806060)

<p align="center">
<img src="readme_content/img/logo_VARAN.png" alt="MarineGEO circle logo" style="height: 400px; width:400px;"/>
</p>

## Index
- [VARAN](#varan)
  - [Overview](#overview)
  - [Quickstart](#quickstart-using-docker-container)
  - [Installation Procedure](#installation-procedure)
    - [Docker (recommended)](#docker-recommended)
    - [Local](#local)
  - [Quickstart](#quickstart)
    - [Block One: STUDY CREATION](#block-one-study-creation)
      - [1. Preparing Input](#1-preparing-input)
      - [2. Launch Varan](#2-launch-varan)
      - [3. Output](#3-output)
    - [Block Two: STUDY MANIPULATION](#block-two-study-manipulation)
      - [1. Preparing Input](#1-preparing-input-1)
      - [2. Launch Varan main](#2-launch-varan-main)
      - [3. Output](#3-output-1)

## Overview

<p align="justify">
Varan is a Python-based application that provides a pipeline to automatically prepare and manipulate cancer genomics data in the specific format supported by the <a href="https://www.cbioportal.org/">cBioPortal</a>.

### Features

* <ins>Study Creation</ins>
<br>This section provides, starting from raw vcf files, a well-structured and validated study folder ready to be uploaded into the local instance of cBioPortal.

* <ins>Study Manipulation</ins>
<br> This section allows users to work on already existing studies. In particular it is possible to merge two studies or to modify a study by extracting/removing samples.

## Quickstart using Docker container

Quick deploy Varan from raw data to launching analysis pipelines in under 30 minutes using a Docker container.


### Requirements
Software Requirements:

* <a href="https://www.docker.com/">Docker</a>

Hardware Requirements:
* ~3 GB storage space for the container image.
* ~25 GB storage space for the VEP cache and the reference genome.

### STEP 1: Clone Varan repository

a) Open a terminal

b) Clone the repository folder
```
git clone https://github.com/bioinformatics-policlinicogemelli/Varan-Pub.git
```
c) Build docker file
```
cd Varan-Pub
docker build -t varan .
```
d) Test installation
```
docker run --rm -it varan -h
```
⚠️ for Windows users: some issues with Git Bash have been reported. It is recommended to launch the docker command through <a href="https://learn.microsoft.com/en-us/powershell/scripting/overview?view=powershell-7.4">Powershell</a>).

### STEP 2: Download the example data



### STEP 3: 







## Installation Procedure

### Docker (recommended)

<details open>
  <summary><b>Prerequisites</b></summary>
  
* [Docker](https://www.docker.com/)
* <a href="https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache"> Vep Cache</a> 
</br> ⚠️ For test data download <a href="https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh37.tar.gz">homo_sapiens_vep_111_GRCh37.tar.gz</a> 

* A reference genome (which you can find <a href="http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache"> here</a>) and its indexed version (<a href="http://www.htslib.org/doc/samtools-faidx.html">how to index the genome</a> with Samtools).
</br>⚠️ You must use the same genome reference version that was used to generate the VCF files. 
</details>

<details open>
  <summary><b>Procedure</b></summary>

1. Open a terminal
2. Clone the repository folder:
```
git clone https://github.com/bioinformatics-policlinicogemelli/Varan-Pub.git
```
3. Build docker file
<br> ⚠️ This step can take about 30 minutes 
```
cd <Varan_folder_path>/Varan-Pub
docker build -t varan .
```
4. Run Varan to test the installation
```
docker run --rm -it varan -h
```

⚠️ for Windows users: some issues with Git Bash (Git for Windows) have been reported. It is recommended to launch the docker command through [Powershell](https://learn.microsoft.com/en-us/powershell/scripting/overview?view=powershell-7.4)
</details>

### Local

<details open>
  <summary><b>Prerequisites</b></summary>

* **Variant Effect Predictor (VEP)**<p align="justify">The Variant Effect Predictor <a href="https://www.ensembl.org/info/docs/tools/vep/index.html">VEP</a> is a tool used to determine the effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions. <br>The steps to install VEP can be found <a href="https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html"> here</a>, while DB and FASTA files can be downloaded <a href="http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache"> here</a>

* **vcf2maf**<p align="justify"><a href="https://github.com/mskcc/vcf2maf/tree/main">vcf2maf</a> is the tool required to convert vcf files into maf format. 
<br>All the installation info can be found <a href="https://github.com/mskcc/vcf2maf/tree/main">here</a>

* **Samtools**<p align="justify"><a href="https://www.htslib.org/">Samtools</a> is a suite of programs for interacting with high-throughput sequencing data. 
<br>All the installation info can be found <a href="https://www.htslib.org/download/">here</a>
</details>

<details open>
  <summary><b>Procedure</b></summary>

To correctly install and use Varan:
1. Open a terminal
2. Type the following command to clone the repository folder: 
```
git clone https://github.com/bioinformatics-policlinicogemelli/Varan-Pub.git
```
3.  Install all the packages required (cbioimporter and oncokb-annotator)
```
cd <varan_folder_path>/Varan-Pub
bash installer.sh
```
 
<p align="justify">
⚠️ <i>Depending on the Python version, pip3 may be required instead of pip</i><br><br>
To test the installation and check if everything works, launch the main script <b>varan.py</b>: 

```
cd <varan_folder_path>/Varan-Pub
python varan.py -h
```

<p align="justify">
⚠️ <i>If any error is printed while launching varan.py, check if step 3 completed without errors</i>
</details>

## Quickstart

<details open>
  <summary><b>Configuration file</b></summary>
The first step to start using Varan is to correctly set the configuration file <i>conf.ini</i>. <br><br>

⚠️ In the next subparagraph, for each field the type of variable requested will be insert between angle brackets <>. For string two possible entry can be find: < 'string' > and < string >. In the first case it is request to insert text inside quotation marks (*i.e. DESCRIPTION='this is the description'*), while on the other one quotation marks are not requested (*i.e. PROJECT_NAME=study*).

This file is divided in 12 subsessions:
<details close>
<summary><ins>Paths</ins></summary>
In this section it's possible to specify the paths for Vep (VEP_PATH) and its cache (VEP_DATA) and fasta (REF_FASTA), vcf2maf (VCF2MAF), and ClinVar (CLINV). In CACHE it is possible to set Ensembl version (i.e. 111)

```
[Annotations]
ANNOTATIONS = < '' >
```
⚠️ clinvar database can be downloaded <a href="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/">here</a>

</details>
<details close>
<summary><ins>Multiple</ins></summary>

In this section it's possible to specify the paths in case of multiple SNV, CNV and/or Combined Variant Output files analysis.   
```
[Multiple]
SNV= < string >
CNV= < string >
COMBOUT= < string >
```
⚠️ This section has to be filled only in case of input by file.

</details>
<details close>
<summary><ins>Zip</ins></summary>

In this section you can decide how to manage MAF and SNV_FILTERED files.
```
[Zip]
ZIP_MAF = < boolean >
ZIP_SNV_FILTERED = < boolean >
COPY_MAF = < boolean >
```
</details>
<details close>
<summary><ins>OncoKB</ins></summary>

In this section is possible to insert the personal <a href="https://www.oncokb.org/">oncoKB</a> key. This key is mandatory to execute the oncoKB annotation.

```
[OncoKB] ONCOKB
ONCOKB= < string >
```
⚠️ The request for the oncoKB key can be done <a href="https://www.oncokb.org/account/register">here</a>


</details>
<details close>
<summary><ins>Project</ins></summary>

In this section is possible to specify project info like study name, ID, description and profile. These info will be insert in meta files.
```
[Project]
PROJECT_ID = < string >
PROJECT_NAME = < string >
DESCRIPTION = < 'string' >
PROFILE_MUT = < 'string' >
PROFILE_CNA = < 'string' >
PROFILE_CNA_HG19 = < 'string' >
PROFILE_SV = < 'string' >
```

</details>
<details close>
<summary><ins>Filters</ins></summary>

Here it is possible to specify the filters' threshold to apply to SNV maf files (for more info about filters threshold setting see [Ex 7)](#ex-7-filter-vcfmaf).
```
[Filters]
CLIN_SIG = < ['string','string',...] >
CONSEQUENCES = < ['string','string',...] >
ONCOKB_FILTER = < ['string','string',...] >
t_VAF_min = < float >
t_VAF_min_novel = < float >
t_VAF_max = < float >
AF = < string > 
POLYPHEN = < ['string','string',...] >
IMPACT = < ['string','string',...] >
SIFT = < ['string','string',...] >
drop_NA_AF = < boolean >
```
⚠️ For AF the field can be populate with </>/<=/>=val (i.e. AF = <0.003)

</details close>
<details>
<summary><ins>Cna</ins></summary>

In this section user can insert CNV genotypes of interest and ploidy. The latter will be used to evaluate copy number discretization using <a href="https://cnvkit.readthedocs.io/en/stable/pipeline.html">cnvkit formula</a>.
```
[Cna]
HEADER_CNV = < ['string','string',...] >
PLOIDY = < int >
CNVkit = < boolean >
```
</details close>
<details>
<summary><ins>TMB</ins></summary>

Here TMB thresholds can be specified.

```
[TMB]
THRESHOLD_TMB = < {'string':'string', 'string':'string', ...} > 
```
i.e. THRESHOLD = {'Low':'<=5','Medium':'<10','High':'>=10'} where the string before : is the label assign to TMB value, while the other is the specific threshold (i.e. for a sample with TMB=15 a label 'High' will be reported in the data clinical sample).

</details close>
<details>
<summary><ins>MSI</ins></summary>

Here MSI thresholds for sites and values can be specified.
```
[MSI]
THRESHOLD_SITES = < string >
THRESHOLD_MSI = < string >
```

⚠️ THRESHOLD_SITES value will be used only if MSI informations are extracted from Combined Variant Output files. If MSI is directly reported as a value inside the input tsv, only THRESHOLD_MSI will be applied.

⚠️ Both THRESHOLD_SITES and THRESHOLD_MSI can be populate with </>/<=/>=val (i.e. THRESHOLD = <20).

</details close>
<details>
<summary><ins>FUSION</ins></summary>

Here Fusions thresholds can be specified.

```
[FUSION]
THRESHOLD_FUSION = < string >
```
⚠️ THRESHOLD can be populate with </>/<=/>=val (i.e. THRESHOLD = >=15)
</details>
<details close>
<summary><ins>ClinicalSample</ins></summary>

Users can customize column names and data types for the data_clinical_sample.txt file.
```
[ClinicalSample]
HEADER_SAMPLE_SHORT = < ['string','string',...] >
HEADER_SAMPLE_LONG = < ['string','string',...] >
HEADER_SAMPLE_TYPE = < ['string','string',...] >
```
⚠️ HEADER_SAMPLE_TYPE accept only STRING, NUMBER, BOOLEAN. If different type is inserted, an error will be raised by Varan.

⚠️ If these fields are left empty, a default Header will be produced.

</details>
<details close>
<summary><ins>ClinicalPatient</ins></summary>

Users can customize column names and data types for the data_clinical_patient.txt file.
```
[ClinicalPatient]
HEADER_PATIENT_SHORT = < ['string','string',...] >
HEADER_PATIENT_LONG = < ['string','string',...] >
HEADER_PATIENT_TYPE = < ['string','string',...] >
```
⚠️ HEADER_PATIENT_TYPE accept only STRING, NUMBER, BOOLEAN. If different type is inserted, an error will be raised by Varan.

⚠️ If these fields are left empty, a default Header will be produced.

</details>
<details close>
<summary><ins>Annotations</ins></summary>
In this section it's possible to insert manual notes that will appear in report_VARAN.html file.

```
[Annotations]
ANNOTATIONS = < ['string','string',...] >
```

</details>
<details close>
<summary><ins>Validation</ins></summary>

If the user has a working cbioportal instance active on his computer, the location (http://localhost:8080) can be insert here. This value will be used for the validation of the output study and will produce a html report with the results.
```
[Validation]
PORT = < string >
```
⚠️ If this field is left empty an offline validation will be conducted.
</details>

⚠️ For Docker a partially compiled configuration file (<i>docker.ini</i>) is available with vep and clinvar path already set.
</details>

<br>
<details open>
  <summary><b>Workflow</b></summary>
<br>

<p align="center">
  <img width="500" height="250" src="readme_content/img/workflow.png">
</p>

<p align="justify"><br>Varan application can be divided in two distinct main blocks that require different inputs and provide different actions. 
<ul>
  <li>The first block contains the functions to create a new study folder ex-novo.</li>
  <li>The second one contains the functions to merge studies and to modify (Update/Extract/Remove samples) an existing study folder.</li> 
</ul>
To keep track of all operations performed, a versioning system and a report file are provided.
</details>

### Block One: STUDY CREATION

<br>

<p align="center">
  <img width="600" height="300" src="readme_content/img/block1.png">
</p>

#### 1. Preparing Input

<p align="justify">
To create a new study folder, Varan requires vcf files as input. Varan can process two types of input: A) <b>Folder</b>; B) <b>File(s)</b>

<details close>
  <summary><ins>Folder</ins></summary>

User must organize an input folder containing all of the vcf and tsv files requested following the structure reported below:

```
input_folder/
├── CNV
│   ├── 001.vcf
│   ├── 002.vcf
│   └── 003.vcf
├── SNV
│   ├── 001.vcf
│   ├── 002.vcf
│   └── 003.vcf
├── CombinedOutput
│   ├── 001_CombinedVariantOutput.tsv
│   ├── 002_CombinedVariantOutput.tsv
│   └── 003_CombinedVariantOutput.tsv
├── FUSIONS
│   └── Fusions.tsv
├── sample.tsv   
└── patient.tsv
```

Where:

* <b>SNV</b> folder contains all vcf files related to single nucleotide variants.
* <b>CNV</b> folder contains all vcf files related to copy number variants.
* <b>CombinedVariantOutput</b> folder contains all the Combined Variant Output file in tsv format (this kind of files contains info about TMB, MSI and Fusions).
* <b>FUSIONS</b> folder contains the template fusion.tsv, which can filled with fusion data if the Combined Variant Output files are unavailable.
* <b>sample.tsv</b> file is a template for clinical data of the samples.
* <b>patient.tsv</b> file is a template for clinical data of the patients.
</details>
<details close>
<summary><ins>File(s)</ins></summary>
User must compile several input file (by filling in specific templates):

* <b>sample.tsv</b> file, a template to fill with samples' clinical info.
* <b>fusions.tsv</b> file, a template to fill with fusions' information.
* <b>patient.tsv</b> file, a template to fill with patients' clinical info.

⚠️ For this input, only the <i>sample.tsv</i> file is mandatory.
</details>
<br>
<details open>
  <summary><b>Templates</b></summary>

For both input by folder and by file, template filling by user is requested.
Below will be briefly explained the structure of these templates:

<details close>
  <summary><i>sample.tsv</i></summary>
This template must be filled by user with all disposable samples' clinical info and will be used to create the <i>data_clinical_sample.txt</i>.

⚠️ This file is mandatory for Varan analysis!

|SAMPLE_ID | PATIENT_ID | MSI | TMB| MSI_THR | TMB_THR| ONCOTREE_CODE| snv_path| cnv_path| comb_path|... |
|:---:|:---:|:---: |:---:   |:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|0000000_DNA| 00000000  | 1| 12| ||BOWEL |path_to_snv|path_to_cnv|path_to_combined_output|...|
|0000001_DNA| 00000001  | 8.0| 8.0| ||UTERUS |path_to_snv|path_to_cnv|path_to_combined_output|...|
|0000003_DNA| 00000003  | 222| 127| ||BOWEL |path_to_snv|path_to_cnv|path_to_combined_output|...|
|...| ...  | ...| ...| ||... |...|...|...|...|

The obligatory fields to keep are:

* <ins>SAMPLE_ID</ins>: IDs of all samples of interest.
* <ins>PATIENT_ID</ins>: IDs of all patient of interest.
* <ins>MSI</ins>: MSI value for each of the samples.
* <ins>TMB</ins>: TMB value for each of the samples.
* <ins>MSI_THR</ins>: MSI categorization based on the threshold set in the <i>conf.ini</i> file. This field has to be left empty and will be filled by Varan.
* <ins>TMB_THR</ins>: TMB categorization based on the threshold set in the <i>conf.ini</i> file. This field has to be left empty and will be filled by Varan.
* <ins>ONCOTREE_CODE</ins>: code to associate for the oncokb annotation. Check <a href="https://oncotree.mskcc.org/#/home">here</a> for more info.
* <ins>snv_path</ins>: path to the SNV vcf file. This column has to be filled in case of input by file.
* <ins>cnv_path</ins>: path to the CNV vcf file. This column has to be filled in case of input by file.
* <ins>comb_path</ins>: path to the Combined Variant Output file. This column has to be filled in case of input by file.

⚠️ The user can add new columns starting from the last one. Modifying or deleting the default ones (even only by changing names) can lead to errors and is strongly not recommended.<br>
⚠️ Both SAMPLE_ID and PATIENT_ID must be filled.<br>
⚠️ The MSI and TMB columns will be read from this file only if no Combined Variant Output file is provided. If the Combined Variant Output files are available, those values will be prioritized.<br>
⚠️ ONCOTREECODE column is mandatory to fill for the oncoKB annotation

</details>
<details close>
  <summary><i>patient.tsv</i></summary>

This template must be filled by user with all disposable patients' clinical info and will be used to create the <i>data_clinical_patient.txt</i>.

⚠️ This file is optional. If missing, a default <i>data_clinical_patient.txt</i> will be create.


|PATIENT_ID | AGE | SEX | ...|
|:---:|:---:|:---:|:---:|
|00000000| 45  | F| ...| 
|00000001| 77  | F| ...| 
|00000003| 23  | M| ...| 
|...| ...  | ...|...| 

The obligatory fields to keep are:
* <ins>PATIENT_ID</ins>: IDs of all patient of interest

⚠️ The user can add new columns starting from PATIENT_ID. Modifying or deleting the default one (even only by changing name) can lead to errors and is strongly not recommended.<br>

</details>

<details close>
  <summary><i>fusions.tsv</i></summary>

This template must be filled by user with all disposable fusion info and will be used to create the <i>data_sv.txt</i>.

⚠️ This file is optional.

|Sample_Id | SV_Status | Site1_Hugo_Symbol |Site2_Hugo_Symbol|...|
|:---:|:---:|:---:|:---:|:---:|
|0000000_DNA| SOMATIC  | APC| BRCA1|...|
|0000001_DNA| SOMATIC  | TP53| BRAF|...|
|0000003_DNA| SOMATIC  | ALK| BRCA2|...|
|...| ... | ...|...|...|

The obligatory fields to keep are:
* <ins>Sample_Id</ins>: IDs of all sample of interest
* <ins>SV_Status</ins>: fusion type
* <ins>Site1_Hugo_Symbol</ins>: first gene involved in fusion 
* <ins>Site2_Hugo_Symbol</ins>: second gene involved in fusion 

⚠️ The user can add new columns starting from Site2_Hugo_Symbol. Modifying or deleting the default ones (even only by changing names) can lead to errors and is strongly not recommended.<br>
⚠️ All mandatory fields must be filled and cannot be left empty.<br>
</details>
</details>

#### 2. Launch Varan

The possible option to launch varan main for block 1 are:

| Options | Description | Type | Required
|:---|:---|:---:|:---:|
|<div style="width:130px">-i --input</div>| <div style="width:220px"><p align="justify">Add this option to insert the path of the input (folder or file(s) </div>| string list | Yes
|-o --output_folder| <p align="justify">Add this option to insert the path where to save the output folder| string | Yes
|-c --cancer| <p align="justify">Add this option to specify the cancer type| string | Yes
|-f --filter| <p align="justify">Add this option to filter out variants from vcf/maf | string | No
|-k --onocoKB| <p align="justify">Add this option to annotate with oncoKB | boolean | No
|-t --analysis_type|<p align="justify">Add this option to specify the type of file (snv/cnv/fus/tab) to analyze. If not specified, all analysis will be done |string|No
|-w --overWrite| <p align="justify">Add this option to overwrite an already exisitng output folder|boolean| No
|-R --resume| <p align="justify">Add this option to resume an already started analysis.<br><br>⚠️This option must be used with caution, because it assumes that the previous VCF to MAF conversion step was successful.| boolean | No
|-m --multiple| <p align="justify">Add this option to specify that the input is a multi-sample vcf file (a single VCF containing information from multiple patients) | boolean | No

<br>
<details open>
  <summary><i>Examples</i></summary>

<details open>
  <summary><i>Docker version</i></summary>

To launch Varan docker version is mandatory to mount sevaral volumes (-v) for granting a correct functioning.
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan <commands>
```



Ex 1) <ins>Launch Varan base analysis with input folder</ins>:

Launch this command to process the contents of the input folder 

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta varan <commands> -v <conf.ini>:/conf.ini varan -i <path_input_folder> -o /output/<output_name> -c <type_of_cancer>
```
Ex 2) <ins>Launch Varan base analysis with input file</ins>:

Launch one of these commands to process the contents of the input file(s):
<li>If you have sample.tsv and patient.tsv as input:</li>

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i sample.tsv patient.tsv -o /output/<output_name> -c mixed 
```
<li>If you have sample.tsv and fusion.tsv as input:</li>

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i sample.tsv "" fusion.tsv -o /output/<output_name> -c mixed
```
<li>If you have sample.tsv, patient.tsv and fusion.tsv as input:</li>

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i sample.tsv patient.tsv fusion.tsv -o /output/<output_name> -c mixed
```
Ex 3) <ins>Multiple vcf analysis</ins>: 

Launch this command to specify that your input is a multi-vcf file or folder:
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -m
```
Ex 4) <ins>Overwrite analysis</ins>:

Launch this command to overwrite the output folder:
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -w
```
Ex 5) <ins>Resume analysis</ins>:

Launch this command to resume an already started analysis:
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -R
```
Ex 6) <ins>Specify analysis</ins>:

Launch one of these commands to specify the analysis' type:

* snv -> only snv analysis.
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -t snv
```
* cnv -> only cnv analysis.
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -t cnv
```
* fus -> only fusion analysis.
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -t fus
```
* tab -> only meta file creation.
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -t tab
```

##### Ex 7) <ins>Filter vcf/maf</ins>:

Launch this command to specify the type of filter to use:

* <code style="color : cyan">d</code> -> filter out from snv mutations with ALT="." and FILTER ≠"PASS".
* <code style="color : cyan">p</code> -> filter out from MAF mutations with FILTER ≠"PASS".
* <code style="color : cyan">v</code>-> filter out from MAF mutations with vaf (t_VF column) values <ins>not in</ins> the ranges [t_VAF_min; t_VAF_max] specified in <i>conf.ini</i>.
* <code style="color : cyan">n</code> -> apply a specific VAF filter to novel mutations (dbSNP_RS = "novel") filtering out from MAF novel mutations with vaf <ins>not in</ins> the range specified in <i>conf.ini</i> t_VAF_min_novel.
* <code style="color : cyan">a</code> -> filter out from MAF mutations with AF values <ins>not in</ins> the range specified in <i>conf.ini</i>. You can choose to drop or keep mutations with no value for this column via 'drop_NA_AF' setting.
* <code style="color : cyan">o</code>-> filter out from MAF mutations with ONCOGENIC values <ins>different</ins> to the ones specified in <i>conf.ini</i> ONCOKB_FILTER field. This requires the OncoKB annotation (-k option).
* <code style="color : cyan">c</code> -> filter out from MAF mutations with CLIN_SIG values <ins>equals</ins> to the ones specified in <i>conf.ini</i> CLIN_SIG field.
* <code style="color : cyan">q</code> -> filter out from MAF mutations with Consequences values <ins>different</ins> to the ones specified in <i>conf.ini</i> CONSEQUENCES field.
* <code style="color : cyan">y</code>-> filter out from MAF mutations with PolyPhen annotation <ins>different</ins> to the ones specified in <i>conf.ini</i> POLYPHEN field.
* <code style="color : cyan">i</code>-> filter out from MAF mutations with IMPACT annotation <ins>equals</ins> to to the ones specified in <i>conf.ini</i> IMPACT field.
* <code style="color : cyan">s</code>-> filter out from MAF mutations with SIFT annotation <ins>equals</ins> to to the ones specified in <i>conf.ini</i> SIFT field.

A few examples of usage are here provided:
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -f q
```
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -f cav
```
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output -v <vep_cache_path>:/vep_cache -v <genomes_folder>:/ref_fasta -v <conf.ini>:/conf.ini varan -i <input_folder> -o /output/<output_name> -c <type_of_cancer> -k -f divo
```

⚠️ More than one filter can be applied simultaneously<br>
</details>

<details open>
  <summary><i>Local version</i></summary>

Ex 1) <ins>Launch Varan base analysis with input folder</ins>:

Launch this command to process the contents of the input folder 

```
python varan.py -i <path_input_folder> -o <path_output_folder> -c <type_of_cancer>
```
Ex 2) <ins>Launch Varan base analysis with input file</ins>:

Launch one of these commands to process the contents of the input file(s)
<li>If you have sample.tsv and patient.tsv as input:</li>

```
python varan.py -i sample.tsv patient.tsv -o output_folder -c mixed
```
<li>If you have sample.tsv and fusion.tsv as input:</li>

```
python varan.py -i sample.tsv "" fusion.tsv -o output_folder -c mixed
```
<li>If you have sample.tsv, patient.tsv and fusion.tsv as input:</li>

```
python varan.py -i sample.tsv patient.tsv fusion.tsv -o output_folder -c mixed
```

Ex 3) <ins>Multiple vcf analysis</ins>: 

Launch this command to specify that it is a multi-sample file.
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -m
```
Ex 4) <ins>Overwrite analysis</ins>:

Launch this command to overwrite the output folder
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -w
```
Ex 5) <ins>Resume analysis</ins>:

Launch this command to resume an already started analysis
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -R
```
Ex 6) <ins>Specify analysis</ins>:

Launch one of these commands to specify the analysis 

* snv -> only snv analysis.
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -t snv
```
* cnv -> only cnv analysis.
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -t cnv
```
* fus -> only fusion analysis.
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -t fus
```
* tab -> only meta file creation.
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -t tab
```
##### Ex 7) <ins>Filter vcf/maf</ins>:

Launch this command to specify the type of filter to use:

* <code style="color : cyan">d</code> -> filter out from snv mutations with ALT="." and FILTER ≠"PASS" 
* <code style="color : cyan">p</code> -> filter out from MAF mutations with FILTER ≠"PASS".
* <code style="color : cyan">v</code>-> filter out from MAF mutations with vaf (t_VF column) values <ins>not in</ins> the ranges [t_VAF_min; t_VAF_max] specified in <i>conf.ini</i>.
* <code style="color : cyan">n</code> -> apply a specific VAF filter to novel mutations (dbSNP_RS = "novel") filtering out from MAF novel mutations with vaf <ins>not in</ins> the range specified in <i>conf.ini</i> t_VAF_min_novel.
* <code style="color : cyan">a</code> -> filter out from MAF mutations with AF values <ins>not in</ins> the range specified in <i>conf.ini</i>. You can choose to drop or keep mutations with no value for this column via 'drop_NA_AF' setting.
* <code style="color : cyan">o</code>-> filter out from MAF mutations with ONCOGENIC values <ins>different</ins> to the ones specified in <i>conf.ini</i> ONCOKB_FILTER field. This requires the OncoKB annotation (-k option).
* <code style="color : cyan">c</code> -> filter out from MAF mutations with CLIN_SIG values <ins>equals</ins> to the ones specified in <i>conf.ini</i> CLIN_SIG field.
* <code style="color : cyan">q</code> -> filter out from MAF mutations with Consequences values <ins>different</ins> to the ones specified in <i>conf.ini</i> CONSEQUENCES field.
* <code style="color : cyan">y</code>-> filter out from MAF mutations with PolyPhen annotation <ins>different</ins> to the ones specified in <i>conf.ini</i> POLYPHEN field.
* <code style="color : cyan">i</code>-> filter out from MAF mutations with IMPACT annotation <ins>equals</ins> to to the ones specified in <i>conf.ini</i> IMPACT field.
* <code style="color : cyan">s</code>-> filter out from MAF mutations with SIFT annotation <ins>equals</ins> to to the ones specified in <i>conf.ini</i> SIFT field.

A few examples of usage are here provided:
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -f q
```
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -f cav
```
```
python varan.py -i <path_to_sample_file> -o <path_output_folder> -c <type_of_cancer> -k -f divo
```

⚠️ More than one filter can be applied simultaneously<br>

</details>
</details>

#### 3. Output

<p align="justify"> After varan.py has run successfully, the resulted output folder will have the following organization and content:

```
study_name
├── case_lists
│   ├── cases_cna.txt
│   ├── cases_sequenced.txt
│   └── cases_sv.txt
├── data_clinical_patient.txt
├── data_clinical_sample.txt
├── data_cna_hg19.seg
├── data_cna_hg19.seg.fc.txt
├── data_cna.txt
├── data_mutations_extended.txt
├── data_sv.txt
├── img
├── MAF_Filtered or MAF_Onco_filtered(*)
├── MAF_OncoKB (**)
├── maf or maf.zip (***)
├── meta_clinical_patient.txt
├── meta_clinical_sample.txt
├── meta_cna_hg19_seg.txt
├── meta_cna.txt
├── meta_mutations_extended.txt
├── meta_study.txt
└── meta_sv.txt
```

(\*) Filtered MAF will be stored respectively inside MAF_Onco_filtered or MAF_Filtered if MAFs were annotated with oncokb or not.<br>
(\**) MAF_OncoKB folder will be create only if -k option is set and will contain all the MAF with oncoKB annotations.
(\***) A folder containing maf files. It will zipped if ZIP_MAF is set as True in conf.ini.

<br>

### Block Two: STUDY MANIPULATION

<br>

<p align="center">
  <img width="500" height="330" src="readme_content/img/block2.png">
</p>

#### 1. Preparing Input

<p align="justify">The input for this block is a study folder correctly populated. It can be the output of the first block or an existing study folder downloaded from cBioPortal. 

Based on the requested type of action (Update, Extract, Remove), the other input can be another study folder or a sample list tsv file.

#### 2. Launch Varan main

The possible option to launch Varan main for block 2 are:

| Options | Input | Type | Required
|----------------|----------------| :---:| :---:|
|-u <br> --Update| <p align="justify">Add this option if you want to update an existing study folder or merge two studies.| boolean | One between -u, -e or -r is required
|-e <br> --Extract| <p align="justify">Add this option if you want to extract samples from an existing study folder.| boolean | One between -u, -e or -r is required
|-r <br> --Remove| <p align="justify">Add this option if you want to remove samples from an existing study folder.| boolean | One between -u, -e or -r is required
|-p <br> --Path| <p align="justify">Add this option to specify the path of the existing study folder to update, or from which to remove/extract samples.| string | Yes
|-n <br> --NewPath| <p align="justify">Add this option to specify the path of the study folder containing updated/new information.| string | Only if the -u option is selected
|-s <br> --SampleList| <p align="justify">Add this option to insert the path of the .txt file containing the list of samples to remove/extract from the study folder.| string | Only if the -e or -r option is selected
|-N <br> --Name| <p align="justify">Add this option if you want to customize the study name (studyID).| string | No
|-o <br> --output_folder| <p align="justify">Add this option to specify the path where to save the output folder. If not provided, it will be created a new version of the existing study.| string | Only if the -N option is selected

<br>
<details open>
  <summary><i>Examples</i></summary>

<details open>
  <summary><i>Docker version</i></summary>

To launch Varan docker version is mandatory to mount at least the input and the output folder for granting a correct functioning.
```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output varan <commands>
```
Ex 1) <ins>Update a study folder</ins>: 
<p align="justify">Launch this command to update a study folder

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output varan -u -p /input/<path_to_old_study_folder> -n /input/<path_to_new_study_folder> -o /output/<path_to_output_folder>
```
Ex 2) <ins>Extract a study with a subset of samples</ins>: 
<p align="justify">Launch this command to extract a list of samples from a study folder and create a new study containing only these samples in the output path.

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output  varan -e -p /input/<path_to_study_folder> -s /input/<path_to_sample_list_file> -o /output/<path_to_output_folder>
```

Ex 3) <ins>Remove samples from a study</ins>:
<p align="justify">Launch this command to remove a list of samples from a study folder and save a new study without them in the output path, assigning a customized the study name.

```
docker run --rm -it -v <input_folder>:/input -v <output_folder>:/output varan -r -p /input/<path_to_study_folder> -s /input/<path_to_sample_list_file> -o /output/<path_to_output_folder> -N <new_studyID_in_meta>
```

</details>

<details open>
  <summary><i>Local version</i></summary>

Ex 1) <ins>Update a study folder</ins>: 
<p align="justify">Launch this command to update a study folder and create a new version of the original study folder.

```
python varan.py -u -p <path_to_old_study_folder> -n <path_to_new_study_folder>
```
Ex 2) <ins>Extract a study with a subset of samples</ins>: 
<p align="justify">Launch this command to extract a list of samples from a study folder and create a new study containing only these samples in the output path.

```
python varan.py -e -p <path_to_study_folder> -s <path_to_sample_list_file> -o <path_to_output_folder>
```

Ex 3) <ins>Remove samples from a study</ins>:
<p align="justify">Launch this command to remove a list of samples from a study folder and save a new study without them in the output path, assigning a customized the study name.

```
python varan.py -r -p <path_to_study_folder> -s <path_to_sample_list_file> -o <path_to_output_folder> -N <new_studyID_in_meta>
```
</details>

#### 3. Output

After varan.py has run successfully, the resulting <i>output_folder</i> will have the organization reported in [block 1](#3-output).
