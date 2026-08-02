FROM ensemblorg/ensembl-vep:release_111.0
USER root
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update &&\
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl gcc g++ gnupg unixodbc-dev openssl git &&\
    apt-get install -y software-properties-common ca-certificates &&\
    apt-get install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev libssl-dev libreadline-dev libffi-dev wget libbz2-dev libsqlite3-dev vcftools bcftools samtools && \
    update-ca-certificates && \
    apt-get update && apt-get install -y python-setuptools && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /python && cd /python && \
    wget https://www.python.org/ftp/python/3.10.8/Python-3.10.8.tgz && \
    tar -zxvf Python-3.10.8.tgz && \
    cd Python-3.10.8 && \
    ls -lhR && \
    ./configure --enable-optimizations && \
    make install && \
    rm -rf /python

COPY requirements.txt /requirements.txt

WORKDIR /

#download oncokb
RUN git clone https://github.com/oncokb/oncokb-annotator.git
RUN python3 -m pip install -r requirements.txt
RUN python3 -m pip install -r oncokb-annotator/requirements/pip3.txt
RUN python3 -m pip install -r oncokb-annotator/requirements/common.txt

#download vcf2maf
ENV VCF2MAF_URL https://api.github.com/repos/mskcc/vcf2maf/tarball/v1.6.22
RUN curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL
RUN tar -zxf mskcc-vcf2maf.tar.gz --one-top-level=mskcc-vcf2maf --strip-components 1
RUN rm mskcc-vcf2maf.tar.gz

#download clinvar
ENV CLINV_URL https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37
#uncomment to install version GRCh38
# ENV CLINV_URL https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38
RUN mkdir clinvar
RUN wget -P /clinvar $CLINV_URL/clinvar.vcf.gz
RUN wget -P /clinvar $CLINV_URL/clinvar.vcf.gz.tbi

#download importer for validator
ENV CBIO_URL https://github.com/cBioPortal/cbioportal-core.git
RUN git clone $CBIO_URL
RUN mv cbioportal-core/scripts/importer .
RUN rm -r cbioportal-core

#install R + SigMA (mutational signature / HRD analysis, run_sigma.R),
#only exercised when varan.py is called with -g/--sigma - see
#SIGMA_INTEGRATION_FEASIBILITY.md and sigma_runner.py for the full design.
#System libraries here are what SigMA's own Bioconductor/CRAN dependency
#tree needs to build from source - libuv1-dev/libharfbuzz-dev/
#libfribidi-dev/libfreetype-dev/libtiff5-dev/libjpeg-dev were added after
#an actual build of this exact block failed on the 'fs' and 'textshaping'
#packages (missing uv.h / hb-ft.h respectively), which cascade-failed
#devtools's own dependency tree (usethis/pkgdown/pkgload/roxygen2/
#testthat all transitively need 'fs'; ragg/textshaping feed rmarkdown/
#bslib/shiny) - see SIGMA_INTEGRATION_FEASIBILITY.md's conda-vs-Docker
#note for why the conda path (prebuilt binaries, no header-hunting) is
#recommended over this one where a conda/mamba environment is available.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        r-base r-base-dev \
        libcurl4-openssl-dev libxml2-dev libpng-dev liblzma-dev libbz2-dev \
        libglpk-dev libuv1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype-dev libtiff5-dev libjpeg-dev gfortran && \
    rm -rf /var/lib/apt/lists/*

#Bioconductor packages first (binary/source via BiocManager), then SigMA
#itself from GitHub (it isn't on CRAN/Bioconductor). Only the hg19 BSgenome
#is installed - Varan's own pipeline is hg19/GRCh37 end to end today (see
#conf.ini's VEP_DATA cache and the CNA pipeline's data_cna_hg19.* naming);
#add BSgenome.Hsapiens.UCSC.hg38 here too if Varan ever grows hg38 support.
RUN R -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' && \
    R -e 'BiocManager::install(c( \
        "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "VariantAnnotation", \
        "GenomicRanges", "IRanges", "gbm", "nnls", "reshape2", "Rmisc", \
        "DT", "gridExtra", "ggplot2", "devtools"), update = FALSE, ask = FALSE)' && \
    R -e 'devtools::install_github("parklab/SigMA", dependencies = FALSE, upgrade = "never")' && \
    R -e 'library(SigMA)'

COPY . /

ENTRYPOINT [ "python3", "/varan.py"]
