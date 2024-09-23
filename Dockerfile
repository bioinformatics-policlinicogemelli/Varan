#FROM buschlab/vep:110-GRCh37
FROM ensemblorg/ensembl-vep:release_111.0
USER root
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update &&\
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends curl gcc g++ gnupg unixodbc-dev openssl git &&\
    apt-get install -y software-properties-common ca-certificates &&\
    apt-get install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev libssl-dev libreadline-dev libffi-dev wget libbz2-dev libsqlite3-dev vcftools bcftools && \
    update-ca-certificates && \
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

COPY . /

CMD ["bash"]
#ENTRYPOINT [ "python3", "/varan.py"]
