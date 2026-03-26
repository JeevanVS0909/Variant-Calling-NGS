############################################################
# ULTRA CLINICAL WES PIPELINE
############################################################

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app

############################################################
# SYSTEM (OPTIMIZED)
############################################################

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    build-essential \
    openjdk-17-jdk \
    perl \
    r-base \
    samtools \
    bcftools \
    bwa \
    fastp \
    tabix \
    bedtools \
    vcftools \
    libbz2-1.0 \
    liblzma5 \
    libgomp1 \
    zlib1g \
    cpanminus \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

############################################################
# GATK
############################################################

RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip \
    && unzip gatk-4.5.0.0.zip \
    && mv gatk-4.5.0.0 /opt/gatk \
    && ln -s /opt/gatk/gatk /usr/local/bin/gatk

############################################################
# VEP (INSTALL + CACHE GRCh38)
############################################################

# Install required Perl modules
RUN cpanm DBI JSON LWP::Simple

# Clone and install VEP
RUN git clone https://github.com/Ensembl/ensembl-vep.git /opt/vep \
    && cd /opt/vep \
    && perl INSTALL.pl --AUTO a --NO_TEST

ENV PATH="/opt/vep:${PATH}"

############################################################
# vcf2maf
############################################################

RUN git clone https://github.com/mskcc/vcf2maf.git /opt/vcf2maf
ENV PATH="/opt/vcf2maf:${PATH}"

############################################################
# Install Miniconda + Configure + Install CNVkit 
############################################################

RUN apt-get update && apt-get install -y wget bzip2 ca-certificates \
    && wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm Miniconda3-latest-Linux-x86_64.sh \
    && /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
    && /opt/conda/bin/conda config --system --add channels defaults \
    && /opt/conda/bin/conda config --system --add channels bioconda \
    && /opt/conda/bin/conda config --system --add channels conda-forge \
    && /opt/conda/bin/conda config --system --set channel_priority strict \
    && /opt/conda/bin/conda install -y mamba -n base -c conda-forge \
    && /opt/conda/bin/mamba install -y cnvkit python=3.10 pandas=1.5 numpy scipy matplotlib pysam \
    && /opt/conda/bin/conda clean -afy

ENV PATH="/opt/conda/bin:$PATH"


############################################################
# MSI SENSOR 2
############################################################

RUN apt-get update && apt-get install -y git build-essential \
    && git clone https://github.com/niu-lab/msisensor2.git \
    && cd msisensor2 \
    && chmod +x msisensor2 \
    && mv msisensor2 /usr/local/bin/

############################################################
# CancerVar 
############################################################

# 1. Install core Python dependencies (from default PyPI)
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scikit-learn \
    cython

# 2. Install PyTorch separately (CPU version)
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu

# 3. Clone only latest version (shallow clone → faster)
RUN git clone --depth 1 https://github.com/WGLab/CancerVar.git /opt/CancerVar

# 4. Set working directory
WORKDIR /opt/CancerVar

# 6. Add to PATH
ENV PATH="/opt/CancerVar:/opt/CancerVar/bin:${PATH}"

# 7. Lightweight test
RUN python3 CancerVar.py -h

############################################################
# FIX WORKDIR 
############################################################

WORKDIR /app

############################################################
# Python libs
############################################################

COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Ensure CIViC API support
RUN pip3 install requests

############################################################
# Copy pipeline
############################################################

COPY pipeline.py .

############################################################
# Create mount points
############################################################

RUN mkdir data resources results

############################################################
# ENTRY
############################################################

CMD ["python3","pipeline.py"]
