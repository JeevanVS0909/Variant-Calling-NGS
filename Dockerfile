# Base image: Ubuntu with Python
FROM ubuntu:22.04

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-21-jre-headless \
    python3 python3-pip \
    wget curl unzip git \
    bwa samtools bcftools \
    fastqc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir pandas

# Set DeepVariant version and directory
ENV DEEPVARIANT_VERSION=1.6.0
ENV DV_DIR=/opt/deepvariant
ENV PATH="$DV_DIR/bin:$PATH"

# Download and extract DeepVariant (CPU version)
RUN mkdir -p ${DV_DIR} && \
    curl -L https://github.com/google/deepvariant/releases/download/v${DEEPVARIANT_VERSION}/deepvariant-${DEEPVARIANT_VERSION}+cl-553159968.cpu.run \
    -o ${DV_DIR}/deepvariant.run && \
    chmod +x ${DV_DIR}/deepvariant.run && \
    ${DV_DIR}/deepvariant.run --extract --quiet && \
    rm -f ${DV_DIR}/deepvariant.run && \
    rm -rf ${DV_DIR}/runfiles && \
    strip ${DV_DIR}/bin/* || true

# Set working directory
WORKDIR /workdir
COPY . /workdir

# Default command
CMD ["python3", "pipeline.py"]
