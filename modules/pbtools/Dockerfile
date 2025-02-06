# Use a base image with Python

FROM python:3.8-slim-buster AS build

# Set environment variables for URLs
ENV PBCPG_URL=https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
ENV TRGT_URL=https://github.com/PacificBiosciences/trgt/releases/download/v1.2.0/trgt-v1.2.0-x86_64-unknown-linux-gnu.tar.gz
ENV MINICONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh



RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    gcc \
    g++ \
    git \
    ibssl1.0.0 \
    zlib1g \
    zlib1g-dev \
    default-libmysqlclient-dev \
    libdbd-mysql-perl \ 
    libkrb5-3 \
    procps \
    samtools \
    bcftools \
    && wget $PBCPG_URL -O /usr/local/bin/pb-CpG-tools.tar.gz \
    && tar -xzf /usr/local/bin/pb-CpG-tools.tar.gz -C /usr/local/bin/ \
    && rm /usr/local/bin/pb-CpG-tools.tar.gz \
    && ln -s /usr/local/bin/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores /usr/local/bin/aligned_bam_to_cpg_scores \
    && wget $TRGT_URL -O /usr/local/bin/trgt-v1.2.0-x86_64-unknown-linux-gnu.tar.gz \
    && gunzip /usr/local/bin/trgt-v1.2.0-x86_64-unknown-linux-gnu.tar.gz \ 
    && tar -xvf /usr/local/bin/trgt-v1.2.0-x86_64-unknown-linux-gnu.tar -C /usr/local/bin/ \
    && ln -s /usr/local/bin/trgt-v1.2.0-x86_64-unknown-linux-gnu/trgt /usr/local/bin/trgt \
    && rm -rf /var/lib/apt/lists/*


RUN wget $MINICONDA_URL -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh


# Add conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Install pbmm2 using conda
RUN conda install -c bioconda -y pbmm2  \
    && conda clean -afy

# install hificnv
RUN conda install -y  hificnv \
    && conda clean -afy

# install pbsv 
RUN conda install -y  -c bioconda pbsv hiphase  \
    && conda clean -afy

#pbgzip
RUN conda install bioconda::pbgzip \
    && conda clean -afy