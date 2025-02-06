# Use a base image with Python

FROM python:3.8-slim-buster AS build

# Set environment variables for URLs
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
    && rm -rf /var/lib/apt/lists/*


RUN wget $MINICONDA_URL -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh


# Add conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Install fibertools using conda https://fiberseq.github.io/fibertools/install.html
RUN conda install -y  -c conda-forge -c bioconda fibertools-rs htslib \
    && conda clean -afy