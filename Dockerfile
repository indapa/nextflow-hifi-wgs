FROM mambaorg/micromamba:1.5.10-noble

# Copy the conda environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml

# Install dependencies into the base environment
# 'micromamba install' will now fetch the pre-built whatshap binary
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba clean -a -y

# Set the container to run as root (optional, based on your previous file)
USER root

# Ensure the conda environment is on the PATH
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"