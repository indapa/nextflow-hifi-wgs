FROM python:3.12-slim

# Install procps so Nextflow can track process metrics
RUN apt-get update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/*

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Copy the script directly into a system PATH directory
COPY bin/parse-alignment-stats.py /usr/local/bin/parse-alignment-stats.py

# Make it executable
RUN chmod +x /usr/local/bin/parse-alignment-stats.py

# Optional: Set a working directory for container runtime
WORKDIR /data