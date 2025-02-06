# Use a base image with Python
FROM python:3.8-slim-buster AS build

# Set environment variables for URLs
ENV MOSDEPTH_URL="https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth"

FROM python:3.8-slim-buster

ENV MOSDEPTH_URL="https://github.com/brentp/mosdepth/releases/download/v0.3.9/mosdepth"

# Install only what’s necessary: wget to download mosdepth, and maybe zlib1g if required
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    zlib1g \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Download and install mosdepth
RUN wget $MOSDEPTH_URL -O /usr/local/bin/mosdepth \
    && chmod +x /usr/local/bin/mosdepth \
    && rm -rf /var/lib/apt/lists/*