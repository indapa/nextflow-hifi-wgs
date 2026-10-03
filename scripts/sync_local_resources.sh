#!/usr/bin/env bash
# Copy the pipeline's S3 reference resources to a local directory so the
# pipeline can run on a workstation without AWS credentials.
#
# Run once by someone with read access to s3://pacbio-datasets and
# s3://aindap-pb-resources, then point local_params.yaml at <dest_dir>.
#
# Usage: scripts/sync_local_resources.sh <dest_dir> [--skip-large]
#   --skip-large  skip the GRCh38 FASTA (~3 GB) and pbmm2 .mmi (~5.4 GB)
#
# Layout written under <dest_dir>:
#   reference/human_GRCh38_no_alt_analysis_set.{fasta,fasta.fai,mmi}
#   resources/...   mirrors s3://pacbio-datasets/resources/

set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

dest="$1"
skip_large="${2:-}"

ref_s3="s3://pacbio-hifi-human-wgs-reference/dataset/hifi-wdl-resources-v2.0.0/GRCh38"
mmi_s3="s3://aindap-pb-resources"
res_s3="s3://pacbio-datasets/resources"

mkdir -p "${dest}/reference" "${dest}/resources"

# Copy a single object unless it already exists locally
fetch() {
    local src="$1" out="$2" extra="${3:-}"
    if [[ -s "$out" ]]; then
        echo "exists: $out"
        return
    fi
    mkdir -p "$(dirname "$out")"
    # shellcheck disable=SC2086
    aws s3 cp $extra "$src" "$out"
}

# --- reference genome (public bucket) + pbmm2 index ---
fetch "${ref_s3}/human_GRCh38_no_alt_analysis_set.fasta.fai" "${dest}/reference/human_GRCh38_no_alt_analysis_set.fasta.fai" --no-sign-request
if [[ "$skip_large" != "--skip-large" ]]; then
    fetch "${ref_s3}/human_GRCh38_no_alt_analysis_set.fasta" "${dest}/reference/human_GRCh38_no_alt_analysis_set.fasta" --no-sign-request
    fetch "${mmi_s3}/human_GRCh38_no_alt_analysis_set.mmi"   "${dest}/reference/human_GRCh38_no_alt_analysis_set.mmi"
fi

# --- sawfish / TRGT ---
fetch "${res_s3}/adotto_strchive_20250827.hg38.bed.gz"                           "${dest}/resources/adotto_strchive_20250827.hg38.bed.gz"
fetch "${res_s3}/expected_cn/expected_cn.hg38.XX.bed"                            "${dest}/resources/expected_cn/expected_cn.hg38.XX.bed"
fetch "${res_s3}/expected_cn/expected_cn.hg38.XY.bed"                            "${dest}/resources/expected_cn/expected_cn.hg38.XY.bed"
fetch "${res_s3}/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz"     "${dest}/resources/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz"
fetch "${res_s3}/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz.tbi" "${dest}/resources/cnv_excluded_regions/annotation_and_common_cnv.hg38.bed.gz.tbi"

# --- GLnexus ---
fetch "${res_s3}/glnexus/all_contigs.bed" "${dest}/resources/glnexus/all_contigs.bed"

# --- FastVEP ---
fetch "${res_s3}/fastvep/GFFs/gencode.v50.annotation.gff3.gz" "${dest}/resources/fastvep/GFFs/gencode.v50.annotation.gff3.gz"
aws s3 sync "${res_s3}/fastvep/SA_files/" "${dest}/resources/fastvep/SA_files/"

# --- DeepVariant / DeepTrio 50 Mb interval chunks ---
aws s3 sync "${res_s3}/deepvariant-regions-chunks/" "${dest}/resources/deepvariant-regions-chunks/"

echo
du -sh "${dest}"
echo "Done. In local_params.yaml, replace /data/hifi-wgs-resources with: $(cd "${dest}" && pwd)"
