# Container Manifest

These are the container images named in the `container` directives of the processes in `modules/`, mapped to the `main.nf` entrypoints that run them. `nextflow.config` doesn't override any container. Docker, Apptainer and Singularity all pull the same images.

Entrypoint columns: **WS** = `WGS_SINGLETON`, **PA** = `POST_ALIGNMENT_ONLY`, **WT** = `WGS_TRIO`, **WTA** = `WGS_TRIO_ALIGNED`, **TFD** = `TRIO_FROM_DEEPTRIO`.

## Images used by the entrypoints

| Image | Process(es) | WS | PA | WT | WTA | TFD |
|---|---|:-:|:-:|:-:|:-:|:-:|
| `quay.io/pacbio/pbmm2:26.2.99_build1` | `PBMM2_ALIGN_SPOT_CHUNK` | ✓ | | ✓ | | |
| `quay.io/pacbio/pbtk:3.5.0_build2` | `MAKE_PBI` | ✓ | | ✓ | | |
| `community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c` | `MERGE_SPOT_CHUNKS`, `bam_stats`, `slice_singleton_bam_by_interval`, `slice_trio_bams_by_interval` | ✓ | ✓ | ✓ | ✓ | |
| `quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0` | `samtools_index` | | | ✓ | ✓ | ✓ |
| `community.wave.seqera.io/library/mosdepth:0.3.10--259732f342cfce27` | `mosdepth_run` | ✓ | ✓ | ✓ | ✓ | |
| `indapa/mosdepth-sex:latest` | `infer_sex` (all four), `plot_dist_coverage` (WS, PA only) | ✓ | ✓ | ✓ | ✓ | |
| `google/deepvariant:1.10.0` | `DEEPVARIANT_CHUNK` | ✓ | ✓ | | | |
| `google/deepvariant:deeptrio-1.10.0` | `deeptrio_wgs_by_chrom` | | | ✓ | ✓ | |
| `community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11` | `concat_chrom_chunks_vcf_singleton`, `concat_full_genome_vcf_singleton` (WS, PA); `concat_chrom_chunks_vcf`, `concat_wgs_vcf` (WT, WTA) | ✓ | ✓ | ✓ | ✓ | |
| `quay.io/mlin/glnexus:v1.2.7` | `glnexus_trio_by_chrom` | | | ✓ | ✓ | ✓ |
| `quay.io/biocontainers/bcftools:1.21--h8b25389_0` | `concat_glnexus_vcf` | | | ✓ | ✓ | ✓ |
| `community.wave.seqera.io/library/bcftools_whatshap:16f1800bbb322710` | `WHATSHAP_PHASE_CHROM` | | | ✓ | ✓ | ✓ |
| `quay.io/biocontainers/bcftools:1.17--haef29d1_0` | `CONCAT_PHASED_VCFS`, `SPLIT_TRIO_VCF_BY_SAMPLE` | | | ✓ | ✓ | ✓ |
| `indapa/whatshap-tabix` (no tag, so `latest`) | `WHATSHAP_STATS_HAPLOTAG` | | | ✓ | ✓ | ✓ |
| `quay.io/pacbio/hiphase:1.5.0_build1` | `hiphase_small_variants` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `docker.io/indapa/fastvep:0.2.0` | `FASTVEP_ANNOTATE_SINGLETON_VCF` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `quay.io/pacbio/pb-cpg-tools:3.0.0_build1` | `cpg_methylation_calling` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `quay.io/pacbio/sawfish:2.2.1_build1` | `sawfish_discover`, `sawfish_joint_call` | ✓ | ✓ | ✓ | ✓ | |
| `quay.io/pacbio/trgt:5.1.0_build2` | `trgt` | ✓ | ✓ | ✓ | ✓ | |

Two of these images use floating tags, so a re-pull can change the tool version: `indapa/mosdepth-sex:latest` and `indapa/whatshap-tabix`. To get reproducible runs, pin them to a version tag or digest.

## Images defined but not called by any entrypoint

| Image | Process(es) |
|---|---|
| `quay.io/pacbio/pbmm2:1.17.0_build1` | `pbmm2_align` (also used in the README to build the `.mmi`) |
| `google/deepvariant:1.8.0` | `deepvariant_targeted_region` |
| `google/deepvariant:1.10.0` | `deepvariant_wgs` |
| `google/deepvariant:deeptrio-1.10.0` | `deeptrio_wgs`, `deeptrio_wgs_50mb_chunk` (GPU) |
| `nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1` | `deepvariant_wgs_parabricks` (GPU) |
| `quay.io/biocontainers/bcftools:1.17--haef29d1_0` | `bcftools_deepvariant_norm` |
| `docker.io/indapa/fastvep:0.2.0` | `FASTVEP_ANNOTATE_TRIO_VCF` |
| `indapa/whatshap-tabix` | `whatshap_trio_phase` |
| `community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c` | `downsample` (used by the standalone `downsample-bam.nf`) |
| `ensemblorg/ensembl-vep:latest` | `annotate_vep` |
| `indapa/hifi-wgs-pipeline:latest` | `PARSE_SAMTOOLS_STATS` |

## Pre-pulling images

To pull every image an entrypoint needs before running, for example on an offline workstation:
```bash
grep -rhoE "container +['\"][^'\"]+" modules | sed -E "s/container +['\"]//" | sort -u | xargs -n1 docker pull
```
This pulls every image in `modules/`, including the unused ones above. The Parabricks image is large and needs an NVIDIA GPU, so leave it out if you don't need it.
