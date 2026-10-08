[Nextflow](https://www.nextflow.io/) secondary analysis pipeline for [PacBio HiFi](https://downloads.pacbcloud.com/public/revio/2022Q4/?utm_source=Website&utm_medium=webpage&utm_term=HomoSapiens-GIAB-trio-HG002-4&utm_content=datasets&utm_campaign=0000-Website-Leads) whole-genome sequencing. It runs on singletons (one sample) and on trios (child plus two parents).

## Main Features
- Read alignment with `pbmm2`
- Small variant calling with `DeepVariant` (singletons) or `DeepTrio` (trios)
- Trio joint calling with `GLnexus`
- Phasing: `HiPhase` (read-backed, for singletons and trio parents) and `WhatsHap` (pedigree-aware, for the trio child)
- Variant annotation with `FastVEP`
- Structural variant calling with `sawfish2`
- Tandem repeat genotyping with `TRGT`
- 5mC methylation calling with `pb-CpG-tools`
- Coverage and sex inference with `mosdepth`

## Requirements
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/), [Apptainer](https://apptainer.org/) or [Singularity](https://sylabs.io/docs/) to run the containers
- [Seqera Platform](https://seqera.io/) is recommended. The default config turns on Wave, Fusion and the `nf-tower` plugin.

## Quick Start
```bash
git clone https://github.com/indapa/nextflow-hifi-wgs.git
cd nextflow-hifi-wgs

# List the available entrypoints
nextflow run main.nf --help

# Default: align and analyze singletons
nextflow run main.nf --samplesheet samples.csv
```

To do a stub run on the local executor without Wave or Fusion, add `-profile test -stub`. The `-stub` flag is what makes Nextflow run the `stub:` blocks; the profile only sets resources.

## Entrypoints
Choose an entrypoint with `-entry <NAME>`. If you don't name one, the pipeline runs `WGS_SINGLETON`.

| Entrypoint | Input | What it runs |
|---|---|---|
| `WGS_SINGLETON` (default) | `--samplesheet`: unaligned BAMs | Alignment, then everything in `POST_ALIGNMENT` |
| `POST_ALIGNMENT_ONLY` | `--samplesheet`: aligned BAMs | `POST_ALIGNMENT` only |
| `WGS_TRIO` | `--trio_samplesheet`: unaligned BAMs | Alignment, then everything in `RUN_TRIO_PIPELINE` |
| `WGS_TRIO_ALIGNED` | `--trio_aligned_samplesheet`: aligned BAMs | `RUN_TRIO_PIPELINE` only |
| `TRIO_FROM_DEEPTRIO` | `--trio_gvcf_samplesheet`: aligned BAMs plus existing DeepTrio VCFs and gVCFs | GLnexus, WhatsHap, FastVEP, HiPhase (parents) and CpG. Skips DeepTrio, sawfish and TRGT. |

### Singleton workflow (`POST_ALIGNMENT`)
1. `samtools` BAM stats, `mosdepth` coverage, sex inference and a coverage distribution plot
2. `DEEPVARIANT_SINGLETON_WGS`: DeepVariant on 50 Mb chunks, merged per chromosome and then genome-wide
3. `hiphase_small_variants`: read-backed phasing and haplotagging
4. `FASTVEP_ANNOTATE_WGS`: annotation of the phased VCF
5. `cpg_methylation_calling` on the haplotagged BAM
6. `sawfish` discover and joint-call, using the expected-copy-number BED that matches the inferred sex
7. `TRGT_GENOTYPING`

### Trio workflow (`RUN_TRIO_PIPELINE`)
1. `DEEPTRIO_WGS`: DeepTrio on 50 Mb chunks, merged per chromosome and then genome-wide. Produces a VCF and a gVCF for each sample.
2. `GLNEXUS_TRIO`: joint calling per chromosome, merged into one VCF per family
3. `WHATSHAP_TRIO_PHASE_BY_CHROM`: pedigree-aware phasing per chromosome, concatenation, then stats and haplotagging of the child's BAM
4. `FASTVEP_ANNOTATE_TRIO_INDIVIDUALS`: splits the phased trio VCF into one VCF per sample and annotates each one
5. `HIPHASE_TRIO_PARENTS`: read-backed phasing and haplotagging of each parent, using that parent's own DeepTrio VCF
6. `CPG_TRIO_METHYLATION`: CpG calling on the haplotagged BAMs of the child (from WhatsHap) and the parents (from HiPhase)
7. `mosdepth` and sex inference, then `SAWFISH_TRIO` (per-sample discover, then a family joint-call)
8. `TRGT_GENOTYPING` on the aligned BAMs, which are not haplotagged

`TRIO_FROM_DEEPTRIO` runs steps 2 to 6 on DeepTrio results you already have.

## Samplesheets
All samplesheets are CSV files with a header row.

**Singletons, unaligned** (`WGS_SINGLETON`, `--samplesheet`):
```
sample_id,bam_file
sample1,/path/to/sample1.unaligned.bam
```

**Singletons, aligned** (`POST_ALIGNMENT_ONLY`, `--samplesheet`):
```
sample_id,bam_file,bai_file
sample1,/path/to/sample1.bam,/path/to/sample1.bam.bai
```

**Trio, unaligned** (`WGS_TRIO`, `--trio_samplesheet`). Use one row per sample. `role` must be `child`, `parent1` or `parent2`.
```
family_id,sample_id,role,bam
FAM1,HG002,child,/path/to/HG002.unaligned.bam
FAM1,HG003,parent1,/path/to/HG003.unaligned.bam
FAM1,HG004,parent2,/path/to/HG004.unaligned.bam
```

**Trio, aligned** (`WGS_TRIO_ALIGNED`, `--trio_aligned_samplesheet`):
```
family_id,sample_id,role,aligned_bam,aligned_bai
FAM1,HG002,child,/path/to/HG002.bam,/path/to/HG002.bam.bai
...
```

**Trio from DeepTrio results** (`TRIO_FROM_DEEPTRIO`, `--trio_gvcf_samplesheet`). Here `vcf` is each sample's DeepTrio VCF, which HiPhase needs for the parents, and `gvcf` is the DeepTrio gVCF, which GLnexus needs.
```
family_id,sample_id,role,aligned_bam,aligned_bai,vcf,vcf_tbi,gvcf,gvcf_tbi
FAM1,HG002,child,/path/HG002.bam,/path/HG002.bam.bai,/path/HG002.vcf.gz,/path/HG002.vcf.gz.tbi,/path/HG002.g.vcf.gz,/path/HG002.g.vcf.gz.tbi
...
```

## Configuration
All parameters live in [nextflow.config](nextflow.config) and you can override any of them on the command line, for example `--output_dir s3://bucket/run1`. The main groups are:
- **Reference:** `reference`, `reference_index`, and `mmi`, the pbmm2 index used by `WGS_TRIO`
- **Intervals:** `intervals_dir` (DeepTrio chunk BEDs), `bed_dir` (DeepVariant chunk BEDs), `chromosomes` (chromosomes for GLnexus and WhatsHap), `glnexus_region_bed`
- **Annotation:** `fastvep_gff`, `fastvep_sa_dir`
- **SV and repeats:** `expected_XX_bed`, `expected_XY_bed`, `excluded_bed`, `trgt_repeats_bed`, `trgt_max_depth`, `trgt_min_mapq`
- **Profiles:**
  - `test`: local executor, 1 CPU and 2 GB per task, Wave and Fusion off. Use with `-stub`.
  - `local`: local executor, with each task capped at `--max_cpus` / `--max_memory`. Wave and Fusion off. Combine with a container profile.
  - `docker`, `apptainer`, `singularity`: choose the container engine. Apptainer and Singularity pull the same images as Docker.

## Running Locally
For a single large workstation instead of AWS Batch. No process definitions change; you only choose a container engine with a profile.

1. **Download the reference resources.** All of them are public. Pick a resources directory and run these steps once, from the repository root:
   ```bash
   RES=/data/hifi-wgs-resources
   mkdir -p "$RES"
   ```
   a. **PacBio HiFi-human-WGS reference bundle** ([Zenodo 14908106](https://zenodo.org/records/14908106)): the GRCh38 FASTA and `.fai`, the TRGT repeat catalog, and the sawfish expected-copy-number and CNV-exclusion BEDs. The tar is 9.5 GB. It unpacks to `hifi-wdl-resources-v2.0.0/GRCh38/`. The pipeline doesn't use the gnomAD and CoLoRSdb files (about 5 GB), so you can delete them after extracting.
   ```bash
   curl -L -o "$RES/hifi-wdl-resources-v2.1.0.tar" \
     "https://zenodo.org/records/14908106/files/hifi-wdl-resources-v2.1.0.tar?download=1"
   tar -xf "$RES/hifi-wdl-resources-v2.1.0.tar" -C "$RES" && rm "$RES/hifi-wdl-resources-v2.1.0.tar"
   GRCH38="$RES/hifi-wdl-resources-v2.0.0/GRCh38"
   ```
   b. **FastVEP gene models**: the GENCODE v50 GFF3 from [gencodegenes.org](https://www.gencodegenes.org/human/).
   ```bash
   mkdir -p "$RES/fastvep/SA_files"
   curl -L -o "$RES/fastvep/gencode.v50.annotation.gff3.gz" \
     https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz
   ```
   c. **FastVEP ClinVar supplementary annotation**: build `clinvar.osa2` from the NCBI ClinVar VCF, following the [fastVEP local setup guide](https://github.com/Huang-lab/fastVEP#local-setup-guide). `sa-build` only converts files; it does not download them. A good build is tens of MB. Anything under 1 MB means the input was empty.
   ```bash
   curl -L -o "$RES/fastvep/clinvar.vcf.gz" \
     https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
   docker run --rm -v "$RES/fastvep:/work" -w /work docker.io/indapa/fastvep:0.2.0 \
     fastvep sa-build --source clinvar -i clinvar.vcf.gz -o SA_files/clinvar --assembly GRCh38
   ls -la "$RES/fastvep/SA_files/"
   ```
   d. **Files generated from the FASTA index**: the 50 Mb DeepVariant/DeepTrio chunk BEDs and the whole-chromosome BED for GLnexus.
   ```bash
   python3 Intervals/generate_bed.py "$GRCH38/human_GRCh38_no_alt_analysis_set.fasta.fai" "$RES/intervals"
   mkdir -p "$RES/glnexus"
   awk -v OFS='\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1, 0, $2}' \
     "$GRCH38/human_GRCh38_no_alt_analysis_set.fasta.fai" > "$RES/glnexus/all_contigs.bed"
   ```
   e. **pbmm2 index** (optional). Only the entrypoints that align unaligned BAMs (`WGS_SINGLETON`, `WGS_TRIO`) need it. It takes about 10 minutes and the output is about 5.4 GB.
   ```bash
   mkdir -p "$RES/reference"
   docker run --rm -v "$RES:$RES" quay.io/pacbio/pbmm2:1.17.0_build1 \
     pbmm2 index --preset HIFI "$GRCH38/human_GRCh38_no_alt_analysis_set.fasta" \
     "$RES/reference/human_GRCh38_no_alt_analysis_set.mmi"
   ```
   These public files differ in three ways from the S3 defaults in [nextflow.config](nextflow.config):
   - The TRGT catalog is `adotto_repeats.updated_pathogenic_repeats` instead of `adotto_strchive_20250827`.
   - The CNV exclusion BED is HiFiCNV's `cnv.excluded_regions.common_50` instead of `annotation_and_common_cnv`.
   - The ClinVar release is the current one.

   If you can read the internal S3 buckets, you can mirror the defaults with `scripts/sync_local_resources.sh "$RES"` instead. That script writes a different directory layout, so adjust the paths in `local_params.yaml` to match.
2. **Edit [local_params.yaml](local_params.yaml).** Replace every `/data/hifi-wgs-resources` with your `$RES` directory, choose an output directory, and set `max_cpus` / `max_memory` to match the machine.
3. **Run** with Docker, or with Apptainer if there's no Docker daemon:
   ```bash
   nextflow run main.nf -profile local,docker    -params-file local_params.yaml -entry WGS_TRIO_ALIGNED --trio_aligned_samplesheet trios.csv
   nextflow run main.nf -profile local,apptainer -params-file local_params.yaml -entry WGS_TRIO_ALIGNED --trio_aligned_samplesheet trios.csv
   ```
   Samplesheet paths (BAMs, VCFs) can be local paths.

Notes:
- Apptainer caches the converted images in `.apptainer_cache/` (set with `--apptainer_cache`), so only the first run pulls them.
- The GPU processes (`deeptrio_wgs_50mb_chunk`, `deepvariant_wgs_parabricks`) use Docker-only `--gpus all` and aren't called by any entrypoint.

## Outputs
Outputs are written under `--output_dir`:

| Directory | Contents |
|---|---|
| `aligned-bams/` | pbmm2-aligned BAMs |
| `bamstats/` | samtools stats |
| `mosdepth/` | Coverage summaries and inferred sex |
| `deepvariant/` | DeepVariant and DeepTrio VCFs and gVCFs (trio results under `DV_trio/`), plus FastVEP-annotated VCFs |
| `glnexus/` | Joint-called trio VCFs |
| `whatshap/` | Trio-phased VCFs, phasing stats, and the child's haplotagged BAM |
| `hiphase/` | HiPhase phased VCFs and haplotagged BAMs |
| `cpg/` | pb-CpG-tools BED and bigWig files (combined, hap1, hap2) |
| `sawfish2/` | sawfish SV calls |
| `trgt/` | TRGT repeat genotypes |

## Repository Layout
```
main.nf          entrypoints, RUN_TRIO_PIPELINE and POST_ALIGNMENT
local_params.yaml  params template for local runs
scripts/         sync_local_resources.sh (mirrors the internal S3 resources; needs bucket access)
Intervals/       generate_bed.py (makes the 50 Mb chunk BEDs from a .fai)
modules/         processes: pbtools, deepvariant, glnexus, whatshap, fastvep,
                 samtools, mosdepth, alignment_metrics, ensemblvep
subworkflows/
  reference_alignment/          PBMM2_SPOT_WGS, PBMM2_ALIGN
  deepvariant_singleton_wgs/    DEEPVARIANT_SINGLETON_WGS (+ _PARABRICKS)
  deeptrio_wgs/                 DEEPTRIO_WGS
  glnexus_trio_merge/           GLNEXUS_TRIO
  whatshap_trio_phase_by_chrom/ WHATSHAP_TRIO_PHASE_BY_CHROM
  fastvep/                      FASTVEP_ANNOTATE_WGS
  fastvep_trio_individuals/     FASTVEP_ANNOTATE_TRIO_INDIVIDUALS
  hiphase_trio_parents/         HIPHASE_TRIO_PARENTS
  cpg_trio_methylation/         CPG_TRIO_METHYLATION
  sawfish_trio/                 SAWFISH_TRIO
  trgt/                         TRGT_GENOTYPING
  concat_and_split_wgs/         CONCAT_AND_SPLIT_WGS
```
