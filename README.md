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
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/docs/) to run the containers
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

To run on the local executor without Wave or Fusion, add `-profile local`.

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
- **Profiles:** `local` uses the local executor with 1 CPU and 2 GB per task, and turns off Wave and Fusion

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
