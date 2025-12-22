[Nextflow](https://www.nextflow.io/) secondary analysis pipeline for the analysis of [PacBio HiFi reads](https://downloads.pacbcloud.com/public/revio/2022Q4/?utm_source=Website&utm_medium=webpage&utm_term=HomoSapiens-GIAB-trio-HG002-4&utm_content=datasets&utm_campaign=0000-Website-Leads). 

## Main Features
- Read alignment using `pbmm2` with support for targeted regions.
- Variant calling using `DeepVariant` with support for targeted regions.
- Phasing of small variants using `Hiphase`.
- Variant annotation using `Ensembl VEP`.

## Requirements
- [Nextflow](https://www.nextflow.io/) version 20.10.0 or higher
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/docs/) for containerized execution
- Sufficient computational resources for processing HiFi reads

## Usage
1. Clone the repository:
   ```bash
    git clone https://github.com/indapa/nextflow-hifi-wgs.git
    cd nextflow-hifi-wgs
    ```

2. Prepare your input data and configuration file (`nextflow.config`).
3. Run the pipeline:
   ```bash
   nextflow run main.nf -c nextflow.config
   ```  

There are two entry points available:
- `-entry ALIGN_DEEP_VARIANT_HIPHASE_VEP_SYT1`: For alignment, DeepVariant variant calling in SYT1 region, phasing, and VEP annotation.
- `-entry HIPHASE_VEP_ONLY`: for phasing and VEP annotation only on pre-called variants and BAM files.


The input samplesheet should be formatted as follows:

```
sample_id,bam_path
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
```

For hiphase samplesheet:

```
sample_id,row.vcf,row.vcf_tbi,row.bam,row.bai
sample1,/path/to/sample1.vcf,/path/to/sample1.vcf.tbi,/path/to/sample1.bam,/path/to/sample1.bai
sample2,/path/to/sample2.vcf,/path/to/sample2.vcf.tbi,/path/to/sample2.bam,/path/to/sample2.bai
```

4. Monitor the progress and check the output files in the specified output directory.

## Configuration
The pipeline can be configured using the `nextflow.config` file. You can specify parameters such as:
- Reference genome path
- Reference gtf path
- Targeted regions string (if applicable)
- Computational resources (CPU, memory)
- Output directories


## Modules
The pipeline is modular, with separate modules for each major step:
- `pbtools`: For read alignment using `pbmm2`and phasing using `hiphase`.
- `deepvariant`: For variant calling using `DeepVariant`.
- `samtools`: For BAM file statistics.
- `ensemblvep`: For variant annotation using `Ensembl VEP`.


## Example

Witht the current configuration, the pipeline will align reads to the specified reference genome, call variants in the targeted SYT1 region, phase the variants, and annotate them using VEP.

```
nextflow run main.nf -entry ALIGN_DEEP_VARIANT_HIPHASE_VEP_SYT1 
```

For running only the phasing and VEP annotation steps:

```
nextflow run main.nf -entry HIPHASE_VEP_ONLY --hiphase_samplesheet path/to/hiphase_samplesheet.csv
```
