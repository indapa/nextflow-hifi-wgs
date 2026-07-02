[Nextflow](https://www.nextflow.io/) secondary analysis pipeline for the analysis of [PacBio HiFi reads](https://downloads.pacbcloud.com/public/revio/2022Q4/?utm_source=Website&utm_medium=webpage&utm_term=HomoSapiens-GIAB-trio-HG002-4&utm_content=datasets&utm_campaign=0000-Website-Leads). 

## Main Features
- Read alignment using `pbmm2`.
- Variant calling using `DeepVariant`.
- Read backed phasing of small variants using `Hiphase`.
- SV calling with ```sawfish2```
- 5mC calling with ```pb-cpg-tools```


## Requirements
- [Nextflow](https://www.nextflow.io/) 
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/docs/) for containerized execution
- Recommend running on Seqera Platform

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




The input samplesheet should be formatted as follows:

```
sample_id,bam_file
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
```

The bam_file should point to un-aligned HiFi reads in BAM format.

For hiphase samplesheet:

```
sample_id,vcf,vcf_tbi,bam,bai
sample1,/path/to/sample1.vcf,/path/to/sample1.vcf.tbi,/path/to/sample1.bam,/path/to/sample1.bai
sample2,/path/to/sample2.vcf,/path/to/sample2.vcf.tbi,/path/to/sample2.bam,/path/to/sample2.bai
```
The bam and bai files should point to aligned HiFi reads in BAM format. The vcf and vcf_tbi files should point to variant calls in VCF format from DeepVariant along with their index files.

For trio samplesheet:
Needs to be in wide format
```
child_id,parent1_id,parent2_id,child_bam,parent1_bam,parent2_bam
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
- `pbtools`: PacBio secondary analysis tools.
- `deepvariant`: DNA small variant calling using `DeepVariant`.
- `samtools`: For BAM file processing.



