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



## Modules
The pipeline is modular, with separate modules for each major step:
- `pbtools`: PacBio secondary analysis tools.
- `deepvariant`: DNA small variant calling using `DeepVariant`.
- `samtools`: For BAM file processing.



