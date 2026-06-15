process PARSE_SAMTOOLS_STATS {

    /* parse samtools stats file and extract key stats including total_bases */
  
    cpus 2
    memory '4 GB'

    tag "${sample_id}"

    publishDir "${params.alignment_metrics_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    container 'indapa/hifi-wgs-pipeline:latest' // Update with your actual Docker Hub username

    input:
    tuple val(sample_id), path(stats)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam.stats.csv"), emit: metrics

    script:
    """
    parse-alignment-stats.py --statsfile ${stats} --sample_id ${sample_id} 
    """

    stub:
    """
    touch ${sample_id}.aligned.bam.stats.csv
    """
}