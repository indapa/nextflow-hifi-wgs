




process bam_stats {

    /* generate stats for the bam file */
    label 'low_memory'

    tag "${sample_id}"

    publishDir "${params.bam_stats_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${bam}.stats"), emit: stats



    script:
    """

    samtools stats --threads 4 ${bam} > ${bam}.stats

    """

    stub:
    """

    touch ${bam}.stats

    """
}



