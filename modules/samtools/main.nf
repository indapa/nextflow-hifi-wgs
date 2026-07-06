
process downsample {

    cpus 4
    memory '4 GB'

    tag { "${sample_id}" }

    publishDir { "${params.downsample_output_dir}/${sample_id}" }, mode: 'copy', overwrite: true
    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), val(coverage), val(frac), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.${coverage}x.bam"), emit: downsampled_bam

    script:
    """
    samtools view -h -s 42.${frac} -b -@ ${task.cpus} -o ${sample_id}.${coverage}x.bam  ${bam}
    """

    stub:
    """
    touch ${sample_id}.${coverage}x.bam
    """
}



process bam_stats {

    /* generate stats for the bam file */
  
    cpus 2
    memory '2 GB'

    tag { "${sample_id}" }

    publishDir { "${params.bam_stats_output_dir}/${sample_id}" }, mode: 'copy', overwrite: true

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${bam}.stats"), emit: stats



    script:
    """

    samtools stats --threads ${task.cpus} ${bam} > ${bam}.stats

    """

    stub:
    """

    touch ${bam}.stats

    """
}



process samtools_index{
    tag { "${sample_id}" }
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0" // or any container with samtools

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """

    stub:
    """
    touch ${bam}.bai    
    """
}