
process slice_trio_bams_by_interval {
    tag { "${family_id} - ${interval_bed.baseName}" }
   
    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'
    
    

    input:
    tuple val(family_id), val(c_id), path(c_bam), path(c_bai), \
          val(p1_id), path(p1_bam), path(p1_bai), \
          val(p2_id), path(p2_bam), path(p2_bai), path(interval_bed)

    output:
    tuple val(family_id), val(c_id), path("${interval_bed.baseName}.child.bam"), path("${interval_bed.baseName}.child.bam.bai"), \
          val(p1_id), path("${interval_bed.baseName}.p1.bam"), path("${interval_bed.baseName}.p1.bam.bai"), \
          val(p2_id), path("${interval_bed.baseName}.p2.bam"), path("${interval_bed.baseName}.p2.bam.bai"), \
          path(interval_bed), emit: sliced_trio_package

    script:
    """
    # generation of regional mini-BAMs
    samtools view -@ ${task.cpus} -b -M -L ${interval_bed} ${c_bam} > ${interval_bed.baseName}.child.bam 
    samtools view -@ ${task.cpus} -b -M -L ${interval_bed} ${p1_bam} > ${interval_bed.baseName}.p1.bam 
    samtools view -@ ${task.cpus} -b -M -L ${interval_bed} ${p2_bam} > ${interval_bed.baseName}.p2.bam 
    

    samtools index ${interval_bed.baseName}.child.bam 
    samtools index ${interval_bed.baseName}.p1.bam 
    samtools index ${interval_bed.baseName}.p2.bam 
    """

    stub:
    """
    touch ${interval_bed.baseName}.child.bam
    touch ${interval_bed.baseName}.child.bam.bai
    touch ${interval_bed.baseName}.p1.bam
    touch ${interval_bed.baseName}.p1.bam.bai
    touch ${interval_bed.baseName}.p2.bam
    touch ${interval_bed.baseName}.p2.bam.bai
    """
}


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