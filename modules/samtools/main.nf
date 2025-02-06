process index_bam {

    /* index bam file */

    tag "${sample_id}"

    publishDir params.bam_index_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${bam}.bai", emit: bai

    script:
    """

    samtools index ${bam}

    """
}

// process to extract out chr20 from the bam file


process extract_chr20 {

    /* extract chr20 from the bam file; helpful for test data */

    tag "${sample_id}"

    publishDir params.aligned_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.chr20.aligned.bam", emit: chr20_bam
    path "${sample_id}.chr20.aligned.bam.bai", emit: chr20_bai

    script:
    """

    samtools view -b ${bam} chr20 > ${sample_id}.chr20.aligned.bam
    samtools index ${sample_id}.chr20.aligned.bam

    """
}

process subsample_bam {
    
    /* subsample unaligned  bam file for the first 25k reads
       Useful for making test dataset  
     */

    tag "${sample_id}"

    publishDir params.test_data_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam)

    output:
    
    path "${bam.simpleName}.subsampled.bam", emit: subsampled_bam
    

    script:
    def basename = bam.simpleName  // Gets name without extension
    """

    samtools view -h ${bam} | head -n 250000  | samtools view -bS - > ${basename}.subsampled.bam
    
    """
}