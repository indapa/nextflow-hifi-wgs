
process SPLIT_BAM {
    tag "${sample_id}"
    publishDir "${params.split_bams_dir}/${sample_id}", mode: 'copy', overwrite: true
    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path "chunk_*.bam", emit: chunks
    path "split_stats.txt", emit: stats

    script:
    """
    # Get total number of reads
    TOTAL_READS=\$(samtools view -c ${bam_file})
    echo "Total reads in input BAM: \$TOTAL_READS" > split_stats.txt
    
    # Calculate number of chunks needed
    CHUNK_SIZE=${params.reads_per_chunk}
    NUM_CHUNKS=\$(((\$TOTAL_READS + \$CHUNK_SIZE - 1) / \$CHUNK_SIZE))
    echo "Splitting into \$NUM_CHUNKS chunks of ~\$CHUNK_SIZE reads each" >> split_stats.txt
    
    # Extract header
    samtools view -H ${bam_file} > header.sam
    
    # Split the BAM file
    samtools view ${bam_file} | split -l \$CHUNK_SIZE - temp_chunk_
    
    # Convert each chunk back to BAM format with header
    chunk_num=1
    for chunk in temp_chunk_*; do
        chunk_name=\$(printf "chunk_%04d.bam" \$chunk_num)
        echo "Processing \$chunk_name..." >> split_stats.txt
        
        # Count reads in this chunk
        reads_in_chunk=\$(wc -l < \$chunk)
        echo "  Reads in \$chunk_name: \$reads_in_chunk" >> split_stats.txt
        
        # Create BAM file with header
        cat header.sam \$chunk | samtools view -bS - > \$chunk_name
        
        # Index the BAM file
        samtools index \$chunk_name
        
        chunk_num=\$((chunk_num + 1))
    done
    
    # Clean up temporary files
    rm temp_chunk_* header.sam
    
    echo "BAM splitting completed successfully" >> split_stats.txt
    """

    stub:
    """
    touch chunk_0001.bam
    touch split_stats.txt
    """
}


process bamToFastq {
    tag "${sample_id}"
    publishDir params.unaligned_fastqs, mode: 'copy'
    label 'low_memory'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), val(file_type), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.unaligned.fastq.gz"), emit: fastq

    script:
    """
    samtools view -f 4 -h "${bam}" | \\
    samtools fastq -0 /dev/stdout  - | gzip > "${sample_id}.unaligned.fastq.gz"
    """

    stub:
    """
    touch ${sample_id}.unaligned.fastq.gz
    """

}

process downsample_unaligned_bam {
    
    /* downsample bam file to specified fraction of original reads */

    tag "${sample_id}"

    publishDir params.downsampled_output_dir, mode: 'copy'
    publishDir "${params.downsampled_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), val(coverage), val(fraction)

    output:
    tuple val(sample_id), path("${sample_id}.downsampled.${coverage}.bam"), emit: downsampled_bam

    script:
    """
    samtools view -s ${fraction} -b ${bam} > ${sample_id}.downsampled.${coverage}.bam
    samtools index ${sample_id}.downsampled.${coverage}.bam
    """

    stub:
    """
    touch ${sample_id}.downsampled.${coverage}.bam
    touch ${sample_id}.downsampled.${coverage}.bam.bai
    """
}

process downsample_bam {
    
    /* downsample bam file to specified fraction of original reads */

    tag "${sample_id}"

    publishDir params.downsampled_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai), val(coverage), val(fraction)

    output:
    tuple val(sample_id), path("${sample_id}.downsampled.${coverage}.bam"), path("${sample_id}.downsampled.${coverage}.bam.bai"), emit: downsampled_bam

    script:
    """
    samtools view -s ${fraction} -b ${bam} > ${sample_id}.downsampled.${coverage}.bam
    samtools index ${sample_id}.downsampled.${coverage}.bam
    """

    stub:
    """
    touch ${sample_id}.downsampled.${coverage}.bam
    touch ${sample_id}.downsampled.${coverage}.bam.bai
    """
}
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

/*
process extract_unmapped_reads_fasta {
    
    // extract unmapped reads from bam file and convert to fasta 

    tag "${sample_id}"
    label 'low_memory'

    publishDir params.unmapped_reads_output_dir, mode: 'copy', overwrite: true

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.unmapped_reads.fa.gz", emit: unmapped_fasta

    script:
    """

    samtools view -f 4 "${bam}" | awk '{print ">" $1 "\n" $10}' | gzip > "${sample_id}.unmapped_reads.fa.gz"

    """

    stub:
    """
    touch ${sample_id}.unmapped_reads.fa.gz
    """
}
*/
// calculate read lengths for the bam file


process read_lengths {

    /* calculate read lengths for the bam file */

    tag "${bam}"
    label 'low_memory'

    publishDir params.bam_stats_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple  path(bam), path(bai)

    output:
    path "${bam}.read_lengths.txt", emit: read_lengths

    script:
    """

    samtools view ${bam} | awk '{print length(\$10)}' | sort -n | uniq -c > ${bam}.read_lengths.txt

    """
}
process filter_mapq_primary_alns {

    label 'low_memory'
    tag "${sample_id}"

    publishDir params.aligned_filtered_output_dir, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai"), emit: filtered_bam_bai

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    script:
    """
    samtools view -b -q 60 -F 0x900 "${bam}" > "${sample_id}.filtered.bam"
    samtools index "${sample_id}.filtered.bam"
    """

    stub:
    """
    touch "${sample_id}.filtered.bam"
    touch "${sample_id}.filtered.bam.bai"
    """
}
// generate stats for the bam file




process bam_stats {

    /* generate stats for the bam file */
    label 'low_memory'

    tag "${sample_id}"

    publishDir "${params.bam_stats_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam)

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



process downsample_bam_stats {

    /* generate stats for the bam file */
    label 'low_memory'

    tag "${bam}"

    publishDir params.downsampled_bam_stats_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple path(bam), path(bai)

    output:
    path "${bam}.stats", emit: stats

    script:
    """

    samtools stats --threads 4 ${bam} > ${bam}.stats

    """

    stub:
    """

    touch ${bam}.stats

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

// process to conver aligned BAM to unaligned BAM

process unaligned_bam {

    /* convert aligned BAM to unaligned BAM */

    tag "${sample_id}"

    publishDir params.unaligned_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.unaligned.bam", emit: unaligned_bam

    script:
    """

    samtools view -h input.bam | \
   awk '{if (\$0 ~ /^@/) {print} else {\$2="4"; \$3="*"; \$4="0"; \$5="0"; \$6="*"; \$7="*"; \$8="0"; \$9="0"; print}}' | \
    samtools view -b > ${sample_id}.unaligned.bam

    """
}