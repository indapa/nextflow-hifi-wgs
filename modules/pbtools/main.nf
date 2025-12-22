process pbmm2_align_syt1_region {
    label 'high_memory'
    publishDir "${params.aligned_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pbmm2:1.17.0_build1" // User confirmed samtools is inside

    input:
    path reference
    tuple val(sample_id), path(bam)
    val threads
    val sort_threads

    output:
    // Fixed: Matches the filenames created in the script
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: aligned_bam

    script:
    """
    pbmm2 --version
    
    # 1. Create the BED file for the region of interest
    echo -e "chr12\t78814774\t79502008" > region.bed

    # 2. Pipeline Explanation:
    #    pbmm2:  Aligns to full genome. Writes to STDOUT (Sorting is disabled by pbmm2 here).
    #    view:   Filters stream against BED file. (-u outputs uncompressed bam for speed).
    #    sort:   Sorts the filtered reads (Mandatory because pbmm2 didn't sort).
    
    pbmm2 align \
        --strip \
        -j ${threads} \
        --preset HIFI \
        --sample ${sample_id} \
        --log-level INFO \
        ${reference} \
        ${bam} \
        | samtools view -u -h -L region.bed \
        | samtools sort -@ ${sort_threads} -o ${sample_id}.aligned.bam

    # 3. Index the sorted BAM (Required for the output channel)
    samtools index ${sample_id}.aligned.bam
    """

    stub:
    """
    touch ${sample_id}.aligned.bam
    touch ${sample_id}.aligned.bam.bai
    """
}

process pbmm2_align_region {
    label 'high_memory'
    publishDir "${params.aligned_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pbmm2:1.17.0_build1" // User confirmed samtools is inside

    input:
    path reference
    tuple val(sample_id), path(bam)
    val region_string
    val threads
    val sort_threads

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: aligned_bam

    script:
    """
    pbmm2 --version
    
    # 1. Create the BED file for the region of interest
    region_bed=$(echo "${region_string}" | sed 's/[:-]/\t/g')
    echo -e "${region_bed}" > region.bed

    # 2. Pipeline Explanation:
    #    pbmm2:  Aligns to full genome. Writes to STDOUT (Sorting is disabled by pbmm2 here).
    #    view:   Filters stream against BED file. (-u outputs uncompressed bam for speed).
    #    sort:   Sorts the filtered reads (Mandatory because pbmm2 didn't sort).
    
    pbmm2 align \
        --strip \
        -j ${threads} \
        --preset HIFI \
        --sample ${sample_id} \
        --log-level INFO \
        ${reference} \
        ${bam} \
        | samtools view -u -h -L region.bed \
        | samtools sort -@ ${sort_threads} -o ${sample_id}.aligned.bam

    # 3. Index the sorted BAM (Required for the output channel)
    samtools index ${sample_id}.aligned.bam
    """

    stub:
    """
    touch ${sample_id}.aligned.bam
    touch ${sample_id}.aligned.bam.bai
    """
}


process pbmm2_align {
    label 'high_memory'
    publishDir "${params.aligned_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pbmm2:1.17.0_build1"

    input:
    path reference
    tuple val(sample_id), path(bam)
    val threads
    val sort_threads

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: aligned_bam    //path "${$sample_id}.read_length_and_quality.tsv", emit: bam_rl_qual
    
    script:
    """
    pbmm2 --version
    echo -e "chr12\t78814774\t79502008" > region.bed

    pbmm2 align \
        --sort \
        --strip \
        -j $threads \
        -J $sort_threads \
        --sort-memory 8G \
        --preset HIFI \
        --sample ${sample_id} \
        --log-level INFO \
        $reference \
        $bam | \
        samtools view -b -h -L region.bed -o ${sample_id}.region_only.bam
        

    
    """

    stub:
    """
    # Create mock aligned BAM file
    touch ${sample_id}.aligned.bam
    
    # Create mock BAM index file
    touch ${sample_id}.aligned.bam.bai
    
    """
}

