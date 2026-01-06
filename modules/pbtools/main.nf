

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
    region_bed=\$(echo "${region_string}" | sed 's/[:-]/\t/g')
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
   
    pbmm2 align \
        --sort \
        -j $threads \
        -J $sort_threads \
        --preset HIFI \
        --sample ${sample_id} \
        --log-level INFO \
        --unmapped \
        --bam-index BAI \
        $reference \
        $bam \
        ${sample_id}.aligned.bam
    
    """

    stub:
    """
    # Create mock aligned BAM file
    touch ${sample_id}.aligned.bam
    
    # Create mock BAM index file
    touch ${sample_id}.aligned.bam.bai
    
    """
}

process hiphase_small_variants {
    /* hiphase small variants only */

    label 'high_memory'
    publishDir "${params.hiphase_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    container "quay.io/pacbio/hiphase:1.5.0_build1"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi), path(pbmm2_bam), path(pbmm2_bai)
    path (reference)
    path (reference_index)

    output:
    tuple val(sample_id), path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi"), emit: phased_vcf 
    tuple val(sample_id), path("*.stats.csv"), path("*.blocks.tsv"), path("*.summary.tsv"), emit: stats
    tuple val(sample_id), path("*.haplotagged.bam"), path("*.haplotagged.bam.bai"), emit: haplotagged_bam

    script:
    def basename = vcf.simpleName  // Gets name without extension
    """
    hiphase --version
    hiphase --reference ${reference} \
        --bam ${pbmm2_bam} \
        --vcf ${vcf} \
        --output-vcf ${basename}.phased.vcf.gz \
        --threads ${task.cpus} \
        --min-mapq 20 \
        --disable-global-realignment \
        --ignore-read-groups \
        --output-bam ${basename}.haplotagged.bam \
        --stats-file ${basename}.stats.csv \
        --blocks-file ${basename}.blocks.tsv \
        --summary-file ${basename}.summary.tsv

    bcftools index -f --tbi ${basename}.phased.vcf.gz
    samtools index ${basename}.haplotagged.bam

    """
}

