

process pbmm2_align {
    publishDir "${params.aligned_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pbmm2:1.17.0_build1"

    input:
    path reference
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: aligned_bam
    
    script:
    // Calculate sort threads as a fraction of total CPUs (common pattern)
    def sort_threads = 4
    def sort_memory = task.ext.args ?: '-m 8G'
    
    """
    # Set TMPDIR to use fast local storage if available
    export TMPDIR=\${TMPDIR:-/tmp}
    
    pbmm2 --version
   
    pbmm2 align \\
        --sort \\
        -j ${task.cpus} \\
        -J ${sort_threads} \\
        --preset HIFI \\
        --sample ${sample_id} \\
        --log-level DEBUG \\
        --unmapped \\
        --bam-index BAI \\
        $reference \\
        $bam \\
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




process cpg_methylation_calling {
    /*
     * 5mC CpG methylation calling from HiFi-aligned BAM using pb-CpG-tools.
     *
     * Assumes a non-haplotagged BAM, so only combined (all-reads) outputs are
     * produced. For haplotagged BAMs, hap1/hap2 bed/bw files would also be
     * emitted — add those to the output block if haplotype tracks are needed.
     *
     * Recommended defaults (model pileup + denovo modsites) are used.
     * --min-mapq 20 and --min-coverage 10 are applied as requested.
     */

   
    tag "$sample_id"
    publishDir "${((params.cpg_output_dir ?: '').toString().trim()) ? "${params.cpg_output_dir}/${sample_id}" : error("Missing required parameter: --cpg_output_dir. Set params.cpg_output_dir in nextflow.config or pass --cpg_output_dir on the command line.")}", mode: 'copy', overwrite: true

    container "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"

    input:
    tuple val(sample_id), path(bam), path(bam_index)  // Aligned (non-haplotagged) BAM + BAI
    path ref                                           // Reference genome FASTA (required for CRAM; staged for BAM)
    path ref_index                                     // Reference FASTA index (.fai)

    output:
    tuple val(sample_id),
          path("${sample_id}.combined.bed.gz"),
          path("${sample_id}.combined.bed.gz.tbi"), emit: combined_bed
    tuple val(sample_id),
         path("${sample_id}.hap1.bed.gz"), 
         path("${sample_id}.hap1.bed.gz.tbi"), emit: hap1_bed
    tuple val(sample_id),
            path("${sample_id}.hap2.bed.gz"),
            path("${sample_id}.hap2.bed.gz.tbi"), emit: hap2_bed
    tuple val(sample_id),
            path("${sample_id}.combined.bw"), emit: combined_bw
    tuple val(sample_id),
            path("${sample_id}.hap1.bw"), emit: hap1_bw
    tuple val(sample_id),
            path("${sample_id}.hap2.bw"), emit: hap2_bw

    script:
    """
    aligned_bam_to_cpg_scores \\
        --bam ${bam} \\
        --output-prefix ${sample_id} \\
        --ref ${ref} \\
        --pileup-mode model \\
        --modsites-mode denovo \\
        --min-mapq 20 \\
        --min-coverage 10 \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}.combined.bed.gz
    touch ${sample_id}.combined.bed.gz.tbi
    touch ${sample_id}.combined.bw
    touch ${sample_id}.hap1.bed.gz
    touch ${sample_id}.hap1.bed.gz.tbi
    touch ${sample_id}.hap1.bw
    touch ${sample_id}.hap2.bed.gz
    touch ${sample_id}.hap2.bed.gz.tbi
    touch ${sample_id}.hap2.bw  
    """
}


process sawfish_discover {
    
    tag "$sample_id"
    container "quay.io/pacbio/sawfish:2.2.1_build1"
    
    input:
    tuple val(sample_id), path(bam), path(bam_index), path(expected_bed)
    path excluded_bed
    path reference
    path reference_index
   
       
    output:
    path "${sample_id}_sawfish_discover", type:'dir', emit: discover_dir
        
    script:
    """
    set -euo pipefail

    sawfish --version
    
    sawfish discover \\
        --threads ${task.cpus} \\
        --bam ${bam} \\
        --ref ${reference} \\
        --expected-cn ${expected_bed} \\
        --cnv-excluded-regions ${excluded_bed} \\
        --output-dir ${sample_id}_sawfish_discover
    """

    stub:
    """
    mkdir -p ${sample_id}_sawfish_discover
    """
}

process sawfish_joint_call {
    
    publishDir "${params.sawfish_output_dir}/joint_call", mode: 'copy'
    container "quay.io/pacbio/sawfish:2.2.1_build1"

    input:
    // "path" here will accept a List of paths because of .collect()
    path all_discover_dirs 
   

    output:
    path "sawfish_joint_call_dir", type:'dir', emit: joint_dir

    script:
    // Groovy magic: transform the list [dir1, dir2, dir3] 
    // into the string "--sample dir1 --sample dir2 --sample dir3"
    def sample_args = all_discover_dirs.collect { "--sample $it" }.join(' ')

    """
    set -euo pipefail

    echo "Running joint call on ${all_discover_dirs.size()} samples..."

    sawfish joint-call \
        --threads ${task.cpus} \
        ${sample_args} \
        --output-dir sawfish_joint_call_dir
    """

    stub:
    """
    mkdir -p sawfish_joint_call_dir
    """
}




process hiphase_small_variants {
    /* hiphase small variants only */
    tag "$sample_id"
    
    

    
    publishDir "${params.hiphase_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    container "quay.io/pacbio/hiphase:1.5.0_build1"

    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi), path(pbmm2_bam), path(pbmm2_bai)
    path reference
    path reference_index

    output:
    tuple val(sample_id), path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi"), emit: phased_vcf 
    tuple val(sample_id), path("*.stats.csv"), path("*.blocks.tsv"), path("*.summary.tsv"), emit: stats
    tuple val(sample_id), path("*.haplotagged.bam"), path("*.haplotagged.bam.bai"), emit: haplotagged_bam

    script:
    def basename = vcf.simpleName
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

    stub:
    def basename = vcf.simpleName
    """
    touch ${basename}.phased.vcf.gz
    touch ${basename}.phased.vcf.gz.tbi
    touch ${basename}.haplotagged.bam
    touch ${basename}.haplotagged.bam.bai
    touch ${basename}.stats.csv
    touch ${basename}.blocks.tsv
    touch ${basename}.summary.tsv
    """
}

