process deepvariant {
    label 'high_memory'
    tag "${sample_id}"
    publishDir params.deepvariant_output_dir, mode: 'copy'
    
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        tuple val(sample_id), path(bam), path(bam_index)          // Input BAM file
        val threads        // Number of shards/threads
        
    
    output:
        path "${sample_id}.deepvariant.vcf.gz", emit: vcf
        path "${sample_id}.deepvariant.vcf.gz.tbi", emit: vcf_tbi
        path "${sample_id}.deepvariant.g.vcf.gz", emit: gvcf
        path "${sample_id}.deepvariant.g.vcf.gz.tbi", emit: gvcf_tbi
        
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${sample_id}.deepvariant.vcf.gz \
        --output_gvcf ${sample_id}.deepvariant.g.vcf.gz \
        --num_shards ${threads} 
    """
}

process deepvariant_chr20 {
    label 'high_memory'
    tag "${sample_id}"
    publishDir params.deepvariant_output_dir, mode: 'copy'
    
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        tuple val(sample_id), path(bam), path(bam_index)          // Input BAM file
        val threads        // Number of shards/threads
        
    
    output:
        path "${sample_id}.deepvariant.vcf.gz", emit: vcf_chr20
        path "${sample_id}.deepvariant.vcf.gz.tbi", emit: vcf_tbi_chr20
        path "${sample_id}.deepvariant.g.vcf.gz", emit: gvcf_chr20
        path "${sample_id}.deepvariant.g.vcf.gz.tbi", emit: gvcf_tbi_chr20
        
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${sample_id}.deepvariant.vcf.gz \
        --output_gvcf ${sample_id}.deepvariant.g.vcf.gz \
        --num_shards ${threads} \
        --regions "chr20"
    """
}