process bcftools_deepvariant_norm {
    label 'process_medium'
    tag "${sample_id}"
    publishDir "${params.deepvariant_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"
    
    input:
        path ref                    // Reference genome
        tuple val(sample_id), path(vcf), path(vcf_tbi)    // VCF from DeepVariant

    output:
        tuple val(sample_id), path("${sample_id}.deepvariant.normalized.vcf.gz"), path("${sample_id}.deepvariant.normalized.vcf.gz.tbi"), emit: vcf_tuple
        
        
    script:
    """
    bcftools norm \\
        -f ${ref} \\
        -m -both \\
        -O z \\
        -o ${sample_id}.deepvariant.normalized.vcf.gz \\
        ${vcf}
    
    bcftools index -t ${sample_id}.deepvariant.normalized.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.deepvariant.normalized.vcf.gz
    touch ${sample_id}.deepvariant.normalized.vcf.gz.tbi
    """
}


process deepvariant_targeted_region {
    label 'high_memory'
    tag "${sample_id}"
    publishDir "${params.deepvariant_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        tuple val(sample_id), path(bam), path(bam_index)          // Input BAM file
        val threads        // Number of shards/threads
        val target_region_string // Target e.g. chr20:1000000-2000000
        
    
    output:
        tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf_tuple
        path "${sample_id}.deepvariant.vcf.gz", emit: vcf
        path "${sample_id}.deepvariant.vcf.gz.tbi", emit: vcf_tbi
        path "${sample_id}.deepvariant.g.vcf.gz", emit: gvcf
        path "${sample_id}.deepvariant.g.vcf.gz.tbi", emit: gvcf_tbi
        
    script:
    def regions_flag = target_region_string ? "--regions \"${target_region_string}\"" : ""
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${sample_id}.deepvariant.vcf.gz \
        --output_gvcf ${sample_id}.deepvariant.g.vcf.gz \
        --num_shards ${threads} \
        ${regions_flag}
        
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz
    touch ${sample_id}.deepvariant.vcf.gz.tbi
    touch ${sample_id}.deepvariant.g.vcf.gz
    touch ${sample_id}.deepvariant.g.vcf.gz.tbi
    """

}


process deepvariant {
    label 'high_memory'
    tag "${sample_id}"
    publishDir "${params.deepvariant_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        tuple val(sample_id), path(bam), path(bam_index)          // Input BAM file
        val threads        // Number of shards/threads
        
    
    output:
        tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf_tuple
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

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz
    touch ${sample_id}.deepvariant.vcf.gz.tbi
    touch ${sample_id}.deepvariant.g.vcf.gz
    touch ${sample_id}.deepvariant.g.vcf.gz.tbi
    """

}

