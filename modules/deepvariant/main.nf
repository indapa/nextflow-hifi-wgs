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

process deeptrio_targeted_region {
    label 'high_memory'
    tag "${family_id}"
    publishDir "${params.deepvariant_output_dir}/DV_trio/${family_id}", mode: 'copy', overwrite: true
    
    // DeepTrio is included in the standard DeepVariant container
    container "google/deepvariant:deeptrio-1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // Reference index
        
        // Input Tuple: Family ID + (ID, BAM, BAI) for Child, Parent1, Parent2
        tuple val(family_id), \
              val(child_id), path(child_bam), path(child_bai), \
              val(p1_id),    path(p1_bam),    path(p1_bai), \
              val(p2_id),    path(p2_bam),    path(p2_bai)

        val threads                 // Number of shards
        val target_region_string    // Target e.g. "chr20:1000000-2000000"
        
    output:
       
        // Child Outputs (Start with family_id)
        tuple val(family_id), val(child_id), path("${child_id}.vcf.gz"), path("${child_id}.vcf.gz.tbi"), emit: child_vcf
        tuple val(family_id), val(child_id), path("${child_id}.g.vcf.gz"), path("${child_id}.g.vcf.gz.tbi"), emit: child_gvcf
        
        // Parent 1 Outputs (Start with family_id)
        tuple val(family_id), val(p1_id), path("${p1_id}.vcf.gz"), path("${p1_id}.vcf.gz.tbi"), emit: p1_vcf
        tuple val(family_id), val(p1_id), path("${p1_id}.g.vcf.gz"), path("${p1_id}.g.vcf.gz.tbi"), emit: p1_gvcf

        // Parent 2 Outputs (Start with family_id)
        tuple val(family_id), val(p2_id), path("${p2_id}.vcf.gz"), path("${p2_id}.vcf.gz.tbi"), emit: p2_vcf
        tuple val(family_id), val(p2_id), path("${p2_id}.g.vcf.gz"), path("${p2_id}.g.vcf.gz.tbi"), emit: p2_gvcf

    script:
    def regions_flag = target_region_string ? "--regions \"${target_region_string}\"" : ""
    
    """
    mkdir -p intermediate_results_dir

    /opt/deepvariant/bin/deeptrio/run_deeptrio \
        --model_type PACBIO \
        --ref ${ref} \
        --reads_child ${child_bam} \
        --reads_parent1 ${p1_bam} \
        --reads_parent2 ${p2_bam} \
        --sample_name_child "${child_id}" \
        --sample_name_parent1 "${p1_id}" \
        --sample_name_parent2 "${p2_id}" \
        --output_vcf_child ${child_id}.vcf.gz \
        --output_vcf_parent1 ${p1_id}.vcf.gz \
        --output_vcf_parent2 ${p2_id}.vcf.gz \
        --output_gvcf_child ${child_id}.g.vcf.gz \
        --output_gvcf_parent1 ${p1_id}.g.vcf.gz \
        --output_gvcf_parent2 ${p2_id}.g.vcf.gz \
        --num_shards ${threads} \
        ${regions_flag}
    """

    stub:
    """
    touch ${child_id}.vcf.gz ${child_id}.vcf.gz.tbi ${child_id}.g.vcf.gz ${child_id}.g.vcf.gz.tbi
    touch ${p1_id}.vcf.gz ${p1_id}.vcf.gz.tbi ${p1_id}.g.vcf.gz ${p1_id}.g.vcf.gz.tbi
    touch ${p2_id}.vcf.gz ${p2_id}.vcf.gz.tbi ${p2_id}.g.vcf.gz ${p2_id}.g.vcf.gz.tbi
    """
}

process glnexus_trio_merge {
    tag "${family_id}"
    label 'high_memory' // GLnexus can use significant RAM
    publishDir "${params.deepvariant_output_dir}/DV_trio/${family_id}"
    
    // Using the container you requested (contains glnexus_cli, bcftools, bgzip)
    container "quay.io/mlin/glnexus:v1.2.7"

    input:
    // We expect a tuple containing the Family ID plus all 3 gVCFs + indices
    tuple val(family_id), \
          val(child_id), path(child_gvcf), path(child_tbi), \
          val(p1_id),    path(p1_gvcf),    path(p1_tbi), \
          val(p2_id),    path(p2_gvcf),    path(p2_tbi)

    output:
    tuple val(family_id), path("${family_id}.joint.vcf.gz"), path("${family_id}.joint.vcf.gz.tbi"), emit: joint_vcf

    script:
    """
    # 1. Run GLnexus
    # We pipe stdout -> bcftools -> bgzip directly without switching containers
    # --threads: Uses Nextflow allocated CPUs
    
    glnexus_cli \
        --config DeepVariant_unfiltered \
        --threads ${task.cpus} \
        ${child_gvcf} \
        ${p1_gvcf} \
        ${p2_gvcf} \
        | bcftools view - \
        | bgzip -c > ${family_id}.joint.vcf.gz

    # 2. Index the resulting VCF
    tabix -p vcf ${family_id}.joint.vcf.gz
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

