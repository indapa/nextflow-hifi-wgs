process bcftools_deepvariant_norm {
    label 'process_medium'
    tag { "${sample_id}" }
    publishDir { "${params.deepvariant_output_dir}/${sample_id}" }, mode: 'copy', overwrite: true
    
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
    tag { "${sample_id}" }
    publishDir { "${params.deepvariant_output_dir}/${sample_id}" }, mode: 'copy', overwrite: true
    
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

process deeptrio_wgs {
    tag { "${family_id}" }
    publishDir { "${params.deepvariant_output_dir}/DV_trio/${family_id}" }, mode: 'copy', overwrite: true
    
    // DeepTrio is included in the standard DeepVariant container
    container "google/deepvariant:deeptrio-1.10.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // Reference index
        
        // Input Tuple: Family ID + (ID, BAM, BAI) for Child, Parent1, Parent2
        tuple val(family_id), \
              val(child_id), path(child_bam), path(child_bai), \
              val(p1_id),    path(p1_bam),    path(p1_bai), \
              val(p2_id),    path(p2_bam),    path(p2_bai)

        
        
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
    
    """

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
        --num_shards ${task.cpus} 
        
    """

    stub:
    """
    touch ${child_id}.vcf.gz ${child_id}.vcf.gz.tbi ${child_id}.g.vcf.gz ${child_id}.g.vcf.gz.tbi
    touch ${p1_id}.vcf.gz ${p1_id}.vcf.gz.tbi ${p1_id}.g.vcf.gz ${p1_id}.g.vcf.gz.tbi
    touch ${p2_id}.vcf.gz ${p2_id}.vcf.gz.tbi ${p2_id}.g.vcf.gz ${p2_id}.g.vcf.gz.tbi
    """
}

process deeptrio_wgs_by_chrom {
    tag { "${family_id}_${interval_bed.baseName}" }
    //publishDir { "${params.deepvariant_output_dir}/DV_trio/${family_id}/by_chrom" }, mode: 'copy', overwrite: true

    container "google/deepvariant:deeptrio-1.10.0"
    
    stageInMode 'copy' // avoid virtual file pointers with deep variant

    input:
        path ref
        path ref_index
        tuple val(family_id), \
              val(child_id), path(child_bam), path(child_bai), \
              val(p1_id),    path(p1_bam),    path(p1_bai), \
              val(p2_id),    path(p2_bam),    path(p2_bai), \
              path (interval_bed) // Targeted region for this chromosome

    output:
        
        tuple val(family_id), val(child_id), val(interval_bed.baseName), path("${child_id}.${interval_bed.baseName}.vcf.gz"),   path("${child_id}.${interval_bed.baseName}.vcf.gz.tbi"),   emit: child_vcf
        tuple val(family_id), val(child_id), val(interval_bed.baseName), path("${child_id}.${interval_bed.baseName}.g.vcf.gz"), path("${child_id}.${interval_bed.baseName}.g.vcf.gz.tbi"), emit: child_gvcf
    
        tuple val(family_id), val(p1_id),    val(interval_bed.baseName), path("${p1_id}.${interval_bed.baseName}.vcf.gz"),     path("${p1_id}.${interval_bed.baseName}.vcf.gz.tbi"),     emit: p1_vcf
        tuple val(family_id), val(p1_id),    val(interval_bed.baseName), path("${p1_id}.${interval_bed.baseName}.g.vcf.gz"),   path("${p1_id}.${interval_bed.baseName}.g.vcf.gz.tbi"),   emit: p1_gvcf
    
        tuple val(family_id), val(p2_id),    val(interval_bed.baseName), path("${p2_id}.${interval_bed.baseName}.vcf.gz"),     path("${p2_id}.${interval_bed.baseName}.vcf.gz.tbi"),     emit: p2_vcf
        tuple val(family_id), val(p2_id),    val(interval_bed.baseName), path("${p2_id}.${interval_bed.baseName}.g.vcf.gz"),   path("${p2_id}.${interval_bed.baseName}.g.vcf.gz.tbi"),   emit: p2_gvcf

    script:
    def model_type = task.ext.model_type ?: 'PACBIO'
    """
    /opt/deepvariant/bin/deeptrio/run_deeptrio \
        --model_type ${model_type} \
        --ref ${ref} \
        --reads_child ${child_bam} \
        --reads_parent1 ${p1_bam} \
        --reads_parent2 ${p2_bam} \
        --sample_name_child "${child_id}" \
        --sample_name_parent1 "${p1_id}" \
        --sample_name_parent2 "${p2_id}" \
        --output_vcf_child ${child_id}.${interval_bed.baseName}.vcf.gz \
        --output_vcf_parent1 ${p1_id}.${interval_bed.baseName}.vcf.gz \
        --output_vcf_parent2 ${p2_id}.${interval_bed.baseName}.vcf.gz \
        --output_gvcf_child ${child_id}.${interval_bed.baseName}.g.vcf.gz \
        --output_gvcf_parent1 ${p1_id}.${interval_bed.baseName}.g.vcf.gz \
        --output_gvcf_parent2 ${p2_id}.${interval_bed.baseName}.g.vcf.gz \
        --num_shards ${task.cpus} \
        --regions ${interval_bed}
    """

    stub:
    """
    touch ${child_id}.${interval_bed.baseName}.vcf.gz ${child_id}.${interval_bed.baseName}.vcf.gz.tbi ${child_id}.${interval_bed.baseName}.g.vcf.gz ${child_id}.${interval_bed.baseName}.g.vcf.gz.tbi
    touch ${p1_id}.${interval_bed.baseName}.vcf.gz ${p1_id}.${interval_bed.baseName}.vcf.gz.tbi ${p1_id}.${interval_bed.baseName}.g.vcf.gz ${p1_id}.${interval_bed.baseName}.g.vcf.gz.tbi
    touch ${p2_id}.${interval_bed.baseName}.vcf.gz ${p2_id}.${interval_bed.baseName}.vcf.gz.tbi ${p2_id}.${interval_bed.baseName}.g.vcf.gz ${p2_id}.${interval_bed.baseName}.g.vcf.gz.tbi
    """
}

process concat_chrom_chunks_vcf {
    tag { "${meta[0]} - ${meta[1]} - ${meta[2]} (${meta[3]})" }
    publishDir { "${params.deepvariant_output_dir}/DV_trio/${meta[0]}/by_chrom/${meta[2]}" }, mode: 'copy', overwrite: true

    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"

    input:
    tuple val(meta), path(chunk_files), path(chunk_indices)

    output:
    tuple val(meta), path("${meta[1]}.${meta[2]}.merged.${meta[3]}"), path("${meta[1]}.${meta[2]}.merged.${meta[3]}.tbi"), emit: merged_file

    script:
    def out_file = "${meta[1]}.${meta[2]}.merged.${meta[3]}"
    def sample_id = meta[1] // e.g. HG003
    
    // 1. Sort files natively in Groovy to guarantee genomic order
    def sorted_chunks = chunk_files.sort { a, b -> a.name <=> b.name }
    """
    rm -f clean_file_list.txt
    echo "${sample_id}" > correct_sample.txt
    mkdir -p sanitized/

    # 1. Loop through the exact array Nextflow gave us, completely ignoring file extensions
    for f in ${sorted_chunks.join(' ')}; do
        bcftools reheader -s correct_sample.txt "\$f" -o sanitized/"\$f"
        bcftools index -t sanitized/"\$f"
    done

    # 2. Re-read the exact same sorted array names out of the sanitized directory 
    # to guarantee they stay in genomic order without using standard ls wildcards
    for f in ${sorted_chunks.join(' ')}; do
        echo "sanitized/\$f" >> clean_file_list.txt
    done

    # 3. Concatenate and index safely
    bcftools concat -a -f clean_file_list.txt -O z -o ${out_file}
    bcftools index -t ${out_file}
    """

    stub:
    """
    touch ${meta[1]}.${meta[2]}.merged.${meta[3]}
    touch ${meta[1]}.${meta[2]}.merged.${meta[3]}.tbi
    """
}

process concat_wgs_vcf {
    tag { "${sample_id} Final WGS (${ext})" }
    publishDir { "${params.deepvariant_output_dir}/DV_trio/${family_id}" }, mode: 'copy', overwrite: true

    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"

    input:
        tuple val(family_id), val(sample_id), path(chrom_files), path(tbis)
        val ext // Pass either "vcf.gz" or "g.vcf.gz"

    output:
        tuple val(family_id), val(sample_id), path("${sample_id}.${ext}"), path("${sample_id}.${ext}.tbi"), emit: merged

    script:
    // 1. Create a map of staged filenames to make lookups fast and simple
    // e.g. ["HG003.chr15.merged.g.vcf.gz": "HG003.chr15.merged.g.vcf.gz"]
    def staged_map = chrom_files.collectEntries { f -> [f.name, f] }

    // 2. Define the desired order of chromosomes
    def chrom_order = (1..22).collect { num -> "chr${num}" } + ['chrX', 'chrY']

    // 3. Build a list of files that exist in our input, ordered by chromosome
    def ordered_files = []
    chrom_order.each { chrom ->
        // Dynamically find the staged file for this chromosome
        def matched_file = staged_map.keySet().find { name -> 
            name.contains(".${chrom}.") && name.endsWith(".${ext}") 
        }
        if (matched_file) {
            ordered_files << matched_file
        }
    }
    """
    rm -f wgs_list.txt
    
    # Write the perfectly ordered staged filenames straight to the file list
    ${ordered_files.collect { file -> "echo '$file' >> wgs_list.txt" }.join('\n')}
    # Safe to use --naive here because whole chromosomes do not overlap coordinates
    bcftools concat \
        --naive \
        -f wgs_list.txt \
        -O z \
        -o ${sample_id}.${ext}

    bcftools index -t ${sample_id}.${ext}
    """
    
    stub:
    """
    touch ${sample_id}.${ext}
    touch ${sample_id}.${ext}.tbi
    """
}



process glnexus_trio_merge {
    tag { "${family_id}" }
    publishDir { "${params.deepvariant_output_dir}/DV_trio/${family_id}" }, mode: 'copy', overwrite: true
    
    container "quay.io/mlin/glnexus:v1.2.7"

    input:
    // Streamlined: Only family ID and the file tracks are needed here
    tuple val(family_id), \
          path(child_gvcf), path(child_tbi), \
          path(p1_gvcf),    path(p1_tbi), \
          path(p2_gvcf),    path(p2_tbi)

    output:
    tuple val(family_id), path("${family_id}.joint.vcf.gz"), path("${family_id}.joint.vcf.gz.tbi"), emit: joint_vcf

    script:
    """
    glnexus_cli \
        --config DeepVariant_unfiltered \
        --threads ${task.cpus} \
        ${child_gvcf} \
        ${p1_gvcf} \
        ${p2_gvcf} \
        | bcftools view - \
        | bgzip -c > ${family_id}.joint.vcf.gz

    tabix -p vcf ${family_id}.joint.vcf.gz
    """

    stub:
    """
    touch ${family_id}.joint.vcf.gz
    touch ${family_id}.joint.vcf.gz.tbi
    """
}


process deepvariant_wgs {
    
    tag { "${sample_id}" }
    publishDir { "${params.deepvariant_output_dir}/${sample_id}" }, mode: 'copy', overwrite: true
    
    container "google/deepvariant:1.10.0"

    input:
    path ref                                                          // Reference genome FASTA
    path ref_index                                                    // Reference index (.fai)
    tuple val(sample_id), path(bam), path(bam_index)                 // Aligned BAM + index

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf_tuple
    tuple val(sample_id), path("${sample_id}.deepvariant.g.vcf.gz"), path("${sample_id}.deepvariant.g.vcf.gz.tbi"), emit: gvcf_tuple
    

    script:
    //def args = task.ext.args ?: ''
    def model_type = task.ext.model_type ?: 'PACBIO'
    
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type ${model_type} \\
        --ref ${ref} \\
        --reads ${bam} \\
        --output_vcf ${sample_id}.deepvariant.vcf.gz \\
        --output_gvcf ${sample_id}.deepvariant.g.vcf.gz \\
        --num_shards ${task.cpus} 
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz
    touch ${sample_id}.deepvariant.vcf.gz.tbi
    touch ${sample_id}.deepvariant.g.vcf.gz
    touch ${sample_id}.deepvariant.g.vcf.gz.tbi
    """
}




