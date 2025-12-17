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

process PARSE_BCFTOOLS_STATS {
    tag "${sample_id}"
    publishDir "${params.deepvariant_stats_dir}/${sample_id}", mode: 'copy', overwrite: true

    container "indapa/indapa-data-analysis:latest"
    
    input:
        tuple val(sample_id), path(stats_overall)
        tuple val(sample_id3), path(filter_counts)
        path(wgs_sample_tracker)

    output:
        path "${sample_id}.deepvariant.stats.csv", emit: parsed_stats
        path "${sample_id}.deepvariant.filter_counts.txt", emit: parsed_filter_counts

        
    script:
    """
    python /opt/bin/parse-bcftools-stats-deepvariant.py  ${sample_id} ${stats_overall} ${filter_counts} ${sample_id}.deepvariant.stats.csv ${sample_id}.deepvariant.filter_counts.csv ${wgs_sample_tracker}
    """

    stub:
    """
    touch ${sample_id}.deepvariant.stats.csv
    touch ${sample_id}.deepvariant.filter_counts.txt
    """
}

process BCFTOOLS_STATS {
    label 'low_memory'
    tag "${sample_id}"
    publishDir "${params.deepvariant_stats_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    conda "bioconda::bcftools=1.17"
    
    input:
        tuple val(sample_id), path(vcf), path(vcf_tbi)
        
    output:
        tuple val(sample_id), path("${sample_id}.deepvariant.vcf.stats"), emit: stats_overall
        tuple val(sample_id), path("${sample_id}.deepvariant.per_chrom.*.vcf.stats"), emit: stats_per_chrom
        tuple val(sample_id), path("${sample_id}.deepvariant.filter_counts.txt"), emit: filter_counts
        
    script:
    """
    # Generate overall stats
    bcftools view -H ${vcf} | cut -f7 | sort | uniq -c > ${sample_id}.deepvariant.filter_counts.txt
    
    bcftools stats ${vcf} > ${sample_id}.deepvariant.vcf.stats
    
    # Generate per-chromosome stats for standard human chromosomes
    for chr in {1..22} X Y; do
        bcftools stats -r "chr\${chr}" ${vcf} > ${sample_id}.deepvariant.per_chrom.chr\${chr}.vcf.stats 2>/dev/null || \
        bcftools stats -r "\${chr}" ${vcf} > ${sample_id}.deepvariant.per_chrom.chr\${chr}.vcf.stats 2>/dev/null || \
        touch ${sample_id}.deepvariant.per_chrom.chr\${chr}.vcf.stats
    done

    # Per-chromosome breakdown
    # Generate per-chromosome FILTER counts
    for chr in {1..22} X Y; do
        echo "chr\${chr}" >> "${sample_id}.deepvariant.filter_counts.txt"
        bcftools view -H -r "chr\${chr}" ${vcf} \\
        | cut -f7 \\
        | sort \\
        | uniq -c >> "${sample_id}.deepvariant.filter_counts.txt" 2>/dev/null || \
        bcftools view -H -r "\${chr}" ${vcf} \\
        | cut -f7 \\
        | sort \\
        | uniq -c >> "${sample_id}.deepvariant.filter_counts.txt" 2>/dev/null || \
        echo "  0 NO_DATA" >> "${sample_id}.deepvariant.filter_counts.txt"
    done
    
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.stats
    
    # Create stub files for chromosomes 1-22, X, Y
    for chr in {1..22} X Y; do
        touch ${sample_id}.deepvariant.per_chrom.chr\${chr}.vcf.stats
    done

     # Create stub filter counts file with realistic content
    for chr in {1..22} X Y; do
        echo "chr\${chr}" >> "${sample_id}.deepvariant.filter_counts.txt"
        echo "    100 PASS" >> "${sample_id}.deepvariant.filter_counts.txt"
        echo "      5 LowQual" >> "${sample_id}.deepvariant.filter_counts.txt"
    done

    touch ${sample_id}.deepvariant.filter_counts.txt
    """
}

process deepvariant_chr20{
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
        tuple val(sample_id), path("${sample_id}.deepvariant.chr20.vcf.gz"), path("${sample_id}.deepvariant.chr20.vcf.gz.tbi"), emit: vcf_tuple
        path "${sample_id}.deepvariant.chr20.g.vcf.gz", emit: gvcf_chr20
        path "${sample_id}.deepvariant.chr20.g.vcf.gz.tbi", emit: gvcf_tbi_chr20
        
        
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${sample_id}.deepvariant.chr20.vcf.gz \
        --output_gvcf ${sample_id}.deepvariant.chr20.g.vcf.gz \
        --num_shards ${threads} \
        --regions "chr20"
    """

    stub:
    """
    touch ${sample_id}.deepvariant.chr20.vcf.gz
    touch ${sample_id}.deepvariant.chr20.vcf.gz.tbi
    touch ${sample_id}.deepvariant.chr20.g.vcf.gz
    touch ${sample_id}.deepvariant.chr20.g.vcf.gz.tbi
    """
}

process deepvariant_biomarkerSNVs {
   
    tag "${sample_id}"
    publishDir params.deepvariant_downsampled_output_dir, mode: 'copy', overwrite: true
    
    // https://github.com/google/deepvariant/issues/916
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        tuple val(sample_id), path(bam), path(bam_index)          // Input BAM file
        path bed_regions // Regions to analyze for biomarker SNVs
        val threads        // Number of shards/threads
        
    
    output:
        path "${sample_id}.deepvariant.biomarkerSNVs.vcf.gz", emit: vcf_biomarkerSNVs
        path "${sample_id}.deepvariant.biomarkerSNVs.vcf.gz.tbi", emit: vcf_tbi_biomarkerSNVs
        path "${sample_id}.deepvariant.biomarkerSNVs.g.vcf.gz", emit: gvcf
        path "${sample_id}.deepvariant.biomarkerSNVs.g.vcf.gz.tbi", emit: gvcf_tbi
    
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${sample_id}.deepvariant.biomarkerSNVs.vcf.gz \
        --output_gvcf ${sample_id}.deepvariant.biomarkerSNVs.g.vcf.gz \
        --regions ${bed_regions} \
        --call_variants_extra_args="allow_empty_examples=true" \
        --num_shards ${threads} 
    """

    stub:
    """
    touch ${sample_id}.deepvariant.biomarkerSNVs.vcf.gz
    touch ${sample_id}.deepvariant.biomarkerSNVs.vcf.gz.tbi
    touch ${sample_id}.deepvariant.biomarkerSNVs.g.vcf.gz
    touch ${sample_id}.deepvariant.biomarkerSNVs.g.vcf.gz.tbi
    """
}

