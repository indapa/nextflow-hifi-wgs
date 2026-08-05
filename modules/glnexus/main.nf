process glnexus_trio_by_chrom {
    tag { "${family_id}_${chrom}" }
    
    container "quay.io/mlin/glnexus:v1.2.7"

    input:
    tuple val(family_id), \
          path(child_gvcf), path(child_tbi), \
          path(p1_gvcf),    path(p1_tbi), \
          path(p2_gvcf),    path(p2_tbi), \
          val(chrom)
    path glnexus_bed

    output:
    tuple val(family_id), path("${family_id}.${chrom}.joint.vcf.gz"), path("${family_id}.${chrom}.joint.vcf.gz.tbi"), emit: joint_vcf

    script:
    def avail_mem = (task.memory.toGiga() * 0.85).toInteger()
    """
    set -euo pipefail

    # 1. Filter BED to target chromosome (escaped \$1 for Nextflow string interpolation)
    awk -v c="${chrom}" '\$1 == c' ${glnexus_bed} > region.bed

    # 2. Run GLnexus directly to binary BCF with explicit memory allocation
    glnexus_cli \
        --config DeepVariant_unfiltered \
        --threads ${task.cpus} \
        --mem-gbytes ${avail_mem} \
        --dir "GLnexus_${chrom}.DB" \
        --bed region.bed \
        ${child_gvcf} \
        ${p1_gvcf} \
        ${p2_gvcf} \
        > ${family_id}.${chrom}.joint.bcf

    # 3. Separate conversion step: BCF -> bgzipped VCF
    bcftools view \
        --threads ${task.cpus} \
        -Oz \
        -o ${family_id}.${chrom}.joint.vcf.gz \
        ${family_id}.${chrom}.joint.bcf

    # 4. Separate indexing step
    bcftools index -t --threads ${task.cpus} ${family_id}.${chrom}.joint.vcf.gz

    # 5. Cleanup temporary RocksDB directory and uncompressed BCF
    rm -rf "GLnexus_${chrom}.DB" "${family_id}.${chrom}.joint.bcf" region.bed
    """

    stub:
    """
    touch ${family_id}.${chrom}.joint.vcf.gz
    touch ${family_id}.${chrom}.joint.vcf.gz.tbi
    """
}

process concat_glnexus_vcf {
    tag { "${family_id}" }
    publishDir { "${params.deepvariant_output_dir}/DV_trio/${family_id}" }, mode: 'copy', overwrite: true

    container "quay.io/biocontainers/bcftools:1.21--h8b25389_0"

    input:
    tuple val(family_id), path(vcfs), path(tbis)

    output:
    tuple val(family_id), path("${family_id}.joint.vcf.gz"), path("${family_id}.joint.vcf.gz.tbi"), emit: merged

    script:
    """
    bcftools concat \
        --allow-overlaps \
        --output-type z \
        --output ${family_id}.joint.vcf.gz \
        ${vcfs}

    tabix -p vcf ${family_id}.joint.vcf.gz
    """

    stub:
    """
    touch ${family_id}.joint.vcf.gz
    touch ${family_id}.joint.vcf.gz.tbi
    """
}