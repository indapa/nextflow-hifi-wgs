process whatshap_trio_phase {
    tag { "${family_id}" }
    publishDir { "${params.deepvariant_output_dir}/DV_trio/${family_id}" }, mode: 'copy', overwrite: true
    container "indapa/whatshap-tabix"

    input:
    // 1. Reference Files
    path reference
    path reference_index
    
    // 2. The JOINED input (VCF + BAMs together)
    // Structure: [family_id, vcf, vcf_tbi, child_id, child_bam...]
    tuple val(family_id), path(vcf), path(vcf_tbi), \
          val(child_id), path(child_bam), path(child_bai), \
          val(p1_id),    path(p1_bam),    path(p1_bai), \
          val(p2_id),    path(p2_bam),    path(p2_bai)

    output:
    tuple val(family_id), path("${family_id}.trio_phased.vcf.gz"), emit: phased_vcf
    path "${child_id}.haplotagged.bam", emit: haplotagged_bam
    path "${family_id}.ped", emit: ped_file
    path "${family_id}.blocks.tsv", emit: block_stats 
    path "${family_id}.gtf", emit: block_gtf  

    script:
    """
    # 1. Create PED file
    
    echo "${family_id} ${child_id} ${p1_id} ${p2_id} 0 0" >> ${family_id}.ped

    # 2. Run Whatshap Phase
    whatshap phase \
        --ped ${family_id}.ped \
        --reference ${reference} \
        --output ${family_id}.trio_phased.vcf.gz \
        ${vcf} \
        ${child_bam} ${p1_bam} ${p2_bam}

    # 3. Index

    whatshap stats \
        --block-list=${family_id}.blocks.tsv \
        ${family_id}.trio_phased.vcf.gz
    

    # GTF 
    whatshap stats --gtf=${family_id}.gtf ${family_id}.trio_phased.vcf.gz

    #haplotag
    tabix -p vcf ${family_id}.trio_phased.vcf.gz
    whatshap haplotag -o ${child_id}.haplotagged.bam --reference ${reference} ${family_id}.trio_phased.vcf.gz ${child_bam}
    """

    stub:
    """
    touch ${family_id}.trio_phased.vcf.gz
    touch ${family_id}.trio_phased.vcf.gz.tbi
    touch ${child_id}.haplotagged.bam
    touch ${child_id}.haplotagged.bam.bai
    touch ${family_id}.ped
    touch ${family_id}.blocks.tsv
    touch ${family_id}.gtf
    """
}


process WHATSHAP_PHASE_CHROM {
    tag "${family_id}:${chrom}"
    label 'process_medium'
    container "indapa/whatshap-tabix"

    input:
    tuple val(family_id), val(chrom), path(vcf), path(vcf_tbi), \
          val(child_id), path(child_bam), path(child_bai), \
          val(p1_id),    path(p1_bam),    path(p1_bai), \
          val(p2_id),    path(p2_bam),    path(p2_bai)
    path reference
    path reference_index

    output:
    tuple val(family_id), path("${family_id}.${chrom}.phased.vcf.gz"), emit: phased_vcf

    script:
    """
    echo "${family_id} ${child_id} ${p1_id} ${p2_id} 0 0" > ${family_id}.ped

    whatshap phase \
        --ped ${family_id}.ped \
        --reference ${reference} \
        --chromosome ${chrom} \
        --output ${family_id}.${chrom}.phased.vcf.gz \
        ${vcf} \
        ${child_bam} ${p1_bam} ${p2_bam}
    """

   
    stub:
    """
    touch ${family_id}.${chrom}.phased.vcf.gz
    """ 
}

process CONCAT_PHASED_VCFS {
    tag "${family_id}"
    label 'process_low'
    container "quay.io/biocontainers/bcftools:1.17--haef29d1_0"

    input:
    tuple val(family_id), path(phased_vcfs)

    output:
    tuple val(family_id), path("${family_id}.trio_phased.vcf.gz"), path("${family_id}.trio_phased.vcf.gz.tbi"), emit: merged_vcf

    script:
    """
    ls *.phased.vcf.gz | sort -V > vcf_list.txt
    bcftools concat --file-list vcf_list.txt -Oz -o ${family_id}.trio_phased.vcf.gz
    bcftools index -t ${family_id}.trio_phased.vcf.gz
    """

    stub:
    """
    touch ${family_id}.trio_phased.vcf.gz
    touch ${family_id}.trio_phased.vcf.gz.tbi
    """
}

// modules/local/whatshap_stats_haplotag/main.nf
process WHATSHAP_STATS_HAPLOTAG {
    tag "${family_id}"
    label 'process_medium'
    publishDir "${params.deepvariant_output_dir}/DV_trio/${family_id}", mode: 'copy', overwrite: true
    container "indapa/whatshap-tabix"

    input:
    tuple val(family_id), path(vcf), path(vcf_tbi), \
          val(child_id), path(child_bam), path(child_bai), \
          val(p1_id), val(p2_id)
    path reference
    path reference_index

    output:
   
    path "${child_id}.haplotagged.bam", emit: haplotagged_bam
    path "${family_id}.ped", emit: ped_file
    path "${family_id}.blocks.tsv", emit: block_stats
    path "${family_id}.gtf", emit: block_gtf

    script:
    """
    echo "${family_id} ${child_id} ${p1_id} ${p2_id} 0 0" > ${family_id}.ped

    whatshap stats \
        --block-list=${family_id}.blocks.tsv \
        ${vcf}

    whatshap stats --gtf=${family_id}.gtf ${vcf}

    whatshap haplotag \
        -o ${child_id}.haplotagged.bam \
        --reference ${reference} \
        ${vcf} ${child_bam}
    """

    stub:
    """
    touch ${family_id}.trio_phased.vcf.gz
    touch ${child_id}.haplotagged.bam
    touch ${family_id}.ped
    touch ${family_id}.blocks.tsv
    touch ${family_id}.gtf
    """
}