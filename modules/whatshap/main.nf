process whatshap_trio_phase {
    tag "${family_id}"
    publishDir "${params.deepvariant_output_dir}/DV_trio/${family_id}", mode: 'copy', overwrite: true
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
}