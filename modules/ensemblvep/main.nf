
process annotate_vep_no_phased {

    /* annotate vcf that has not been phased */
    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    publishDir params.vep_output_dir, mode: 'copy', pattern: '*.vcf.gz*'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path pigeon_gtf
    path pigeon_tbi
    path reference

    output:
    tuple val(sample_id), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vep_vcf

    script:

    """
    vep -i ${vcf} -o ${sample_id}.vep.vcf.gz --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --fork 8 --compress_output bgzip
    tabix -p vcf ${sample_id}.vep.vcf.gz
    """
    

}


process annotate_vep {

    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    publishDir params.vep_output_dir, mode: 'copy', pattern: '*.vcf.gz*'

    input:

    tuple val(sample_id), path(phased_vcf), path(phased_tbi)
    path pigeon_gtf
    path pigeon_tbi
    path reference
    

    output:
    path "${sample_id}.deepvariant.phased.vep.vcf.gz", emit: vep_phased_deepvariant
    path "${sample_id}.deepvariant.phased.vep.vcf.gz.tbi", emit: vep_phased_deepvariant_tbi

    script:
    
    """

    vep -i ${phased_vcf} -o ${sample_id}.deepvariant.phased.vep.vcf.gz --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --fork 8 --compress_output bgzip
    tabix -p vcf ${sample_id}.deepvariant.phased.vep.vcf.gz
    """


}