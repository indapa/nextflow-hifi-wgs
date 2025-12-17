
process annotate_sv_vcf {

    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    publishDir "${params.vep_output_dir}", mode: 'copy', pattern: '*.vcf.gz*'

    input:
    tuple val(sample_id), path(pbsv_vcf)
    path pigeon_gtf
    path pigeon_tbi
    path reference

    output:
    path "${sample_id}.pbsv.vep.vcf.gz", emit: vep_pbsv_vcf
    path "${sample_id}.pbsv.vep.vcf.gz.tbi", emit: vep_pbsv_tbi

   

   script:
    """
    tabix -p vcf ${pbsv_vcf}
    vep -i ${pbsv_vcf} -o ${sample_id}.pbsv.vep.vcf.gz --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --overlaps  --fork 8 --compress_output bgzip
    
    tabix -p vcf ${sample_id}.pbsv.vep.vcf.gz
    """
}

process annotate_cnv_vcf {

    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    publishDir params.vep_output_dir, mode: 'copy', pattern: '*.vcf.gz*'

    input:
    tuple val(sample_id), path(cnv_vcf)
    path pigeon_gtf
    path pigeon_tbi
    path reference

    output:
    path "${sample_id}.hificnv.vep.vcf.gz", emit: vep_cnv_vcf
    path "${sample_id}.hificnv.vep.vcf.gz.tbi", emit: vep_cnv_tbi

    script:

    """
    tabix -p vcf ${cnv_vcf}
    vep -i ${cnv_vcf} -o ${sample_id}.hificnv.vep.vcf.gz --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --overlaps  --fork 8 --compress_output bgzip 
    tabix -p vcf ${sample_id}.hificnv.vep.vcf.gz
    """
}

process annotate_vep_no_phased {

    /* annotate vcf that has not been phased */
    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    publishDir "${params.vep_output_dir}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.vcf.gz*'

    input:
    tuple val(sample_id), path(deepvariant_vcf), path(deepvariant_tbi)
    path pigeon_gtf
    path pigeon_tbi
    path reference

    output:
    tuple val(sample_id), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vep_vcf

    script:

    """
    vep -i ${deepvariant_vcf} -o ${sample_id}.deepvariant.vep.vcf.gz --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --fork 8 --compress_output bgzip
    tabix -p vcf ${sample_id}.deepvariant.vep.vcf.gz
    """
    
    stub:
    """
    touch ${sample_id}.deepvariant.vep.vcf.gz
    touch ${sample_id}.deepvariant.vep.vcf.gz.tbi
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