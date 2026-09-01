process FASTVEP_ANNOTATE_SINGLETON_VCF {
    tag "${sample_id}"
    container 'docker.io/indapa/fastvep:0.2.0'
    publishDir "${params.deepvariant_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), path(vcf)
    path gff3
    path fasta
    path("sa_dir/*")

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz"), emit: annotated_vcf

    script:
    """
    fastvep \
        --input ${vcf} \
        --gff ${gff3} \
        --fasta ${fasta} \
        --sa-dir sa_dir \
        --output ${sample_id}.annotated.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.annotated.vcf.gz
    """
}

process FASTVEP_ANNOTATE_TRIO_VCF {
    container 'docker.io/indapa/fastvep:0.2.0'
    publishDir { "${params.deepvariant_output_dir}/${family_id}" }, mode: 'copy', overwrite: true
    cpus 4
    memory '16 GB'

    input:
    tuple val(family_id), path(vcf)
    path gff3
    path fasta
    path("sa_dir/*")

    output:
    tuple val(family_id), path("${family_id}.annotated.vcf.gz")

    script:
    """
    fastvep annotate \
        -i ${vcf} \
        -o ${family_id}.annotated.vcf.gz \
        --gff3 ${gff3} \
        --fasta ${fasta} \
        --sa-dir sa_dir \
        --hgvs
    """

    stub:
    """
    touch ${family_id}.annotated.vcf.gz
    """
}