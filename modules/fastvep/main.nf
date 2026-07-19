process FASTVEP_ANNOTATE {
    container 'docker.io/indapa/fastvep:0.2.0'
    publishDir { "${params.deepvariant_output_dir}/${sample_id}" }, mode: 'copy', overwrite: true
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), path(vcf)
    path gff3
    path fasta
    path sa_dir

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz")

    script:
    """
    fastvep annotate \
        -i ${vcf} \
        -o ${sample_id}.annotated.vcf.gz \
        --gff3 ${gff3} \
        --fasta ${fasta} \
        --sa-dir ${sa_dir} \
        --hgvs
    """
}