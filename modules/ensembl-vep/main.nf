process vep-annotate {

    tag "${sample_id}"
    container 'quay.io/biocontainers/ensembl-vep:112.0--pl5321h2a3209d_0'
    publishDir params.deepvariant_output_dir, mode: 'copy'

    input:

    path phased_vcf
    tuple val(sample_id), path(haplotagged_bam)
    path pigeon_gtf
    path pigeon_tbi
    path reference
    

    output:
    path "${sample_id}.deepvariant.phased.vep.vcf.gz", emit: vep_phased_deepvariant

    script:
    """

    vep --help

    vep -i ${phased_vcf} -o ${sample_id}.deepvariant.phased.vep.vcf.gz --force_overwrite --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --fork 8

    """


}