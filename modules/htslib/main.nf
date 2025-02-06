process tabix {
    tag "${sample_id}"

    publishDir params.tabix_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'

    input:
    tuple val(sample_id), path(vcf)

    output:
    path "${vcf}.tbi", emit: tbi

    script:
    """

    tabix -p vcf ${vcf}

    """
}