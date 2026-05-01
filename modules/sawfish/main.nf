process sawfish_discover {
    /*
     * Structural variant discovery from PacBio HiFi aligned BAMs using sawfish.
     *
     * Sawfish identifies large structural variants (insertions, deletions,
     * duplications, inversions, translocations) from long-read alignments.
     */

    label 'high_memory'
    tag "$sample_id"
    publishDir "${params.sawfish_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    container "quay.io/pacbio/sawfish@sha256:18ba096219fea38d6b32f5706fb794a05cc5d1d6cc16e2a09e3a13d62d8181d4"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_index

    output:
    tuple val(sample_id), path("${sample_id}.sawfish.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${sample_id}.sawfish.vcf.gz.tbi"), emit: vcf_index

    script:
    """
    sawfish discover \
        --ref ${ref} \
        --bam ${bam} \
        --threads ${task.cpus} \
        --output-dir sawfish_output

    # Rename and index the output VCF
    mv sawfish_output/variants.vcf.gz ${sample_id}.sawfish.vcf.gz
    mv sawfish_output/variants.vcf.gz.tbi ${sample_id}.sawfish.vcf.gz.tbi
    """

    stub:
    """
    touch ${sample_id}.sawfish.vcf.gz
    touch ${sample_id}.sawfish.vcf.gz.tbi
    """
}
