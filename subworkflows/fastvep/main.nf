/*
========================================================================================
    MODULE IMPORT
========================================================================================
*/
include { FASTVEP_ANNOTATE_SINGLETON_VCF } from '../../modules/fastvep/'

/*
========================================================================================
    SUBWORKFLOW: FASTVEP_ANNOTATE_WGS
========================================================================================
*/
workflow FASTVEP_ANNOTATE_WGS {

    take:
    vcf_ch   // channel: tuple(sample_id, vcf, tbi) OR tuple(sample_id, vcf)
    gff3     // path(gff3)
    fasta    // path(fasta)
    sa_dir   // channel/path to sa_dir files

    main:
    // 1. Ensure input tuple matches the 2-element tuple expected by process: [sample_id, vcf]
    ch_vep_input = vcf_ch
        .map { item ->
            def sample_id = item[0]
            def vcf       = item[1]
            tuple(sample_id, vcf)
        }

    // 2. Call FastVEP annotation process
    FASTVEP_ANNOTATE_SINGLETON_VCF(
        ch_vep_input,
        gff3,
        fasta,
        sa_dir
    )

    emit:
    FASTVEP_ANNOTATE_SINGLETON_VCF.out.annotated_vcf // tuple(sample_id, annotated_vcf)
}