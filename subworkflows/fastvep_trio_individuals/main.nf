include { SPLIT_TRIO_VCF_BY_SAMPLE }      from '../../modules/whatshap'
include { FASTVEP_ANNOTATE_SINGLETON_VCF } from '../../modules/fastvep'

// Splits the WhatsHap trio-phased VCF back into one VCF per individual
// (proband + both parents), then annotates each individually-phased VCF with
// FastVEP.
workflow FASTVEP_ANNOTATE_TRIO_INDIVIDUALS {

    take:
    phased_vcf_ch // tuple(family_id, vcf, tbi) -- e.g. WHATSHAP_TRIO_PHASE_BY_CHROM.out.phased_vcf
    trio_bam_ch   // tuple(family_id, child_id, child_bam, child_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai) -- e.g. trio_bams_assembled
    gff3          // path(gff3)
    fasta         // path(fasta)
    sa_dir        // channel/path to sa_dir files

    main:
    // Only the sample IDs are needed to split the joint VCF
    ch_split_input = phased_vcf_ch
        .join(trio_bam_ch, by: 0)
        .map { family_id, vcf, tbi, child_id, _cb, _cbi, p1_id, _p1b, _p1bi, p2_id, _p2b, _p2bi ->
            tuple(family_id, vcf, tbi, child_id, p1_id, p2_id)
        }

    SPLIT_TRIO_VCF_BY_SAMPLE(ch_split_input)

    ch_individual_vcfs = channel.empty()
        .mix(
            SPLIT_TRIO_VCF_BY_SAMPLE.out.child_vcf,
            SPLIT_TRIO_VCF_BY_SAMPLE.out.p1_vcf,
            SPLIT_TRIO_VCF_BY_SAMPLE.out.p2_vcf
        )
        // -> tuple(sample_id, vcf)

    FASTVEP_ANNOTATE_SINGLETON_VCF(
        ch_individual_vcfs,
        gff3,
        fasta,
        sa_dir
    )

    emit:
    FASTVEP_ANNOTATE_SINGLETON_VCF.out.annotated_vcf // tuple(sample_id, annotated_vcf)
}
