// subworkflows/local/whatshap_trio_phase_by_chrom.nf
include { WHATSHAP_PHASE_CHROM }     from '../../modules/whatshap'
include { CONCAT_PHASED_VCFS }       from '../../modules/whatshap'
include { WHATSHAP_STATS_HAPLOTAG }  from '../../modules/whatshap'

workflow WHATSHAP_TRIO_PHASE_BY_CHROM {
    take:
    ch_trio_input   // [family_id, vcf, vcf_tbi, child_id, child_bam, child_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai]
    ch_reference
    ch_ref_index
    chromosomes     // channel.fromList(params.chromosomes)

    main:
    // Scatter: combine each trio with every chromosome
    // Result: [family_id, chrom, vcf, tbi, child_id, child_bam, child_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai]
    ch_phase_input = ch_trio_input.combine(chromosomes)
        .map { family_id, vcf, vcf_tbi, child_id, child_bam, child_bai,
               p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai, chrom ->
            tuple(family_id, chrom, vcf, vcf_tbi,
                  child_id, child_bam, child_bai,
                  p1_id, p1_bam, p1_bai,
                  p2_id, p2_bam, p2_bai)
        }

    // Phase each chromosome independently
    WHATSHAP_PHASE_CHROM(ch_phase_input, ch_reference, ch_ref_index)

    // Gather: group per-chrom results by family, then concatenate
    ch_to_merge = WHATSHAP_PHASE_CHROM.out.phased_vcf
        .groupTuple(by: 0)
        // [family_id, [chr1.phased.vcf.gz, chr2.phased.vcf.gz, ...]]

    CONCAT_PHASED_VCFS(ch_to_merge)

    // Post-phase: stats + haplotag on the merged VCF
    ch_child_info = ch_trio_input.map { row ->
        tuple(row[0], row[3], row[4], row[5], row[6], row[9])
        // [family_id, child_id, child_bam, child_bai, p1_id, p2_id]
    }

    ch_stats_input = CONCAT_PHASED_VCFS.out.merged_vcf
        .join(ch_child_info, by: 0)
        // [family_id, merged_vcf, merged_tbi, child_id, child_bam, child_bai, p1_id, p2_id]

    WHATSHAP_STATS_HAPLOTAG(ch_stats_input, ch_reference, ch_ref_index)

    emit:
    phased_vcf      = WHATSHAP_STATS_HAPLOTAG.out.phased_vcf
    haplotagged_bam = WHATSHAP_STATS_HAPLOTAG.out.haplotagged_bam
    block_stats     = WHATSHAP_STATS_HAPLOTAG.out.block_stats
    block_gtf       = WHATSHAP_STATS_HAPLOTAG.out.block_gtf
}