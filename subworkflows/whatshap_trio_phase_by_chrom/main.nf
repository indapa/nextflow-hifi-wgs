// subworkflows/local/whatshap_trio_phase_by_chrom/main.nf
include { WHATSHAP_PHASE_CHROM }    from '../../modules/whatshap'
include { CONCAT_PHASED_VCFS }      from '../../modules/whatshap'
include { WHATSHAP_STATS_HAPLOTAG } from '../../modules/whatshap'

workflow WHATSHAP_TRIO_PHASE_BY_CHROM {
    take:
    ch_trio_input   // [family_id, vcf, vcf_tbi, child_id, child_bam, child_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai]
    ch_reference    // value channel: path to reference fasta
    ch_ref_index    // value channel: path to reference fasta index
    chromosomes     // channel.fromList(params.chromosomes)

    main:
    // Split ch_trio_input into two branches so it isn't consumed twice
    ch_trio_input.multiMap { row ->
        phase: row
        child_info: tuple(row[0], row[3], row[4], row[5], row[6], row[9])
        // [family_id, child_id, child_bam, child_bai, p1_id, p2_id]
    }.set { ch_split }

    // Scatter: combine each trio with every chromosome
    ch_phase_input = ch_split.phase
        .combine(chromosomes)
        .map { family_id, vcf, vcf_tbi, child_id, child_bam, child_bai,
               p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai, chrom ->
            tuple(family_id, chrom, vcf, vcf_tbi,
                  child_id, child_bam, child_bai,
                  p1_id, p1_bam, p1_bai,
                  p2_id, p2_bam, p2_bai)
        }

    // Phase each chromosome independently
    WHATSHAP_PHASE_CHROM(ch_phase_input, ch_reference, ch_ref_index)

    // Gather: group per-chrom results by family_id, then concatenate
    // WHATSHAP_PHASE_CHROM.out.phased_vcf emits [family_id, phased_vcf]
    ch_to_merge = WHATSHAP_PHASE_CHROM.out.phased_vcf
        .groupTuple(by: 0)
        // [family_id, [chr1.vcf.gz, chr2.vcf.gz, ...]]

    CONCAT_PHASED_VCFS(ch_to_merge)

    // Post-phase: stats + haplotag on the merged VCF
    // CONCAT_PHASED_VCFS.out.merged_vcf emits [family_id, merged_vcf, merged_tbi]
    ch_stats_input = CONCAT_PHASED_VCFS.out.merged_vcf
        .join(ch_split.child_info, by: 0)
        // [family_id, merged_vcf, merged_tbi, child_id, child_bam, child_bai, p1_id, p2_id]

    WHATSHAP_STATS_HAPLOTAG(ch_stats_input, ch_reference, ch_ref_index)

    emit:
    phased_vcf      = CONCAT_PHASED_VCFS.out.merged_vcf       // [family_id, vcf, tbi]
    haplotagged_bam = WHATSHAP_STATS_HAPLOTAG.out.haplotagged_bam
    block_stats     = WHATSHAP_STATS_HAPLOTAG.out.block_stats
    block_gtf       = WHATSHAP_STATS_HAPLOTAG.out.block_gtf
}