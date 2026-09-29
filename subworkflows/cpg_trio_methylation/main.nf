include { cpg_methylation_calling } from '../../modules/pbtools'
include { samtools_index }          from '../../modules/samtools'

// Runs CpG methylation calling (pb-CpG-tools) on every haplotagged BAM in a
// trio: the child's haplotagged BAM (from WHATSHAP_TRIO_PHASE_BY_CHROM) and
// each parent's haplotagged BAM (from HIPHASE_TRIO_PARENTS).
workflow CPG_TRIO_METHYLATION {

    take:
    child_haplotagged_bam_ch  // path: WHATSHAP_TRIO_PHASE_BY_CHROM.out.haplotagged_bam (bare path, no index)
    parent_haplotagged_bam_ch // tuple(sample_id, bam, bai): HIPHASE_TRIO_PARENTS.out.haplotagged_bam
    ref_file                  // path(reference_fasta)
    ref_index_file            // path(reference_fai)

    main:
    // The child's haplotagged BAM has no sample_id or index attached yet --
    // recover the sample_id from the filename and index it
    child_cpg_input_ch = child_haplotagged_bam_ch
        .map { bam -> tuple(bam.baseName.replaceAll(/\.haplotagged/, ''), bam) }

    samtools_index(child_cpg_input_ch)

    // Combine child + both parents into one tuple(sample_id, bam, bai) stream
    all_haplotagged_bam_ch = samtools_index.out.mix(parent_haplotagged_bam_ch)

    cpg_methylation_calling(
        all_haplotagged_bam_ch,
        ref_file,
        ref_index_file
    )

    emit:
    combined_bed = cpg_methylation_calling.out.combined_bed // tuple(sample_id, bed, tbi)
    hap1_bed     = cpg_methylation_calling.out.hap1_bed     // tuple(sample_id, bed, tbi)
    hap2_bed     = cpg_methylation_calling.out.hap2_bed     // tuple(sample_id, bed, tbi)
    combined_bw  = cpg_methylation_calling.out.combined_bw  // tuple(sample_id, bw)
    hap1_bw      = cpg_methylation_calling.out.hap1_bw      // tuple(sample_id, bw)
    hap2_bw      = cpg_methylation_calling.out.hap2_bw      // tuple(sample_id, bw)
}
