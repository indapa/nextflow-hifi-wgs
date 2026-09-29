include { hiphase_small_variants } from '../../modules/pbtools'

// Runs the existing per-sample hiphase_small_variants process on just the two
// parents of a trio, using each parent's own DeepTrio-called small-variant VCF
// (read-backed phasing, independent of the child's PED-based WhatsHap phase)
// to emit a phased VCF + haplotagged BAM per parent.
workflow HIPHASE_TRIO_PARENTS {

    take:
    deeptrio_vcf_ch // tuple(family_id, sample_id, vcf, tbi) -- e.g. DEEPTRIO_WGS.out.vcf
    aligned_bam_ch  // tuple(sample_id, bam, bai) -- e.g. individual_aligned_bams
    sample_roles_ch // tuple(sample_id, role)
    ref_file        // path(reference_fasta)
    ref_index_file  // path(reference_fai)

    main:
    // Keep only parent1/parent2 samples and drop family_id/role -- hiphase
    // here phases each parent independently, one sample at a time
    parents_vcf_ch = deeptrio_vcf_ch
        .map { _family_id, sample_id, vcf, tbi -> tuple(sample_id, vcf, tbi) }
        .join(sample_roles_ch, by: 0)
        .filter { _sample_id, _vcf, _tbi, role -> role == 'parent1' || role == 'parent2' }
        .map { sample_id, vcf, tbi, _role -> tuple(sample_id, vcf, tbi) }

    hiphase_input_ch = parents_vcf_ch
        .join(aligned_bam_ch, by: 0)
        // -> tuple(sample_id, vcf, tbi, bam, bai)

    hiphase_small_variants(
        hiphase_input_ch,
        ref_file,
        ref_index_file
    )

    emit:
    phased_vcf      = hiphase_small_variants.out.phased_vcf       // tuple(sample_id, vcf, tbi)
    haplotagged_bam = hiphase_small_variants.out.haplotagged_bam  // tuple(sample_id, bam, bai)
    stats           = hiphase_small_variants.out.stats            // tuple(sample_id, stats, blocks, summary)
}
