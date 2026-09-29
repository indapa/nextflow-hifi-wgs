include { deeptrio_wgs_by_chrom; concat_chrom_chunks_vcf } from '../../modules/deepvariant'
include { concat_wgs_vcf as concat_genome_vcf; concat_wgs_vcf as concat_genome_gvcf } from '../../modules/deepvariant'
include { slice_trio_bams_by_interval } from '../../modules/samtools'


workflow DEEPTRIO_WGS {

    take:
    ref_file       // path(reference_fasta)
    ref_index_file // path(reference_fai)
    trio_bam_ch    // channel: tuple(family_id, child_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai)
    bed_chunks_ch  // channel: path(bed) e.g., Channel.fromPath("beds/*.bed")

    main:
    // 1. Combine trio BAMs with each 50MB BED chunk
    scatter_ch = trio_bam_ch
        .combine(bed_chunks_ch)
        .map { family_id, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai, bed ->
            tuple(family_id, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai, bed)
        }

    // 2. Extract 50MB regional mini-BAMs for child/parent1/parent2 using samtools
    slice_trio_bams_by_interval(scatter_ch)

    // 3. Call DeepTrio across 50MB chunks in parallel
    deeptrio_wgs_by_chrom(
        ref_file,
        ref_index_file,
        slice_trio_bams_by_interval.out.sliced_trio_package
    )

    // 4. Fan the per-member VCF/gVCF outputs into a single channel keyed by
    //    [family_id, sample_id, chrom, file_type]
    ch_chunks = Channel.empty()
        .mix(
            deeptrio_wgs_by_chrom.out.child_vcf,
            deeptrio_wgs_by_chrom.out.child_gvcf,
            deeptrio_wgs_by_chrom.out.p1_vcf,
            deeptrio_wgs_by_chrom.out.p1_gvcf,
            deeptrio_wgs_by_chrom.out.p2_vcf,
            deeptrio_wgs_by_chrom.out.p2_gvcf
        )
        .map { family_id, sample_id, chrom_chunk, vcf, tbi ->
            def chrom     = chrom_chunk.split('_')[0]   // e.g. "chr1_0_50000000" -> "chr1"
            def file_type = vcf.name.endsWith('g.vcf.gz') ? 'g.vcf.gz' : 'vcf.gz'
            def meta = [family_id, sample_id, chrom, file_type]
            tuple(meta, vcf, tbi)
        }
        .groupTuple(by: 0)

    // 5. Merge 50MB chunks into single per-chromosome VCFs and gVCFs
    concat_chrom_chunks_vcf(ch_chunks)

    // 6. Group per-chromosome files by [family_id, sample_id, file_type] for genome-wide merging
    ch_genome_inputs = concat_chrom_chunks_vcf.out.merged_file
        .map { meta, vcf, tbi ->
            tuple(meta[0], meta[1], meta[3], vcf, tbi)
        }
        .groupTuple(by: [0, 1, 2])

    // 7. Split into VCF / gVCF branches so each can be concatenated with the right extension
    ch_genome_inputs
        .branch { _family_id, _sample_id, file_type, _vcf, _tbi ->
            gvcf: file_type == 'g.vcf.gz'
            vcf:  file_type == 'vcf.gz'
        }
        .set { by_type }

    // 8. Concatenate per-chromosome files into final genome-wide VCF/gVCF per trio member
    ch_vcf_out = concat_genome_vcf(
        by_type.vcf.map { family_id, sample_id, _file_type, vcf, tbi -> tuple(family_id, sample_id, vcf, tbi) },
        'vcf.gz'
    )
    ch_gvcf_out = concat_genome_gvcf(
        by_type.gvcf.map { family_id, sample_id, _file_type, gvcf, gtbi -> tuple(family_id, sample_id, gvcf, gtbi) },
        'g.vcf.gz'
    )

    emit:
    vcf  = ch_vcf_out  // tuple(family_id, sample_id, vcf, tbi)
    gvcf = ch_gvcf_out // tuple(family_id, sample_id, gvcf, gtbi)
}
