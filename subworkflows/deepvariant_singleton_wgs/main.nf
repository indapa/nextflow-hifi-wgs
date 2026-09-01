include { DEEPVARIANT_CHUNK; concat_full_genome_vcf_singleton; concat_chrom_chunks_vcf_singleton} from '../../modules/deepvariant'
include {slice_singleton_bam_by_interval} from '../../modules/samtools'


workflow DEEPVARIANT_SINGLETON_WGS {

    take:
    ref_file       // path(reference_fasta)
    ref_index_file // path(reference_fai)
    aligned_bam_ch // channel: tuple(sample_id, bam, bai)
    bed_chunks_ch  // channel: path(bed) e.g., Channel.fromPath("beds/*.bed")

    main:
    // 1. Combine sample BAM with each 50MB BED chunk
    scatter_ch = aligned_bam_ch
        .combine(bed_chunks_ch)
        .map { sample_id, bam, bai, bed ->
            tuple(sample_id, bam, bai, bed)
        }

    // 2. Extract 50MB regional mini-BAMs using samtools
    slice_singleton_bam_by_interval(scatter_ch)

    // 3. Call variants across 50MB chunks in parallel
    DEEPVARIANT_CHUNK(
        ref_file,
        ref_index_file,
        slice_singleton_bam_by_interval.out.sliced_singleton_package
    )

    // 4. Prepare VCF chunks channel: tuple([sample_id, chrom, "vcf.gz"], vcf, tbi)
    ch_vcf_prepared = DEEPVARIANT_CHUNK.out.vcf
        .map { sample_id, vcf, tbi ->
            def chrom = vcf.name.split('_')[0].replaceAll("^${sample_id}\\.", "") // Extracts "chr1"
            def meta  = [sample_id, chrom, "vcf.gz"]
            tuple(meta, vcf, tbi)
        }
        .groupTuple(by: 0)

    // 5. Prepare gVCF chunks channel: tuple([sample_id, chrom, "g.vcf.gz"], gvcf, gtbi)
    ch_gvcf_prepared = DEEPVARIANT_CHUNK.out.gvcf
        .map { sample_id, gvcf, gtbi ->
            def chrom = gvcf.name.split('_')[0].replaceAll("^${sample_id}\\.", "") // Extracts "chr1"
            def meta  = [sample_id, chrom, "g.vcf.gz"]
            tuple(meta, gvcf, gtbi)
        }
        .groupTuple(by: 0)

    // 6. Merge 50MB chunks into single per-chromosome VCFs and gVCFs
    ch_all_chrom_chunks = ch_vcf_prepared.mix(ch_gvcf_prepared)
    concat_chrom_chunks_vcf_singleton(ch_all_chrom_chunks)

    // 7. Group per-chromosome files by sample_id & file_type for genome-wide merging
    ch_genome_inputs = concat_chrom_chunks_vcf_singleton.out.chrom_merged_file
        .groupTuple(by: [0, 1]) // Groups by [sample_id, file_type]

    // 8. Call the concatenation process on the grouped channel
    concat_full_genome_vcf_singleton(ch_genome_inputs)

    // 9. Separate output channel into VCF and gVCF branches
    concat_full_genome_vcf_singleton.out.genome_vcf
        .branch { _sample_id, vcf, _tbi ->
            gvcf: vcf.name.endsWith("g.vcf.gz")
            vcf:  !vcf.name.endsWith("g.vcf.gz")
        }
        .set { genome_files }

    emit:
    vcf  = genome_files.vcf   // tuple(sample_id, vcf, tbi)
    gvcf = genome_files.gvcf  // tuple(sample_id, gvcf, gtbi)
}