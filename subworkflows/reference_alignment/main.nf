include { pbmm2_align; PBMM2_ALIGN_SPOT_CHUNK;  MERGE_SPOT_CHUNKS; MAKE_PBI} from '../../modules/pbtools'

workflow PBMM2_ALIGN {
    take:
    ref_file      // Path to reference genome file
    alignment_ch  // Channel: tuple(sample_id, bam)

    main:
    pbmm2_align(ref_file, alignment_ch)

    emit:
    pbmm2_align.out.aligned_bam  // Adjust process output name if needed
}

workflow PBMM2_SPOT_WGS {

    take:
    ref_mmi          // path(ref.mmi)
    unaligned_bam_ch // channel: tuple(sample_id, path(bam))

    main:
    // 1. Branch channel to check if .pbi exists alongside the BAM file
    unaligned_bam_ch
        .branch { _sample_id, bam ->
            has_pbi:   file("${bam}.pbi").exists()
            needs_pbi: !file("${bam}.pbi").exists()
        }
        .set { bams_branched }

    // 2. Generate .pbi for BAMs that lack one
    MAKE_PBI(bams_branched.needs_pbi)

    // 3. Attach existing .pbi files for BAMs that already have one
    ch_existing_pbi = bams_branched.has_pbi
        .map { sample_id, bam ->
            tuple(sample_id, bam, file("${bam}.pbi"))
        }

    // 4. Combine both branches back into a unified tuple: [sample_id, bam, pbi]
    ch_ready_for_scatter = MAKE_PBI.out.bam_pbi.mix(ch_existing_pbi)

    // 5. Scale chunk count dynamically (1 chunk per 4-5 GB of input BAM, Min 12, Max 24)
    ch_scatter = ch_ready_for_scatter
        .flatMap { sample_id, bam, pbi ->
            def bam_gb = bam.size() / (1024 * 1024 * 1024)
            def num_chunks = Math.max(12, Math.min(24, (bam_gb / 4.5).intValue()))
            
            (1..num_chunks).collect { chunk_id ->
                tuple(sample_id, bam, pbi, chunk_id, num_chunks)
            }
        }

    // 6. Align scattered chunks in parallel on AWS Batch Spot nodes
    PBMM2_ALIGN_SPOT_CHUNK(ref_mmi, ch_scatter)

    // 7. Group output chunks per sample ID
    ch_grouped_chunks = PBMM2_ALIGN_SPOT_CHUNK.out.chunk_bam
        .groupTuple(by: 0)

    // 8. Merge into final genome BAM & index
    MERGE_SPOT_CHUNKS(ch_grouped_chunks)

    emit:
    MERGE_SPOT_CHUNKS.out.aligned_bam // tuple(sample_id, aligned_bam, bai)
}