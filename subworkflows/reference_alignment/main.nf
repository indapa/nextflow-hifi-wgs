include { pbmm2_align } from '../../modules/pbtools'

workflow PBMM2_ALIGN {
    take:
    ref_file      // Path to reference genome file
    alignment_ch  // Channel: tuple(sample_id, bam)

    main:
    pbmm2_align(ref_file, alignment_ch)

    emit:
    pbmm2_align.out  // Adjust process output name if needed
}