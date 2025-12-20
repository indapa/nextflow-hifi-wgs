#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align; hiphase_small_variants } from './modules/pbtools'
include { deepvariant; BCFTOOLS_STATS; bcftools_deepvariant_norm; deepvariant_targeted_region} from './modules/deepvariant'
include { bam_stats } from './modules/samtools'

// Remove the top-level checkSamplesheet call too!
// Only check samplesheet in workflows that need it

workflow ALIGN_DEEP_VARIANT_BCFTOOLS_STATS_SYT1 {
    // ✅ Check samplesheet only in this workflow
    if (!params.samplesheet) {
        error "Parameter 'samplesheet' is required for this workflow!"
    }
    
    if (!file(params.samplesheet).exists()) {
        exit 1, "Samplesheet file not found: ${params.samplesheet}"
    }
    
    // ✅ Create channel only in this workflow
    def input_bams_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def bam_file = file(row.bam_file)
            return tuple(sample_id, bam_file)
        }
    
    /* read alignment */
    pbmm2_align(
        file(params.reference),
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )

  

    /* deepvariant  */
    deepvariant(
        file(params.reference), 
        file(params.reference_index), 
        pbmm2_align.out.aligned_bam, 
        params.deepvariant_threads 
        
    )
    // Transform [id, bam, bai] -> [id, bam]
    ch_for_stats = pbmm2_align.out.aligned_bam
        .map { sample_id, bam, bai -> tuple(sample_id, bam) }
        
    bam_stats(ch_for_stats)

    /* bcftools normalization 
    bcftools_deepvariant_norm(
        file(params.reference), 
        deepvariant.out.vcf_tuple
    )
 */
   /*
    hiphase_input_ch = deepvariant.out.vcf_tuple
        .join( pbmm2_align.out.aligned_bam )


    // Run HiPhase
    
    hiphase_small_variants(
        hiphase_input_ch,
        params.reference,
        params.reference_index
    )
    */
}

