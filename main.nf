#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align; hiphase_small_variants } from './modules/pbtools'
include { deepvariant; BCFTOOLS_STATS; bcftools_deepvariant_norm; deepvariant_targeted_region} from './modules/deepvariant'

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

    /* deepvariant - targeted region */
    deepvariant_targeted_region(
        file(params.reference), 
        file(params.reference_index), 
        pbmm2_align.out.aligned_bam, 
        params.deepvariant_threads, 
        params.syt1_region
    )

    /* bcftools normalization */
    bcftools_deepvariant_norm(
        file(params.reference), 
        deepvariant_targeted_region.out.vcf_tuple
    )
}

workflow HIPHASE_PHASING {
    // ✅ Check hiphase_samplesheet only in this workflow
    if (!params.hiphase_samplesheet) {
        error "Parameter 'hiphase_samplesheet' is required for this workflow!"
    }
    
    if (!file(params.hiphase_samplesheet).exists()) {
        exit 1, "Samplesheet file not found: ${params.hiphase_samplesheet}"
    }
    
    // ✅ Create channel only in this workflow
    def hiphase_input_ch = channel.fromPath(params.hiphase_samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def vcf = file(row.vcf)
            def vcf_tbi = file(row.vcf_tbi)
            def bam = file(row.bam)
            def bai = file(row.bai)
            
            return tuple(sample_id, vcf, vcf_tbi, bam, bai)
        }
    
    log.info "Running HiPhase phasing workflow..."
    log.info "Reference: ${params.reference}"
    log.info "Samplesheet: ${params.hiphase_samplesheet}"
    
    hiphase_small_variants(
        hiphase_input_ch,
        file(params.reference),
        file(params.reference_index)
    )
}