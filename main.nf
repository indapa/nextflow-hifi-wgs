#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_methylation_calling } from './modules/pbtools'
include {  glnexus_trio_merge; deepvariant_wgs } from './modules/deepvariant'
include { bam_stats } from './modules/samtools'
include { annotate_vep } from './modules/ensemblvep'
include { whatshap_trio_phase } from './modules/whatshap'
include { mosdepth_run; infer_sex; plot_dist_coverage } from './modules/mosdepth'



// =========================================================================
//  WORKFLOW READ ALIGNMENT + POST ALIGNMENT 
// =========================================================================




workflow  {
    // ✅ Check samplesheet only in this workflow

    if (params.help) {
        println """
        Available workflows:
        1. DEFAULT: nextflow run main.nf --samplesheet samples.csv performs read alignment and post-alignment analyses (e.g., bam stats)
        2. RUN_DEEPTRIO: nextflow run main.nf -entry RUN_DEEPTRIO --trio_samplesheet trios.csv performs DeepTrio and WhatsHap phasing on
        trios in the samplesheet with aligned bams
        
        """.stripIndent()
        exit 0
    }
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
        
    )

    /* post alignment */
    POST_ALIGNMENT(
        pbmm2_align.out.aligned_bam
    )

}

workflow POST_ALIGNMENT {
    // This can include processes like variant calling, phasing, annotation, etc.

    take:
    aligned_bam_ch

    main:
    bam_stats(
        aligned_bam_ch
    )

    
    deepvariant_wgs(
        file(params.reference),
        file(params.reference_index),
        aligned_bam_ch
    )

    cpg_methylation_calling(
        aligned_bam_ch,           // tuple(sample_id, bam, bai) from pbmm2_align
        file(params.reference),
        file(params.reference_index)
    )

     mosdepth_run(aligned_bam_ch)

    infer_sex(mosdepth_run.out.summary)
    
    plot_dist_coverage(mosdepth_run.out.global_dist)
}










