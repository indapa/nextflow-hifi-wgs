#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_methylation_calling; sawfish_discover; sawfish_joint_call } from './modules/pbtools'
include { glnexus_trio_merge; deepvariant_wgs } from './modules/deepvariant'
include { bam_stats } from './modules/samtools'
include { annotate_vep } from './modules/ensemblvep'
include { whatshap_trio_phase } from './modules/whatshap'
include { mosdepth_run; infer_sex; plot_dist_coverage } from './modules/mosdepth'


// =========================================================================
//  WORKFLOW: READ ALIGNMENT + POST ALIGNMENT
// =========================================================================

workflow {
    if (params.help) {
        println """
        Available workflows:
        1. DEFAULT: nextflow run main.nf --samplesheet samples.csv performs read alignment and post-alignment analyses (e.g., bam stats)
        2. POST_ALIGNMENT_ONLY: nextflow run main.nf -entry POST_ALIGNMENT_ONLY --samplesheet samples.csv runs post-alignment only on pre-aligned BAMs
        3. RUN_DEEPTRIO: nextflow run main.nf -entry RUN_DEEPTRIO --trio_samplesheet trios.csv performs DeepTrio and WhatsHap phasing on
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


// =========================================================================
//  ENTRY POINT: POST-ALIGNMENT ONLY (skip alignment)
// =========================================================================

workflow POST_ALIGNMENT_ONLY {
    if (!params.samplesheet) {
        error "Parameter 'samplesheet' is required! CSV must have columns: sample_id, bam_file, bai_file"
    }

    if (!file(params.samplesheet).exists()) {
        exit 1, "Samplesheet file not found: ${params.samplesheet}"
    }

    def aligned_bam_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def bam = file(row.bam_file)
            def bai = file(row.bai_file)
            return tuple(sample_id, bam, bai)
        }

    POST_ALIGNMENT(aligned_bam_ch)
}


// =========================================================================
//  SUB-WORKFLOW: POST ALIGNMENT
// =========================================================================

workflow POST_ALIGNMENT {
    take:
    aligned_bam_ch

    main:
    /*
    bam_stats(aligned_bam_ch)

    deepvariant_wgs(
        file(params.reference),
        file(params.reference_index),
        aligned_bam_ch
    )

    cpg_methylation_calling(
        aligned_bam_ch,
        file(params.reference),
        file(params.reference_index)
    )
    */
    mosdepth_run(aligned_bam_ch)
    

    infer_sex(mosdepth_run.out.summary)

    plot_dist_coverage(mosdepth_run.out.global_dist)

    // Determine expected BED file based on inferred sex
    def expected_bed_ch = infer_sex.out.sex.map { sample_id, sex_csv ->
        def lines = sex_csv.readLines()
        def sex = lines.size() > 1 ? lines[1].split(',')[3].trim() : 'UNKNOWN'

        def expected_bed
        if (sex == 'FEMALE') {
            expected_bed = file(params.expected_XX_bed)
        } else if (sex == 'MALE') {
            expected_bed = file(params.expected_XY_bed)
        } else {
            throw new Exception("Error: Invalid or missing sex '${sex}' inferred for sample ${sample_id}. Expected 'FEMALE' or 'MALE'.")
        }

        return tuple(sample_id, expected_bed)
    }

    // Join BAM channel with the correct BED file channel by sample_id
    // Creates: [sample_id, bam, bai, expected_bed]
    def sawfish_in_ch = aligned_bam_ch.join(expected_bed_ch, by: 0)

    sawfish_discover(
        sawfish_in_ch,
        file(params.excluded_bed),
        file(params.reference),
        file(params.reference_index)
    )

    sawfish_joint_call(
        sawfish_discover.out.discover_dir.collect(),
        file(params.reference),
        file(params.reference_index)
    )
}










