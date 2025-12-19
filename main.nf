#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { pbmm2_align; hiphase_small_variants } from './modules/pbtools'
include { deepvariant; BCFTOOLS_STATS; bcftools_deepvariant_norm; deepvariant_targeted_region} from './modules/deepvariant'

def required_params = ['reference', 'samplesheet'  ]
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}

def checkSamplesheet(samplesheet_file) {
    if (!file(samplesheet_file).exists()) {
        exit 1, "Samplesheet file not found: ${samplesheet_file}"
    }
    return file(samplesheet_file)
}

ss_status = checkSamplesheet(params.samplesheet)

// Create channels for input BAM files (DSL2 style)
def input_bams_ch = channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample_id
        def bam_file = file(row.bam_file)
        
        return tuple(sample_id, bam_file)
    }



workflow ALIGN_DEEP_VARIANT_BCFTOOLS_STATS_SYT1 {
    /* read alignment */
    pbmm2_align(
        file(params.reference),
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )

    /* deepvariant - NOTE: using deepvariant_targeted_region */
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


