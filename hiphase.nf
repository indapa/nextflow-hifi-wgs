#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { hiphase_small_variants } from './modules/pbtools'

// Parameters
params.samplesheet        = ''
params.reference          = ''
params.reference_index    = ''
params.hiphase_output_dir = ''

// Parameter validation
def required_params = ['samplesheet', 'reference', 'reference_index', 'hiphase_output_dir']
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}

// Channels
ref_ch     = Channel.value(file(params.reference, checkIfExists: true))
ref_idx_ch = Channel.value(file(params.reference_index, checkIfExists: true))

input_ch = Channel.fromPath(file(params.samplesheet, checkIfExists: true))
    .splitCsv(header: true)
    .map { row ->
        tuple(row.sample_id, file(row.vcf), file(row.vcf_tbi), file(row.pbmm2_bam), file(row.pbmm2_bai))
    }

workflow {
    // Call the imported process exactly like before
    hiphase_small_variants(input_ch, ref_ch, ref_idx_ch)
}