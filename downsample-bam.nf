#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { downsample } from './modules/samtools'

// Parameters
params.input        = ''
params.output = ''

// Parameter validation
def required_params = ['samplesheet', 'output']
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}

//uple val(sample_id), val(coverage), val(frac), path(bam)
// Channels
input_ch = Channel.fromPath(file(params.samplesheet, checkIfExists: true))
    .splitCsv(header: true)
    .map { row ->
        tuple(row.sample_id, row.coverage, row.frac, file(row.bam, checkIfExists: true))
    }

workflow {
    // Call the imported process exactly like before
    downsample(input_ch)
}