#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { downsample } from './modules/samtools'

// Parameters
params.samplesheet = ''
params.downsample_output_dir = ''

// Parameter validation
def required_params = ['samplesheet', 'downsample_output_dir']
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}



// tuple val(sample_id), val(coverage), val(frac), path(bam)
input_ch = Channel.fromPath(file(params.samplesheet, checkIfExists: true))
    .splitCsv(header: true)
    .map { row ->
        tuple(row.sample_id, row.coverage, row.frac, file(row.bam, checkIfExists: true))
    }

workflow {
    // Call the imported process exactly like before
    downsample(input_ch)
}