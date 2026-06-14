#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

// Include modules
include { bam_stats } from './modules/samtools'

// Initialize parameters
params.samplesheet = ''
params.bam_stats_output_dir = ''

// Validate required parameters
def required_params = ['samplesheet', 'bam_stats_output_dir']
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}

// Safely instantiate the samplesheet file (checks existence implicitly if requested)
def csv_file = file(params.samplesheet, checkIfExists: true)

// Create channels for input BAM files
Channel.fromPath(csv_file)
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample_id
        def bam = file(row.bam_file, checkIfExists: true) // Automatically errors if missing
        return tuple(sample_id, bam)
    }
    .set { input_bam_ch }

workflow {
    // Corrected view syntax (uncomment to debug)
    // input_bam_ch.view { sample_id, bam -> "ID: $sample_id BAM: $bam" }

    // Invoke the process
    bam_stats_ch = bam_stats(input_bam_ch)
}