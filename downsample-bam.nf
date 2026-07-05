include { downsample } from './modules/samtools'

// Parameters
params.samplesheet = ''
params.downsample_output_dir = ''

// Parameter validation


workflow {
    def required_params = ['samplesheet', 'downsample_output_dir']
    required_params.each { param ->
        if (!params[param]) {
            error "Parameter '$param' is required!"
        }
    }
    input_ch = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, row.coverage, row.frac, file(row.bam, checkIfExists: true))
        }

    downsample(input_ch)
}

