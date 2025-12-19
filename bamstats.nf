#!/usr/local/bin/nextflow

nextflow.enable.dsl=2



include {bam_stats} from './modules/samtools'



def required_params = ['samplesheet', 'aligned_output_dir']
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}


// Replace the input_bams channel definition with this:
def checkSamplesheet(samplesheet_file) {
    if (!file(samplesheet_file).exists()) {
        exit 1, "Samplesheet file not found: ${samplesheet_file}"
    }
    return file(samplesheet_file)
}

ss_status = checkSamplesheet(params.samplesheet)

// Create channels for input BAM files
Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample_id
        def bam = file(row.bam)
       
        if (!bam.exists()) {
            error "bam file not found: ${bam}"
        }
        return tuple(sample_id, bam)
    }
    .set { input_bam_ch }





workflow {
    
    //input_bam_ch.view {  bam, bam_index -> "BAM: $bam BAI: $bam_index" }

    bam_stats_ch = bam_stats (input_bam_ch)

    
 



}