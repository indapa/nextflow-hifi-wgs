#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { tabix} from './modules/htslib'




def required_params = ['samplesheet', 'tabix_output_dir']
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
        def vcf = file(row.small_variant_vcf)
        if (!vcf.exists()) {
            error "VCF file not found: ${vcf}"
        }
        return tuple(sample_id, vcf)
    }
    .set { input_vcfs_ch }





workflow {
    
    input_vcfs_ch.view { sample_id, vcf -> "Sample ID: $sample_id, VCF: $vcf" }

    
    tabix {
        input_vcfs_ch
    }






}


