#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_pileup; hificnv; trgt; pb_discover; pb_call; hiphase } from './modules/pbtools'
include { mosdepth } from './modules/mosdepth'
include { deepvariant; deepvariant_chr20 } from './modules/deepvariant'
include {fibertools_extract} from './modules/fiberseq'
include {annotate_vep} from './modules/ensemblvep'




def required_params = ['reference', 'samplesheet', 'cpgmodel', 'karyotype', 'sv_output_dir']
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
        def bam_file = file(row.bam_file)
        if (!bam_file.exists()) {
            error "BAM file not found: ${bam_file}"
        }
        return tuple(sample_id, bam_file)
    }
    .set { input_bams_ch }



def REGIONS = [
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'chrX'
   
]

// Create a channel from the fixed regions
Channel
    .fromList(REGIONS)
    .set { regions_ch }

workflow {
    
    input_bams_ch.view { sample_id, bam -> "Sample ID: $sample_id, BAM: $bam" }

    regions_ch.view { region -> "BED region: $region" } 

/* read alignment */
    pbmm2_align(
        file(params.reference),
        //input_bams,
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )
    // print the pbmm2_align channel

    pbmm2_align.out.aligned_bam.view { "Aligned BAM: $it" }

 /* cpg calling  */
    cpg_pileup(
        pbmm2_align.out.aligned_bam,
        file(params.reference),
        file(params.reference_index),
        file(params.cpgmodel)
    )

/* read depth analysis  */
    mosdepth(pbmm2_align.out.aligned_bam)





/* aligned bam channel used for cnv, tandem repeat and sv analysis */
    bam_bai_ch = pbmm2_align.out.aligned_bam.map { bam, bai -> 
    def sample_id = bam.baseName.replaceFirst(/\..*$/, '')
    tuple(sample_id, bam, bai)
}

    bam_bai_ch.view { sample_id, bam, bai -> "Input: $sample_id, BAM: $bam, BAI: $bai" }

/* cnv analysis */
    hificnv(
        bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.exclude_bed),
        file(params.expected_bed),
        params.cpu
    )


/*tandem repeat analysis  */
    trgt(   bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.trgt_repeats),
        params.karyotype,
        params.cpu
    )




/* SV analysis
        regions
        discover
        call 
*/
  
       bam_bai_ch
        .combine(regions_ch)
        .map { sample_id, bam, bai, region -> 
            // For each region, create the tuple: [sample_id, region, bam, bai]
            [sample_id, region, bam, bai]
        }
        .set { bam_regions_ch }

        
        bam_regions_ch.view { sample_id, region, bam, bai -> "sample: $sample_id, bed: $region, BAM: $bam, BAI: $bai"}

        pb_discover_results=pb_discover(bam_regions_ch, params.trf_bed )

         // Group svsig files by sample_id
        svsig_files_by_sample = pb_discover_results.groupTuple()

        // Create a channel for the reference genome
        reference_ch = channel.fromPath(params.reference)

        // Run pb_call process
        pb_call(svsig_files_by_sample, reference_ch)

        

        //deepvariant

        deepvariant(params.reference, params.reference_index, bam_bai_ch, params.deepvariant_threads )
        //deepvariant_chr20(params.reference, params.reference_index, bam_bai_ch, params.deepvariant_threads )

    


        //hiphase

        hiphase(
            deepvariant.out.vcf,
            deepvariant.out.vcf_tbi,
            pb_call.out.pb_call,
            trgt.out.repeat_vcf,
            pbmm2_align.out.aligned_bam,
            file(params.reference)
        )


        //fibertools extract 

        fibertools_extract (
            bam_bai_ch
        )

        //vep annotate

        annotate_vep(
            hiphase.out.phased_deepvariant,
            file(params.pigeon_gtf),
            file(params.pigeon_tbi),
            file(params.reference)
        )

}


