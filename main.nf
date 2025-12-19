#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { pbmm2_align; hiphase_small_variants } from './modules/pbtools'
include { deepvariant; deepvariant_chr20; BCFTOOLS_STATS; bcftools_deepvariant_norm; deepvariant_targeted_region} from './modules/deepvariant'

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
def input_bams_ch = Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample_id
        def bam_file = file(row.bam_file)
        if (!bam_file.exists()) {
            error "BAM file not found: ${bam_file}"
        }
        return tuple(sample_id, bam_file)
    }

Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> row.sample_id }
    .set { sample_ids_ch }

def REGIONS = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
    'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
]

// Create a channel from the fixed regions
Channel
    .fromList(REGIONS)
    .set { regions_ch }

// Named workflow for post-alignment analysis
workflow POST_ALIGNMENT_ANALYSIS {
    take:
    aligned_bam_ch  // Channel: tuple(bam, bai)
    sample_ids_ch   // Channel: val(sample_id)

    main:
   

    /* aligned bam channel used for cnv, tandem repeat and sv analysis */
    bam_bai_ch = aligned_bam_ch.map { bam, bai -> 
        def sample_id = bam.baseName.replaceFirst(/\..*$/, '')
        tuple(sample_id, bam, bai)
    }

     /* read depth analysis */
    mosdepth(bam_bai_ch)

    /* cnv analysis */
    hificnv(
        bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.exclude_bed),
        file(params.expected_bed),
        params.cpu
    )

    /* tandem repeat analysis */
    trgt(
        bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.trgt_repeats),
        params.karyotype,
        params.cpu
    )

    /* SV analysis - regions, discover, call */
    bam_bai_ch
        .combine(regions_ch)
        .map { sample_id, bam, bai, region -> 
            [sample_id, region, bam, bai]
        }
        .set { bam_regions_ch }

    pb_discover_results = pb_discover(bam_regions_ch, params.trf_bed)

    // Group svsig files by sample_id
    svsig_files_by_sample = pb_discover_results.groupTuple()

    // Run pb_call process
    pb_call(svsig_files_by_sample, file(params.reference))

    /* deepvariant */
    deepvariant(params.reference, params.reference_index, bam_bai_ch, params.deepvariant_threads)

    /* bcftools normalization */
    bcftools_deepvariant_norm(params.reference, deepvariant.out.vcf_tuple)

    // Run bcftools stats
    BCFTOOLS_STATS(bcftools_deepvariant_norm.out.vcf_tuple)

    /* run cpg pileup */

    cpg_pileup(aligned_bam_ch)
    split_bed_by_chrom(cpg_pileup.out.pileup_beds)

    // Extract the cpg results to parquet format
    cpg_tools_extract_to_parquet(cpg_pileup.out.pileup_beds)

   

    

    
    /* vep annotate */
    annotate_vep_no_phased(
        bcftools_deepvariant_norm.out.vcf_tuple,
        file(params.pigeon_gtf_bgzip),
        file(params.pigeon_gtf_tbi),
        file(params.reference)
    )

    bam_stats(bam_bai_ch)

    //PARSE_SAMTOOLS_BAM_STATS(bam_stats.out.stats)
}

// Main workflow - full pipeline
workflow {
    /* read alignment */
    pbmm2_align(
        file(params.reference),
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )

    // Run post-alignment analysis
    POST_ALIGNMENT_ANALYSIS(pbmm2_align.out.aligned_bam, sample_ids_ch)
}

// Named workflow for starting from aligned BAMs
workflow FROM_ALIGNED_BAMS {
    // Create channel for pre-aligned BAM files
    // Assumes samplesheet has aligned BAM files instead of raw reads
    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def bam_file = file(row.bam_file)
            def bai_file = file(row.bai_file) // Assumes BAI file is in samplesheet
            if (!bam_file.exists()) {
                error "BAM file not found: ${bam_file}"
            }
            if (!bai_file.exists()) {
                error "BAI file not found: ${bai_file}"
            }
            return tuple(bam_file, bai_file)
        }
        .set { aligned_bam_ch }

    // Run post-alignment analysis
    POST_ALIGNMENT_ANALYSIS(aligned_bam_ch, sample_ids_ch)
}

// Entry point 3: Alignment only (pbmm2_align only)
workflow ALIGN_ONLY {
    /* read alignment */
    pbmm2_align(
        file(params.reference),
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )

     // create bam_bai_ch for deepvariant
    bam_bai_ch_align_only = pbmm2_align.out.aligned_bam.map { bam, bai -> 
        def sample_id = bam.baseName.replaceFirst(/\..*$/, '')
        tuple(sample_id, bam, bai)
    }

    /* call bam stats */
    bam_stats(bam_bai_ch_align_only)

      /* read depth analysis */
    mosdepth(bam_bai_ch_align_only)

    //PARSE_SAMTOOLS_BAM_STATS(bam_stats.out.stats)

}




//entry point 4: DeepVariant only from aligned BAMs and BCFtools stats
// start with aligned BAMs from samplesheet
workflow ALIGN_DEEP_VARIANT_BCFTOOLS_STATS {
    /* read alignment */
    pbmm2_align(
        file(params.reference),
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )

   

    /* deepvariant */
    deepvariant(params.reference, params.reference_index, bam_bai_ch, params.deepvariant_threads)

    /* bcftools normalization */
    bcftools_deepvariant_norm(params.reference, deepvariant.out.vcf_tuple)
     
    // Run bcftools stats
    BCFTOOLS_STATS(bcftools_deepvariant_norm.out.vcf_tuple)
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
        params.reference, 
        params.reference_index, 
        pbmm2_align.out.aligned_bam, 
        params.deepvariant_threads, 
        params.syt1_region
    )

    /* bcftools normalization */
    bcftools_deepvariant_norm(
        params.reference, 
        deepvariant_targeted_region.out.vcf_tuple
    )
}


workflow DEEPVARIANT_ONLY {
 
      Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def bam_file = file(row.bam_file)
            def bai_file = file(row.bai_file)
            if (!bam_file.exists()) {
                error "BAM file not found: ${bam_file}"
            }
            if (!bai_file.exists()) {
                error "BAI file not found: ${bai_file}"
            }
            return tuple(sample_id, bam_file, bai_file)
        }
        .set { aligned_bam_ch }

    /* deepvariant starting from samplesheet of aligned bams */
    deepvariant(params.reference, params.reference_index, aligned_bam_ch, params.deepvariant_threads)

    /* bcftools normalization */
    bcftools_deepvariant_norm(params.reference, deepvariant.out.vcf_tuple)

    // Run bcftools stats
    BCFTOOLS_STATS(bcftools_deepvariant_norm.out.vcf_tuple)
}


