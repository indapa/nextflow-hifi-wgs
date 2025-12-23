#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align; pbmm2_align_syt1_region;  hiphase_small_variants } from './modules/pbtools'
include {  deepvariant_targeted_region; deeptrio_targeted_region } from './modules/deepvariant'
include { bam_stats } from './modules/samtools'
include { annotate_vep } from './modules/ensemblvep'

// =========================================================================
//  WORKFLOW 1: ALIGNMENT + DEEPVARIANT TARGETED REGION + HIPHASE + VEP (Entry Point)
// =========================================================================


workflow ALIGN_DEEP_VARIANT_HIPHASE_VEP_SYT1 {
    // ✅ Check samplesheet only in this workflow
    if (!params.samplesheet) {
        error "Parameter 'samplesheet' is required for this workflow!"
    }
    
    if (!file(params.samplesheet).exists()) {
        exit 1, "Samplesheet file not found: ${params.samplesheet}"
    }
    
    // ✅ Create channel only in this workflow
    def input_bams_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample_id
            def bam_file = file(row.bam_file)
            return tuple(sample_id, bam_file)
        }
    
    /* read alignment */
    pbmm2_align_region(
        file(params.reference),
        input_bams_ch,
        params.syt1_region,  // Add the region parameter for SYT1 
        params.cpu,
        params.sort_threads
    )

  

    /* deepvariant targeted region */
    deepvariant_targeted_region(
        file(params.reference), 
        file(params.reference_index), 
        pbmm2_align_syt1_region.out.aligned_bam, 
        params.deepvariant_threads,
        params.syt1_region
        
    )
    // Transform [id, bam, bai] -> [id, bam]
    ch_for_stats = pbmm2_align_syt1_region.out.aligned_bam
        .map { sample_id, bam, bai -> tuple(sample_id, bam) }
        
    bam_stats(ch_for_stats)

 
   
    hiphase_input_ch = deepvariant_targeted_region.out.vcf_tuple
        .join( pbmm2_align_syt1_region.out.aligned_bam )


    // Run HiPhase
    
    hiphase_small_variants(
        hiphase_input_ch,
        params.reference,
        params.reference_index
    )

    annotate_vep (
        hiphase_small_variants.out.phased_vcf,
        params.pigeon_gtf,
        params.pigeon_gtf_tbi,
        params.reference,
        params.reference_index
    )
}
    
// =========================================================================
//  WORKFLOW 2: DEEPTRIO TARGETED REGION: run DeepTrio on trios in a samplesheet with aligned bams
// ========================================================================



workflow RUN_DEEPTRIO {
    // Safety Checks
    if (!params.trio_samplesheet) {
        error "Parameter 'trio_samplesheet' is required for the RUN_DEEPTRIO workflow!"
    }
    
    // 1. Parse CSV and map to [family_id, meta_map]
    ch_samples = channel.fromPath(params.trio_samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.family_id,
                sample_id: row.sample_id,
                bam: file(row.bam),
                bai: file(row.bai),
                role: row.role.toLowerCase()  // normalize role to lowercase
            ]
            return tuple(row.family_id, meta)
        }

    // 2. Group by Family ID -> Returns [family_id, [meta1, meta2, meta3]]
    ch_trios = ch_samples
        .groupTuple(by: 0) 
        .map { family_id, members ->
            
            // Check we have exactly 3 members (optional safety check)
            if (members.size() != 3) {
                log.warn "Family ${family_id} does not have exactly 3 members. Skipping."
                return null 
            }

            // Sort members into variables based on 'role' column
            def child = members.find { it.role == 'child' || it.role == 'proband' }
            def father = members.find { it.role == 'father' || it.role == 'dad' }
            def mother = members.find { it.role == 'mother' || it.role == 'mom' }

            // Ensure all roles were found
            if (!child || !father || !mother) {
                error "Family ${family_id} is missing a role! Found: ${members.role}"
            }

            // Return the specific Tuple structure required by the process
            return tuple(
                family_id,
                child.id, child.bam, child.bai,
                father.id, father.bam, father.bai,
                mother.id, mother.bam, mother.bai
            )
        }

        // view the trio channel

        ch_trios.view()
    

        // 3. Run DeepTrio
    deeptrio_targeted_region(
        file(params.reference),
        file(params.reference_index),
        ch_trios,
        params.deepvariant_threads,
        params.syt1_region  
    )
}



// =========================================================================
//  WORKFLOW 2: HIPHASE + VEP ONLY (Entry Point)
// =========================================================================
workflow HIPHASE_VEP_ONLY {
    
    // Safety Checks
    if (!params.hiphase_samplesheet) {
        error "Parameter 'hiphase_samplesheet' is required for the HIPHASE_ONLY workflow!"
    }

    // Parse HiPhase Samplesheet 
    // EXPECTED COLUMNS: sample_id, vcf, vcf_tbi, bam, bai
    def hiphase_ch = channel.fromPath(params.hiphase_samplesheet)
        .splitCsv(header: true)
        .map { row ->
            return tuple(
                row.sample_id,
                file(row.vcf),
                file(row.vcf_tbi),
                file(row.bam),
                file(row.bai)
            )
        }

    // Run Process
    hiphase_small_variants(
        hiphase_ch,
        params.reference,
        params.reference_index
    )

    annotate_vep (
        hiphase_small_variants.out.phased_vcf,
        params.pigeon_gtf,
        params.pigeon_gtf_tbi,
        params.reference,
        params.reference_index
    )
}

// Default Workflow (runs if you don't specify -entry)
workflow {
    ALIGN_DEEP_VARIANT_HIPHASE_VEP_SYT1()
}

