#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align } from './modules/pbtools'
include {  deepvariant_targeted_region; deeptrio_targeted_region; glnexus_trio_merge } from './modules/deepvariant'
include { bam_stats } from './modules/samtools'
include { annotate_vep } from './modules/ensemblvep'
include { whatshap_trio_phase } from './modules/whatshap'



// =========================================================================
//  WORKFLOW READ ALIGNMENT + POST ALIGNMENT 
// =========================================================================




workflow  {
    // ✅ Check samplesheet only in this workflow

    if (params.help) {
        println """
        Available workflows:
        1. DEFAULT: nextflow run main.nf --samplesheet samples.csv performs read alignment and post-alignment analyses (e.g., bam stats)
        2. RUN_DEEPTRIO: nextflow run main.nf -entry RUN_DEEPTRIO --trio_samplesheet trios.csv performs DeepTrio and WhatsHap phasing on
        trios in the samplesheet with aligned bams
        
        """.stripIndent()
        exit 0
    }
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
    pbmm2_align(
        file(params.reference),
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )

    /* post alignment */
    POST_ALIGNMENT(
        pbmm2_align.out.aligned_bam
    )

}

workflow POST_ALIGNMENT {
    // This can include processes like variant calling, phasing, annotation, etc.

    take:
    aligned_bam_ch

    main:
    bam_stats(
        aligned_bam_ch
    )

}


/* =========================================================================
WORKFLOW 2: DEEPTRIO TARGETED REGION: run DeepTrio on trios in a samplesheet with aligned bams

1. Parse samplesheet to group by family_id
2. For each family, run deeptrio_targeted_region process
3. Merge gVCFs with GLnexus
4. Phase variants with Whatshap using trio information


========================================================================

*/ 

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
                id: row.sample_id,      // <--- Make sure this is row.sample_id, NOT row.family_id
                bam: file(row.bam),
                bai: file(row.bai),
                role: row.role.toLowerCase()
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

        //ch_trios.view()
    

        // 3. Run DeepTrio
    
    deeptrio_targeted_region(
        file(params.reference),
        file(params.reference_index),
        ch_trios,
        params.deepvariant_threads,
        params.syt1_region  
    )

    ch_glnexus_input = deeptrio_targeted_region.out.child_gvcf
        .join(deeptrio_targeted_region.out.p1_gvcf)
        .join(deeptrio_targeted_region.out.p2_gvcf)

    // view ch_glnexus_input
    //ch_glnexus_input.view()

    glnexus_trio_merge(ch_glnexus_input)

    ch_phasing_input = glnexus_trio_merge.out.joint_vcf
        .join(ch_trios)
    ch_phasing_input.view()

    whatshap_trio_phase(
        file(params.reference),
        file(params.reference_index),
        ch_phasing_input
    )


    
}







