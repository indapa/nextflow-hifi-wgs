#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_methylation_calling; sawfish_discover; sawfish_joint_call; hiphase_small_variants } from './modules/pbtools'
include { glnexus_trio_merge; deeptrio_wgs } from './modules/deepvariant'
include { bam_stats } from './modules/samtools'
include { annotate_vep } from './modules/ensemblvep'
include { whatshap_trio_phase } from './modules/whatshap'
include { mosdepth_run; infer_sex; plot_dist_coverage } from './modules/mosdepth'
include { samtools_index } from './modules/samtools'


// =========================================================================
//  WORKFLOW: READ ALIGNMENT + POST ALIGNMENT
// =========================================================================

workflow {
    if (params.help) {
        println """
        Available workflows:
        1. DEFAULT: nextflow run main.nf --samplesheet samples.csv performs read alignment and post-alignment analyses (e.g., bam stats)
        2. POST_ALIGNMENT_ONLY: nextflow run main.nf -entry POST_ALIGNMENT_ONLY --samplesheet samples.csv runs post-alignment only on pre-aligned BAMs
        3. RUN_DEEPTRIO: nextflow run main.nf -entry RUN_DEEPTRIO --trio_samplesheet trios.csv performs DeepTrio and WhatsHap phasing on
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
    )

    /* post alignment */
    POST_ALIGNMENT(
        pbmm2_align.out.aligned_bam
    )
}

// =========================================================================
//  WORKFLOW: TRIO ANALYSIS (DeepTrio + WhatsHap phasing); 
//  cpg calling
// sv calling with sawfish2
// =========================================================================

workflow WGS_TRIO {
    if (!params.trio_samplesheet) {
        error "Parameter 'trio_samplesheet' is required for this workflow!"
    }

    if (!file(params.trio_samplesheet).exists()) {
        exit 1, "Trio samplesheet file not found: ${params.trio_samplesheet}"
    }

    def trio_bams_ch = channel.fromPath(params.trio_samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def family_id = row.family_id
            def child_id = row.child_id
            def p1_id = row.parent1_id
            def p2_id = row.parent2_id
            def child_bam = file(row.child_bam)
            def p1_bam = file(row.parent1_bam)
            def p2_bam = file(row.parent2_bam)
            return tuple(family_id, child_id, child_bam, p1_id, p1_bam, p2_id, p2_bam)
        }

         // take the id and uBAMs and make channel for pbmm2_align
        flat_bams_ch = trio_bams_ch
        .flatMap { family_id, child_id, child_bam, p1_id, p1_bam, p2_id, p2_bam ->
            [
                tuple(child_id, child_bam),
                tuple(p1_id, p1_bam),
                tuple(p2_id, p2_bam)
            ]
        }
        .unique()


   

    // 3. Run Alignment
    pbmm2_align(file(params.reference), flat_bams_ch)

    // pbmm2_align.out.aligned_bam emits: [ sample_id, aligned_bam, aligned_bai ]
    aligned_bams_ch = pbmm2_align.out.aligned_bam

    // 4. Reconstruct the Trio structure using cross-referencing joins
    // We join the original trio channel sequentially on child, parent1, and parent2 IDs
    deeptrio_input_ch = trio_bams_ch
        // Join for Child: matches trio_bams_ch(child_id) with aligned_bams_ch(sample_id)
        .map { fam, c_id, c_b, p1_id, p1_b, p2_id, p2_b -> tuple(c_id, fam, p1_id, p1_b, p2_id, p2_b) }
        .join(aligned_bams_ch) // yields: [ c_id, fam, p1_id, p1_b, p2_id, p2_b, c_bam, c_bai ]
        
        // Join for Parent 1
        .map { c_id, fam, p1_id, p1_b, p2_id, p2_b, c_bam, c_bai -> tuple(p1_id, fam, c_id, c_bam, c_bai, p2_id, p2_b) }
        .join(aligned_bams_ch) // yields: [ p1_id, fam, c_id, c_bam, c_bai, p2_id, p2_b, p1_bam, p1_bai ]
        
        // Join for Parent 2
        .map { p1_id, fam, c_id, c_bam, c_bai, p2_id, p2_b, p1_bam, p1_bai -> tuple(p2_id, fam, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai) }
        .join(aligned_bams_ch) // yields: [ p2_id, fam, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_bam, p2_bai ]
        
        // Final map to restore your exact process input order
        .map { p2_id, fam, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_bam, p2_bai ->
            tuple(fam, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai)
        }

    // 5. Run DeepTrio
    deeptrio_wgs(params.reference, params.reference_index, deeptrio_input_ch)

    // GLNexus merge
    // Reconstruct the merged input tuple grouped by family_id
    glnexus_input_ch = deeptrio_wgs.out.child_gvcf
        // Join Child with Parent 1 on family_id (the first element)
        .join(deeptrio_wgs.out.p1_gvcf) 
        // Current state: [fam, child_id, child_gvcf, child_tbi, p1_id, p1_gvcf, p1_tbi]
        
        // Join the result with Parent 2 on family_id
        .join(deeptrio_wgs.out.p2_gvcf)
        // Current state: [fam, child_id, child_gvcf, child_tbi, p1_id, p1_gvcf, p1_tbi, p2_id, p2_gvcf, p2_tbi]

    // Run GLnexus
    glnexus_trio_merge(glnexus_input_ch)

    



    // Assemble Whatshap Input
    // glnexus_trio_merge.out.joint_vcf is: [family_id, joint_vcf, joint_vcf_tbi]
    // deeptrio_input_ch is: [family_id, child_id, child_bam, child_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai]
    
    whatshap_input_ch = glnexus_trio_merge.out.joint_vcf
        .join(deeptrio_input_ch)
        // This yields exactly the required structure:
        // [family_id, vcf, vcf_tbi, child_id, child_bam, child_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai]

    // Run Whatshap Phase
    whatshap_trio_phase(
        params.reference, 
        params.reference_index, 
        whatshap_input_ch
    )

    // 1. Gather all parent VCFs into a single channel and format to [sample_id, vcf, vcf_tbi]
    // deeptrio_wgs.out.p1_vcf / p2_vcf are: [family_id, parent_id, vcf, vcf_tbi]
    parent_vcfs_ch = deeptrio_wgs.out.p1_vcf
        .mix(deeptrio_wgs.out.p2_vcf)
        .map { family_id, parent_id, vcf, vcf_tbi -> 
            tuple(parent_id, vcf, vcf_tbi) 
        }

    // 2. Join the parent VCFs with their corresponding individual aligned BAMs
    // aligned_bams_ch is: [sample_id, bam, bai]
    hiphase_input_ch = parent_vcfs_ch
        .join(aligned_bams_ch)
        // Yields: [parent_id, vcf, vcf_tbi, bam, bai]

    // 3. Run HiPhase for parents (Arguments must match input declaration order)
    hiphase_small_variants(
        hiphase_input_ch,        // 1st input: tuple
        file(params.reference),        // 2nd input: path reference
        file(params.reference_index)   // 3rd input: path reference_index
    )

    // --- Prepare Child Haplotagged BAM ---
    // whatshap_trio_phase outputs: child_id via deeptrio_input_ch
    // We map whatshap's output back to its child_id using a join or map.
    // The easiest way is to extract child_id from trio_bams_ch or deeptrio_input_ch:
    child_id_ch = trio_bams_ch.map { fam, child_id, c_b, p1_id, p1_b, p2_id, p2_b -> tuple(fam, child_id) }
    
    child_bam_ch = whatshap_trio_phase.out.haplotagged_bam
        // whatshap outputs a single file, so we pair it with its family ID via channel matching or join
        // Since whatshap_trio_phase.out.phased_vcf has [family_id, vcf], we can use that to track family context:
        .map { bam -> tuple(bam.baseName.replaceAll(/\.haplotagged/, ''), bam) } // maps filename to child_id
    
    // Run samtools index on the child bam
    samtools_index(child_bam_ch)
    child_cpg_input = samtools_index.out
    // --- Prepare Parents Haplotagged BAM ---
    // hiphase_small_variants.out.haplotagged_bam already matches: [sample_id, bam, bai]

    // --- Combine All Family Members Together ---
    all_haplotagged_bams_ch = child_cpg_input.mix(parent_cpg_input)

    // --- Run CpG Methylation Calling ---
    cpg_methylation_calling(
        all_haplotagged_bams_ch,
        file(params.reference),
        file(params.reference_index)
    )

    mosdepth_run(all_haplotagged_bams_ch)
    infer_sex(mosdepth_run.out.summary)
    plot_dist_coverage(mosdepth_run.out.global_dist)

    expected_bed_ch = infer_sex.out.sex.map { sample_id, sex_csv ->
        def lines = sex_csv.readLines()
        def sex = lines.size() > 1 ? lines[1].split(',')[3].trim() : 'UNKNOWN'

        def expected_bed
        if (sex == 'FEMALE') {
            expected_bed = file(params.expected_XX_bed)
        } else if (sex == 'MALE') {
            expected_bed = file(params.expected_XY_bed)
        } else {
            throw new Exception("Error: Invalid or missing sex '${sex}' inferred for sample ${sample_id}. Expected 'FEMALE' or 'MALE'.")
        }

        return tuple(sample_id, expected_bed)
    }

    sawfish_in_ch = aligned_bams_ch.join(expected_bed_ch, by: 0)

    awfish_discover(
        sawfish_in_ch,
        file(params.excluded_bed),
        file(params.reference),
        file(params.reference_index)
    )

    sawfish_joint_call(
        sawfish_discover.out.discover_dir.collect()
        
    )
    
}


// =========================================================================
//  ENTRY POINT: POST-ALIGNMENT ONLY (skip alignment)
// =========================================================================

workflow POST_ALIGNMENT_ONLY {
    if (!params.samplesheet) {
        error "Parameter 'samplesheet' is required! CSV must have columns: sample_id, bam_file, bai_file"
    }

    if (!file(params.samplesheet).exists()) {
        exit 1, "Samplesheet file not found: ${params.samplesheet}"
    }

    def aligned_bam_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def bam = file(row.bam_file)
            def bai = file(row.bai_file)
            return tuple(sample_id, bam, bai)
        }

    POST_ALIGNMENT(aligned_bam_ch)
}


// =========================================================================
//  SUB-WORKFLOW: POST ALIGNMENT
// =========================================================================

workflow POST_ALIGNMENT {
    take:
    aligned_bam_ch

    main:
    
    bam_stats(aligned_bam_ch)

    deepvariant_wgs(
        file(params.reference),
        file(params.reference_index),
        aligned_bam_ch
    )

    //join aligned_bam_ch with the output of deepvariant_wgs by sample_id

    aligned_bam_ch.join(deepvariant_wgs.out.vcf_tuple, by: 0)
        .map { sample_id, bam, bai, vcf, vcf_tbi ->
            tuple(sample_id, vcf, vcf_tbi, bam, bai)
        }
        .set { aligned_bam_with_vcf_ch }


    hiphase_small_variants(
        aligned_bam_with_vcf_ch,
        file(params.reference),
        file(params.reference_index)
    )    

    cpg_methylation_calling(
        hiphase_small_variants.out.haplotagged_bam,
        file(params.reference),
        file(params.reference_index)
    )
    
    mosdepth_run(aligned_bam_ch)
    

    infer_sex(mosdepth_run.out.summary)

    plot_dist_coverage(mosdepth_run.out.global_dist)

    // Determine expected BED file based on inferred sex
    expected_bed_ch = infer_sex.out.sex.map { sample_id, sex_csv ->
        def lines = sex_csv.readLines()
        def sex = lines.size() > 1 ? lines[1].split(',')[3].trim() : 'UNKNOWN'

        def expected_bed
        if (sex == 'FEMALE') {
            expected_bed = file(params.expected_XX_bed)
        } else if (sex == 'MALE') {
            expected_bed = file(params.expected_XY_bed)
        } else {
            throw new Exception("Error: Invalid or missing sex '${sex}' inferred for sample ${sample_id}. Expected 'FEMALE' or 'MALE'.")
        }

        return tuple(sample_id, expected_bed)
    }

    // Join BAM channel with the correct BED file channel by sample_id
    // Creates: [sample_id, bam, bai, expected_bed]
     sawfish_in_ch = aligned_bam_ch.join(expected_bed_ch, by: 0)

    sawfish_discover(
        sawfish_in_ch,
        file(params.excluded_bed),
        file(params.reference),
        file(params.reference_index)
    )

    sawfish_joint_call(
        sawfish_discover.out.discover_dir.collect(),
        
    )
}










