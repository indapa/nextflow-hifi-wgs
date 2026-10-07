#!/usr/bin/env nextflow



include { cpg_methylation_calling; sawfish_discover; sawfish_joint_call; hiphase_small_variants } from './modules/pbtools'

include { bam_stats; samtools_index } from './modules/samtools'

include { WHATSHAP_TRIO_PHASE_BY_CHROM } from './subworkflows/whatshap_trio_phase_by_chrom'
include { CONCAT_AND_SPLIT_WGS } from './subworkflows/concat_and_split_wgs'
include { GLNEXUS_TRIO }         from './subworkflows/glnexus_trio_merge'
include { PBMM2_SPOT_WGS } from './subworkflows/reference_alignment'
include {FASTVEP_ANNOTATE_WGS} from './subworkflows/fastvep'
include {DEEPVARIANT_SINGLETON_WGS; DEEPVARIANT_SINGLETON_WGS_PARABRICKS} from './subworkflows/deepvariant_singleton_wgs'
include { DEEPTRIO_WGS } from './subworkflows/deeptrio_wgs'
include { SAWFISH_TRIO } from './subworkflows/sawfish_trio'
include { HIPHASE_TRIO_PARENTS } from './subworkflows/hiphase_trio_parents'
include { CPG_TRIO_METHYLATION } from './subworkflows/cpg_trio_methylation'
include { FASTVEP_ANNOTATE_TRIO_INDIVIDUALS } from './subworkflows/fastvep_trio_individuals'
include { TRGT_GENOTYPING } from './subworkflows/trgt'

include { mosdepth_run; infer_sex; plot_dist_coverage } from './modules/mosdepth'




// =========================================================================
//  WORKFLOW: READ ALIGNMENT + POST ALIGNMENT (SINGLETONS)
// =========================================================================

workflow {
    if (params.help) {
        println """
        Available workflows:
        1. DEFAULT: nextflow run main.nf --samplesheet samples.csv 
           Performs read alignment and post-alignment analyses on singletons.
        2. POST_ALIGNMENT_ONLY: nextflow run main.nf -entry POST_ALIGNMENT_ONLY --samplesheet samples.csv 
           Runs post-alignment singletons analyses on pre-aligned BAMs.
        3. WGS_TRIO: nextflow run main.nf -entry WGS_TRIO --trio_samplesheet trios.csv 
           Performs alignment, DeepTrio, Phasing, CpG, and SV calls on unaligned trio inputs.
        4. WGS_TRIO_ALIGNED: nextflow run main.nf -entry WGS_TRIO_ALIGNED --trio_aligned_samplesheet trios.csv
           Performs DeepTrio and downstream pipelines on pre-aligned trio inputs.
        5. TRIO_FROM_DEEPTRIO: nextflow run main.nf -entry TRIO_FROM_DEEPTRIO --trio_gvcf_samplesheet trios.csv
           Runs GLnexus joint calling, WhatsHap trio phasing, FastVEP annotation, HiPhase parent phasing,
           and CpG methylation on existing DeepTrio VCFs/gVCFs + pre-aligned BAMs.
        """.stripIndent()
        exit 0
    }

    if (params.entry == 'WGS_TRIO_ALIGNED') {
        WGS_TRIO_ALIGNED()
    
    }
    
    else if ( params.entry == 'WGS_TRIO') {
        WGS_TRIO()
    }
    else if (params.entry == 'WGS_SINGLETON') {
        WGS_SINGLETON()
    }
    else if (params.entry == 'TRIO_FROM_DEEPTRIO') {
        TRIO_FROM_DEEPTRIO()
    }
    else if (params.entry == 'POST_ALIGNMENT_ONLY') {
        POST_ALIGNMENT_ONLY()

    }
    else {
        WGS_SINGLETON()
    }
}
// =========================================================================
//  WORKFLOW: TRIO ANALYSIS ENTRYPOINTS
// =========================================================================

// --- Entrypoint 1: Starts from Raw Unaligned BAMs ---

workflow WGS_SINGLETON {
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
    PBMM2_SPOT_WGS(
        file(params.reference),
        input_bams_ch
    )

    /* post alignment */
    POST_ALIGNMENT(
        PBMM2_SPOT_WGS.out
    )

}

workflow WGS_TRIO {
    

    if (!file(params.trio_samplesheet).exists()) {
        exit 1, "Trio samplesheet file not found: ${params.trio_samplesheet}"
    }

    raw_samples_ch = channel.fromPath(params.trio_samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            tuple(row.family_id, row.sample_id, row.role, file(row.bam)) 
        }

    align_input_ch = raw_samples_ch.map { _fam, sample_id, _role, bam -> tuple(sample_id, bam) }
    PBMM2_SPOT_WGS(
        file(params.mmi),
        align_input_ch
    )
    trio_bams_assembled = raw_samples_ch
        .map { fam, sample_id, role, _raw_bam -> tuple(sample_id, fam, role) }
        .join(PBMM2_SPOT_WGS.out)
        .map { sample_id, fam, role, bam, bai -> 
            tuple(fam, [role: role, id: sample_id, bam: bam, bai: bai]) 
        }
        .groupTuple(by: 0)
        .map { fam, members ->
            def c  = members.find { m -> m.role == 'child' }
            def p1 = members.find { m -> m.role == 'parent1' }
            def p2 = members.find { m -> m.role == 'parent2' }

        return tuple(fam, c.id, c.bam, c.bai, p1.id, p1.bam, p1.bai, p2.id, p2.bam, p2.bai)
    }


    sample_roles_ch = raw_samples_ch
        .map { _fam, sample_id, role, _bam -> tuple(sample_id, role) }

    sample_to_family_ch = raw_samples_ch
        .map { fam, sample_id, _role, _bam -> tuple(sample_id, fam) }

    // Isolate single aligned BAM trackers for downstream tools (Sawfish/HiPhase)
    individual_aligned_bams = PBMM2_SPOT_WGS.out

    RUN_TRIO_PIPELINE(trio_bams_assembled, individual_aligned_bams, sample_roles_ch, sample_to_family_ch)
}

// --- Entrypoint 2: Starts from Pre-Aligned BAMs ---
workflow WGS_TRIO_ALIGNED {
    

    if (!file(params.trio_aligned_samplesheet).exists()) {
        exit 1, "Aligned samplesheet file not found: ${params.trio_aligned_samplesheet}"
    }

    trio_bams_assembled = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.family_id, [
                role: row.role, 
                id: row.sample_id, 
                bam: file(row.aligned_bam), 
                bai: file(row.aligned_bai)
            ])
        }
        .groupTuple(by: 0)
        .map { fam, members ->
            def c  = members.find { m -> m.role == 'child' }
            def p1 = members.find { m -> m.role == 'parent1' }
            def p2 = members.find { m -> m.role == 'parent2' }

        
            return tuple(fam, c.id, c.bam, c.bai, p1.id, p1.bam, p1.bai, p2.id, p2.bam, p2.bai)
        }

    sample_roles_ch = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.role) }

    sample_to_family_ch = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.family_id) }

    // Reconstruct flat stream of individual aligned BAMs for downstream hooks
    individual_aligned_bams = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.aligned_bam), file(row.aligned_bai)) }

    RUN_TRIO_PIPELINE(trio_bams_assembled, individual_aligned_bams, sample_roles_ch, sample_to_family_ch)
}

// --- Entrypoint 3: Starts from existing DeepTrio VCFs/gVCFs + pre-aligned BAMs ---
// Runs GLnexus joint calling -> WhatsHap trio phasing -> FastVEP annotation,
// plus HiPhase read-backed phasing of the parents -> CpG methylation (child + parents).
// Samplesheet columns: family_id, sample_id, role, aligned_bam, aligned_bai, vcf, vcf_tbi, gvcf, gvcf_tbi
workflow TRIO_FROM_DEEPTRIO {

    if (!params.trio_gvcf_samplesheet || !file(params.trio_gvcf_samplesheet).exists()) {
        exit 1, "Trio gVCF samplesheet file not found: ${params.trio_gvcf_samplesheet}"
    }

    samples_ch = channel.fromPath(params.trio_gvcf_samplesheet)
        .splitCsv(header: true)

    trio_bams_assembled = samples_ch
        .map { row ->
            tuple(row.family_id, [
                role: row.role,
                id: row.sample_id,
                bam: file(row.aligned_bam),
                bai: file(row.aligned_bai)
            ])
        }
        .groupTuple(by: 0)
        .map { fam, members ->
            def c  = members.find { m -> m.role == 'child' }
            def p1 = members.find { m -> m.role == 'parent1' }
            def p2 = members.find { m -> m.role == 'parent2' }

            return tuple(fam, c.id, c.bam, c.bai, p1.id, p1.bam, p1.bai, p2.id, p2.bam, p2.bai)
        }

    sample_roles_ch = samples_ch
        .map { row -> tuple(row.sample_id, row.role) }

    individual_aligned_bams = samples_ch
        .map { row -> tuple(row.sample_id, file(row.aligned_bam), file(row.aligned_bai)) }

    deeptrio_vcf_ch = samples_ch
        .map { row -> tuple(row.family_id, row.sample_id, file(row.vcf), file(row.vcf_tbi)) }

    deeptrio_gvcf_ch = samples_ch
        .map { row -> tuple(row.family_id, row.sample_id, file(row.gvcf), file(row.gvcf_tbi)) }

    GLNEXUS_TRIO(
        deeptrio_gvcf_ch,
        sample_roles_ch,
        channel.fromList(params.chromosomes),
        file(params.glnexus_region_bed)
    )

    WHATSHAP_TRIO_PHASE_BY_CHROM(
        GLNEXUS_TRIO.out,
        trio_bams_assembled,
        file(params.reference),
        file(params.reference_index),
        channel.fromList(params.chromosomes)
    )

    FASTVEP_ANNOTATE_TRIO_INDIVIDUALS(
        WHATSHAP_TRIO_PHASE_BY_CHROM.out.phased_vcf,
        trio_bams_assembled,
        file(params.fastvep_gff),
        file(params.reference),
        channel.fromPath("${params.fastvep_sa_dir}/*").collect()
    )

    // Read-backed phasing + haplotagging for each parent independently
    HIPHASE_TRIO_PARENTS(
        deeptrio_vcf_ch,
        individual_aligned_bams,
        sample_roles_ch,
        file(params.reference),
        file(params.reference_index)
    )

    // CpG methylation calling on the haplotagged child + parent BAMs
    CPG_TRIO_METHYLATION(
        WHATSHAP_TRIO_PHASE_BY_CHROM.out.haplotagged_bam,
        HIPHASE_TRIO_PARENTS.out.haplotagged_bam,
        file(params.reference),
        file(params.reference_index)
    )
}

// =========================================================================
//  SUB-WORKFLOW: SHARED TRIO DOWNSTREAM ENGINE
// =========================================================================

workflow RUN_TRIO_PIPELINE {
    take:
    trio_bams_assembled
    individual_aligned_bams
    sample_roles_ch
    sample_to_family_ch

    main:

    // 1. Load all interval BED files from the specified directory
    raw_intervals_ch = channel.fromPath("${params.intervals_dir}/*.bed")

    // 2. Filter to autosomes + X/Y chunks
    intervals_ch = raw_intervals_ch.filter { f -> f.baseName =~ /^chr([1-9]|1[0-9]|2[0-2]|[XY])_/ }

    // =========================================================================
    // Pre-calculate how many chunks exist per chromosome so groupTuple can
    // emit eagerly via groupKey without waiting for the entire channel to close.
    // =========================================================================
   

    // Run DeepTrio end-to-end: 50MB scatter -> per-chromosome merge -> genome-wide merge
    DEEPTRIO_WGS(
        file(params.reference),
        file(params.reference_index),
        trio_bams_assembled,
        intervals_ch
    )

    deeptrio_vcf_ch  = DEEPTRIO_WGS.out.vcf   // tuple(family_id, sample_id, vcf, tbi)
    deeptrio_gvcf_ch = DEEPTRIO_WGS.out.gvcf  // tuple(family_id, sample_id, gvcf, gtbi)

    // Joint-call each family per chromosome, then merge into one VCF per family
    ch_chroms_glnexus = channel.fromList(params.chromosomes)

    GLNEXUS_TRIO(
        deeptrio_gvcf_ch,
        sample_roles_ch,
        ch_chroms_glnexus,
        file(params.glnexus_region_bed)
    )

    // Phase each family's joint-called VCF by chromosome, then concatenate back together
    ch_chroms_whatshap = channel.fromList(params.chromosomes)

    WHATSHAP_TRIO_PHASE_BY_CHROM(
        GLNEXUS_TRIO.out,
        trio_bams_assembled,
        file(params.reference),
        file(params.reference_index),
        ch_chroms_whatshap
    )

    // Split the trio-phased VCF back into per-individual VCFs and annotate each with FastVEP
    FASTVEP_ANNOTATE_TRIO_INDIVIDUALS(
        WHATSHAP_TRIO_PHASE_BY_CHROM.out.phased_vcf,
        trio_bams_assembled,
        file(params.fastvep_gff),
        file(params.reference),
        channel.fromPath("${params.fastvep_sa_dir}/*").collect()
    )

    // Read-backed phasing + haplotagging for each parent independently
    HIPHASE_TRIO_PARENTS(
        deeptrio_vcf_ch,
        individual_aligned_bams,
        sample_roles_ch,
        file(params.reference),
        file(params.reference_index)
    )

    // CpG methylation calling on the haplotagged child + parent BAMs
    CPG_TRIO_METHYLATION(
        WHATSHAP_TRIO_PHASE_BY_CHROM.out.haplotagged_bam,
        HIPHASE_TRIO_PARENTS.out.haplotagged_bam,
        file(params.reference),
        file(params.reference_index)
    )

    // Structural variants: per-sample discovery + joint genotyping per family (sawfish)
    mosdepth_run(individual_aligned_bams)
    infer_sex(mosdepth_run.out.summary)

    SAWFISH_TRIO(
        individual_aligned_bams,
        sample_to_family_ch,
        infer_sex.out.sex,
        file(params.expected_XY_bed),
        file(params.expected_XX_bed),
        file(params.excluded_bed),
        file(params.reference),
        file(params.reference_index)
    )

    // Tandem repeat genotyping on the non-haplotagged aligned BAMs
    TRGT_GENOTYPING(
        individual_aligned_bams,
        file(params.reference),
        file(params.reference_index),
        file(params.trgt_repeats_bed),
        file(params.expected_XY_bed),
        file(params.expected_XX_bed),
        infer_sex.out.sex
    )
}


// =========================================================================
//  ENTRY POINT: SINGLETON POST-ALIGNMENT ONLY
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
//  SUB-WORKFLOW: SINGLETON POST ALIGNMENT
// =========================================================================

workflow POST_ALIGNMENT {
    take:
    aligned_bam_ch // tuple(sample_id, bam, bai)

    main:
   
   
    bam_stats(aligned_bam_ch)

    mosdepth_run(aligned_bam_ch)
    infer_sex(mosdepth_run.out.summary)
    
    plot_dist_coverage(mosdepth_run.out.global_dist)
    
 
    
    // Call singletons variant calling subworkflow (50MB shards + concat)
    
    DEEPVARIANT_SINGLETON_WGS(
        file(params.reference),
        file(params.reference_index),
        aligned_bam_ch,
        channel.fromPath("${params.bed_dir}/*.bed")
    )
    
    
    
    aligned_bam_ch
        .join(DEEPVARIANT_SINGLETON_WGS.out.vcf, by: 0)
        .map { sample_id, bam, bai, vcf, vcf_tbi ->
            // Reorder to match: tuple(sample_id, vcf, vcf_tbi, bam, bai)
            tuple(sample_id, vcf, vcf_tbi, bam, bai)
        }
        .set { aligned_bam_with_vcf_ch }

    hiphase_small_variants(
        aligned_bam_with_vcf_ch,
        file(params.reference),
        file(params.reference_index)
    )    

    FASTVEP_ANNOTATE_WGS(
        hiphase_small_variants.out.phased_vcf,          // Directly passes phased VCF channel
        file(params.fastvep_gff),
        file(params.reference),
        channel.fromPath("${params.fastvep_sa_dir}/*").collect() // Staged into sa_dir/*
    )

    cpg_methylation_calling(
        hiphase_small_variants.out.haplotagged_bam,
        file(params.reference),
        file(params.reference_index)
    )
    
    

    expected_bed_ch = infer_sex.out.sex.map { sample_id, sex_csv ->
        def lines = sex_csv.readLines()
        def sex = lines.size() > 1 ? lines[1].split(',')[3].trim() : 'UNKNOWN'

        def expected_bed = (sex == 'FEMALE') ? file(params.expected_XX_bed) :
                           (sex == 'MALE')   ? file(params.expected_XY_bed) : null
        
        if (!expected_bed) {
            throw new Exception("Error: Invalid sex '${sex}' inferred for sample ${sample_id}.")
        }
        return tuple(sample_id, expected_bed)
    }
    
    sawfish_in_ch = aligned_bam_ch.join(expected_bed_ch, by: 0)

    sawfish_discover(
        sawfish_in_ch,
        file(params.excluded_bed),
        file(params.reference),
        file(params.reference_index)
    )

    sawfish_joint_call(
        sawfish_discover.out.discover_dir.collect()
    )
   
       TRGT_GENOTYPING(
        aligned_bam_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.trgt_repeats_bed),
        file(params.expected_XY_bed),
        file(params.expected_XX_bed),
        infer_sex.out.sex
    )
     
 
}











