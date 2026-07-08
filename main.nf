#!/usr/bin/env nextflow



include { pbmm2_align; cpg_methylation_calling; sawfish_discover; sawfish_joint_call; hiphase_small_variants } from './modules/pbtools'
include { 
    glnexus_trio_merge; 
    deeptrio_wgs; 
    deeptrio_wgs_by_chrom; 
    deepvariant_wgs;
    bcftools_concat as concat_p1_vcf;
    bcftools_concat as concat_p2_vcf;
    bcftools_concat as concat_child_gvcf;
    bcftools_concat as concat_p1_gvcfs;
    bcftools_concat as concat_p2_gvcfs 
} from './modules/deepvariant'
include { bam_stats } from './modules/samtools'
include { whatshap_trio_phase } from './modules/whatshap'
include { mosdepth_run; infer_sex; plot_dist_coverage } from './modules/mosdepth'
include { samtools_index } from './modules/samtools'


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
        """.stripIndent()
        exit 0
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
//  WORKFLOW: TRIO ANALYSIS ENTRYPOINTS
// =========================================================================

// --- Entrypoint 1: Starts from Raw Unaligned BAMs ---
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
    pbmm2_align(file(params.reference), align_input_ch)

    trio_bams_assembled = raw_samples_ch
        .map { fam, sample_id, role, _raw_bam -> tuple(sample_id, fam, role) }
        .join(pbmm2_align.out.aligned_bam)
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

    // Isolate single aligned BAM trackers for downstream tools (Sawfish/HiPhase)
    individual_aligned_bams = pbmm2_align.out.aligned_bam

    RUN_TRIO_PIPELINE(trio_bams_assembled, individual_aligned_bams)
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

    // Reconstruct flat stream of individual aligned BAMs for downstream hooks
    individual_aligned_bams = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.aligned_bam), file(row.aligned_bai)) }

    RUN_TRIO_PIPELINE(trio_bams_assembled, individual_aligned_bams)
}

// =========================================================================
//  SUB-WORKFLOW: SHARED TRIO DOWNSTREAM ENGINE
// =========================================================================

workflow RUN_TRIO_PIPELINE {
    take:
    trio_bams_assembled
    individual_aligned_bams

    main:
    chromosomes_ch = channel.of(
        'chr1','chr2','chr3','chr4','chr5','chr6',
        'chr7','chr8','chr9','chr10','chr11','chr12',
        'chr13','chr14','chr15','chr16','chr17','chr18',
        'chr19','chr20','chr21','chr22','chrX','chrY'
    )

    // Run DeepTrio (Scatter)
    deeptrio_wgs_by_chrom(
        file(params.reference),
        file(params.reference_index),
        trio_bams_assembled,
        chromosomes_ch
    )

    // Group scattered outputs per sample/role to prepare for concat
    //child_vcfs_ch  = deeptrio_wgs_by_chrom.out.child_vcf.groupTuple(by: [0, 1])
    child_gvcfs_ch = deeptrio_wgs_by_chrom.out.child_gvcf.groupTuple(by: [0, 1])

    p1_vcfs_ch     = deeptrio_wgs_by_chrom.out.p1_vcf.groupTuple(by: [0, 1])
    p1_gvcfs_ch    = deeptrio_wgs_by_chrom.out.p1_gvcf.groupTuple(by: [0, 1])

    p2_vcfs_ch     = deeptrio_wgs_by_chrom.out.p2_vcf.groupTuple(by: [0, 1])
    p2_gvcfs_ch    = deeptrio_wgs_by_chrom.out.p2_gvcf.groupTuple(by: [0, 1])

    // Concat scattered chromosomes via the generic reusable process
    //child_vcf_merged  = bcftools_concat(child_vcfs_ch, 'vcf.gz').merged
    p1_vcf_merged     = concat_p1_vcf(p1_vcfs_ch, 'vcf.gz').merged
    p2_vcf_merged     = concat_p2_vcf(p2_vcfs_ch, 'vcf.gz').merged

    child_gvcf_merged = concat_child_gvcf(child_gvcfs_ch, 'g.vcf.gz').merged
    p1_gvcf_merged    = concat_p1_gvcfs(p1_gvcfs_ch, 'g.vcf.gz').merged
    p2_gvcf_merged    = concat_p2_gvcfs(p2_gvcfs_ch, 'g.vcf.gz').merged

    // Normalize and join gVCF channels strictly by family_id for GLnexus
    glnexus_input_ch = child_gvcf_merged.map { fam, _id, v, t -> tuple(fam, v, t) }
        .join( p1_gvcf_merged.map { fam, _id, v, t -> tuple(fam, v, t) } )
        .join( p2_gvcf_merged.map { fam, _id, v, t -> tuple(fam, v, t) } )

    glnexus_trio_merge(glnexus_input_ch)

    // Assemble WhatsHap Input: Join joint VCF with our complete family structure channel
    whatshap_input_ch = glnexus_trio_merge.out.joint_vcf.join(trio_bams_assembled)

    whatshap_trio_phase(
        file(params.reference), 
        file(params.reference_index), 
        whatshap_input_ch
    )

    // Prepare HiPhase input tracking for Parents 
    parent_vcfs_ch = p1_vcf_merged.mix(p2_vcf_merged)
        .map { _fam, sample_id, vcf, tbi -> tuple(sample_id, vcf, tbi) }

    hiphase_input_ch = parent_vcfs_ch.join(individual_aligned_bams)

    hiphase_small_variants(
        hiphase_input_ch,        
        file(params.reference),        
        file(params.reference_index)   
    )

    // Reconstruct Child ID track mapping to align filename context for Samtools Index
    child_bam_ch = whatshap_trio_phase.out.haplotagged_bam
        .map { bam -> tuple(bam.baseName.replaceAll(/\.haplotagged/, ''), bam) }
    
    samtools_index(child_bam_ch)
    child_cpg_input = samtools_index.out

    // Mix parent and child haplotagged streams cleanly
    parent_cpg_input = hiphase_small_variants.out.haplotagged_bam
    all_haplotagged_bams_ch = child_cpg_input.mix(parent_cpg_input)

    cpg_methylation_calling(
        all_haplotagged_bams_ch,
        file(params.reference),
        file(params.reference_index)
    )

    // QC, Sex inference, and Sawfish Structural Variant Pipeline
    mosdepth_run(individual_aligned_bams)
    infer_sex(mosdepth_run.out.summary)
    plot_dist_coverage(mosdepth_run.out.global_dist)

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

    sawfish_in_ch = individual_aligned_bams.join(expected_bed_ch, by: 0)

    sawfish_discover(
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
    aligned_bam_ch

    main:
    bam_stats(aligned_bam_ch)

    deepvariant_wgs(
        file(params.reference),
        file(params.reference_index),
        aligned_bam_ch
    )

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
}










