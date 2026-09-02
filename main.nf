#!/usr/bin/env nextflow



include { pbmm2_align; cpg_methylation_calling; sawfish_discover; sawfish_joint_call; hiphase_small_variants } from './modules/pbtools'
include { 
    deeptrio_wgs; 
    deeptrio_wgs_by_chrom; 
    deepvariant_wgs;
    concat_chrom_chunks_vcf;
    concat_wgs_vcf

} from './modules/deepvariant'
include { bam_stats; slice_trio_bams_by_interval; samtools_index } from './modules/samtools'
include { FASTVEP_ANNOTATE_TRIO_VCF; FASTVEP_ANNOTATE_SINGLETON_VCF } from './modules/fastvep'

include { WHATSHAP_TRIO_PHASE_BY_CHROM } from './subworkflows/whatshap_trio_phase_by_chrom'
include { CONCAT_AND_SPLIT_WGS } from './subworkflows/concat_and_split_wgs'
include { GLNEXUS_TRIO }         from './subworkflows/glnexus_trio_merge'
include { PBMM2_ALIGN; PBMM2_SPOT_WGS} from './subworkflows/reference_alignment'
include {FASTVEP_ANNOTATE_WGS} from './subworkflows/fastvep'
include {DEEPVARIANT_SINGLETON_WGS} from './subworkflows/deepvariant_singleton_wgs'
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
        """.stripIndent()
        exit 0
    }

    if (params.entry == 'WGS_TRIO_ALIGNED') {
        WGS_TRIO_ALIGNED()
    }   else{
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
        PBMM2_ALIGN.out
    )

    }
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


    sample_roles_ch = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.role) }

    sample_to_family_ch = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.family_id) }

    // Isolate single aligned BAM trackers for downstream tools (Sawfish/HiPhase)
    individual_aligned_bams = pbmm2_align.out.aligned_bam

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


    sample_roles_for_hiphase = channel.fromPath(params.trio_aligned_samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, row.role) }


    // 1. Load all interval BED files from the specified directory
    raw_intervals_ch = channel.fromPath("${params.intervals_dir}/*.bed")

    // 2. Filter to autosomes + X/Y chunks
    intervals_ch = raw_intervals_ch.filter { f -> f.baseName =~ /^chr([1-9]|1[0-9]|2[0-2]|[XY])_/ }

    // =========================================================================
    // Pre-calculate how many chunks exist per chromosome so groupTuple can
    // emit eagerly via groupKey without waiting for the entire channel to close.
    // =========================================================================
    def counts_by_chrom = [:]
    file(params.intervals_dir).list().each { name ->
        if (name.endsWith(".bed") && name.startsWith("chr")) {
            def chrom = name.split('_')[0]
            counts_by_chrom[chrom] = (counts_by_chrom[chrom] ?: 0) + 1
        }
    }

    // Combine trio tuple with each interval BED file to create a scatter input channel
    slicing_matrix_ch = trio_bams_assembled.combine(intervals_ch)

    //log.info "DEBUG: trio interval slice channel contents:"
    //slicing_matrix_ch.view { fam, c_id, c_bam, c_bai, p1_id, p1_bam, p1_bai, p2_id, p2_bam, p2_bai, interval_bed ->
    //    return "Family: ${fam}, Child: ${c_id}, ${c_bam}, ${c_bai}, ${p1_id}, ${p1_bam}, ${p1_bai}, ${p2_id}, ${p2_bam}, ${p2_bai}, Interval: ${interval_bed.baseName}"
    //}

    slice_trio_bams_by_interval(slicing_matrix_ch)

    // Run DeepTrio (Scatter)
    deeptrio_wgs_by_chrom(
        file(params.reference),
        file(params.reference_index),
        slice_trio_bams_by_interval.out.sliced_trio_package
    )

    all_chunks_ch = channel.empty()
        .mix(
            deeptrio_wgs_by_chrom.out.child_gvcf.map { fam, id, bed, f, t -> [ [fam, id, file(bed).baseName.split('_')[0], 'g.vcf.gz'], f, t ] },
            deeptrio_wgs_by_chrom.out.p1_gvcf.map    { fam, id, bed, f, t -> [ [fam, id, file(bed).baseName.split('_')[0], 'g.vcf.gz'], f, t ] },
            deeptrio_wgs_by_chrom.out.p2_gvcf.map    { fam, id, bed, f, t -> [ [fam, id, file(bed).baseName.split('_')[0], 'g.vcf.gz'], f, t ] },
            deeptrio_wgs_by_chrom.out.p1_vcf.map     { fam, id, bed, f, t -> [ [fam, id, file(bed).baseName.split('_')[0], 'vcf.gz'],   f, t ] },
            deeptrio_wgs_by_chrom.out.p2_vcf.map     { fam, id, bed, f, t -> [ [fam, id, file(bed).baseName.split('_')[0], 'vcf.gz'],   f, t ] },
            deeptrio_wgs_by_chrom.out.child_vcf.map  { fam, id, bed, f, t -> [ [fam, id, file(bed).baseName.split('_')[0], 'vcf.gz'],   f, t ] }
        )

    // Group chunks safely using the exact unique metadata profile
    grouped_chrom_chunks = all_chunks_ch
        .map { meta, f, t ->
            // Pass the literal meta array [fam, id, chrom, type] directly
            tuple(meta, f, t)
        }
        // Group by index 0 (the meta array) so samples never cross-contaminate
        .groupTuple(by: 0)

    // Concatenate per-chromosome chunks safely
    concat_chrom_chunks_vcf(grouped_chrom_chunks)

    CONCAT_AND_SPLIT_WGS(
        concat_chrom_chunks_vcf.out.merged_file,
        sample_roles_ch
    )
    ch_chroms_glnexus = channel.fromList(params.chromosomes)
    ch_chroms_whatshap = channel.fromList(params.chromosomes)
    GLNEXUS_TRIO(
        CONCAT_AND_SPLIT_WGS.out.glnexus_input,
        ch_chroms_glnexus,
        file(params.glnexus_region_bed)
    )
    

    
    whatshap_input_ch = GLNEXUS_TRIO.out
    .join(trio_bams_assembled, by: 0)


    WHATSHAP_TRIO_PHASE_BY_CHROM(
        whatshap_input_ch,
        params.reference,
        params.reference_index,
        ch_chroms_whatshap
    )


    // annotate the phased VCF with FASTVEP
    ch_gff3     = channel.fromPath(params.fastvep_gff)
    ch_fasta    = channel.fromPath(params.reference)
    ch_sa_files = channel.fromPath("${params.fastvep_sa_dir}/**").collect()

    FASTVEP_ANNOTATE_TRIO_VCF(
        WHATSHAP_TRIO_PHASE_BY_CHROM.out.phased_vcf,
        ch_gff3,
        ch_fasta,
        ch_sa_files
    )

 

 
   hiphase_parents_input_ch = CONCAT_AND_SPLIT_WGS.out.vcf_merged
    .map { _family_id, sample_id, file, tbi ->
        // Key by sample_id to join with individual BAMs (4 elements in signature)
        tuple(sample_id.toString(), file, tbi)
    }
    // Join with your individual BAMs by sample_id (index 0)
    .join(individual_aligned_bams, by: 0) // -> [sample_id, vcf, tbi, bam, bai]
    // Join with roles to find parent1 and parent2
    .join(sample_roles_for_hiphase, by: 0)         // -> [sample_id, vcf, tbi, bam, bai, role]
    // Keep only the parents
    .filter { _sample_id, _vcf, _tbi, _bam, _bai, role -> 
        role == 'parent1' || role == 'parent2' 
    }
    // Drop the role to match HiPhase's exact input tuple
    .map { sample_id, vcf, tbi, bam, bai, _role -> 
        tuple(sample_id, vcf, tbi, bam, bai) 
    }

   hiphase_parents_input_ch.view { sample_id, vcf, tbi, bam, bai ->
       return "DEBUG: HiPhase Input - Parent Sample: ${sample_id} | VCF: ${vcf.name}  TBI: ${tbi.name} | BAM: ${bam.name} | BAI : ${bai.name}"
   }


     
    hiphase_small_variants(
        hiphase_parents_input_ch,        
        file(params.reference),        
        file(params.reference_index)   
    )
    
    // Reconstruct Child ID track mapping to align filename context for Samtools Index
    child_bam_ch = WHATSHAP_TRIO_PHASE_BY_CHROM.out.haplotagged_bam
        .map { bam -> tuple(bam.baseName.replaceAll(/\.haplotagged/, ''), bam) }
    
    samtools_index(child_bam_ch)
    child_cpg_input = samtools_index.out

    // Mix parent and child haplotagged streams cleanly
    parent_cpg_input = hiphase_small_variants.out.haplotagged_bam
    all_haplotagged_bams_ch = child_cpg_input.mix(parent_cpg_input)
    all_haplotagged_bams_ch.view() { sample_id, bam, bai ->
        return "DEBUG: CpG Input - Sample: ${sample_id} | BAM: ${bam.name} | BAI: ${bai.name}"
    }
    

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

    expected_bed_ch.view { sample_id, bed ->
        return "DEBUG: Expected BED for Sawfish - Sample: ${sample_id} | BED: ${bed.name}"
    }
    
    sawfish_in_ch = individual_aligned_bams.join(expected_bed_ch, by: 0)
    sawfish_in_ch.view { sample_id, bam, bai, expected_bed ->
        return "DEBUG: Sawfish Input - Sample: ${sample_id} | BAM: ${bam.name} | BAI: ${bai.name} | Expected BED: ${expected_bed.name}"
    }
   
    sawfish_discover(
        sawfish_in_ch,
        file(params.excluded_bed),
        file(params.reference),
        file(params.reference_index)
    )

    sawfish_discover.out.discover_dir // Emits: [sample_id, discover_dir] (Ensure the process emits sample_id in its tuple)
    // Join with sample_roles_ch to get the family context
    // Assuming you have/can make a channel that maps: [sample_id, family_id]
    .join(sample_to_family_ch, by: 0) // -> [sample_id, discover_dir, family_id]
    .map { _sample_id, discover_dir, family_id ->
        tuple(family_id, discover_dir)
    }
    .groupTuple(by: 0) // Groups all 3 members under the family_id -> [family_id, [dir_child, dir_p1, dir_p2]]
    .set { sawfish_joint_input_ch }

    sawfish_joint_call(
        sawfish_joint_input_ch
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
}











