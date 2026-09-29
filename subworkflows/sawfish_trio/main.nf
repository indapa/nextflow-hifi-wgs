include { sawfish_discover; sawfish_joint_call } from '../../modules/pbtools'

// Runs sawfish's two-step trio/pedigree workflow (per PacBio's user guide:
// https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md):
//   1. `discover` independently per sample
//   2. `joint-call` across all discover directories for a family
// sawfish tracks each sample's expected copy number independently through
// joint-call, so mixed-sex trios (e.g. an XX mother/father and XY child) are
// handled correctly without any special-casing here.
workflow SAWFISH_TRIO {

    take:
    aligned_bam_ch      // channel: tuple(sample_id, bam, bai)
    sample_to_family_ch // channel: tuple(sample_id, family_id)
    infer_sex_ch        // channel: output of modules/mosdepth/infer_sex, tuple(sample_id, sex_csv)
    expected_XY_bed     // path: expected-cn BED for XY samples
    expected_XX_bed     // path: expected-cn BED for XX samples
    excluded_bed        // path: CNV-excluded regions BED
    ref_file            // path(reference_fasta)
    ref_index_file      // path(reference_fai)

    main:
    // 1. Resolve each sample's expected-cn BED from its inferred sex
    expected_bed_ch = infer_sex_ch.map { sample_id, sex_csv ->
        def lines = sex_csv.readLines()
        if (lines.size() <= 1) {
            throw new Exception("Error: Missing inferred sex information for sample ${sample_id}.")
        }
        def sex = lines[1].split(',')[3].trim()
        def expected_bed = (sex == 'FEMALE') ? expected_XX_bed :
                           (sex == 'MALE')   ? expected_XY_bed : null
        if (!expected_bed) {
            throw new Exception("Error: Invalid sex '${sex}' inferred for sample ${sample_id}. Expected 'FEMALE' or 'MALE'.")
        }
        tuple(sample_id, expected_bed)
    }

    // 2. Per-sample discovery -- independent per sample, runs in parallel
    sawfish_in_ch = aligned_bam_ch.join(expected_bed_ch, by: 0)

    sawfish_discover(
        sawfish_in_ch,
        excluded_bed,
        ref_file,
        ref_index_file
    )

    // 3. Regroup each family's discover directories for joint genotyping
    sawfish_joint_input_ch = sawfish_discover.out.discover_dir
        .join(sample_to_family_ch, by: 0)
        .map { _sample_id, discover_dir, family_id -> tuple(family_id, discover_dir) }
        .groupTuple(by: 0)

    // 4. Multi-sample joint-call per family
    sawfish_joint_call(sawfish_joint_input_ch)

    emit:
    discover_dir = sawfish_discover.out.discover_dir  // tuple(sample_id, discover_dir)
    joint_dir    = sawfish_joint_call.out.joint_dir    // tuple(family_id, joint_call_dir)
}
