include { glnexus_trio_by_chrom; concat_glnexus_vcf } from '../../modules/glnexus'

workflow GLNEXUS_TRIO {

    take:
    glnexus_input_ch   // tuple( sample_id, family_id, gvcf, tbi, role )
    ch_chroms          // channel of chromosome names (e.g. 'chr1', 'chr2', ...)
    region_bed         // path to full BED file

    main:
    // Regroup by family, extracting role-tagged file maps
    glnexus_prepared_ch = glnexus_input_ch
        .map { sample_id, family_id, vcf, tbi, role ->
            tuple(family_id, [role: role, vcf: vcf, tbi: tbi])
        }
        .groupTuple(by: 0)
        .map { family_id, members ->
            def child   = members.find { it.role == 'child' }
            def parent1 = members.find { it.role == 'parent1' }
            def parent2 = members.find { it.role == 'parent2' }

            tuple(
                family_id,
                child.vcf,   child.tbi,
                parent1.vcf, parent1.tbi,
                parent2.vcf, parent2.tbi
            )
        }

    // Scatter: combine each prepared trio with every chromosome
    scattered_ch = glnexus_prepared_ch
        .combine(ch_chroms)
        // -> tuple(family_id, child_gvcf, child_tbi, p1_gvcf, p1_tbi, p2_gvcf, p2_tbi, chrom)

    // Run GLnexus per family per chromosome
    glnexus_trio_by_chrom(scattered_ch, region_bed)

    // Gather: group all per-chrom results by family_id, then concat
    concat_input_ch = glnexus_trio_by_chrom.out.joint_vcf
        .groupTuple(by: 0)
        // -> tuple(family_id, [vcf1, vcf2, ...], [tbi1, tbi2, ...])

    concat_glnexus_vcf(concat_input_ch)

    emit:
    joint_vcf = concat_glnexus_vcf.out.merged  // tuple( family_id, joint.vcf.gz, joint.vcf.gz.tbi )
}