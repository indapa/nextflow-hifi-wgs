include { glnexus_trio_by_chrom; concat_glnexus_vcf } from '../../modules/glnexus'

workflow GLNEXUS_TRIO {

    take:
    deeptrio_gvcf_ch  // tuple(family_id, sample_id, gvcf, tbi) -- e.g. DEEPTRIO_WGS.out.gvcf
    sample_roles_ch   // tuple(sample_id, role) -- role is 'child', 'parent1', or 'parent2'
    ch_chroms         // channel of chromosome names (e.g. 'chr1', 'chr2', ...)
    region_bed        // path to full BED file

    main:
    // Attach role to each per-sample gVCF, then regroup by family into a single
    // child/parent1/parent2 tuple
    glnexus_prepared_ch = deeptrio_gvcf_ch
        .map { family_id, sample_id, gvcf, tbi -> tuple(sample_id, family_id, gvcf, tbi) }
        .join(sample_roles_ch, by: 0)
        .map { _sample_id, family_id, gvcf, tbi, role ->
            tuple(family_id, [role: role, gvcf: gvcf, tbi: tbi])
        }
        .groupTuple(by: 0)
        .map { family_id, members ->
            def child   = members.find { member -> member.role == 'child' }
            def parent1 = members.find { member -> member.role == 'parent1' }
            def parent2 = members.find { member -> member.role == 'parent2' }

            tuple(
                family_id,
                child.gvcf,   child.tbi,
                parent1.gvcf, parent1.tbi,
                parent2.gvcf, parent2.tbi
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
    concat_glnexus_vcf.out.merged  // tuple( family_id, joint.vcf.gz, joint.vcf.gz.tbi )
}
