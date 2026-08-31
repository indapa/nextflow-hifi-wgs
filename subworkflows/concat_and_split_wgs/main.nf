include { concat_wgs_vcf } from '../../modules/deepvariant'

workflow CONCAT_AND_SPLIT_WGS {

    take:
    ch_chrom_chunks    // tuple( [family_id, sample_id, chrom, ext], vcf, tbi )
    sample_roles_ch    // tuple( sample_id_string, role )

    main:
    // Regroup all chromosomes together per sample/type
    wgs_input_ch = ch_chrom_chunks
        .map { meta, vcf, tbi ->
            // meta is: [family_id, sample_id, chrom, ext]
            // Regroup by: [ family_id, sample_id, ext ]
            tuple( [meta[0], meta[1], meta[3]], vcf, tbi )
        }
        .groupTuple(by: 0)
        .map { group_key, vcfs, tbis ->
            // group_key is: [family_id, sample_id, ext]
            // Split to match concat_wgs_vcf input signature:
            // tuple(family_id, sample_id, vcfs, tbis) and the separate ext
            tuple( tuple(group_key[0], group_key[1], vcfs, tbis), group_key[2] )
        }

    concat_wgs_vcf(
        wgs_input_ch.map { tup, _ext -> tup },
        wgs_input_ch.map { _tup, ext -> ext }
    )

    // Branch merged output by file extension
    concat_wgs_vcf.out.merged
        .branch { _family_id, _sample_id, file, _tbi ->
            vcf:  file.name.endsWith('vcf.gz') && !file.name.endsWith('g.vcf.gz')
            gvcf: file.name.endsWith('g.vcf.gz')
        }
        .set { split_concat_ch }

    // Prepare gVCF stream with String-typed join key
    prep_concat_ch = split_concat_ch.gvcf
        .map { family_id, sample_id, file, tbi ->
            tuple(sample_id.toString(), family_id, file, tbi)
        }

    // Join gVCFs with sample roles on sample_id
    glnexus_input_ch = prep_concat_ch
        .join(sample_roles_ch, by: 0)

    emit:
    vcf_merged    = split_concat_ch.vcf       // tuple( family_id, sample_id, vcf, tbi )
    glnexus_input = glnexus_input_ch          // tuple( sample_id, family_id, gvcf, tbi, role )
}