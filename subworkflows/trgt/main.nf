/*
========================================================================================
    MODULE IMPORT
========================================================================================
*/

include { trgt } from '../../modules/pbtools/'


/*
========================================================================================
    SUBWORKFLOW: TRGT gentoyping
========================================================================================
*/

workflow TRGT_GENOTYPING {
    
    take:
    aligned_bam_ch    // channel: tuple(sample_id, bam, bai)
    fasta             // path(fasta)
    fasta_idx         // path (fasta idx)
    repeat_bed        // path(trgt repeats bed)
    expected_XY_bed   // path expected XY bed
    expected_XX_bed   // path expected XX bed
    infer_sex_ch     // output of modules/mosdepth/infer_sex process

    

    main:

      sex_val_ch = infer_sex_ch.map { sample_id, sex_csv ->
        def lines = sex_csv.readLines()

        if (lines.size() <= 1) {
            throw new Exception("Error: Missing inferred sex information for sample ${sample_id}. Expected 'FEMALE' or 'MALE'.")
        }

        def fields = lines[1].split(',')
        if (fields.size() <= 3) {
            throw new Exception("Error: Malformed inferred sex CSV for sample ${sample_id}. Could not read sex from line 2, column 4.")
        }

        def sex = fields[3].trim()
        if (sex != 'FEMALE' && sex != 'MALE') {
            throw new Exception("Error: Invalid sex '${sex}' inferred for sample ${sample_id}. Expected 'FEMALE' or 'MALE'.")
        }
        return tuple(sample_id, sex)
      }

    trgt_input_ch = aligned_bam_ch.join(sex_val_ch)
        .map { sample_id, bam, bai, sex ->
               tuple(sample_id, bam, bai, sex)
        }
    trgt(
        trgt_input_ch,
        fasta,
        fasta_idx,
        repeat_bed,
        expected_XY_bed,
        expected_XX_bed
    )

    emit:
      
    repeat_vcf  = trgt.out.repeat_vcf   // tuple(sample_id, vcf, tbi)
    repeat_bam = trgt.out.spanning_reads // tuple(sample_id, bam, bai)
    dropout_file = trgt.out.dropouts // tuple(sample_id, dropout_file)
    stats_file = trgt.out.stats // tuple(sample_id, stats_file)
    
}