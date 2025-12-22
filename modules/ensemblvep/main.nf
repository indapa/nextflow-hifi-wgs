



process annotate_vep {
    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    containerOptions "--user \$(id -u):\$(id -g)"
    publishDir "${params.vep_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(phased_vcf), path(phased_tbi)
    path pigeon_gtf
    path pigeon_tbi
    path reference, stageAs: 'reference.fasta'
    path reference_fai, stageAs: 'reference.fasta.fai'
    
    output:
    tuple val(sample_id), path("${outfile}"), path("${outfile}.tbi"), emit: vep_vcf
    
    script:
    outfile = "${sample_id}.phased.vep.vcf.gz"
    warning_file = "${sample_id}.vep_warnings.txt"
   
    """
    vep -i ${phased_vcf} \
        -o ${outfile} \
        --format vcf \
        --vcf \
        --compress_output bgzip \
        --fork ${task.cpus} \
        --gtf ${pigeon_gtf} \
        --fasta reference.fasta \
        --force_overwrite \
        --warning_file ${warning_file} \
        --no_stats

    tabix -p vcf ${outfile}
    """
}