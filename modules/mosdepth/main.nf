process mosdepth_run {

    tag "$sample_id"
    publishDir "${params.mosdepth_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    container "community.wave.seqera.io/library/mosdepth:0.3.10--259732f342cfce27"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}.mosdepth.region.dist.txt"), emit: region_dist
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"),path("${sample_id}.regions.bed.gz.csi"), emit: bed

    script:
    def args = task.cpus > 1 ? "--threads ${task.cpus - 1}" : ""
    def prefix = "${sample_id}"
    """
    mosdepth --version

    mosdepth \
        ${args} \
        -n \
        --fast-mode \
        --by 10000 \
        ${prefix} \
        ${bam}

   
    """
}

process plot_dist_coverage {
    tag "$sample_id"
    container "indapa/mosdepth-sex:latest" 
    publishDir "${params.mosdepth_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(global_summary)

    output:
    tuple val(sample_id), path("${sample_id}.dist.html"), emit: plot
    

    script:
    """
    python3 /opt/bin/plot-dist.py -o ${sample_id}.dist.html  ${global_summary}
    """ 
}


process infer_sex {
    tag "$sample_id"
    container "indapa/mosdepth-sex:latest" 
    publishDir "${params.mosdepth_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(summary)

    output:
    tuple val(sample_id), path("${sample_id}_inferred_sex.csv"), emit: sex
    tuple val(sample_id), path("${sample_id}_mean_depth.txt", optional: true), emit: mean

    script:
    """
    python3 /opt/bin/infer-sex.py --summary ${summary} --sample_id ${sample_id}
    """
}