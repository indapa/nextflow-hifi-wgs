process mosdepth {
    publishDir "${params.bam_stats_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "community.wave.seqera.io/library/mosdepth:0.3.10--259732f342cfce27"
    
    input:
    tuple val(sample_id),path(bam), path(bam_index)

    output:
    path "*.mosdepth.global.dist.txt", emit: global_dist
    path "*.mosdepth.summary.txt", emit: summary
    path "*.regions.bed.gz", emit: region_bed

    script:
    prefix = bam.baseName
    """
    mosdepth --version

    mosdepth \
        --threads ${task.cpus - 1} \
        -n \
        --fast-mode \
        --by 10000 \
        ${prefix} \
        ${bam}
    """

    stub:
    def prefix = bam.baseName
    """
    
    
    # Create mock mosdepth output files
    echo -e "depth\tcount" > ${prefix}.mosdepth.global.dist.txt
    echo -e "0\t1000" >> ${prefix}.mosdepth.global.dist.txt
    echo -e "1\t5000" >> ${prefix}.mosdepth.global.dist.txt
    echo -e "10\t15000" >> ${prefix}.mosdepth.global.dist.txt
    
    echo -e "chrom\tlength\tbases\tmean\tmin\tmax" > ${prefix}.mosdepth.summary.txt
    echo -e "chr1\t248956422\t240000000\t30.5\t0\t150" >> ${prefix}.mosdepth.summary.txt
    echo -e "total\t3088269832\t2900000000\t28.2\t0\t150" >> ${prefix}.mosdepth.summary.txt
    
    echo -e "chr1\t0\t10000\t25.5" | gzip > ${prefix}.regions.bed.gz
    echo -e "chr1\t10000\t20000\t32.1" | gzip >> ${prefix}.regions.bed.gz
    """

}