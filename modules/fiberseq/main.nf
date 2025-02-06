process fibertools_extract  {

    publishDir params.fibertools_output_dir, mode: 'copy'
    tag "$sample_id"


    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    path "${sample_id}.m6a.bed.gz", emit: m6a_bed 
    path "${sample_id}.m6a.bed.gz.tbi", emit: m6a_bed_tbi

    script:

    """

    ft --version
    ft extract ${bam} --m6a ${sample_id}.m6a.bed.gz
    tabix -p bed ${sample_id}.m6a.bed.gz

    """


    

}