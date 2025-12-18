process pbmm2_align {
    label 'high_memory'
    publishDir "${params.aligned_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pbmm2:1.17.0_build1"

    input:
    path reference
    tuple val(sample_id), path(bam)
    val threads
    val sort_threads

    output:
    tuple path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: aligned_bam
    //path "${$sample_id}.read_length_and_quality.tsv", emit: bam_rl_qual
    
    script:
    """
    pbmm2 --version
    pbmm2 align \\
        --sort \\
        -j $threads \\
        -J $sort_threads \\
        --preset HIFI \\
        --sample ${sample_id} \\
        --log-level INFO \\
        --unmapped \\
        --bam-index BAI \\
        $reference \\
        $bam \\
        ${sample_id}.aligned.bam

    
    """

    stub:
    """
    # Create mock aligned BAM file
    touch ${sample_id}.aligned.bam
    
    # Create mock BAM index file
    touch ${sample_id}.aligned.bam.bai
    
    """
    


}

process cpg_pileup_haplotagged {
    publishDir "${params.cpg_output_dir}", mode: 'copy'
    tag "$bam"
    container "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"

    input:
    tuple path(bam), path(bam_index)
    

    output:
    path "${bam.baseName}.hap1.bed.gz", emit: hap1_bed
    path "${bam.baseName}.hap2.bed.gz", emit: hap2_bed
    path "${bam.baseName}.hap1.bw", emit: hap1_bigwig
    path "${bam.baseName}.hap2.bw", emit: hap2_bigwig
    
    path "${bam.baseName}.combined.bed.gz", emit: combined_bed
    path "${bam.baseName}.combined.bw", emit: combined_bigwig


    script:
    """
    aligned_bam_to_cpg_scores \\
        --threads 4 \\
        --bam ${bam} \\
        --output-prefix ${bam.baseName} \\
        --min-mapq 20 \\
        --min-coverage 10 
        
    """
}

process cpg_pileup_downsampled {

    publishDir "${params.cpg_downsampled_output_dir}", mode: 'copy', overwrite: true
    tag "$bam"
    container "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"
    errorStrategy 'ignore'  // Add this line to skip failed samples

    input:
    tuple path(bam), path(bam_index)
    

    output:
    path "${bam.baseName}.*.bed.gz", optional: true, emit: pileup_beds
    path "${bam.baseName}.*.bw", optional: true, emit: pileup_bigwigs

    script:
    """
    aligned_bam_to_cpg_scores \\
        --threads 4 \\
        --bam ${bam} \\
        --output-prefix ${bam.baseName} \\
        --min-mapq 20 \\
        --min-coverage 4 

    """


     stub:
    """
    echo -e "chr1\t1000\t1001\t0.85" | gzip > test.aligned.combined.bed.gz
    touch test.aligned.combined.bw
    """

}


process cpg_pileup {
    
    publishDir "${params.cpg_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"
    errorStrategy 'ignore'
    

    input:
    tuple val(sample_id), path(bam), path(bam_index)

    output:
    tuple val(sample_id), path("${bam.baseName}.*.bed.gz"), optional: true, emit: pileup_beds
    tuple val(sample_id), path("${bam.baseName}.*.bw"), optional: true, emit: pileup_bigwigs
    tuple val(sample_id), path("${bam.baseName}.cpg_params.yml"), emit: params_yml

    script:
    """
    # Create parameters YAML file
    cat > ${bam.baseName}.cpg_params.yml <<EOF
    tool: aligned_bam_to_cpg_scores
    container: quay.io/pacbio/pb-cpg-tools:3.0.0_build1
    parameters:
      threads: 4
      bam: ${bam}
      output_prefix: ${bam.baseName}
      min_mapq: ${params.cpg_min_mapq}
      min_coverage: ${params.cpg_min_coverage}
    sample_id: ${sample_id}
    timestamp: \$(date -Iseconds)
    EOF

    aligned_bam_to_cpg_scores \\
        --threads 4 \\
        --bam ${bam} \\
        --output-prefix ${bam.baseName} \\
        --min-mapq ${params.cpg_min_mapq} \\
        --min-coverage ${params.cpg_min_coverage}
    """

    stub:
    """
    echo -e "chr1\t1000\t1001\t0.85" | gzip > ${bam.baseName}.aligned.combined.bed.gz
    touch ${bam.baseName}.aligned.combined.bw
    touch ${bam.baseName}.cpg_params.yml
    """
}

process cpg_pileup_filtered {
    
    publishDir "${params.cpg_filtered_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/pb-cpg-tools:3.0.0_build1"
    errorStrategy 'ignore'
    

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    

    output:
    tuple val(sample_id), path("${bam.baseName}.*.bed.gz"), optional: true, emit: pileup_beds
    tuple val(sample_id), path("${bam.baseName}.*.bw"), optional: true, emit: pileup_bigwigs

    script:
    """
    aligned_bam_to_cpg_scores \\
        --threads 4 \\
        --bam ${bam} \\
        --output-prefix ${bam.baseName} \\
        --min-mapq 60 \\
        --min-coverage 4
        
        
    """

    stub:
    """
    # Create mock BED files (compressed)
    echo -e "chr1\t1000\t1001\t0.85" | gzip > ${bam.baseName}.aligned.combined.bed.gz
    
    # Create mock BigWig files
    touch ${bam.baseName}.aligned.combined.bw
    """
}


process hificnv {
    publishDir "${params.cnv_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/hificnv:1.0.1_build1"
    errorStrategy 'ignore'  // Add this line to skip failed samples

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path reference
    path reference_index
    path exclude_bed
    path expected_bed
    val(cpus)

    output:
    tuple val(sample_id), path("*.vcf.gz"), optional: true, emit: hifi_vcf
    tuple val(sample_id), path("*.bw"), optional: true, emit: hifi_bigwig
    tuple val(sample_id), path("*.bedgraph"), optional: true, emit: hifi_bedgraph
    tuple val(sample_id), path("*.log"), emit: hifi_log


    script:
    """
    set -euo pipefail

    hificnv --version
    hificnv --bam ${bam} \\
    --ref ${reference} \\
    --exclude  ${exclude_bed} \\
    --expected-cn ${expected_bed} \\
    --threads ${cpus} \\
    --output-prefix hificnv.${sample_id}
    """

    stub:
    """
    # Create mock VCF file (compressed)
    echo -e "##fileformat=VCFv4.2" > hificnv.${sample_id}.vcf
    echo -e "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" >> hificnv.${sample_id}.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> hificnv.${sample_id}.vcf
    echo -e "chr1\t1000000\t.\tN\t<DUP>\t60\tPASS\tSVTYPE=DUP" >> hificnv.${sample_id}.vcf
    gzip hificnv.${sample_id}.vcf

    # Create mock BigWig file
    touch hificnv.${sample_id}.bw
    echo "Mock BigWig CNV data" > hificnv.${sample_id}.depth.bw

    # Create mock BedGraph file
    echo -e "chr1\t0\t1000000\t2.0" > hificnv.${sample_id}.bedgraph
    echo -e "chr1\t1000000\t2000000\t3.0" >> hificnv.${sample_id}.copynum.bedgraph

    # Create mock log file
    echo "HiFiCNV analysis started" > hificnv.${sample_id}.log
    echo "Processing sample: ${sample_id}" >> hificnv.${sample_id}.log
    echo "Analysis completed successfully" >> hificnv.${sample_id}.log
    """



}


process trgt {
    publishDir "${params.trgt_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    container "quay.io/pacbio/trgt:3.0.0_build1"
    
    input:
        tuple val(sample_id), path(bam), path(bam_index)
        path reference
        path reference_index
        path tandem_repeat_bed
        val(karyotype)
        val(cpus)

    output:
        tuple val(sample_id), path("${sample_id}.trgt.spanning.sorted.bam"), path("${sample_id}.trgt.spanning.sorted.bam.bai"), optional: true, emit: spanning_reads
        tuple val(sample_id), path("${sample_id}.trgt.sorted.vcf.gz"), path("${sample_id}.trgt.sorted.vcf.gz.tbi"), optional: true, emit: repeat_vcf
    
    script:
    """
    set -euo pipefail

    trgt --version
    trgt genotype \\
        --threads ${cpus} \\
        --karyotype ${karyotype} \\
        --genome ${reference} \\
        --repeats ${tandem_repeat_bed} \\
        --reads ${bam} \\
        --output-prefix ${sample_id}.trgt

    bcftools --version
    bcftools sort \\
        --output-type z \\
        --output-file ${sample_id}.trgt.sorted.vcf.gz \\
        ${sample_id}.trgt.vcf.gz

    bcftools index \\
        --threads ${cpus} \\
        --tbi \\
        ${sample_id}.trgt.sorted.vcf.gz

    samtools --version
    samtools sort \\
        -@ ${cpus} \\
        -o ${sample_id}.trgt.spanning.sorted.bam \\
        ${sample_id}.trgt.spanning.bam

    samtools index \\
        -@ ${cpus} \\
        ${sample_id}.trgt.spanning.sorted.bam
    """

    stub:
    """
    # Create mock sorted BAM file for spanning reads
    echo "Mock BAM header" > ${sample_id}.trgt.spanning.sorted.bam
    echo "Mock spanning reads data" >> ${sample_id}.trgt.spanning.sorted.bam
    
    # Create mock BAM index file
    touch ${sample_id}.trgt.spanning.sorted.bam.bai
    
    
    touch ${sample_id}.trgt.sorted.vcf.gz
    
    # Create mock tabix index file
    touch ${sample_id}.trgt.sorted.vcf.gz.tbi
    """


}

/* next set of processes deal with SVs */


process pb_discover {
    //publishDir params.sv_output_dir, mode: 'copy'
    tag "${sample_id}:${region}"
    container "quay.io/pacbio/pbsv:2.11.0_build1"
    input:
    tuple val(sample_id), val(region), path(bam), path(bam_index)
    path trf_bed
       
    output:
    tuple val(sample_id), path("${sample_id}.${region}.svsig.gz"), emit: pb_discover
        
    script:
    """
    set -euo pipefail

    pbsv --version
    pbsv discover --hifi --tandem-repeats ${trf_bed} --region ${region} ${bam} ${sample_id}.${region}.svsig.gz
    """

    stub:
    """
    # Create mock structural variant signature file (compressed)
    echo -e "# PBSV signature file" > ${sample_id}.${region}.svsig
    echo -e "# Sample: ${sample_id}" >> ${sample_id}.${region}.svsig
    echo -e "# Region: ${region}" >> ${sample_id}.${region}.svsig
    echo -e "chr1\t1000000\t1000500\tDEL\t500\t60" >> ${sample_id}.${region}.svsig
    echo -e "chr1\t2000000\t2001000\tINS\t1000\t55" >> ${sample_id}.${region}.svsig
    echo -e "chr1\t3000000\t3002000\tDUP\t2000\t65" >> ${sample_id}.${region}.svsig
    gzip ${sample_id}.${region}.svsig
    """

}




process pb_call {
    label 'high_memory_spot'
    container "quay.io/pacbio/pbsv:2.11.0_build1"
    publishDir "${params.sv_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tag "$sample_id"
    input:
    tuple val(sample_id), path(svsig_files)
    path reference
       
    output:
    tuple val(sample_id), path("${sample_id}.pbsv.vcf.gz"), path("${sample_id}.pbsv.vcf.gz.tbi"), emit: pb_call
        
    script:
    """
    set -euo pipefail

    pbsv --version
    pbsv call -j 8 ${reference} ${svsig_files.join(' ')} ${sample_id}.pbsv.vcf

    bgzip ${sample_id}.pbsv.vcf
    bcftools index --tbi ${sample_id}.pbsv.vcf.gz
    

    """

    stub:
    """
    # Create mock PBSV VCF file with structural variants
    echo -e "##fileformat=VCFv4.2" > ${sample_id}.pbsv.vcf
    echo -e "##source=pbsv-2.11.0" >> ${sample_id}.pbsv.vcf
    echo -e "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" >> ${sample_id}.pbsv.vcf
    echo -e "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" >> ${sample_id}.pbsv.vcf
    echo -e "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" >> ${sample_id}.pbsv.vcf
    echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> ${sample_id}.pbsv.vcf
    echo -e "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">" >> ${sample_id}.pbsv.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_id}" >> ${sample_id}.pbsv.vcf
    echo -e "chr1\t1000000\tpbsv.DEL.1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-500;END=1000500\tGT:AD\t0/1:15,10" >> ${sample_id}.pbsv.vcf
    echo -e "chr1\t2000000\tpbsv.INS.1\tN\t<INS>\t55\tPASS\tSVTYPE=INS;SVLEN=1000;END=2000001\tGT:AD\t0/1:12,8" >> ${sample_id}.pbsv.vcf
    echo -e "chr1\t3000000\tpbsv.DUP.1\tN\t<DUP>\t65\tPASS\tSVTYPE=DUP;SVLEN=2000;END=3002000\tGT:AD\t0/1:18,12" >> ${sample_id}.pbsv.vcf
    echo -e "chr2\t5000000\tpbsv.INV.1\tN\t<INV>\t58\tPASS\tSVTYPE=INV;SVLEN=0;END=5001500\tGT:AD\t0/1:14,9" >> ${sample_id}.pbsv.vcf
    
    # Compress the VCF file
    bgzip ${sample_id}.pbsv.vcf
    
    # Create mock tabix index file
    touch ${sample_id}.pbsv.vcf.gz.tbi
    """
    
}


process hiphase_small_variants {
    /* hiphase small variants only */

    label 'high_memory'
    publishDir "${params.hiphase_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    container "quay.io/pacbio/hiphase:1.5.0_build1"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf), path(vcf_tbi), path(pbmm2_bam), path(pbmm2_bai)
    path (reference)

    output:
    tuple val(sample_id), path("*.phased.vcf.gz"), path("*.phased.vcf.gz.tbi"), emit: phased_vcf 
    tuple val(sample_id), path("*.stats.csv"), path("*.blocks.tsv"), path("*.summary.tsv"), emit: stats
    tuple val(sample_id), path("*.haplotagged.bam"), path("*.haplotagged.bam.bai"), emit: haplotagged_bam

    script:
    def basename = vcf.simpleName  // Gets name without extension
    """
    hiphase --version

    hiphase --reference ${reference} \
            --bam ${pbmm2_bam} \
            --vcf ${vcf} \
            --output-vcf ${basename}.phased.vcf.gz \
            --threads ${task.cpus} \
            --min-mapq 20 \
            --disable-global-realignment \
            --output-bam ${sample_id}.haplotagged.bam \
            --stats-file ${basename}.stats.csv \
            --blocks-file ${basename}.blocks.tsv \
            --summary-file ${basename}.summary.tsv

    bcftools index -f  --tbi ${basename}.phased.vcf.gz
    samtools index ${sample_id}.haplotagged.bam
    """
}

process hiphase {
    label 'high_memory'
    publishDir "${params.hiphase_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    container "quay.io/pacbio/hiphase:1.5.0_build1"
    tag "$sample_id"
    
    input:
    path deepvariant_vcf
    path deepvariant_tbi
    tuple val(sample_id), path(pbsv_vcf), path(pbsv_tbi)
    tuple val(sample_id), path(trgt_vcf), path(trgt_tbi)
    tuple path(pbmm2_bam), path(pbmm2_bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.phased.vcf.gz"), path("${sample_id}.deepvariant.phased.vcf.gz.tbi"), emit: phased_deepvariant 
    tuple val(sample_id), path("${sample_id}.pbsv.phased.vcf.gz"), path("${sample_id}.pbsv.phased.vcf.gz.tbi"), emit: phased_pbsv 
    tuple val(sample_id), path("${sample_id}.trgt.phased.vcf.gz"), path("${sample_id}.trgt.phased.vcf.gz.tbi"), emit: phased_trgt 
    

    script:
    """
    hiphase --version

    hiphase --reference ${reference} \
            --bam ${pbmm2_bam} \
            --vcf ${deepvariant_vcf} \
            --output-vcf ${sample_id}.deepvariant.phased.vcf.gz \
            --vcf ${pbsv_vcf} \
            --output-vcf ${sample_id}.pbsv.phased.vcf.gz \
            --vcf ${trgt_vcf} \
            --output-vcf ${sample_id}.trgt.phased.vcf.gz 
            
            

   

    bcftools index --force  --tbi ${sample_id}.deepvariant.phased.vcf.gz
    bcftools index --force --tbi ${sample_id}.pbsv.phased.vcf.gz
    bcftools index  --force --tbi ${sample_id}.trgt.phased.vcf.gz
    """
}