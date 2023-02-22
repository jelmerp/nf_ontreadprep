
process CONCAT_FQ {
    publishDir "${params.outdir}/guppy", mode: 'copy'

    input:
    path fq_files
    
    output:
    path "concat.fastq.gz", emit: fq
    
    script:
    """
    cat ${fq_files} > concat.fastq.gz
    """
}
