
process CONCAT_FQ {
    publishDir "${params.outdir}/guppy", mode: 'copy'
    label "process_low"
    
    input:
    path fq_files
    val run_id
    
    output:
    path "*.fastq.gz", emit: fq
    
    script:
    """
    cat ${fq_files} > ${run_id}.fastq.gz
    """
}
