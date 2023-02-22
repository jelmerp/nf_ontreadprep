
process RM_ORGANEL {
    publishDir "${params.outdir}/rm_organel", mode: 'copy'

    input:
    path fq_in
    path ref_assembly
    val seqids

    output:
    path "*fastq.gz", emit: fq
    path "logs/*", emit: log
    
    script:
    """
    rm_organel.sh \
        --fq_in ${fq_in} \
        --fq_out concat_rmorganel.fastq.gz \
        --ref ${ref_assembly} \
        --seqids ${organel_contigs}

    cp -v .command.log logs/slurm.log
    """
}
