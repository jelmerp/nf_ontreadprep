
process READFILTER_CONTIGBLACKLIST {
    publishDir "${params.outdir}/readfilter_contigblacklist", mode: 'copy'
    label "process_high"

    input:
    path fq_in              // FASTQ file
    path ref_assembly       // Ref. assembly nucleotide FASTA file
    val seq_ids             // Contig/scaffold IDs: reads mapping to these will be removed

    output:
    path "*fastq.gz", emit: fq
    path "logs/*", emit: log
    
    script:
    """
    rm_organel.sh \
        --fq_in ${fq_in} \
        --fq_out concat_filtered.fastq.gz \
        --ref ${ref_assembly} \
        --seqids ${seq_ids}

    cp -v .command.log logs/slurm.log
    """
}
