
process PORECHOP {
    publishDir "${params.outdir}/porechop", mode: 'copy', pattern: "logs/*"
    label 'process_high_memory'
    
    input:
    path fastq
    
    output:
    path "out/*fastq.gz", emit: fq
    path "out/logs/*", emit: log
    
    script:
    """
    porechop.sh \
        --infile ${fastq} \
        --outdir out

    cp -v .command.log out/logs/slurm.log
    """
}
