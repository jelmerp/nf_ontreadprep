
process GUPPY {
    publishDir "${params.outdir}/guppy", mode: 'copy', pattern: "logs/*"

    input:
    tuple val(file_id), path(fast5)
    path guppy_config
    
    output:
    path "*fastq.gz", emit: fq
    path "*sequencing_summary.txt", emit: seqsum
    path "logs/*", emit: log
    
    script:
    """
    guppy_gpu.sh \
        --infile ${fast5} \
        --config ${guppy_config} \
        --outdir .

    mv -v pass/*fastq.gz ${file_id}.fastq.gz
    mv -v sequencing_summary.txt ${file_id}_sequencing_summary.txt

    mv -v  guppy*log logs/guppy_${file_id}.log
    cp -v .command.log logs/slurm_${file_id}.log
    """
}
