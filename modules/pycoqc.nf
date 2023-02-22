process PYCOQC {
    publishDir "${params.outdir}/pycoqc", mode: 'copy'

    input:
    path seqsum
    val min_pass_len
    val min_pass_qual
    
    output:
    path "*html", emit: report
    path "*log", emit: log
    
    script:
    len_arg = min_pass_len ? "--min_pass_len ${min_pass_len}" : ""
    outfile = min_pass_len ? "pycoqc_minlen.html" : "pycoqc.html" 
    """
    pycoQC \
        -f ${seqsum} \
        -o ${outfile} \
        --min_pass_qual ${min_pass_qual} \
        ${len_arg}

    pycoQC --version > version.log
    cp .command.log command.log
    """
}
