process PYCOQC {
    publishDir "${params.outdir}/pycoqc", mode: 'copy'

    //TODO Quality argument

    input:
    path seqsum
    val min_pass_len
    
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
        ${len_arg}

    pycoQC --version > version.log
    cp .command.log command.log
    """
}
