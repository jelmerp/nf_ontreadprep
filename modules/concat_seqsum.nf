
process CONCAT_SEQSUM {
    publishDir "${params.outdir}/guppy", mode: 'copy'
    label "process_low"
    
    input:
    path seqsum_files
    
    output:
    path "seqsum.txt", emit: seqsum
    
    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' ${seqsum_files} > seqsum.txt
    """
}
