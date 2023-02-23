#!/usr/bin/env nextflow

// Nextflow workflow to prep ONT sequences for assembly

//TODO - Don't use symlinks for scripts! Won't be included in repo
//TODO - Get base config to work!

// Constants
WF_VERSION=0.1

// Help function
def helpMessage() {
    log.info"""
    ============================================================================
                            A S S E M B L Y   P R E P   W O R K F L O W
    ============================================================================
    REQUIRED OPTIONS:
      --fast5_dir       <dir>   Dir with input FAST5 files
      --guppy_config    <file>  Guppy config file
                                  - The appropriate config file depends on the flowcell + kit combination.
                                  - Modify this file to e.g. set a different qual-score threshold. 
    
    OPTIONAL DATA I/O OPTIONS:
      --outdir          <dir>   Final output dir for workflow results                           [default: 'results/nf_ontreadprep']
      --ref_assembly    <file>  Reference genome assembly (for organel/contig removal)          [default: none]
      --contig_blacklist <str>  Comma-separated string with contig names from '--ref_assembly'. [default: none]
                                  - Reads that map to these contigs will be removed, but only
                                    when both '--ref_assembly' and '--organel_contigs' are used.
      --pyco_minq       <int>   Min. quality score for PycoQC to consider a read 'passed'.      [default: 10]
      --pyco_minlen     <int>   Min. read length in bp for PycoQC to consider a read 'passed'.  [default: 1000]
                                  - NOTE: PycoQC is run twice: with and without this threshold.
    
    UTILITY OPTIONS
      --help                      Print this help message and exit
      --version                   Print the workflow version and exit
    """.stripIndent()
    exit 0
}

def versionMessage() {
    nf_ontreadprep workflow version ${WF_VERSION}
}

// =============================================================================
//                           GENERAL SETUP
// =============================================================================
// Print help if required
if (params.help) helpMessage()
if (params.version) versionMessage()

// Include modules
include { GUPPY } from './modules/guppy'
include { CONCAT_FQ } from './modules/concat_fq'
include { CONCAT_SEQSUM } from './modules/concat_seqsum'
include { PYCOQC as PYCOQC_ALL } from './modules/pycoqc'
include { PYCOQC as PYCOQC_1K } from './modules/pycoqc'
include { READFILTER_CONTIGBLACKLIST } from './modules/readfilter_contigblacklist'

// Check parameters
if (!params.fast5_dir) { exit 1, '\nERROR: Input dir not specified! Use "--fast5_dir"\n' }
if (!params.guppy_config) { exit 1, '\nERROR: Guppy config not specified! Use "--guppy_config"\n' }

def check_path_list = [ params.fast5_dir, params.ref_assembly, params.guppy_config ]
for (par in check_path_list) { if (par) { file(par, checkIfExists: true) } }

// =============================================================================
//                            DETERMINE WHAT TO RUN
// =============================================================================
// Process parameters
fast5_files = params.fast5_dir + "/" + "*fast5"

// Determine what to run
filter_reads = params.ref_assembly && params.contig_blacklist ? true : false

// =============================================================================
//                              REPORT
// =============================================================================
// Run info
log.info ""
log.info "======================================================================"
log.info "          A S S E M B L Y   P R E P   W O R K F L O W"
log.info "======================================================================"
log.info "Parameters:"
params.each { k, v -> if (v) { println "${k}: ${v}" } }
log.info ""
println "Remove reads mapping to blacklisted contigs?        ${filter_reads}"
log.info "====================================================================\n"

// =============================================================================
//                            SET-UP INPUT CHANNELS
// =============================================================================
// Reads
ch_reads = Channel
    .fromPath( fast5_files, checkIfExists: true )
    .map { file -> tuple(file.baseName, file) }
println "Number of input files:"; ch_reads.count().view()

// Guppy config
ch_guppy_config = Channel.fromPath(params.guppy_config, checkIfExists: true)
ch_guppy_config = ch_guppy_config.collect()

// Reference genomes
ch_ref = Channel.fromPath(params.ref_assembly, checkIfExists: true)

// =============================================================================
//                            THE WORKFLOW
// =============================================================================
// Workflow
workflow {

    // Guppy basecalling
    GUPPY(ch_reads, ch_guppy_config)
        
    // Concatenate basecalled FASTQs
    ch_fq_separate = GUPPY.out.fq.collect()
    ch_fq = CONCAT_FQ(ch_fq_separate).fq

    // Concatenate seqsums
    ch_seqsum_separate = GUPPY.out.seqsum.collect()
    ch_seqsum = CONCAT_SEQSUM(ch_seqsum_separate).seqsum

    // Read QC
    PYCOQC_ALL(ch_seqsum, 0, params.pyco_minq)
    PYCOQC_1K(ch_seqsum, params.pyco_minlen, params.pyco_minq)

    if (filter_reads) {
        READFILTER_CONTIGBLACKLIST(ch_fq, ch_ref, params.contig_blacklist)
    }
}

// Report completion/failure of workflow
workflow.onComplete {
    println ( workflow.success ? """
        ========================================================================
                            PIPELINE SUCCESSFULLY COMPLETED
        ========================================================================
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        ========================================================================
                                    PIPELINE FAILED
        ========================================================================
        """
        .stripIndent()
    )
}
