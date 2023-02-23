
# ONT Read Prep Workflow

This small, simple Nextflow workflow will take as its main input a directory of
Oxford Nanopore Technologies (ONT) `FAST5` files, and will:

- Run `Guppy` in GPU mode to base-call and produce `FASTQ` files
  (one file per input `FAST5`)

- Concatenate all `FASTQ` files into 1 large, gzipped file

- Concatenates all `sequencing_summary` files produced by Guppy

- Run `PycoQC` for read QC using the concatenated `sequencing_summary` file.
  (Run twice, once with a minimum read length, to get stats for longer reads.)

- Optionally, remove reads from the concatenated `FASTQ` that map to specific contigs
  (in practice, I am using this to remove organelle-derived reads prior
  to assembling a plant genome with these reads).
  To do so, you'll have to input a reference genome nucleotide FASTA file,
  and a list of "blacklisted" contig/scaffold IDs.

## Workflow usage

After installing Nextflow,
run `nextflow run https://github.com/jelmerp/nf_ontreadprep` with some of the
following options:

```
============================================================================
                  O N T    R E A D   P R E P   W O R K F L O W
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
```

## Shell wrapper script

If you're at the Ohio Supercomputer Center,
you can also run the shell wrapper script `nf_ontreadprep.sh`,
using `sbatch -A <PROJECT-NAME> nf_ontreadprep.sh` and some of the following
options:

```
        nf_ontreadprep.sh (v. 1.0): Run nf_ontreadprep
        ==============================================
DESCRIPTION:
  Nextflow workflow to prepare ONT FAST5 reads (e.g. basecalling, QC)

USAGE:
  sbatch nf_ontreadprep.sh -i <input-file> -o <output-dir> [...]
  bash nf_ontreadprep.sh -h

REQUIRED OPTIONS:
  -i/--fast5_dir  <dir>   Input dir with FAST5 files
  --guppy_config  <file>  Guppy config file

OTHER KEY OPTIONS:
  --ref_assembly  <file>  Reference genome assembly nucleotide FASTA file
  --contig_blacklist<str> Comma-separated list of contigs to remove reads for
  -o/--outdir     <dir>   Output dir (will be created if needed)
  --more_args     <str>   Quoted string with more argument(s) for nf_ontreadprep

NEXTFLOW-RELATED OPTIONS:
  -profile        <str>  Profile(s) to use from config files          [default: 'singularity']
  -work-dir       <dir>  Scratch (work) dir for the workflow          [default: '/fs/scratch/<OSC-proj>/<user>/ont_readprep']
  -c/-config      <file  Additional config file                       [default: none]
                           - Settings in this file will override default settings
                           - Note that the mcic-scripts OSC config file will always be included, too
                             (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)

UTILITY OPTIONS:
  -h                      Print this help message and exit
  --help                  Print the help for nf_ontreadprep and exit
  -v                      Print the version of this script and exit
  -v/--version            Print the version of nf_ontreadprep and exit

EXAMPLE COMMANDS:
  sbatch nf_ontreadprep.sh -i data/fast5 --guppy_config config/dna_r10.4.1_e8.2_260bps_sup.cfg
```
