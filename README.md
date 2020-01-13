[![pipeline status](https://gitlab.com/schmeierlab/workflows/crisprseek-grna-id/badges/master/pipeline.svg)](https://gitlab.com/schmeierlab/workflows/crisprseek-grna-id/commits/master)

# PROJECT: workflow-crispr-grna-identification

- AUTHOR: Sebastian Schmeier (s.schmeier@pm.me)
- DATE: 2019 
- VERSION: 0.2.0

## Overview

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow to identify gRNAs and potential offtargets using [CRISPRseek](https://www.bioconductor.org/packages/release/bioc/html/CRISPRseek.html).

The workflow uses user defined transcript annotations in form of [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4)-files and extracts for each transcript potential gRNAs for CRISPR experiments (e.g. CRISPRi).

The workflow also requires a genome file in [FASTA](http://genetics.bwh.harvard.edu/pph/FASTA.html)-format.
Currently the workflow supports the following genomes: hg38, hg19, mm10 and mm9.


## Workflow details

The following steps will be performed by the workflow:

1. Compile a set of unique transcripts from all supplied GTF-files.
2. Get the 5'-coordinate (TSS) of each transcript.
3. Get the user-defined upstream and downstream region around the TSS (e.g. -1000, +500). This is the region that will be searched for gRNAs.
4. Extract sequence in fasta for each coordinate from step 3.
5. Run [CRISPRseek](https://www.bioconductor.org/packages/release/bioc/html/CRISPRseek.html) on each sequence.
6. Extract user-defined number of top gRNA per transcript.
7. Run offtarget analysis for each gRNA found in step 6.
8. Select user-defined top gRNAs per transcript based on efficiancy and offtarget analysis.


## Installation

### Install Miniconda

```bash
# Install miniconda
# LINUX:
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# MACOSX:
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

### Make a Snakemake environment

```bash
# Install snakemake
conda create -n snakemake snakemake>5.8.0
conda activate snakemake
```

### Clone workflow

```bash
git clone https://gitlab.com/schmeierlab/workflows/workflow-crispr-grna-identification.git
cd workflow-crispr-grna-identification
```


## Running the workflow

### Create a file-sheet

Create a files.txt with path to GTF-files. 
See the example for details on the file.

### Adjust config.yaml

Currently it is required to submit config parameters via `--configfile config.yaml`.
Change `config.yaml` accordingly. In particular:

1. Add the correct path to the sheet in config.
2. Add the correct path to the genome FASTA-file 
3. Add correct genome version. 

### Execute workflow

```bash
# Do a dryrun of the workflow, show rules, order, and commands
snakemake -np --configfile config.yaml

# Prefered way for reproducibility: use Singularity container.
# Singularity needs to be installed system wide.
#
# Bind dir outside your home (e.g. here: /mnt/data/) with:
# --singularity-args "--bind /mnt/data/"
#
snakemake -p --use-singularity --configfile config.yaml --jobs 32 2> run.log

# Or just use conda without Singularity
snakemake -p --use-conda --configfile config.yaml --jobs 32 2> run.log

# To delete all created result files use
snakemake -p clean --configfile config.yaml
```


## Testing the workflow

### Singularity-mode

Needs a working [Singularity](https://sylabs.io/singularity/) installation on the system.

If necessary, bind direcories outside your working dir (here: `/mnt/data`) with `--singularity-args "--bind /mnt/data/"`.

```bash
snakemake -p --configfile tests/test_config.yaml --use-singularity 2> test.log
```


### Conda-mode

```bash
snakemake -p --configfile tests/test_config.yaml --use-conda 2> test.log
```
