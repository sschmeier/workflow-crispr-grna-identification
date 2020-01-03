# PROJECT: crisprseek-grna-id

- AUTHOR: Sebastian Schmeier (s.schmeier@pm.me)
- DATE: 2019 
- VERSION: 0.1.0

## Overview

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow to identify gRNAs and potential offtargets using [CRISPRseek](https://www.bioconductor.org/packages/release/bioc/html/CRISPRseek.html).

The workflow uses user defined transcript annotations in form of [GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4)-files and extracts for each tx potential gRNA for CRISPR experiments (e.g. CRISPRi).

The workflow also requires a genome file in [FASTA](http://genetics.bwh.harvard.edu/pph/FASTA.html)-format. Currently the workflow supports the following genomes: hg38, hg19, mm10 and mm9.


## Installation


```bash
# Install miniconda
# LINUX:
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# MACOSX:
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh

# Install snakemake
conda create -n snakemake snakemake>5.8.0
conda activate snakemake
```

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

# Prefered way for reproducibility: use singularity container
# Singularity needs to be installed system wide
snakemake -p --use-singularity --singularity-args "--bind /mnt/data/" --configfile config.yaml --jobs 16 2> run.log

# Or just use conda without singularity
snakemake -p --use-conda --configfile config.yaml --jobs 16 2> run.log

# show a detailed summary of the produced files and used commands
snakemake -D

# To delete all created result files use
snakemake -p clean --configfile config.yaml
```
