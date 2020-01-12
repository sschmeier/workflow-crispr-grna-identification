## =============================================================================
## WORKFLOW PROJECT: get-crispr-grna
## INIT DATE: 2019
import pandas as pd
import glob, os, os.path, datetime, sys, csv
from os.path import join, abspath
from snakemake.utils import validate, min_version

## =============================================================================
## set minimum snakemake version #####
min_version("5.8.0")

## =============================================================================
## SETUP
## =============================================================================

## SET SOME DEFAULT PATH
DIR_SCRIPTS = abspath("scripts")
DIR_REPORT  = abspath("report")
DIR_ENVS    = abspath("envs")
DIR_SCHEMAS = abspath("schemas")

## LOAD VARIABLES FROM CONFIGFILE
## submit on command-line via --configfile
if config=={}:
    print('Please submit config-file with "--configfile <file>". Exit.')
    sys.exit(1)

sys.stderr.write("********** Submitted config: **********\n")
for k,v in config.items():
    sys.stderr.write("{}: {}\n".format(k,v))
sys.stderr.write("***************************************\n")

## Validate configfile with yaml-schema
validate(config, schema=join(DIR_SCHEMAS, "config.schema.yaml"))

## define global Singularity image for reproducibility
## USE: "--use-singularity" to run all jobs in container
singularity: "shub://sschmeier/container-crisprseek:0.0.2"
#singularity: "docker://continuumio/miniconda3"

## Setup result dirs
DIR_BASE       = abspath(config["resultdir"])
DIR_LOGS       = join(DIR_BASE, "logs")
DIR_BENCHMARKS = join(DIR_BASE, "benchmarks")
DIR_RES        = join(DIR_BASE, "results")

## Workflow specific setup
CFG_GENOME     = abspath(config["ref"]["genome"])
CFG_GENOME_VER = config["ref"]["version"]
if CFG_GENOME_VER not in ["hg38", "hg19", "mm10", "mm9"]:
    print("Error: Currently only hg19, hg38, mm10, and mm9 are supported. EXIT.")
    sys.exit(1)

CFG_UP   = config["gtf"]["upstream"]
CFG_DN   = config["gtf"]["dnstream"]
CFG_NUM  = config["crispr"]["max_num_test"]
CFG_MIN  = config["crispr"]["min_efficacy"]
CFG_NUM_FINAL = config["crispr"]["num_select"]
CFG_MM = config["crispr"]["num_mm"]
CFG_SM = config["crispr"]["scoring_method"]
if CFG_SM not in ["Hsu-Zhang", "CFDscore"]:
    print("Error: Scoring method needs to be one of Hsu-Zhang or CFDscore. EXIT.")
    sys.exit(1)
    
## =============================================================================
## SAMPLES
files = pd.read_csv(config["gtf"]["files"], sep=",").set_index("gtf", drop=False)
validate(files, schema=join(DIR_SCHEMAS, "files.schema.yaml"))

## reading samplename from samplesheet
sys.stderr.write('Reading files from sheet: "{}" ...\n'.format(config["gtf"]["files"]))

## test if sample in dir
for fname in files["gtf"]:
    if not os.path.isfile(fname):
        sys.stderr.write("File '{}' from sheet can not be found. Make sure the file exists. Exit\n".format(fname))
        sys.exit()


## =============================================================================
## SETUP FINAL TARGETS
## =============================================================================
TARGETS = join(DIR_RES, "stats.txt")

## =============================================================================
## FUNCTIONS
## =============================================================================

## =============================================================================
## RULES
## =============================================================================

## Pseudo-rule to state the final targets, so that the whole runs
rule all:
    input:
        TARGETS


## 1. Extract from gtfs a set of unique tx based on 5' positions
rule extract_unique_5prime:
    input:
        files=files["gtf"]
    output:
        join(DIR_RES, "unique_tx_pos.bed")
    log:
        join(DIR_LOGS, "extract_unique_5prime.log")
    benchmark:
        join(DIR_BENCHMARKS, "extract_unique_5prime.txt")
    params:
        script=join(DIR_SCRIPTS, "get_unique_5prime.py")
    shell:
         "cat {input.files} | python {params.script} - > {output} 2> {log}"
    

# 2. extract genome sizes for respective genome
rule extract_genomesizes:
    input:
        CFG_GENOME
    output:
        join(DIR_RES, "genome.sizes.txt")
    log:
        join(DIR_LOGS, "extract_genomesizes.log")
    benchmark:
        join(DIR_BENCHMARKS, "extract_genomesizes.txt")
    conda:
        join(DIR_ENVS, "biopy.yaml")
    params:
        script=join(DIR_SCRIPTS, "get_genome_sizes.py"),
        extra=""
    shell:
        "python {params.script} {params.extra} {input} > {output} 2> {log}"
    

## 3. Get upstream and dn-stream around TSS
rule extract_coordinates:
    input:
        bed=join(DIR_RES, "unique_tx_pos.bed"),
        sizes=join(DIR_RES, "genome.sizes.txt")
    output:
        join(DIR_RES, "unique_tx_pos_updn.bed")
    log:
        join(DIR_LOGS, "extract_coordinates.log")
    benchmark:
        join(DIR_BENCHMARKS, "extract_coordinates.txt")
    conda:
        join(DIR_ENVS, "gff.yaml")
    shell:
         "cat {input.bed} | bedtools slop -r {CFG_DN} -l {CFG_UP}"
         " -s -i stdin -g {input.sizes} > {output} 2> {log}"
    

## 4. Extract seqs
rule extract_fasta:
    input:
        join(DIR_RES, "unique_tx_pos_updn.bed")
    output:
        join(DIR_RES, "unique_tx_pos_updn.fa")
    log:
        join(DIR_LOGS, "extract_fasta.log")
    benchmark:
        join(DIR_BENCHMARKS, "extract_fasta.txt")
    conda:
        join(DIR_ENVS, "gff.yaml")
    shell:
         "bedtools getfasta -nameOnly -fi {CFG_GENOME} -bed {input} > {output} 2> {log}"


## 5. split fasta
checkpoint split_fasta:
    input:
        join(DIR_RES, "unique_tx_pos_updn.fa")
    output:
        dir=directory(join(DIR_RES, "00_tmp")),
        info=join(DIR_RES, "location_info.tsv")
    log:
        join(DIR_LOGS, "split_fasta.log")
    benchmark:
        join(DIR_BENCHMARKS, "split_fasta.txt")
    conda:
        join(DIR_ENVS, "biopy.yaml")
    params:
        script=join(DIR_SCRIPTS, "split_seqs.py"),
        extra=""
    shell:
        "python {params.script} {params.extra} {input} {output.dir} {output.info} 2> {log}"


rule crispr:
    # Find gRNAs in the regions around the TSS
    input:
        join(DIR_RES, "00_tmp/{i}.fa")
    output:
        directory(join(DIR_RES, "01_crispr/{i}"))
    log:
        join(DIR_LOGS, "crispr/{i}.log")
    benchmark:
        join(DIR_BENCHMARKS, "crispr/{i}.txt")
    priority: 10
    conda:
        join(DIR_ENVS, "r.yaml")
    threads: 1  # quick with one seq ~ under 1 min
    params:
        genome_version=CFG_GENOME_VER
    script:
        join(DIR_SCRIPTS, "crisprseek_snake.R")


rule select_bestEff_grna:
    # From all gRNAs select the ones with best efficacy
    input:
        join(DIR_RES, "01_crispr/{i}")
    output:
        eff=join(DIR_RES, "02_crispr_bestEff/{i}.tsv"),
        fa=join(DIR_RES, "02_crispr_bestEff/{i}.fa")
    log:
        join(DIR_LOGS, "select_bestEff_grna/{i}.log")
    benchmark:
        join(DIR_BENCHMARKS, "select_bestEff_grna/{i}.txt")
    priority: 9
    threads: 1  # quick with one seq ~ under 1 min
    conda:
        join(DIR_ENVS, "pandas.yaml")
    params:
        in_eff=join(DIR_RES, "01_crispr/{i}/gRNAefficacy.xls"),
        num=CFG_NUM,
        min=CFG_MIN
    script:
        join(DIR_SCRIPTS, "select_top_grna.py")


rule crispr_offtargets:
    # For the selected gRNAs in previous step, run offtarget analysis
    input:
        join(DIR_RES, "02_crispr_bestEff/{i}.fa")
    output:
        directory(join(DIR_RES, "03_crispr_offtarget/{i}"))
    log:
        join(DIR_LOGS, "crispr_offtargets/{i}.log")
    benchmark:
        join(DIR_BENCHMARKS, "crispr_offtargets/{i}.txt")
    priority: 8
    conda:
        join(DIR_ENVS, "r.yaml")
    threads: 1  # quick with one seq ~ under 1 min
    params:
        genome_version=CFG_GENOME_VER,
        missmatches=CFG_MM,
        scoring_method=CFG_SM
    script:
        join(DIR_SCRIPTS, "crisprseek_offtarget_snake.R")


rule select_gRNA_min_offtargets:
    # 1. Select all gRNA withput offtargets, sort by best efficacy
    # 2. Select from those the top number
    # 3. If number cannot be filled as here are not enough gRNAs 
    #    withoput offtargets, do
    #    
    #    - Calculate number of offtargets for the ones with offtargets
    #    - Select the ones with minimum offtargets 
    #    - stop if we have the requred number selected or running out of gRNAs
    input:
        offt_dir=join(DIR_RES, "03_crispr_offtarget/{i}"),
        grna_file=join(DIR_RES, "02_crispr_bestEff/{i}.tsv")
    output:
        join(DIR_RES, "04_crispr_final_gRNA/{i}.tsv")
    log:
        join(DIR_LOGS, "select_gRNA_min_offtargets/{i}.log")
    benchmark:
        join(DIR_BENCHMARKS, "select_gRNA_min_offtargets/{i}.txt")
    priority: 7
    threads: 1  # quick with one seq ~ under 1 min
    conda:
        join(DIR_ENVS, "pandas.yaml")
    params:
        infile_summary=join(DIR_RES, "03_crispr_offtarget/{i}/Summary.xls"),
        num=CFG_NUM_FINAL
    script:
        join(DIR_SCRIPTS, "select_top_grna_offtargets.py")


def aggregate_input(wildcards):
    '''
    aggregate the file names of the fa-files
    generated at the split step
    '''
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    #return expand(join(DIR_RES, "03_crispr_offtarget/{i}"),     
    return expand(join(DIR_RES, "04_crispr_final_gRNA/{i}.tsv"),
                  i=glob_wildcards(join(checkpoint_output, '{i}.fa')).i)


# rule collect:
#     input:
#         aggregate_input,
#     output:
#         touch(join(DIR_RES, "collect.done.txt"))

rule stats:
    input:
        infiles=aggregate_input,
        info=join(DIR_RES, "location_info.tsv"),
        bed=join(DIR_RES, "unique_tx_pos_updn.bed")
    output:
        join(DIR_RES, "stats.txt")
    log:
        join(DIR_LOGS, "stats.log")
    benchmark:
        join(DIR_BENCHMARKS, "stats.txt")
    threads: 1
    conda:
        join(DIR_ENVS, "pandas.yaml")
    params:
        script=join(DIR_SCRIPTS, "stats_locations.py"),
        extra=""
    shell:
        "python {params.script} {params.extra} -o {output}"
        " {input.bed} {input.info} {input.infiles} 2> {log}"
            

rule clean:
    shell:
        "rm -rf {DIR_BASE}/*"
