##CRISPRseek Scenario 5: Target and off-target analysis for user specified gRNAs to identify off-target analysis:

##set working directory
setwd("/Volumes/scratch/diermeierlab/Kat/Automate")

##load libraries
library(CRISPRseek)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

## setting args from bash
args <- commandArgs()
top25var <- args[6]

##create output directory
dir.create(file.path(getwd(),"script3out"), showWarnings = FALSE)

##set output directory
outputDir <- file.path("script3out",top25var)

## input file path
gRNAFilePath <- paste(file.path("script2out", "top25fa/"), top25var, sep="") 

## Scenario 5: Target and off-target analysis for user specified gRNAs
results <- offTargetAnalysis(inputFilePath = gRNAFilePath, enable.multicore = FALSE, n.cores.max = 1, annotateExon = FALSE, findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE, findgRNAs = FALSE, BSgenomeName = Hsapiens, chromToSearch = "all", txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, orgAnn = org.Hs.egSYMBOL, max.mismatch = 0, outputDir = outputDir, overwrite = TRUE)

##If need OffTarget.RDS for each file:
##outfile <- gsub(".fa", ".rds", top25var)
##system(paste('mv offTargets.RDS ',outfile))
