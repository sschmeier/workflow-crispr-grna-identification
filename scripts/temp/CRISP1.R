##CRISPRseek gRNA design for all gRNA from https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/CRISPR/doc/CRISPRdemo.html
#Scenario 7. Quick gRNA finding with gRNA efficacy prediction.

##set working directory
setwd("/Volumes/scratch/diermeierlab/Kat/Automate")

##load libraries
library(CRISPRseek)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

##grab & set bash arguments
args <-	commandArgs()
splitFileVar <-args[6]

##create and set output directory
dir.create(file.path(getwd(),"script1out"), showWarnings = FALSE)
outputDir <- paste(file.path(getwd(),"script1out/"),splitFileVar, sep="")

## input file path
inputFilePath <- paste(file.path(getwd(),"fa/"),splitFileVar, sep="")

## Scenario 7. Quick gRNA finding with gRNA efficacy prediction
results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = FALSE, enable.multicore = TRUE,
                             n.cores.max = 6, annotateExon = FALSE, findPairedgRNAOnly = FALSE, chromToSearch = "",
                             max.mismatch = 0, exportAllgRNAs = "fasta", annotatePaired = FALSE, BSgenomeName = Hsapiens, 
                             outputDir = outputDir, overwrite = TRUE)
                             
