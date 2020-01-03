## Load dplyr package
library(dplyr)

## set working directory
setwd("/Volumes/scratch/diermeierlab/Kat/Automate")

## setting args from bash
args <- commandArgs()
script1var <- args[6]

## creating folders
dir.create(file.path(getwd(),"script2out"), showWarnings = FALSE)
dir.create(file.path(getwd(),"script2out","top25fa"), showWarnings = FALSE)

## renaming gRNAefficacy xls to a txt file
file.rename(paste(getwd(),"script1out", script1var, "gRNAefficacy.xls", sep="/"), paste(getwd(), "script1out", script1var, "gRNAefficacy.txt", sep="/"))

## setting var for xls
gRNAvarCSV <- paste(getwd(), "script1out", script1var, "gRNAefficacy.txt", sep="/")

## Read in data file
geneDat <- read.table(gRNAvarCSV, sep='\t', header=TRUE, stringsAsFactors = FALSE)

## Convert gRNefficiency to numeric (sets text values to NAs)
geneDat$gRNAefficacy <- as.numeric(geneDat$gRNAefficacy)

## Extract gene names from "name" variable
geneDat$geneName <- strsplit(as.vector(geneDat$name),"_") %>% lapply(., function(x) x[1]) %>% unlist()

## Extract the unique gene names
genes <- unique(geneDat$geneName)

## Function to extract the row relating to the most efficient gRNA for a given gene
## getMaxEff <- function(x){
##  geneDat %>% dplyr::filter(geneDat$geneName==x) %>% 
##    dplyr::filter(gRNAefficacy==max(gRNAefficacy, na.rm=TRUE))
## }

## Test of function - gets the data for the first gene
## getMaxEff(genes[1])

## Get the data for most efficient gRNA for each gene
## geneMax <- lapply(genes, function(x) getMaxEff(x))

## Turn the above data into a matrix
## geneMaxMatrix <- geneMax[[1]]
## for(i in 2:length(geneMax)) geneMaxMatrix <- rbind(geneMaxMatrix, geneMax[[i]])

## Write it out to a csv file
## write.csv(geneMaxMatrix, file = paste("script2out/gRNA_max_efficacy_", script1var, ".csv", sep=""), row.names=FALSE)

## Write the sequences to a fasta file
## This is pretty ugly, but it works...

## write.table(paste(">", geneMaxMatrix$name[1]), 
##            file=paste("script2out/gRNA_max_e", script1var, ".fa", sep=""), append=FALSE,
##            col.names=FALSE, row.names=FALSE, quote=FALSE)
## write.table(geneMaxMatrix$gRNAplusPAM[1],
##            file=paste("script2out/gRNA_max_e", script1var, ".fa", sep=""), append=TRUE,
##            col.names=FALSE, row.names=FALSE, quote=FALSE)
## write.table(" ",
##            file=paste("script2out/gRNA_max_e", script1var, ".fa", sep=""), append=TRUE,
##            col.names=FALSE, row.names=FALSE, quote=FALSE)


## for(i in 2:nrow(geneMaxMatrix)){
##  write.table(paste(">", geneMaxMatrix$name[i]), 
##              file=paste("script2out/gRNA_max_e", script1var, ".fa", sep=""), append=TRUE,
##              col.names=FALSE, row.names=FALSE, quote=FALSE)
##  write.table(geneMaxMatrix$gRNAplusPAM[i],
##              file=paste("script2out/gRNA_max_e", script1var, ".fa", sep=""), append=TRUE,
##              col.names=FALSE, row.names=FALSE, quote=FALSE)
##  write.table(" ",
##              file=paste("script2out/gRNA_max_e", script1var, ".fa", sep=""), append=TRUE,
##              col.names=FALSE, row.names=FALSE, quote=FALSE)
  
## }


## Function to get the 25 most efficient gRNAs per gene
getSortEff <- function(x){
  geneDat %>% dplyr::filter(geneDat$geneName==x) %>% 
    dplyr::arrange(., desc(gRNAefficacy)) %>% 
    head(., 25)
}

## Test of function - gets the data for the first gene
getSortEff(genes[1])

## Get the 25 most efficient gRNAs for each gene
geneSort <- lapply(genes, function(x) getSortEff(x))

## Create a matrix of the above data
geneSortMatrix <- geneSort[[1]]
for(i in 2:length(geneSort)) geneSortMatrix <- rbind(geneSortMatrix, geneSort[[i]])

## Write the matrix out to a csv file
write.csv(geneSortMatrix, file = paste("script2out/top25fa/gRNA_top25_", script1var, ".csv", sep=""), row.names=FALSE)

## Write the sequences out to fastq

write.table(paste(">", geneSortMatrix$name[1]), 
            file=paste("script2out/top25fa/gRNA_top25_", script1var, sep=""), append=FALSE,
            col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(geneSortMatrix$gRNAplusPAM[1],
            file=paste("script2out/top25fa/gRNA_top25_", script1var, sep=""), append=TRUE,
            col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(" ",
            file=paste("script2out/top25fa/gRNA_top25_", script1var, sep=""), append=TRUE,
            col.names=FALSE, row.names=FALSE, quote=FALSE)


for(i in 2:nrow(geneSortMatrix)){
  write.table(paste(">", geneSortMatrix$name[i]), 
              file=paste("script2out/top25fa/gRNA_top25_", script1var, sep=""), append=TRUE,
              col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(geneSortMatrix$gRNAplusPAM[i],
              file=paste("script2out/top25fa/gRNA_top25_", script1var, sep=""), append=TRUE,
              col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(" ",
              file=paste("script2out/top25fa/gRNA_top25_", script1var, sep=""), append=TRUE,
              col.names=FALSE, row.names=FALSE, quote=FALSE)
  
}
