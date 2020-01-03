## CRISPRseek gRNA design for all gRNA from:
## https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/CRISPR/doc/CRISPRdemo.html
## Scenario 8. Quick gRNA finding with gRNA efficacy prediction.

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

get_grna <- function(data_path, out_path, threads, genome_version) {
    # R code
    library(CRISPRseek)
    if ( genome_version == "hg38" ) {
        library(BSgenome.Hsapiens.UCSC.hg38)
        ## Scenario 8. Quick gRNA finding with gRNA efficacy prediction
        results <- offTargetAnalysis(data_path,
                                    findgRNAsWithREcutOnly = FALSE,
                                    enable.multicore = TRUE,
                                    n.cores.max = threads,
                                    annotateExon = FALSE,
                                    findPairedgRNAOnly = FALSE,
                                    chromToSearch = "",
                                    max.mismatch = 0,
                                    exportAllgRNAs = "fasta",
                                    annotatePaired = FALSE,
                                    BSgenomeName = Hsapiens,
                                    outputDir = out_path,
                                    overwrite = TRUE)
    } else if ( genome_version == "hg19" ) {
        library(BSgenome.Hsapiens.UCSC.hg19)
        ## Scenario 8. Quick gRNA finding with gRNA efficacy prediction
        results <- offTargetAnalysis(data_path,
                                    findgRNAsWithREcutOnly = FALSE,
                                    enable.multicore = TRUE,
                                    n.cores.max = threads,
                                    annotateExon = FALSE,
                                    findPairedgRNAOnly = FALSE,
                                    chromToSearch = "",
                                    max.mismatch = 0,
                                    exportAllgRNAs = "fasta",
                                    annotatePaired = FALSE,
                                    BSgenomeName = Hsapiens,
                                    outputDir = out_path,
                                    overwrite = TRUE)
    } else if ( genome_version == "mm10" ) {
        library(BSgenome.Mmusculus.UCSC.mm10)
        ## Scenario 8. Quick gRNA finding with gRNA efficacy prediction
        results <- offTargetAnalysis(data_path,
                                    findgRNAsWithREcutOnly = FALSE,
                                    enable.multicore = TRUE,
                                    n.cores.max = threads,
                                    annotateExon = FALSE,
                                    findPairedgRNAOnly = FALSE,
                                    chromToSearch = "",
                                    max.mismatch = 0,
                                    exportAllgRNAs = "fasta",
                                    annotatePaired = FALSE,
                                    BSgenomeName = Mmusculus,
                                    outputDir = out_path,
                                    overwrite = TRUE)
    } else if ( genome_version == "mm9" ) {
        library(BSgenome.Mmusculus.UCSC.mm9)
        ## Scenario 8. Quick gRNA finding with gRNA efficacy prediction
        results <- offTargetAnalysis(data_path,
                                    findgRNAsWithREcutOnly = FALSE,
                                    enable.multicore = TRUE,
                                    n.cores.max = threads,
                                    annotateExon = FALSE,
                                    findPairedgRNAOnly = FALSE,
                                    chromToSearch = "",
                                    max.mismatch = 0,
                                    exportAllgRNAs = "fasta",
                                    annotatePaired = FALSE,
                                    BSgenomeName = Mmusculus,
                                    outputDir = out_path,
                                    overwrite = TRUE)
    } else {
        print("Genome version not yet supported, please contact us.")
        stop("EXIT") 
    }
}

get_grna(snakemake@input[[1]], 
         snakemake@output[[1]], 
         snakemake@threads, 
         snakemake@params[["genome_version"]])




