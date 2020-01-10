## CRISPRseek gRNA design for all gRNA from:
## https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/CRISPR/doc/CRISPRdemo.html
## Scenario 5: Target and off-target analysis for user specified gRNAs

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
# for the paranoid
set.seed(42)

offtarget_run <- function(data_path, out_path, threads, genome_version, missmatches, smethod) {
    # R code
    library(CRISPRseek)

    if ( genome_version == "hg38" ) {
        library(BSgenome.Hsapiens.UCSC.hg38)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        library(org.Hs.eg.db)
        results <- offTargetAnalysis(inputFilePath = data_path, 
                                    enable.multicore = FALSE, 
                                    n.cores.max = 1, 
                                    annotateExon = FALSE, 
                                    findgRNAs = FALSE, # this will only use the specifie gRNAs and skip the search
                                    findgRNAsWithREcutOnly = FALSE, 
                                    findPairedgRNAOnly = FALSE, 
                                    chromToSearch = "all", 
                                    BSgenomeName = Hsapiens, # not sure how to do this with a var
                                    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene, # not sure how to do this with a var
                                    orgAnn = org.Hs.egSYMBOL, # not sure how to do this with a var
                                    max.mismatch = missmatches,
				    scoring.method = smethod,
                                    outputDir = out_path, 
                                    overwrite = TRUE)

    } else if ( genome_version == "hg19" ) {
        library(BSgenome.Hsapiens.UCSC.hg19)
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        library(org.Hs.eg.db)
        results <- offTargetAnalysis(inputFilePath = data_path, 
                                        enable.multicore = FALSE, 
                                        n.cores.max = 1, 
                                        annotateExon = FALSE, 
                                        findgRNAs = FALSE, # this will only use the specifie gRNAs and skip the search
                                        findgRNAsWithREcutOnly = FALSE, 
                                        findPairedgRNAOnly = FALSE, 
                                        chromToSearch = "all", 
                                        BSgenomeName = Hsapiens, # not sure how to do this with a var
                                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, # not sure how to do this with a var
                                        orgAnn = org.Hs.egSYMBOL, # not sure how to do this with a var
                                        max.mismatch = missmatches,
     					scoring.method = smethod,					
                                        outputDir = out_path, 
                                        overwrite = TRUE)
     } else if ( genome_version == "mm9" ) {
        library(BSgenome.Mmusculus.UCSC.mm9)
        library(TxDb.Mmusculus.UCSC.mm9.knownGene)
        library(org.Mmu.eg.db)
        results <- offTargetAnalysis(inputFilePath = data_path, 
                                        enable.multicore = FALSE, 
                                        n.cores.max = 1, 
                                        annotateExon = FALSE, 
                                        findgRNAs = FALSE, # this will only use the specifie gRNAs and skip the search
                                        findgRNAsWithREcutOnly = FALSE, 
                                        findPairedgRNAOnly = FALSE, 
                                        chromToSearch = "all", 
                                        BSgenomeName = Mmusculus, 
                                        txdb = TxDb.Mmusculus.UCSC.mm9.knownGene, 
                                        orgAnn = org.Mmu.egSYMBOL, 
                                        max.mismatch = missmatches,
     					scoring.method = smethod,
                                        outputDir = out_path, 
                                        overwrite = TRUE)
    } else if ( genome_version == "mm10" ) {
        library(BSgenome.Mmusculus.UCSC.mm10)
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
        library(org.Mmu.eg.db)
        results <- offTargetAnalysis(inputFilePath = data_path, 
                                        enable.multicore = FALSE, 
                                        n.cores.max = 1, 
                                        annotateExon = FALSE, 
                                        findgRNAs = FALSE, # this will only use the specifie gRNAs and skip the search
                                        findgRNAsWithREcutOnly = FALSE, 
                                        findPairedgRNAOnly = FALSE, 
                                        chromToSearch = "all", 
                                        BSgenomeName = Mmusculus, 
                                        txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                        orgAnn = org.Mmu.egSYMBOL, 
                                        max.mismatch = missmatches,
					scoring.method = smethod,
                                        outputDir = out_path, 
                                        overwrite = TRUE)
    } else {
        print("Genome version not yet supported, please contact us.")
        stop("EXIT") 
    }
}

offtarget_run(snakemake@input[[1]], 
              snakemake@output[[1]], 
              snakemake@threads, 
              snakemake@params[["genome_version"]],
              snakemake@params[["missmatches"]],
	      snakemake@params[["scoring_method"]])


