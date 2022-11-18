library("optparse")
library("rGREAT")

#Parse command line parameters
option_list <- list()

option_list$bed <- make_option('--bed', type='character', help="The path to the peaks bed file to test for enrichment")
option_list$topN <- make_option('--topN', type='integer', help="Number of top peaks to select for enrichment testing")


opt <- parse_args(OptionParser(option_list=option_list,description="Get enriched pathways and phenotypes from genomic peaks"))

bed <- read.delim(opt$bed)
bed <- bed[!grepl("chrGL|chrKI",bed$seqnames),]

#get the enrichment using either the increased or decreased peaks
getEnrichment <- function(allPeaks,dir,topN=5000){
  
  peaks <- allPeaks[ allPeaks$FDR <= 0.01,]
  
  if(dir == "up") {
    peaks <- peaks[ order(peaks$Fold,decreasing = TRUE ),1:3][1:topN,]
    
  } else {
    peaks <- peaks[ order(peaks$Fold,decreasing = FALSE ),1:3][1:topN,]
  }
  
  job <- submitGreatJob(peaks,allPeaks,species = 'hg38',rule="oneClosest")
  
  mousePheno <- getEnrichmentTables(job,ontology="Mouse Phenotype",download_by="tsv")[[1]]
  mousePheno <- mousePheno[ order(mousePheno$HyperP,decreasing = FALSE),]
  write.table(mousePheno,file=sprintf("mousePhenotypes_%s.txt",dir), col.names = T,row.names = F, sep = "\t", quote = F)

  
  humanPheno <- getEnrichmentTables(job,ontology="Human Phenotype",download_by="tsv")[[1]]
  humanPheno <- humanPheno[ order(humanPheno$HyperP,decreasing = FALSE),]
  write.table(humanPheno,file=sprintf("humanPhenotypes_%s.txt",dir), col.names = T,row.names = F, sep = "\t", quote = F)
  
  
  GO <- getEnrichmentTables(job,ontology="GO Biological Process",download_by="tsv")[[1]]
  GO <- GO[ order(GO$HyperP,decreasing = FALSE),]
  write.table(GO,file=sprintf("GO_%s.txt",dir), col.names = T,row.names = F, sep = "\t", quote = F)
  
}

getEnrichment(bed,"up",opt$topN)
getEnrichment(bed,"down",opt$topN)















