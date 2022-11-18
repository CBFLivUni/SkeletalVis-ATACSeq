library("optparse")
library("tidyverse")
library("DiffBind")
library("BiocParallel")

#Parse command line parameters
option_list <- list()

option_list$peaks <- make_option('--peaks', type='character' ,help="The path to the narrow peaks from MACS2")
option_list$bams <- make_option('--bams', type='character', help="The path to the bam files to merge per condition")
option_list$sampleTable <- make_option('--sampleTable', type='character', help="sampleTable mapping sample to condition")
option_list$diffPeakTable <- make_option('--diffPeakTable', type='character', help="Results table of peaks")
option_list$cpuCores <- make_option('--cpuCores', type='integer', help="Number of cpu cores to use")
option_list$comparison <- make_option('--comparison', type='character', help="Factor to use in contrast")
option_list$numerator <- make_option('--numerator', type='character', help="Factor to use in numerator")
option_list$denominator<- make_option('--denominator', type='character', help="Factor to use in denominator")
option_list$useConsensus<- make_option('--useConsensus', type='logical', default=FALSE,help="Should the peaks be filtered so they are present in at least 2 of the tissue types?")


opt <- parse_args(OptionParser(option_list=option_list,description="Identify differential peaks"))

############
#run Diffbind
############


getTable <- function(fileList,fileType){
  fileList <- unlist(strsplit(fileList, ' '))
  names(fileList) <- tools::file_path_sans_ext(basename(fileList))
  fileTable <- stack(fileList)
  colnames(fileTable) <- c(fileType,"SampleID")
  fileTable$SampleID <- gsub("_peaks|_noMT_noDups_unique","",  fileTable$SampleID)
  return(fileTable)
}

bams <- getTable(opt$bams,"bamReads")

#exclude the indexes if present 
bams <- bams[ !grepl(".bai",bams$bamReads),]
peaks <- getTable(opt$peaks,"Peaks")

print(bams)
print(peaks)

sampleTable <- read.delim(opt$sampleTable)
colnames(sampleTable)[1] <- "SampleID"

print(sampleTable)
#create the design formula
factors <- colnames(sampleTable)[-1:-2]
designFormula <- as.formula(paste("~", paste(factors, 
 collapse = " + ")))


sampleTable <- merge(sampleTable,bams,by="SampleID")
sampleTable <- merge(sampleTable,peaks,by="SampleID")

print(sampleTable)
sampleTable$PeakCaller <- "bed"



sampleSheet <- dba(sampleSheet = sampleTable, minOverlap = 2,config=data.frame(AnalysisMethod=DBA_DESEQ2,th=1,
 DataType=DBA_DATA_GRANGES, RunParallel=FALSE,
 cores=opt$cpuCores,
 minQCth=15, fragmentSize=125,
 bCorPlot=FALSE, reportInit="DBA",
 bUsePval=FALSE, design=TRUE,
 doBlacklist=TRUE, doGreylist=TRUE))

png("allPeaksCorrelation.png")
correlationPlot <- plot(sampleSheet)
dev.off()

if(opt$useConsensus) {
  consensus <- dba.peakset(sampleSheet, consensus=DBA_TISSUE, minOverlap=2)

  png("consensusOverlap.png")
  cons.ol <- dba.plotVenn(consensus, consensus$masks$Consensus)
  dev.off()

  consensus <-  dba(consensus, mask=consensus$masks$Consensus,config=data.frame(AnalysisMethod=DBA_DESEQ2,th=0.05,
    DataType=DBA_DATA_GRANGES, RunParallel=FALSE,
    cores=opt$cpuCores,
    minQCth=15, fragmentSize=125,
    bCorPlot=FALSE, reportInit="DBA",
    bUsePval=FALSE, design=TRUE,
    doBlacklist=TRUE, doGreylist=TRUE))
  consensusPeaks <- dba.peakset(consensus, bRetrieve=TRUE)
  saveRDS(consensusPeaks,file="consensusPeaks.RDS")


  
  df <- data.frame(seqnames=seqnames(consensusPeaks),
   starts=start(consensusPeaks)-1,
   ends=end(consensusPeaks),
   names=c(rep(".", length(consensusPeaks))),
   scores=c(rep(".", length(consensusPeaks))),
   strands=strand(consensusPeaks),stringsAsFactors = F)

  df$seqnames <- as.character(df$seqnames)
  df[ !grepl("KI|GL",df$seqnames),"seqnames"] <- paste0("chr",df[ !grepl("KI|GL",df$seqnames),"seqnames"])
  write.table(df, file="consensusPeaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


  sampleSheet2 <- dba.count(sampleSheet, peaks=consensusPeaks,minOverlap=1,summits=FALSE)
  saveRDS(sampleSheet2,file="dbaCountConsensus.RDS")



  
  } else{

    #get peaks for diffBind later
  consensusPeaks <- dba.peakset(sampleSheet, bRetrieve=TRUE)
 
  
  df <- data.frame(seqnames=seqnames(consensusPeaks),
   starts=start(consensusPeaks)-1,
   ends=end(consensusPeaks),
   names=c(rep(".", length(consensusPeaks))),
   scores=c(rep(".", length(consensusPeaks))),
   strands=strand(consensusPeaks),stringsAsFactors = F)

  df$seqnames <- as.character(df$seqnames)
  df[ !grepl("KI|GL",df$seqnames),"seqnames"] <- paste0("chr",df[ !grepl("KI|GL",df$seqnames),"seqnames"])
  write.table(df, file="consensusPeaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


    sampleSheet2 <- dba.count(sampleSheet, minOverlap = 2)
    peaks <- dba.peakset(sampleSheet2, bRetrieve=TRUE)
    saveRDS(peaks,file="consensusPeaks.RDS")
    saveRDS(sampleSheet2,file="dbaCount.RDS")
  }




  png("affinityCorrelation.png")
  correlationPlot2 <- plot(sampleSheet2)
  dev.off()

  sampleSheet2 <- dba.contrast(sampleSheet2, design=designFormula,contrast=c(opt$comparison,opt$numerator,opt$denominator),minMembers=2)
  sampleSheet2 <- dba.analyze(sampleSheet2, method = DBA_DESEQ2)
  plot(sampleSheet2, contrast = 1)

  saveRDS(sampleSheet2,file="dbaAnalyse.RDS")

  sigPeaks <- as.data.frame(dba.report(sampleSheet2))
  sigPeaks$seqnames = paste("chr",
    sigPeaks$seqnames,
    sep = "")
  if(opt$useConsensus) {
    write.table(sigPeaks, file = paste0("consensus_",opt$diffPeakTable), sep = "\t", row.names = FALSE, quote = FALSE)
    } else {
      write.table(sigPeaks, file = opt$diffPeakTable, sep = "\t", row.names = FALSE, quote = FALSE)
    }

    png("pca.png")
    dba.plotPCA(sampleSheet2, DBA_TREATMENT, label = DBA_ID, contrast = 1,
      labelSize = 0.6)
    dev.off()

    png("MA.png")
    dba.plotMA(sampleSheet2, bNormalized = TRUE, 
     th = 0.01, fold = 1)
    dev.off()

    png("volcano.png")
    dba.plotVolcano(sampleSheet2, dotSize = 0.8, 
      th = 0.01,
      fold = 1)
    dev.off()

    png("boxplot.png")
    dba.plotBox(sampleSheet2)
    dev.off()

    png("heatmap.png")
    dba.plotHeatmap(sampleSheet2, contrast=1, correlations=FALSE, scale = "row")
    dev.off()
