library("optparse")

#Parse command line parameters
option_list <- list()

option_list$bigWigs <- make_option('--bigWigs', type='character', help="The path to the bigwig files to merge per condition")
option_list$chrSizes <- make_option('--chrSizes', type='character', help="chromosome sizes")
option_list$sampleTable <- make_option('--sampleTable', type='character', help="sampleTable mapping sample to condition")
option_list$tmpDir <- make_option('--tmpDir', type='character', default = "/tmp/", help="The path to the tmp directory to be used to sort bedfiles")


opt <- parse_args(OptionParser(option_list=option_list,description="Merge bigwig files using by Treatment"))



old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/opt/miniconda3/envs/diffbind/bin", sep = ":"))
Sys.setenv(TMPDIR=opt$tmpDir)


getTable <- function(fileList,fileType){
  fileList <- unlist(strsplit(fileList, ' '))
  names(fileList) <- tools::file_path_sans_ext(basename(fileList))
  fileTable <- stack(fileList)
  colnames(fileTable) <- c(fileType,"SampleID")
  fileTable$SampleID <- gsub("_peaks|_noMT_noDups_unique","",  fileTable$SampleID)
  return(fileTable)
}

bigwigs <- getTable(opt$bigWigs,"bigWigs")

sampleTable <- read.delim(opt$sampleTable)
colnames(sampleTable)[1] <- "SampleID"



sampleTable <- merge(sampleTable,bigwigs,by="SampleID")


#Merge all the bigwig files per treatment
mergeBigWigs <- function(condition,sampleTable,chrSizes) {
  bigWigs <- paste(sampleTable[sampleTable$combined==condition,"bigWigs"],collapse = " ")
  cmd <- sprintf("bigWigMerge %s ATAC_%s.bedGraph",bigWigs,condition)
  print(cmd)
  system(cmd)
  cmd <- sprintf("LC_COLLATE=C sort -k1,1 -k2,2n ATAC_%s.bedGraph > ATAC_%s.bedGraph.sorted",condition,condition)
  print(cmd)
  system(cmd)
  cmd <- sprintf("bedGraphToBigWig ATAC_%s.bedGraph.sorted %s ATAC_%s.bw",condition,chrSizes,condition)
  print(cmd)
  system(cmd)
}


#make a combined factor if multiple
factors <- sampleTable[,!grepl("SampleID|File|bigWigs",colnames(sampleTable)),drop=FALSE]
sampleTable$combined <- apply(factors,1,function(x) paste(x,collapse = "_"))


#loop through the treatments and merge together all the bigwig files for that treatment
sapply(unique(sampleTable$combined),mergeBigWigs,sampleTable,opt$chrSizes)
