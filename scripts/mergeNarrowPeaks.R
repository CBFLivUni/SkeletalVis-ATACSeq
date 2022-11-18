library("optparse")
library("tidyverse")

#Parse command line parameters
option_list <- list()

option_list$narrowPeaks <- make_option('--narrowPeaks', type='character', help="The path to the narrowPeak files to merge per condition")
option_list$sampleTable <- make_option('--sampleTable', type='character', help="sampleTable mapping sample to condition")
option_list$tmpDir <- make_option('--tmpDir', type='character', default = "/tmp/", help="The path to the tmp directory to be used to sort bedfiles")


opt <- parse_args(OptionParser(option_list=option_list,description="Merge narrowPeaks files using by Treatment"))


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

narrowPeaks <- getTable(opt$narrowPeaks,"narrowPeaks")

sampleTable <- read.delim(opt$sampleTable)
colnames(sampleTable)[1] <- "SampleID"



sampleTable <- merge(sampleTable,narrowPeaks,by="SampleID")


#Merge all the bigwig files per treatment
mergeNarrowPeaks <- function(condition,sampleTable) {
  
  narrowPeaks <- paste(sampleTable[sampleTable$combined==condition,"narrowPeaks"],collapse = " ")
  
  #cat all the files togther
  cmd <- sprintf("cat %s > %s_combined.narrowPeaks",narrowPeaks,condition)
  print(cmd)
  system(cmd)
  
  #sort by coordinate
  cmd <- sprintf("sort -k1,1 -k2,2n -k3,3n %s_combined.narrowPeaks > %s_sorted.narrowPeaks",condition,condition)
  print(cmd)
  system(cmd)
  
  #extract the relevant info
  cmd <- sprintf("cat %s_sorted.narrowPeaks | awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$7}' >  %s_sorted.bed",condition,condition)
  print(cmd)
  system(cmd)
  
  #merge peaks within 10bp
  cmd <- sprintf("bedtools merge -d 10 -c 4 -o mean -i %s_sorted.bed > %s_mergedpeaks.bed",condition,condition)
  print(cmd)
  system(cmd)
  
  #add a column  with the number of replicates that each peak was found in.
  cmd <- sprintf("bedtools intersect -wa -wb -filenames -a %s_mergedpeaks.bed -b %s -sorted -F 1.0 > %s_replicates.bed",condition,narrowPeaks,condition)
  print(cmd)
  system(cmd)
  
  bed <- read.delim(sprintf("%s_replicates.bed",condition),header=F)
  bed$peakID <- paste(bed[,1],bed[,2],bed[,3],sep="_")
  bed <- bed %>% group_by(peakID) %>% mutate(num_replicates=n_distinct(V5)) %>%
    slice_head(n=1) %>% ungroup() %>% select(V1,V2,V3,V4,num_replicates) %>% arrange(V1,V2)
  
  write.table(bed,file = sprintf("%s_replicates.bed",condition),col.names=FALSE,row.names = FALSE,quote = FALSE,sep="\t")
  
  percentageReplicated <- sum(table(bed$num_replicates)[-1])/sum(table(bed$num_replicates))*100
  
  return(table(bed$num_replicates))

}


#make a combined factor if multiple
factors <- sampleTable[,!grepl("SampleID|File|narrowPeaks",colnames(sampleTable)),drop=FALSE]
sampleTable$combined <- apply(factors,1,function(x) paste(x,collapse = "_"))


#loop through the treatments and merge together all the narrowpeak files for that treatment
sapply(unique(sampleTable$combined),mergeNarrowPeaks,sampleTable)
