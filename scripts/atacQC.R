library("optparse")

#Parse command line parameters
option_list <- list()

option_list$bamFile <- make_option('--bamFile', type='character' ,help="bam file to calculate fragment length read density")
option_list$fragmentPlot <- make_option('--fragmentPlot', type='character' ,help="Plot of fragment length")
option_list$fragmentTable <- make_option('--fragmentTable', type='character' ,help="Table of plot data")

opt <- parse_args(OptionParser(option_list=option_list,description="Calculate fragment lengths for atac-seq data"))

library("ATACseqQC")
library("cowplot")

png(opt$fragmentPlot)
fragSize <- fragSizeDist(opt$bamFile,basename(opt$bamFile))
dev.off()

fragSizeTable <- stack(fragSize[[1]])[,2:1]
colnames(fragSizeTable) <- c("ReadDensity","FragLength")
write.table(fragSizeTable, file=opt$fragmentTable, sep = "\t",row.names=F,quote = F,col.names = T)

