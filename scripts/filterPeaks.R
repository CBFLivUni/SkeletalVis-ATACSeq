library("optparse")
library("tidyverse")

#Parse command line parameters
option_list <- list()

option_list$bed <- make_option('--bed', type='character', help="The path to the peaks bed file to use for homer")
option_list$topN <- make_option('--topN', type='integer', help="Number of top peaks to select")


opt <- parse_args(OptionParser(option_list=option_list,description="Get the top up and down peaks for homer input"))


sig_peaks <- read.table(opt$bed,
                        sep = "\t", 
                        header = TRUE, 
                        quote = "", 
                        stringsAsFactors = FALSE)


# split into up and down reg peaks
up_peaks <- sig_peaks %>% filter(Fold > 0,FDR<=0.01) %>% slice_max(Fold,n=5000)
down_peaks <- sig_peaks %>% filter(Fold < 0,FDR<=0.01) %>% slice_min(Fold,n=5000)

# write output
write.table(up_peaks, file = "SigPeaksUP.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(down_peaks, file = "SigPeaksDOWN.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)





