library("optparse")
library("tidyverse")
library("bedr")

#Parse command line parameters
option_list <- list()

option_list$diffPeaks <- make_option('--diffPeaks', type='character', help="The path to the bigwig files to merge per condition")


opt <- parse_args(OptionParser(option_list=option_list,description="get sequences of sig atac-seq peaks in fasta format for meme-chip"))



old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/opt/miniconda3/envs/diffbind/bin", sep = ":"))


## Load significant peaks ----

sig_peaks <- read.table(opt$diffPeaks,
                        sep = "\t", 
                        header = TRUE, 
                        quote = "", 
                        stringsAsFactors = FALSE)

sig_peaks$seqnames <- gsub("chr","",sig_peaks$seqnames)

# split into up and down reg peaks
up_peaks <- sig_peaks %>% filter(Fold > 2,FDR<=0.01) %>% 
            .[,c(1:3)] %>% 
            bedr.sort.region(check.chr = FALSE)
down_peaks <- sig_peaks %>% filter(Fold < -2,FDR<=0.01) %>% 
            .[,c(1:3)] %>% 
            bedr.sort.region(check.chr = FALSE)

# need indexed fasta file
up_fasta <- bedr::get.fasta(x = up_peaks,
                            fasta = opt$fasta,check.chr = FALSE) %>%
            mutate(index = paste(">", index, sep = ""))

down_fasta <- bedr::get.fasta(x = down_peaks,
                            fasta = opt$fasta,check.chr = FALSE) %>%
            mutate(index = paste(">", index, sep = ""))

# format fasta
up_fasta2 <- vector(mode = "list")
for(i in 1:length(rownames(up_fasta))){
    print(i)
    up_fasta2[[i]] <- rbind(up_fasta$index[i], up_fasta$sequence[i])
}
up_fasta2 <- do.call(rbind, up_fasta2)

down_fasta2 <- vector(mode = "list")
for(i in 1:length(rownames(down_fasta))){
    print(i)
    down_fasta2[[i]] <- rbind(down_fasta$index[i], down_fasta$sequence[i])
}
down_fasta2 <- do.call(rbind, down_fasta2)

# write
write.table(up_fasta2, file = "SigPeaksUP.fasta",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(down_fasta2, file = "SigPeaksDOWN.fasta",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)





