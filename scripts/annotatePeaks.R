library("optparse")
library("tidyverse")


old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/opt/miniconda3/envs/diffbind/bin", sep = ":"))

#Parse command line parameters
option_list <- list()

option_list$diffPeaks <- make_option('--diffPeaks', type='character' ,help="Annotate conserved peaks bedfile from diffbind")
option_list$annotated <- make_option('--annotated', type='character' ,help="Name of the output annotated bedfile")


opt <- parse_args(OptionParser(option_list=option_list,description="Annotate peaks"))

sigPeaks <- read.delim(opt$diffPeaks)

# format for homer annotatePeaks.pl
annotateBED <- sigPeaks %>%
  mutate(PeakID = paste("peak_", 1:length(rownames(.)), sep = ""),
         strand = "+") %>%
  .[,c(1:3, 12, 4:5)]

sigPeaks <- sigPeaks %>%
  mutate(PeakID = paste("peak_", 1:length(rownames(.)), sep = ""))

write.table(annotateBED, file = "tmp.bed",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


annotatePeaks <- paste("/opt/miniconda3/envs/diffbind/bin/annotatePeaks.pl", "tmp.bed",
                       "hg38", ">", opt$annotated,
                       sep = " ")

system(annotatePeaks)

annotatedPeaks <- read.table(opt$annotated,
                             sep = "\t",
                             quote = "",
                             header = TRUE,
                             fill = TRUE)

colnames(annotatedPeaks)[1] <- "PeakID"

annotatedPeaks <- annotatedPeaks %>%
  .[,c(1, 8:19)] %>%
  left_join(sigPeaks, .,
            by = c("PeakID"))

write.table(annotatedPeaks, file = opt$annotated,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
