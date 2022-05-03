
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
inbred_name = args[1]

insertion_df = read.csv(paste("~/goodman_insertion_df/",args[1],"_insertion_df.csv", sep = ""), header = T)
filtered_dat = insertion_df[which(nchar(insertion_df$tsds) == 9),]
write.csv(filtered_dat, file = args[2])

