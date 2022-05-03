#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# combine simple df of each inbred

library(plyr)
library(dplyr)                                                  
library(readr)  

files_all <- list.files(path = "/Users/gcampidilli/Documents/research/mu_tir/step_4/blast_out_flanking_filtered",     # Identify all csv files in folder
                       pattern = "*.csv", full.names = TRUE) 
  
data_all = lapply(1:length(files_all), function(x){
  read.csv(files_all[x], header = T, row.names = 1)})

inbred_list = read.delim("/Users/gcampidilli/Documents/research/mu_tir/step_4/inbred_list.txt")

#inbred1 = read.csv(files_all[[1]], header = T, row.names = 1)
total_df = data.frame("inbred"=NA,"chrom"=NA,"exactpos"=NA, "id"=NA)
x = 1
for(x in 1:length(files_all)) {
  current_df = data_all[[x]]
  temp_df = data.frame(chrom = data_all[[x]]$sname, exactpos = data_all[[x]]$s_start)
  temp_df = temp_df[!duplicated(temp_df),]
  temp_df$chrom = strtoi(substring(temp_df$chrom,first = 4))
  temp_df$id = paste(temp_df$chrom, "_",temp_df$exactpos, sep = "")
  n_reps = nrow(temp_df)
  return_df = data.frame(inbred = rep(inbred_list[x,1], n_reps), chrom = temp_df$chrom, exactpos = temp_df$exactpos, id = temp_df$id)
  total_df = rbind(total_df, return_df)
}

write.csv(total_df, file = "v2_total_prematrix_df.csv")




