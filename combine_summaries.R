#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# combine simple df of each inbred

library(plyr)
library(dplyr)                                                  
library(readr)  

files_all <- list.files(path = "/Users/gcampidilli/Documents/research/mu_tir/step_5/generegion_summary",     # Identify all csv files in folder
                        pattern = "*.csv", full.names = TRUE) 

data_all = lapply(1:length(files_all), function(x){
  read.csv(files_all[x], header = T, row.names = 1)})

total_df = data.frame("inbred_name_vector"=NA,	"region"=NA,	"bp_count"=NA,	"bp_proportion"=NA,	"insertion_count"=NA,
                      "insertion_proportion"=NA,	"insertions_per_1kb"=NA)

x = 1
for(x in 1:length(files_all)) {
  current_df = data_all[[x]]
  temp_df = data.frame(inbred_name_vector = data_all[[x]]$inbred_name_vector,	region= data_all[[x]]$region,	bp_count= data_all[[x]]$bp_count,
                       bp_proportion= data_all[[x]]$bp_proportion,insertion_count= data_all[[x]]$insertion_count,
                       insertion_proportion= data_all[[x]]$insertion_proportion,	insertions_per_1kb= data_all[[x]]$insertions_per_1kb)
  temp_df = temp_df[!duplicated(temp_df),]
  total_df = rbind(total_df, temp_df)
}

write.csv(total_df, file = args[1])

