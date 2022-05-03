#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# replace w22 with inbred name 
inbred_name = args[1]

# wd = "/Users/gcampidilli/Documents/research/mu_tir/step_4"
# setwd(wd)
# inbred_name = "B57"

library(plyr) 
library(dplyr)                                                  
library(readr) 

# filter flanking seq blast outputs 
inbred_flanking = read.delim(paste("blast_out_flanking/", inbred_name, "_blast_out_flanking.txt", sep = ""), header = F)
colnames(inbred_flanking) = c("qname", "sname", "percent_identity", "alignment_length", "mismatches", 
                              "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

duplicates = unique(inbred_flanking[duplicated(inbred_flanking$qname),]$qname)

remove_duplicates = function(x) {
  df = inbred_flanking[inbred_flanking$qname==duplicates[x],]
  num_duplicates = nrow(df)
  df$duplicate_num = 1:num_duplicates
  # we care both about percent identity
  df_sorted_percent = df[order(-df$percent_identity),]
  return(row.names(df_sorted_percent[(2:num_duplicates),]))
}

if(length(duplicates) > 0){
  # remove duplicates from inbred_flanking
  rows.remove = unlist(lapply(1:length(duplicates), remove_duplicates))
}

inbred_flanking_filtered = inbred_flanking[!(row.names(inbred_flanking) %in% rows.remove),]

print("filtering phase 1 done")
head(inbred_flanking_filtered)

# work backwards, otherwise this code would be go over the same insertion multiple times
end = length(inbred_flanking_filtered$s_start)
while(end > 1){
  duplicate_index_vector = vector()
  duplicate_index_vector = append(duplicate_index_vector,end)
  x = 1
  for(x in 1:(end-1)){
    if(x != 0 && end != 0) {
      # tests to see whether insertions are within 18 bp of each other
      # if they are, then the index are added to vector
      if(abs(inbred_flanking_filtered$s_start[end]-inbred_flanking_filtered$s_start[x])<19){
        duplicate_index_vector = append(duplicate_index_vector, x)
      }
    }
  }
  
  # use the duplicate vector to identify which duplicate(s) to remove
  # first look at percent_identity, next alignment length
  if(length(duplicate_index_vector)>1){
    remove_vector = vector()
    if(length(which(inbred_flanking_filtered$percent_identity[duplicate_index_vector]==max(inbred_flanking_filtered$percent_identity[duplicate_index_vector]))) == 1) { 
      remove_vector = which(inbred_flanking_filtered$percent_identity[duplicate_index_vector]!=max(inbred_flanking_filtered$percent_identity[duplicate_index_vector]))
    } else{
      if(length(which(inbred_flanking_filtered$alignment_length[duplicate_index_vector]==max(inbred_flanking_filtered$alignment_length[duplicate_index_vector]))) == 1){
        remove_vector = which(inbred_flanking_filtered$alignment_length[duplicate_index_vector]!=max(inbred_flanking_filtered$alignment_length[duplicate_index_vector]))} 
      else {
        remove_vector = duplicate_index_vector[-c(1)]
      }
    }
    inbred_flanking_filtered = inbred_flanking_filtered[-remove_vector,]
    end = end - (length(remove_vector)+1)
  }else{
    end = end - 1
  } 
} 
head("filtering phase 2 done")
head(inbred_flanking_filtered)

library(plyr)
inbred_flanking_filtered$rounded_hundred = round_any(inbred_flanking_filtered$s_start, 100)
inbred_flanking_filtered$rounded_thousand = round_any(inbred_flanking_filtered$s_start, 1000)

inbred_flanking_filtered$id_exact = paste(inbred_flanking_filtered$sname, "_",inbred_flanking_filtered$s_start, sep = "")
inbred_flanking_filtered$id_hundred = paste(inbred_flanking_filtered$sname, "_",inbred_flanking_filtered$rounded_hundred, sep = "")
inbred_flanking_filtered$id_thousand = paste(inbred_flanking_filtered$sname, "_",inbred_flanking_filtered$rounded_thousand, sep = "")

write.csv(inbred_flanking_filtered, file = args[2])

simple_prematrix_df_exact = data.frame(inbred = inbred_name, insertion = unique(inbred_flanking_filtered$id_exact))
write.csv(simple_prematrix_df_exact, file = args[3])

simple_prematrix_df_hundred = data.frame(inbred = inbred_name, insertion = unique(inbred_flanking_filtered$id_hundred))
write.csv(simple_prematrix_df_hundred, file = args[4])

simple_prematrix_df_thousand = data.frame(inbred = inbred_name, insertion = unique(inbred_flanking_filtered$id_thousand))
write.csv(simple_prematrix_df_thousand, file = args[5])


##### COMMENT OUT WHEN USING SCRIPT ON TERMINAL ###
# write.csv(inbred_flanking_filtered, file = "blast_out_flanking_filtered/38-11Goodman-Buckler_flanking_filtered.csv")
# 
# simple_prematrix_df_exact = data.frame(inbred = inbred_name, insertion = unique(inbred_flanking_filtered$id_exact))
# write.csv(simple_prematrix_df_exact, file = "prematrix_df/exact/38-11Goodman-Buckler_prematrix_df_exact.csv")
# 
# simple_prematrix_df_hundred = data.frame(inbred = inbred_name, insertion = unique(inbred_flanking_filtered$id_hundred))
# write.csv(simple_prematrix_df_hundred, file = "prematrix_df/hundred/38-11Goodman-Buckler_prematrix_df_hundred.csv")
# 
# simple_prematrix_df_thousand = data.frame(inbred = inbred_name, insertion = unique(inbred_flanking_filtered$id_thousand))
# write.csv(simple_prematrix_df_thousand, file ="prematrix_df/thousand/38-11Goodman-Buckler_prematrix_df_thousand.csv")

# Next: combine simple prematrix_df for all inbreds, then to get total count of TEs for each inbred:

