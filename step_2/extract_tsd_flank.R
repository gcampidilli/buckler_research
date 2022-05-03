#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

library(dplyr)
library(seqinr)

inbred_name = args[1]
print(paste("arg_2:", args[2],sep=""))

print("loading data")
bam_df = read.csv(paste("~/goodman_bam_df/",inbred_name,"_bam_df.csv", sep = ""), header = T)
print(paste("bam_df dimensions:", dim(bam_df),sep = ""))

inbred_insertions_9TSD = read.csv(paste("~/nov19_goodman_insertion_df/",inbred_name,"_insertion_df_9TSD.csv", sep = ""), header = T)
print(paste("inbred_insertion_9TSD dimensions:", dim(inbred_insertions_9TSD), sep = ""))

print("making subset")
bam_subset = subset(bam_df, bam_df$fastaname %in% inbred_insertions_9TSD$qname)
print(paste("bam subset dimensions:", dim(bam_subset), sep = ""))


print("taking substring of flanking + TSD")
inbred_seq_array = array()
inbred_seq_array = unlist(lapply(1:length(bam_subset$strand), function(x){
  if(bam_subset$strand[x] == "+"){
    mu_start = inbred_insertions_9TSD$q_start[x]
    append(inbred_seq_array, substring(bam_subset$seq[x], 1, mu_start))
  } else {
    mu_start = inbred_insertions_9TSD$q_end[x]
    append(inbred_seq_array, substring(bam_subset$seq[x], mu_start, 150))}
}))

inbred_seq_array = inbred_seq_array[!is.na(inbred_seq_array)]

inbred_subset_fastaname = paste(bam_subset$fastaname, "_flanking", sep = "")

print("writing to fasta")

write.fasta(sequences = inbred_seq_array[1], names = inbred_subset_fastaname[1], file.out = args[2])
lapply(2:length(inbred_seq_array), function(x)write.fasta(sequences = inbred_seq_array[x],names = inbred_subset_fastaname[x], open = "a", file.out =args[2]))
