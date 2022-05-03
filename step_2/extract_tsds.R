#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# extract TSDs
# wd = "/Users/gcampidilli/Documents/research/mu_tir"
# setwd(wd)

library(GenomicRanges)
library(IRanges)
library(dplyr)
library(stringr)
library(rtracklayer)
library(xml2)
library(GenomicFeatures)
library(seqinr)
library(Rsamtools)

# replace w22 with inbred name 
inbred_name = args[1]

blast_results = read.delim(paste("~/goodman_blast_out/", inbred_name, "_mu_blast_out.txt", sep = ""), header = T)
colnames(blast_results) = c("qname", "sname", "percent_identity", "alignment_length", "mismatches", 
                                "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

bam_df = read.csv(paste("~/goodman_bam_df/",inbred_name,"_bam_df.csv",sep = ""), header = T )
bam_df = bam_df[,-c(1,14)]
# now examine all hits that has s_start = 1, this will ensure that we are only looking at + strand hits for mutirs
blast_results_start1 = subset(blast_results, blast_results$s_start == 1)

# make sure there are no repeats
duplicates = blast_results_start1[duplicated(blast_results_start1$qname),]$qname

if(length(duplicates) > 0){
  # print duplicates
  blast_results_start1[blast_results_start1$qname==duplicates,]
  
  # remove duplicates from blast_results_start1
  x = 1
  for(x in 1:length(duplicates)){
    df = blast_results_start1[blast_results_start1$qname==duplicates[x],]
    num_duplicates = nrow(df)
    df_sorted = df[order(-df$percent_identity),]
    rows_remove = row.names(df_sorted[(2:num_duplicates),])
    blast_results_start1 = subset(blast_results_start1, !(row.names(blast_results_start1) %in% rows_remove))
  }
}

extract_tsd = function(x){
  current_hit = blast_results_start1[x,]
  hit_name = current_hit$qname
  # extract reference bam query
  hit_bam_ref = bam_df[grepl(paste(hit_name,"$",sep = ""), bam_df$fastaname),]
  if(hit_bam_ref$strand == "-") {
    qend = current_hit$q_end
    # extract bam query sequence
    seq = as.character(hit_bam_ref$seq)
    # extract TSD
    substr(seq, (qend+1), (qend+9))
  } else{
    hit_bam_ref = subset(hit_bam_ref, hit_bam_ref$strand == "+")
    qstart = current_hit$q_start
    # extract bam query sequence
    seq = as.character(hit_bam_ref$seq)
    # extract TSD
    substr(seq, (qstart-9), (qstart-1))
  }
}

tsds_total = unlist(lapply(1:nrow(blast_results_start1), function(x)extract_tsd(x)))

# create fastaname vector
strands_total = subset(bam_df,bam_df$fastaname %in% blast_results_start1$qname)

insertion_df = data.frame(blast_results_start1, tsds= tsds_total, strand = strands_total$strand)

write.csv(insertion_df, file = args[2])

