#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

# combine bam files into fasta for any given inbred
# make sure wd is where .bam files and .bai files are located
# 
# wd = "/Users/gcampidilli/Documents/research/mu_tir/all_bams"
# setwd(wd)

inbred_name = args[1]
bam_list = args[3:length(args)]

print(bam_list)
num_bams = length(bam_list)

CPU=24 ## adjust as needed
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
# inbred_name = "W22"
# load all bams
# bam_list = c("W22_HGCJCBGXX_GCCAAT_srt_dedup.sortedfiltered.MuTIR.bam",
#              "W22_HVMFTCCXX_L8.clean_srt_dedup.sortedfiltered.MuTIR.bam",
#              "W22_HVMFTCCXX_L7.clean_srt_dedup.sortedfiltered.MuTIR.bam")
combined_seq = vector()
combined_names = vector()
combined_df = data.frame("qname"=NA,"flag"=NA,"rname"=NA,"strand"=NA,"pos"=NA,"qwidth"=NA,
                         "mapq"=NA,"cigar"=NA,"mrnm"=NA, "mpos" = NA, "isize" = NA, "seq" = NA,
                         "qual"=NA, "fastaname"=NA)
x = 1
for(x in 1:num_bams) {
  bamname = bam_list[x]
 print(bamname)
  samplename = gsub("_srt_dedup.sortedfiltered.MuTIR.bam", "",bamname)
  samplename = gsub(".clean", "", samplename)
  samplename = gsub("./", "", samplename)
  bam = scanBam(bamname, index = paste(bamname,".bai", sep =""))
  # transform bam into a usable data frame (expand following section of code)
  #names of the BAM fields
  names(bam[[1]])
  # [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"
  # [9] "mrnm"   "mpos"   "isize"  "seq"    "qual"
  #distribution of BAM flags
  table(bam[[1]]$flag)
  #function for collapsing the list of lists into a single list
  #as per the Rsamtools vignette
  .unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  #store names of BAM fields
  bam_field <- names(bam[[1]])
  #go through each BAM field and unlist
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  #store as data frame
  bam_df <- data.frame(list)
  names(bam_df) <- bam_field
  bam_df$qual = NA
  
  # remember for our blast results, w22_bam is our QUERY, and mutir groups are our SUBJECT
  # we want our SUBJECT to start at 1
  bam_df$fastaname = paste0(samplename, '_read', 1:nrow(bam_df))
  
  combined_seq = c(combined_seq, as.character(bam_df$seq))
  combined_names = c(combined_names, bam_df$fastaname)
 
  combined_df = rbind(combined_df, bam_df)
 
  print(samplename)
}

combined_seq_DNAString = lapply(1:length(combined_seq), function(x)DNAString(combined_seq[x]))

#write.fasta(sequences = combined_seq_DNAString,names = combined_names,file.out = args[2])

write.csv(combined_df[-c(1),], file = args[2])


