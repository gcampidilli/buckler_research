# convert v5 insertion matricies to v4 so they are compatible with Merritt's data

wd = '/Users/gcampidilli/Documents/research/mu_tir/step_7'
setwd(wd)


library(data.table)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(stringr)
library(rtracklayer)
library(xml2)
library(GenomicFeatures)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")

# upload v5 insertion matrix
v5_ins_mat = read.csv('/Users/gcampidilli/Documents/research/mu_tir/step_4/jan21_insertion_matrix_hundred.csv', header = T, row.names = 1)
v5_colnames = colnames(v5_ins_mat)

# remove scaf entries - keeping them causes NAs
scaf_cols = grep('scaf', v5_colnames)

# create GRanges object
chr_list = unlist(lapply(1:length(v5_colnames), function(x){
  str_split(v5_colnames[x], "_")[[1]][1]
}))
chr_list = chr_list[-c(scaf_cols)]

bp_list = as.numeric(unlist(lapply(1:length(v5_colnames), function(x){
  str_split(v5_colnames[x], "_")[[1]][2]
})))
bp_list = bp_list[-c(scaf_cols)]

v5_colnames_gff = GRanges(seqnames = chr_list, ranges = IRanges(start = bp_list, end = bp_list), strand = '+')

# save as GFF
export(v5_colnames_gff, 'v5_colnames_hundred.gff', format = 'GFF')

# run GFF at https://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter?db=core
# upload v4
v4_colnames_df = readGFF('/Users/gcampidilli/Documents/research/mu_tir/step_7/v4_colnames_hundred.gff')

v4_colnames = paste(v4_colnames_df$group,v4_colnames_df$start,sep='_')
length(v4_colnames) # 1354 PROBLEM!!!! 
length(v5_colnames_gff) # 1359 - we lost 5 entries somewhere

# figure out which entries were lost
# compare counts by chr number
table(v4_colnames_df$group)
table(sort(chr_list)) 
# lost entries:
# chr1 - 1
# chr5 - 1
# chr6 - 1
# chr9 - 1
# chr10 - 1

# not sure how to go about this - ask Michelle


v4_ins_mat = v5_ins_mat
v4_ins_mat = v4_ins_mat[,-c(scaf_cols)]
# dim(ins_mat_v4)
colnames(v4_ins_mat) = v4_colnames

write.csv(v4_ins_mat, file = "insertion_matrix_hundred_v4.csv")



######## 
# do same thing on ins_mat_thousand
# remove all variables from environment
rm(list=ls()) 

# upload v5 insertion matrix
v5_ins_mat = read.csv('/Users/gcampidilli/Documents/research/mu_tir/step_4/jan21_insertion_matrix_thousand.csv', header = T, row.names = 1)
v5_colnames = colnames(v5_ins_mat)

# remove scaf entries - keeping them causes NAs
scaf_cols = grep('scaf', v5_colnames)

# create GRanges object
chr_list = unlist(lapply(1:length(v5_colnames), function(x){
  str_split(v5_colnames[x], "_")[[1]][1]
}))
chr_list = chr_list[-c(scaf_cols)]

bp_list = as.numeric(unlist(lapply(1:length(v5_colnames), function(x){
  str_split(v5_colnames[x], "_")[[1]][2]
})))
bp_list = bp_list[-c(scaf_cols)]

v5_colnames_gff = GRanges(seqnames = chr_list, ranges = IRanges(start = bp_list, end = bp_list), strand = '+')

# save as GFF
export(v5_colnames_gff, 'v5_colnames_thousand.gff', format = 'GFF')

# run GFF at https://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter?db=core
# upload v4
v4_colnames_df = readGFF('/Users/gcampidilli/Documents/research/mu_tir/step_7/v4_colnames_thousand.gff')

v4_colnames = paste(v4_colnames_df$group,v4_colnames_df$start,sep='_')

v4_ins_mat = v5_ins_mat
v4_ins_mat = v4_ins_mat[,-c(scaf_cols)]
# dim(ins_mat_v4)
colnames(v4_ins_mat) = v4_colnames

write.csv(v4_ins_mat, file = "insertion_matrix_thousand_v4.csv")



