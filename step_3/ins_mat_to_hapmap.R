
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

ins_mat = read.csv('insertion_matrix_exact_v4.csv', row.names = 1, header = T)

# rs#    		alleles    	chrom    pos    strand    assembly#    center    protLSID    assayLSID    panelLSID    QCcode  Sample columns: ex. Parental_1  Parental_2
# Instead of Parental, we will add column for each row in ins mat i.e.  33-16  38-11   etc

# assume homozygous bc inbreds
# 'AA' for TE absent
# 'TT' for TE present

rs_num = colnames(ins_mat)

alleles = rep('A/T', length(rs_num))
  
chrom = unlist(lapply(1:length(rs_num), function(x){
  as.numeric(sub('chr','',str_split(rs_num[x], "_")[[1]][1]))
}))

pos = as.numeric(unlist(lapply(1:length(rs_num), function(x){
  as.numeric(str_split(rs_num[x], "_")[[1]][2])
})))

hm_partial = data.frame(rs = rs_num, alleles = "A/T", chrom = chrom, pos = pos, strand = '+',
                      assembly = NA, center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA)

samples = t(ins_mat)
samples[samples == 0] = 'AA'
samples[samples == 1] = 'TT'

hm_total = cbind(hm_partial, samples)

write.table(hm_total, "step_7_ins_mat_test.hmp.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


# next steps

# merritt pheno - GWAS run on all SNPs for each inbred in goodman, associate with traits
# current counts from Merritt's genotype can act as benchmark for what SNPs/genotypes can explain phenotypic variances


# if there are TE associations that have zero SNP associations, those would be the ones we are interested in 
# which would indicate TE genotyping is telling us something that SNP genotyping isn't

# bc mutator insert into promoter, next step would be
# generate null model that allows random sampling of 1kb window - ask whether there exists greater or fewer flowering time associations with genes that mu is inserted into

# plot # association overlap along chromosomes - we might expect that TEs may tag some of the known QTL associated with traits
