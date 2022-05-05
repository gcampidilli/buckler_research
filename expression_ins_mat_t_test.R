

# SUBSET INSERTION MATRIX TO HAVE SAME INBREDS AS EXPRESSION DATASET

load('ins_mat_v4.rda')
load('kremling_expression/GShoot_kremling_formatted_v4_hapmapids.rda')
load('correct_inbred_spellings.rda')

subset_ins_mat_by_inbred_for = function(expression_data){
  kremling_inbred_spellings = expression_data$Taxa
  kremling_inbred_spellings = sub('282set_', '', kremling_inbred_spellings)

  included_inbreds = kremling_inbred_spellings[kremling_inbred_spellings %in% correct_inbred_spellings]
  ins_mat_updated = ins_mat_v4[ins_mat_v4$X %in% included_inbreds,]
  colnames(ins_mat_updated)[1] = 'Taxa'
  return(ins_mat_updated)
}

# ins_mat_for_shoot_exp = subset_ins_mat_by_inbred_for(GShoot_kremling_formatted_v4_hapmapids)
# save(ins_mat_for_shoot_exp, file='ins_mat_for_shoot_exp.rda')


# CONSTRUCT INSERTION MATRIX GRANGES OBJECT
load('ins_mat_for_shoot_exp.rda')
df_to_granges = function(ins_mat){
  # prepare vectors for granges object, location, chromosome, bp
  locs = colnames(ins_mat)[-c(1)]
  ins_mat_chrom = sub("_[^_]+$", "", locs)
  ins_mat_chrom = unlist(lapply(1:length(ins_mat_chrom), function(x){
    as.numeric(substr(ins_mat_chrom[x], 4, nchar(ins_mat_chrom[x])))
  }))
  ins_mat_bp = as.numeric(sub(".*_", "", locs))
  # make Granges object from data frame
  ins_mat_granges_df = data.frame(chr = ins_mat_chrom, start = ins_mat_bp, end = ins_mat_bp, strand = rep('+', length(ins_mat_chrom)))
  ins_mat_granges = makeGRangesFromDataFrame(ins_mat_granges_df)
  return(ins_mat_granges)
}

# ins_mat_for_shoot_exp_granges = df_to_granges(ins_mat_for_shoot_exp)
# save(ins_mat_for_shoot_exp_granges, file = 'ins_mat_for_shoot_exp_granges.rda')


# CONVERT INSERTION MATRIX COLUMN NAMES FROM LOCI TO GENE NAMES, SOME COLUMNS WILL BE LOST

library(data.table)
library(GenomicRanges)
library(ggplot2)
library(IRanges)
library(dplyr)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)

# load B73 4.0 reference genes
load('ref_genes.rda')
# load insertion matrix granges object
load('ins_mat_for_shoot_exp_granges.rda')

ins_mat_columns_to_genes = function(ins_mat, ins_mat_granges){
  # find overlaps between insertion matrix and reference genes
  overlaps = findOverlaps(ins_mat_granges, ref_genes, type = 'within')
  # create blank insertion matrix only using gene names
  ins_mat_genes_colnames = data.frame(matrix(0, nrow = nrow(ins_mat), ncol = length(overlaps)))
  # add insertion matrix rows to ins_mat_genes_colnames from ins_mat
  x = 1
  for(x in 1:length(overlaps)){
    ins_mat_genes_colnames[,x] = ins_mat[,queryHits(overlaps)[x]]
    colnames(ins_mat_genes_colnames)[x] = ref_genes[subjectHits(overlaps)[x]]$gene_id
  }
  rownames(ins_mat_genes_colnames) = ins_mat$Taxa
  return(ins_mat_genes_colnames)
}
indv_gene_ins_mat_for_shoot_exp = ins_mat_columns_to_genes(ins_mat_for_shoot_exp, ins_mat_for_shoot_exp_granges)
save(indv_gene_ins_mat_for_shoot_exp, file = 'indv_gene_ins_mat_for_shoot_exp.rda')

# COMBINE COLUMNS THAT ARE THE SAME GENE


unique_genes = unique(colnames(gene_ins_mat_for_kremling_shoot))

unique_gene_ins_mat_for_kremling_shoot = data.frame(matrix(0, nrow = nrow(gene_ins_mat_for_kremling_shoot), ncol = length(unique_genes)))

for(x in 1:length(unique_genes)){
  current_gene = unique_genes[x]
  gene_cols = grep(current_gene, colnames(gene_ins_mat_for_kremling_shoot))

  if(length(gene_cols) > 1){
    gene_cols_subset = gene_ins_mat_for_kremling_shoot[,gene_cols]
    gene_ins_mat_combined = unlist(lapply(1:nrow(gene_cols_subset), function(x){
      max(gene_cols_subset[x,])
    }))
    unique_gene_ins_mat_for_kremling_shoot[,x] = gene_ins_mat_combined
    colnames(unique_gene_ins_mat_for_kremling_shoot)[x] = current_gene
  }else{
    unique_gene_ins_mat_for_kremling_shoot[,x] = gene_ins_mat_for_kremling_shoot[,gene_cols]
    colnames(unique_gene_ins_mat_for_kremling_shoot)[x] = current_gene
  }
}


save(combined_gene_ins_mat_for_shoot_exp, file = 'combined_gene_ins_mat_for_shoot_exp.rda')












# CONVERT EXPRESSION DATA TO LONG DF

library(tidyverse)
library(rstatix)
library(ggpubr)
load('kremling_expression/GShoot_kremling_formatted_v4_hapmapids.rda')
load('correct_inbred_spellings.rda')

exp_to_long_df = function(exp_dat, included_taxa){
  exp_dat$Taxa = sub('282set_', '',exp_dat$Taxa)
  if(length(!duplicated(exp_dat$Taxa)) > 0){
    exp_dat = exp_dat[!duplicated(exp_dat$Taxa),]
  }
  exp_dat = exp_dat[exp_dat$Taxa %in% included_taxa,]
  # use pivot longer to create long df with columns 'Taxa', 'gene', 'expression'
  dat_long = dat %>%
    pivot_longer(-Taxa, names_to = 'gene')
  colnames(dat_long)[3] = 'expression'
  return(dat_long)
}

# exp_shoot_long_df = exp_to_long_df(GShoot_kremling_formatted_v4_hapmapids, ins_mat_for_shoot_exp$Taxa)
# save(exp_shoot_long_df, file = 'exp_shoot_long_df.rda')

# CALCULATE T-TESTS AND RETURN PVALUE FOR EXPRESSION ~ INSERTIONS FOR EACH GENE





















# CONSTRUCTION OF REF_GENES OBJECT
#load B73 4.0 reference genome
gff_db =import.gff3("/Users/gcampidilli/Documents/research/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz", format = "gff3")
# only examine genes
ref_genes = subset(gff_db, gff_db$type == 'gene')
save(ref_genes, file = 'ref_genes.rda')

