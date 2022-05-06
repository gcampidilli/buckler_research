load('ref_genes.rda')

subset_ins_mat_by_inbred_for = function(expression_data){
  kremling_inbred_spellings = expression_data$Taxa
  kremling_inbred_spellings = sub('282set_', '', kremling_inbred_spellings)

  included_inbreds = kremling_inbred_spellings[kremling_inbred_spellings %in% correct_inbred_spellings]
  ins_mat_updated = ins_mat_v4[ins_mat_v4$X %in% included_inbreds,]
  colnames(ins_mat_updated)[1] = 'Taxa'
  return(ins_mat_updated)
}


# CONSTRUCT INSERTION MATRIX GRANGES OBJECT
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

# CONVERT INSERTION MATRIX COLUMN NAMES FROM LOCI TO GENE NAMES, SOME COLUMNS WILL BE LOST

library(data.table)
library(GenomicRanges)
library(ggplot2)
library(IRanges)
library(dplyr)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)

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


# COMBINE COLUMNS THAT ARE THE SAME GENE

combine_gene_cols = function(indv_gene_ins_mat){
  unique_genes = unique(colnames(indv_gene_ins_mat))
  combined_gene_ins_mat = data.frame(matrix(0, nrow = nrow(indv_gene_ins_mat), ncol = length(unique_genes)))

  for(x in 1:length(unique_genes)){
    current_gene = unique_genes[x]
    gene_cols = grep(current_gene, colnames(indv_gene_ins_mat))

    if(length(gene_cols) > 1){
      gene_cols_subset = indv_gene_ins_mat[,gene_cols]
      gene_ins_mat_combined = unlist(lapply(1:nrow(gene_cols_subset), function(x){
        max(gene_cols_subset[x,])
      }))
      combined_gene_ins_mat[,x] = gene_ins_mat_combined
      colnames(combined_gene_ins_mat)[x] = current_gene
    }else{
      combined_gene_ins_mat[,x] = indv_gene_ins_mat[,gene_cols]
      colnames(combined_gene_ins_mat)[x] = current_gene
    }
  }

  return(combined_gene_ins_mat)
}

# CONVERT EXPRESSION DATA TO LONG DF

library(tidyverse)
library(rstatix)
library(ggpubr)

exp_to_long_df = function(exp_dat, included_taxa){
  exp_dat$Taxa = sub('282set_', '',exp_dat$Taxa)
  if(length(!duplicated(exp_dat$Taxa)) > 0){
    exp_dat = exp_dat[!duplicated(exp_dat$Taxa),]
  }
  exp_dat = exp_dat[exp_dat$Taxa %in% included_taxa,]
  # use pivot longer to create long df with columns 'Taxa', 'gene', 'expression'
  dat_long = exp_dat %>%
    pivot_longer(-Taxa, names_to = 'gene')
  colnames(dat_long)[3] = 'expression'
  return(dat_long)
}


# CALCULATE T-TESTS AND RETURN PVALUE FOR EXPRESSION ~ INSERTIONS FOR EACH GENE

pval = function(gene_name, combined_gene_ins_mat, dat_long){
  if(length(grep(gene_name, dat_long$gene)) > 0){
    gene_exp_subset = dat_long[grep(gene_name, dat_long$gene),]
    insertions = combined_gene_ins_mat[,gene_name]
    tmp_df = cbind(expression = gene_exp_subset$expression, insertions)
    result =  t.test(expression~insertions, tmp_df)
    return(result$p.value)
  } else{
    return(NA)
  }
}

entire_df_pval = function(combined_gene_ins_mat, dat_long){
  ins_count = unlist(lapply(1:ncol(combined_gene_ins_mat), function(x){
    sum(combined_gene_ins_mat[,x])
  }))
  candidate_genes = colnames(combined_gene_ins_mat)[which(ins_count >= 4)]
  candidate_genes_pval = unlist(lapply(1:length(candidate_genes), function(x){
    pval(candidate_genes[x], combined_gene_ins_mat, dat_long)
  }))
  return_df = data.frame(genes = candidate_genes, pval = candidate_genes_pval)
  return(return_df)
}




# CONSTRUCTION OF REF_GENES OBJECT
#load B73 4.0 reference genome
# gff_db =import.gff3("/Users/gcampidilli/Documents/research/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz", format = "gff3")
# # only examine genes
# ref_genes = subset(gff_db, gff_db$type == 'gene')
# save(ref_genes, file = 'ref_genes.rda')

