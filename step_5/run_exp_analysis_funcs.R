# SOURCE SCRIPT WITH FUNCTIONS
source('expression_analysis_funcs.R')

# LOAD NECESSARY DATA OBJECTS
load('ins_mat_v4.rda')
load('correct_inbred_spellings.rda')

# load all expression data
load('kremling_expression/GShoot_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/GRoot_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/Kern_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/L3_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/L3Base_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/L3Tip_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/LMAD_kremling_formatted_v4_hapmapids.rda')
load('kremling_expression/LMAN_kremling_formatted_v4_hapmapids.rda')

# SUBSET INSERTION MATRIX TO HAVE SAME INBREDS AS EXPRESSION DATASET

ins_mat_for_shoot_exp = subset_ins_mat_by_inbred_for(GShoot_kremling_formatted_v4_hapmapids)
save(ins_mat_for_shoot_exp, file='exp_ins_mat/orig_subset/ins_mat_for_shoot_exp.rda')

ins_mat_for_Root_exp = subset_ins_mat_by_inbred_for(GRoot_kremling_formatted_v4_hapmapids)
save(ins_mat_for_Root_exp, file='exp_ins_mat/orig_subset/ins_mat_for_Root_exp.rda')

ins_mat_for_Kern_exp = subset_ins_mat_by_inbred_for(Kern_kremling_formatted_v4_hapmapids)
save(ins_mat_for_Kern_exp, file='exp_ins_mat/orig_subset/ins_mat_for_Kern_exp.rda')

ins_mat_for_L3_exp = subset_ins_mat_by_inbred_for(L3_kremling_formatted_v4_hapmapids)
save(ins_mat_for_L3_exp, file='exp_ins_mat/orig_subset/ins_mat_for_L3_exp.rda')

ins_mat_for_L3Base_exp = subset_ins_mat_by_inbred_for(L3Base_kremling_formatted_v4_hapmapids)
save(ins_mat_for_L3Base_exp, file='exp_ins_mat/orig_subset/ins_mat_for_L3Base_exp.rda')

ins_mat_for_L3Tip_exp = subset_ins_mat_by_inbred_for(L3Tip_kremling_formatted_v4_hapmapids)
save(ins_mat_for_L3Tip_exp, file='exp_ins_mat/orig_subset/ins_mat_for_L3Tip_exp.rda')

ins_mat_for_LMAD_exp = subset_ins_mat_by_inbred_for(LMAD_kremling_formatted_v4_hapmapids)
save(ins_mat_for_LMAD_exp, file='exp_ins_mat/orig_subset/ins_mat_for_LMAD_exp.rda')

ins_mat_for_LMAN_exp = subset_ins_mat_by_inbred_for(LMAN_kremling_formatted_v4_hapmapids)
save(ins_mat_for_LMAN_exp, file='exp_ins_mat/orig_subset/ins_mat_for_LMAN_exp.rda')

# CONSTRUCT INSERTION MATRIX GRANGES OBJECT

ins_mat_for_shoot_exp_granges = df_to_granges(ins_mat_for_shoot_exp)
save(ins_mat_for_shoot_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_shoot_exp_granges.rda')

ins_mat_for_Root_exp_granges = df_to_granges(ins_mat_for_Root_exp)
save(ins_mat_for_Root_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_Root_exp_granges.rda')

ins_mat_for_Kern_exp_granges = df_to_granges(ins_mat_for_Kern_exp)
save(ins_mat_for_Kern_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_Kern_exp_granges.rda')

ins_mat_for_L3_exp_granges = df_to_granges(ins_mat_for_L3_exp)
save(ins_mat_for_L3_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_L3_exp_granges.rda')

ins_mat_for_L3Base_exp_granges = df_to_granges(ins_mat_for_L3Base_exp)
save(ins_mat_for_L3Base_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_L3Base_exp_granges.rda')

ins_mat_for_L3Tip_exp_granges = df_to_granges(ins_mat_for_L3Tip_exp)
save(ins_mat_for_L3Tip_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_L3Tip_exp_granges.rda')

ins_mat_for_LMAD_exp_granges = df_to_granges(ins_mat_for_LMAD_exp)
save(ins_mat_for_LMAD_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_LMAD_exp_granges.rda')

ins_mat_for_LMAN_exp_granges = df_to_granges(ins_mat_for_LMAN_exp)
save(ins_mat_for_LMAN_exp_granges, file = 'exp_ins_mat/granges/ins_mat_for_LMAN_exp_granges.rda')


# CONVERT INSERTION MATRIX COLUMN NAMES FROM LOCI TO GENE NAMES, SOME COLUMNS WILL BE LOST

indv_gene_ins_mat_for_shoot_exp = ins_mat_columns_to_genes(ins_mat_for_shoot_exp, ins_mat_for_shoot_exp_granges)
save(indv_gene_ins_mat_for_shoot_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_shoot_exp.rda')

indv_gene_ins_mat_for_Root_exp = ins_mat_columns_to_genes(ins_mat_for_Root_exp, ins_mat_for_Root_exp_granges)
save(indv_gene_ins_mat_for_Root_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_Root_exp.rda')

indv_gene_ins_mat_for_Kern_exp = ins_mat_columns_to_genes(ins_mat_for_Kern_exp, ins_mat_for_Kern_exp_granges)
save(indv_gene_ins_mat_for_Kern_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_Kern_exp.rda')

indv_gene_ins_mat_for_L3_exp = ins_mat_columns_to_genes(ins_mat_for_L3_exp, ins_mat_for_L3_exp_granges)
save(indv_gene_ins_mat_for_L3_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_L3_exp.rda')

indv_gene_ins_mat_for_L3Base_exp = ins_mat_columns_to_genes(ins_mat_for_L3Base_exp, ins_mat_for_L3Base_exp_granges)
save(indv_gene_ins_mat_for_L3Base_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_L3Base_exp.rda')

indv_gene_ins_mat_for_L3Tip_exp = ins_mat_columns_to_genes(ins_mat_for_L3Tip_exp, ins_mat_for_L3Tip_exp_granges)
save(indv_gene_ins_mat_for_L3Tip_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_L3Tip_exp.rda')

indv_gene_ins_mat_for_LMAD_exp = ins_mat_columns_to_genes(ins_mat_for_LMAD_exp, ins_mat_for_LMAD_exp_granges)
save(indv_gene_ins_mat_for_LMAD_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_LMAD_exp.rda')

indv_gene_ins_mat_for_LMAN_exp = ins_mat_columns_to_genes(ins_mat_for_LMAN_exp, ins_mat_for_LMAN_exp_granges)
save(indv_gene_ins_mat_for_LMAN_exp, file = 'exp_ins_mat/indv_gene/indv_gene_ins_mat_for_LMAN_exp.rda')

# COMBINE COLUMNS THAT ARE THE SAME GENE

combined_gene_ins_mat_for_shoot_exp = combine_gene_cols(indv_gene_ins_mat_for_shoot_exp)
save(combined_gene_ins_mat_for_shoot_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_shoot_exp.rda')

combined_gene_ins_mat_for_Root_exp = combine_gene_cols(indv_gene_ins_mat_for_Root_exp)
save(combined_gene_ins_mat_for_Root_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_Root_exp.rda')

combined_gene_ins_mat_for_Kern_exp = combine_gene_cols(indv_gene_ins_mat_for_Kern_exp)
save(combined_gene_ins_mat_for_Kern_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_Kern_exp.rda')

combined_gene_ins_mat_for_L3_exp = combine_gene_cols(indv_gene_ins_mat_for_L3_exp)
save(combined_gene_ins_mat_for_L3_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_L3_exp.rda')

combined_gene_ins_mat_for_L3Base_exp = combine_gene_cols(indv_gene_ins_mat_for_L3Base_exp)
save(combined_gene_ins_mat_for_L3Base_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_L3Base_exp.rda')

combined_gene_ins_mat_for_L3Tip_exp = combine_gene_cols(indv_gene_ins_mat_for_L3Tip_exp)
save(combined_gene_ins_mat_for_L3Tip_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_L3Tip_exp.rda')

combined_gene_ins_mat_for_LMAD_exp = combine_gene_cols(indv_gene_ins_mat_for_LMAD_exp)
save(combined_gene_ins_mat_for_LMAD_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_LMAD_exp.rda')

combined_gene_ins_mat_for_LMAN_exp = combine_gene_cols(indv_gene_ins_mat_for_LMAN_exp)
save(combined_gene_ins_mat_for_LMAN_exp, file = 'exp_ins_mat/combined_gene/combined_gene_ins_mat_for_LMAN_exp.rda')

# CONVERT EXPRESSION DATA TO LONG DF
library(tidyverse)
library(rstatix)
library(ggpubr)

exp_shoot_long_df = exp_to_long_df(GShoot_kremling_formatted_v4_hapmapids, ins_mat_for_shoot_exp$Taxa)
save(exp_shoot_long_df, file = 'exp_ins_mat/exp_dat_long/exp_shoot_long_df.rda')

exp_Root_long_df = exp_to_long_df(GRoot_kremling_formatted_v4_hapmapids, ins_mat_for_Root_exp$Taxa)
save(exp_Root_long_df, file = 'exp_ins_mat/exp_dat_long/exp_Root_long_df.rda')

exp_Kern_long_df = exp_to_long_df(Kern_kremling_formatted_v4_hapmapids, ins_mat_for_Kern_exp$Taxa)
save(exp_Kern_long_df, file = 'exp_ins_mat/exp_dat_long/exp_Kern_long_df.rda')

exp_L3_long_df = exp_to_long_df(L3_kremling_formatted_v4_hapmapids, ins_mat_for_L3_exp$Taxa)
save(exp_L3_long_df, file = 'exp_ins_mat/exp_dat_long/exp_L3_long_df.rda')

exp_L3Base_long_df = exp_to_long_df(L3Base_kremling_formatted_v4_hapmapids, ins_mat_for_L3Base_exp$Taxa)
save(exp_L3Base_long_df, file = 'exp_ins_mat/exp_dat_long/exp_L3Base_long_df.rda')

exp_L3Tip_long_df = exp_to_long_df(L3Tip_kremling_formatted_v4_hapmapids, ins_mat_for_L3Tip_exp$Taxa)
save(exp_L3Tip_long_df, file = 'exp_ins_mat/exp_dat_long/exp_L3Tip_long_df.rda')

exp_LMAD_long_df = exp_to_long_df(LMAD_kremling_formatted_v4_hapmapids, ins_mat_for_LMAD_exp$Taxa)
save(exp_LMAD_long_df, file = 'exp_ins_mat/exp_dat_long/exp_LMAD_long_df.rda')

exp_LMAN_long_df = exp_to_long_df(LMAN_kremling_formatted_v4_hapmapids, ins_mat_for_LMAN_exp$Taxa)
save(exp_LMAN_long_df, file = 'exp_ins_mat/exp_dat_long/exp_LMAN_long_df.rda')


# CALCULATE T-TESTS AND RETURN PVALUE FOR EXPRESSION ~ INSERTIONS FOR EACH GENE

pval_exp_ins_mat_for_shoot = entire_df_pval(combined_gene_ins_mat_for_shoot_exp, exp_shoot_long_df)
save(pval_exp_ins_mat_for_shoot , file = 'pval_exp_ins_mat_for_shoot.rda')
save(pval_exp_ins_mat_for_shoot, file = 'pval_exp_ins_mat_for_shoot.csv')

pval_exp_ins_mat_for_Root = entire_df_pval(combined_gene_ins_mat_for_Root_exp, exp_Root_long_df)
save(pval_exp_ins_mat_for_Root , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_Root.rda')
write.csv(pval_exp_ins_mat_for_Root, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_Root.csv')

pval_exp_ins_mat_for_Kern = entire_df_pval(combined_gene_ins_mat_for_Kern_exp, exp_Kern_long_df)
save(pval_exp_ins_mat_for_Kern , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_Kern.rda')
write.csv(pval_exp_ins_mat_for_Kern, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_Kern.csv')

pval_exp_ins_mat_for_L3 = entire_df_pval(combined_gene_ins_mat_for_L3_exp, exp_L3_long_df)
save(pval_exp_ins_mat_for_L3 , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3.rda')
write.csv(pval_exp_ins_mat_for_L3, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3.csv')

pval_exp_ins_mat_for_L3Base = entire_df_pval(combined_gene_ins_mat_for_L3Base_exp, exp_L3Base_long_df)
save(pval_exp_ins_mat_for_L3Base , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3Base.rda')
write.csv(pval_exp_ins_mat_for_L3Base, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3Base.csv')

pval_exp_ins_mat_for_L3Tip = entire_df_pval(combined_gene_ins_mat_for_L3Tip_exp, exp_L3Tip_long_df)
save(pval_exp_ins_mat_for_L3Tip , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3Tip.rda')
write.csv(pval_exp_ins_mat_for_L3Tip, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3Tip.csv')

pval_exp_ins_mat_for_LMAD = entire_df_pval(combined_gene_ins_mat_for_LMAD_exp, exp_LMAD_long_df)
save(pval_exp_ins_mat_for_LMAD , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_LMAD.rda')
write.csv(pval_exp_ins_mat_for_LMAD, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_LMAD.csv')

pval_exp_ins_mat_for_LMAN = entire_df_pval(combined_gene_ins_mat_for_LMAN_exp, exp_LMAN_long_df)
save(pval_exp_ins_mat_for_LMAN , file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_LMAN.rda')
write.csv(pval_exp_ins_mat_for_LMAN, file = 'exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_LMAN.csv')
