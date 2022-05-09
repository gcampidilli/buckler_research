# combine all pvalue tables into one

load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_shoot.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_Root.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_Kern.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3Base.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_L3Tip.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_LMAD.rda')
load('exp_ins_mat/pval_t_test_df/pval_exp_ins_mat_for_LMAN.rda')

all_possible_genes = unique(c(pval_exp_ins_mat_for_Kern$genes, pval_exp_ins_mat_for_shoot$genes, pval_exp_ins_mat_for_Root$genes,
                              pval_exp_ins_mat_for_L3$genes, pval_exp_ins_mat_for_L3Base$genes, pval_exp_ins_mat_for_L3Tip$genes,
                              pval_exp_ins_mat_for_LMAD$genes, pval_exp_ins_mat_for_LMAN$genes))

pval_df = data.frame(matrix(NA, nrow = length(all_possible_genes), ncol = 8))
rownames(pval_df) = all_possible_genes
colnames(pval_df) = c('Root', 'Shoot', 'Kern','L3','L3Base','L3Tip','LMAD','LMAN')

shoot_pval = data.frame(shoot_pval = pval_exp_ins_mat_for_shoot$pval)
rownames(shoot_pval) = pval_exp_ins_mat_for_shoot$genes

root_pval = data.frame(root_pval = pval_exp_ins_mat_for_Root$pval)
rownames(root_pval) = pval_exp_ins_mat_for_Root$genes

kern_pval = data.frame(kern_pval = pval_exp_ins_mat_for_Kern$pval)
rownames(kern_pval) = pval_exp_ins_mat_for_Kern$genes

l3_pval = data.frame(l3_pval = pval_exp_ins_mat_for_L3$pval)
rownames(l3_pval) = pval_exp_ins_mat_for_L3$genes

l3tip_pval = data.frame(l3tip_pval = pval_exp_ins_mat_for_L3Tip$pval)
rownames(l3tip_pval) = pval_exp_ins_mat_for_L3Tip$genes

l3base_pval = data.frame(l3base_pval = pval_exp_ins_mat_for_L3Base$pval)
rownames(l3base_pval) = pval_exp_ins_mat_for_L3Base$genes

lmad_pval = data.frame(lmad_pval = pval_exp_ins_mat_for_LMAD$pval)
rownames(lmad_pval) = pval_exp_ins_mat_for_LMAD$genes

lman_pval = data.frame(lman_pval = pval_exp_ins_mat_for_LMAN$pval)
rownames(lman_pval) = pval_exp_ins_mat_for_LMAN$genes




df1 = transform(merge(pval_df, shoot_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df2 = transform(merge(df1, root_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df3 = transform(merge(df2, kern_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df4 = transform(merge(df3, l3_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df5 = transform(merge(df4, l3tip_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df6 = transform(merge(df5, l3base_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df7 = transform(merge(df6, lmad_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)
df8 = transform(merge(df7, lman_pval, by = 'row.names'),row.names=Row.names, Row.names=NULL)


combined_pval_t_test_df = df8[,-c(1:8)]

save(combined_pval_t_test_df, file = 'exp_ins_mat/combined_pval_t_test_df.rda')
write.csv(combined_pval_t_test_df, file = 'exp_ins_mat/combined_pval_t_test_df.csv')



