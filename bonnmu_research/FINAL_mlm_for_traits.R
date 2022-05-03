
# extending linear model from chlor_pheno ~ copy_number + insertion_into_gene to dividing copy number into category
# do the same for germ_pheno
# example: chlor ~ 5utr_copy_number + intron_copy_number + cds_copy_number + 3utr_copy_number + promoter_copy_number + intergenic_copy_number + insertion_into_gene

library(data.table)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(rtracklayer)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")

wd = "/Users/gcampidilli/Documents/research"
setwd(wd)

# this value will be used to determine the proportion of genes selected for downstream analyses
CUTOFF = 0.05

# go to FINAL_insertion_categorization.R for construction of cleaned_bonnMu_rm_double.RData
load("cleaned_bonnMu_rm_double.RData")
head(cleaned_bonnMu_rm_double)

# cbm_cat_updated created in FINAL_insertion_categorization.R
load("cbm_cat_updated.RData")
head(cbm_cat_updated)


##############  not pulling in genes from literature due to low significance in previous analyses ########

# construct df that will be used for linear model chlorosis ~ copy # + 0/1 insertion in associated gene for each category

load("logit_phenoDF.RData")
fam_names = unique(sort(cleaned_bonnMu_rm_double$Stock))
copy_number = data.frame(cleaned_bonnMu_rm_double) %>% group_by(Stock) %>% dplyr::summarize(count=n())
insertions = matrix(nrow = 1671, ncol = 6); colnames(insertions) = c("utr5", "cds","intron","utr3","promoter","upstream")

chlor_lm_df = data.frame(copy_num = copy_number$count, chlor_pheno = logit_phenoDF$logit_chlor, insertions)
rownames(chlor_lm_df) = fam_names

# identify which families have insertions in top 5% [Current CUTOFF value = 5%] of logit_glm for chlor trait
# using pvals from ["logit_fiveplus_glm"]: GLM analysis conducted on genes that are inserted into in at least 5 
# [Current FREQ value = 5] different F2 families

logit_glm = read.csv("glm_stats_FREQ_5PLUS.csv", header = T)
logit_glm$Trait = str_replace_all(logit_glm$Trait, "logit_chlor", "chlor")
chlor_glm= logit_glm[grepl("chlor", logit_glm$Trait),]

chlor_log10_p = -log10(chlor_glm$p)
chlor_glm = cbind(chlor_glm, chlor_log10_p)
chlor_glm = chlor_glm[order(-chlor_glm$chlor_log10_p),]

#save(chlor_glm, file = "chlor_glm_FREQ_5PLUS_TOTAL.RData")
#load("chlor_glm_FREQ_5PLUS_TOTAL.RData")

chlor_genes = data.frame(chlor_glm$Marker)
# subset_cutoff = nrow(chlor_genes)*CUTOFF
# geneset = data.frame(chlor_genes[1:subset_cutoff,])
geneset = 
colnames(geneset) = "genes"
colnames(chlor_lm_df)[8] = "upstream_of_promoter"



# assign has_insertion to 1 for families that have insertion in a gene in union_geneset
get_insertion_categories_chlor = function(x){
        # for current bonnMu family, create df that includes all insertions and insertion info from cbm_cat for the current family
        fam_insertions = cbm_cat_updated[grepl(rownames(chlor_lm_df[x,]), cbm_cat_updated$Stock),]
        # if any of the insertions insert into a gene of interest, record which category it is inserting into
        fam_chlor_insertions=fam_insertions$insertion_category[fam_insertions$GM %in% geneset$genes]
        return(unique(fam_chlor_insertions))
}

x = 1
for(x in 1:nrow(chlor_lm_df)){
        ins_cats=get_insertion_categories_chlor(x)
        ins_cats=ins_cats[!is.na(ins_cats)]
        if(length(ins_cats)>0){
                chlor_lm_df[x,ins_cats]=1
        }
}

chlor_lm_df[is.na(chlor_lm_df)] = 0

chlor_lm_df$chlor_pheno = as.numeric(chlor_lm_df$chlor_pheno)

#save(chlor_lm_df, file = "chlor_lm_df_FREQ_5PLUS_TOTAL.RData")

chlor_lm = lm(chlor_pheno ~ copy_num + utr5 + cds + intron + utr3 + promoter + upstream_of_promoter + intergenic, data = chlor_lm_df)
summary(chlor_lm)

sink("/Users/gcampidilli/Documents/research/lm_outputs/chlor_lm_summary_FREQ_5PLUS_TOTAL.txt")
print(summary(chlor_lm))
sink()
########################    DO SAME FOR GERM TRAIT     ################################################

germ_lm_df = data.frame(copy_num = copy_number$count, germ_pheno = logit_phenoDF$logit_germ, insertions)
rownames(germ_lm_df) = fam_names

# identify which families have insertions in top 5% [Current CUTOFF value = 5%] of logit_glm for germ trait
# using pvals from ["logit_fiveplus_glm"]: GLM analysis conducted on genes that are inserted into in at least 5 
# [Current FREQ value = 5] different F2 families

logit_glm = read.csv("glm_stats_FREQ_5PLUS.csv", header = T)
logit_glm$Trait = str_replace_all(logit_glm$Trait, "logit_germ", "germ")
germ_glm= logit_glm[grepl("germ", logit_glm$Trait),]

germ_log10_p = -log10(germ_glm$p)
germ_glm = cbind(germ_glm, germ_log10_p)
germ_glm = germ_glm[order(-germ_glm$germ_log10_p),]

# save(germ_glm, file = "germ_glm_FIVEPLUS_INS_TOTAL.RData")
# load("germ_glm_FIVEPLUS_INS_TOTAL.RData")

germ_genes = data.frame(germ_glm$Marker)
subset_cutoff = nrow(germ_genes)*CUTOFF
germ_geneset = data.frame(germ_genes[1:subset_cutoff,])
colnames(germ_geneset) = "genes"
colnames(germ_lm_df)[8] = "upstream_of_promoter"


# assign has_insertion to 1 for families that have insertion in a gene in germ_geneset
get_insertion_categories = function(x){
        ## all insertions in that bonnmu family
        fam_insertions = cbm_cat_updated[grepl(rownames(germ_lm_df[x,]), cbm_cat_updated$Stock),]
        fam_germ_insertions=fam_insertions$insertion_category[fam_insertions$GM %in% germ_geneset$genes]
        return(unique(fam_germ_insertions))
}

x = 1
for(x in 1:nrow(germ_lm_df)){
        ins_cats=get_insertion_categories(x)
        ins_cats=ins_cats[!is.na(ins_cats)]
        if(length(ins_cats)>0){
                germ_lm_df[x,ins_cats]=1
        }
}

germ_lm_df[is.na(germ_lm_df)] = 0


#save(germ_lm_df, file = "germ_lm_df_FREQ_5PLUS_TOTAL.RData")


germ_lm = lm(germ_pheno ~ copy_num + utr5 + cds + intron + utr3 + promoter + upstream_of_promoter + intergenic, data = germ_lm_df)
summary(germ_lm)

sink("/Users/gcampidilli/Documents/research/lm_outputs/germ_lm_summary_FREQ_5PLUS_TOTAL.txt")
print(summary(germ_lm))
sink()
rm(CUTOFF)

### make summary pheno file with copynum
logit_phenoDF_w_copynum = data.frame(chlor_lm_df[,1:2], germ_lm_df[,2])
colnames(logit_phenoDF_w_copynum)[3] = "germ_pheno"
save(logit_phenoDF_w_copynum, file = "logit_phenoDF_w_copynum.RData")


