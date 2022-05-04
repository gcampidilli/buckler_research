########################################################################################
#     LOAD PACKAGES + DATA, FILTER, AND FORMAT

library(data.table)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(GenomicFeatures)

wd = "/Users/gcampidilli/Documents/research/mu_tir/step_6"
setwd(wd)

# inbred_te_count is from "find . -name '*.csv' | xargs wc -l > num_ins_step5.txt" command in terminal
# in insertion_dat_categorized directory
# inbred_te_count txt file formatting at bottom of script

inbred_te_count = read.table("../step_5/te_count_step5.txt", head = T)

# the TE count we see here will likely be greater than the number of entries each inbred has in the insertion_matrix
ins_mat_te_count = read.csv("../step_4/ins_matrix_te_count.csv", head = F)

# ex. we have 75 insertions for 33-16 according to te_count_step5.txt, but there are only 57 insertions for 33-16 in the insertion matrix (0/1)
# if we look at 33-16 insertion_dat_categorized, we see that 10 insertions are unassigned to a gene region, leaving us with 65 insertions
# though 65 and 57 are not super close, 65 is the more precise number, as the insertion_matrix construction rounds to nearest 1000 bp
# use te_count from step_5 for visualization

goodman_panel_info = read.csv("../flint_garcia_panel_info.csv", head = T)

setdiff(inbred_te_count$inbred, goodman_panel_info$Inbred)
setdiff(goodman_panel_info$Inbred, inbred_te_count$inbred)


# random note - search flint=garcia in finder to find rscripts that discuss coverage

# go to te_count text document CHANGE FLINT GARCIA PANEL INFO, remove "Goodman_Buckler
# include flint-garcia info file in doc with updated names

gpi_reduced = goodman_panel_info[which(goodman_panel_info$Inbred %in% inbred_te_count$inbred),]
gpi_reduced = gpi_reduced[order(gpi_reduced$Inbred),]

# check to make sure all inbred names in same order
table(gpi_reduced$Inbred == inbred_te_count$inbred)

gpi_reduced$te_count = inbred_te_count$te_count

# load te count summary so we can subtract tes that weren't assigned in each inbred

gr_summary_total = read.csv("../step_5/total_generegion_summary.csv", head = T)
gr_summary_total = gr_summary_total[-c(1),-c(1)]

not_assigned_te = gr_summary_total[which(gr_summary_total$region == "Not assigned"),]
not_assigned_te = not_assigned_te[order(not_assigned_te$inbred_name_vector),]

# are we missing any inbreds in inbred_te_count - no
setdiff(not_assigned_te$inbred_name_vector, gpi_reduced$Inbred) 
gpi_reduced$te_count_assigned = gpi_reduced$te_count - not_assigned_te$insertion_count

# add columns to gpi for proportion of insertions in each gene region 
gpi_reduced$promoter_ins_prop = gr_summary_total[which(gr_summary_total$region == "Promoter"),]$insertion_proportion
gpi_reduced$utr5_ins_prop = gr_summary_total[which(gr_summary_total$region == "UTR5"),]$insertion_proportion
gpi_reduced$exon_ins_prop = gr_summary_total[which(gr_summary_total$region == "Exon"),]$insertion_proportion
gpi_reduced$intron_ins_prop = gr_summary_total[which(gr_summary_total$region == "Intron"),]$insertion_proportion
gpi_reduced$utr3_ins_prop = gr_summary_total[which(gr_summary_total$region == "UTR3"),]$insertion_proportion
gpi_reduced$intergenic_ins_prop = gr_summary_total[which(gr_summary_total$region == "Intergenic"),]$insertion_proportion

gpi_reduced$promoter_ins_per_1kb = gr_summary_total[which(gr_summary_total$region == "Promoter"),]$insertions_per_1kb
gpi_reduced$utr5_ins_per_1kb = gr_summary_total[which(gr_summary_total$region == "UTR5"),]$insertions_per_1kb
gpi_reduced$exon_ins_per_1kb = gr_summary_total[which(gr_summary_total$region == "Exon"),]$insertions_per_1kb
gpi_reduced$intron_ins_per_1kb = gr_summary_total[which(gr_summary_total$region == "Intron"),]$insertions_per_1kb
gpi_reduced$utr3_ins_per_1kb = gr_summary_total[which(gr_summary_total$region == "UTR3"),]$insertions_per_1kb
gpi_reduced$intergenic_ins_per_1kb = gr_summary_total[which(gr_summary_total$region == "Intergenic"),]$insertions_per_1kb

gpi_gene_region_summary = gpi_reduced
save(gpi_gene_region_summary, file = "gpi_gene_region_summary.RData")



########################################################################################
# TE Count (assigned), histogram

load("gpi_gene_region_summary.RData")
gpi_gene_region_summary[which(gpi_gene_region_summary$te_count_assigned > 300),]$Inbred
# * "A6"  "B57"  "VaW6" excluded from barchart
adj_gene_region_summary = gpi_gene_region_summary[which(gpi_gene_region_summary$te_count_assigned < 300),]
ggplot(data = adj_gene_region_summary, aes(x=te_count_assigned)) + geom_histogram(bins = 30)+
  labs(title = "Inbred Mutator TE Count Histogram*") + ylab("Frequency") + xlab("TE Count") +
  theme_gray(base_size = 14) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+theme(axis.title.x = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(size = 15)) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))  


# TE Count by population, histogram
   
ggplot(data = adj_gene_region_summary, aes(x=te_count_assigned)) + geom_histogram(bins = 30)+
  labs(title = "Inbred Mutator TE Count Histogram*") + ylab("Frequency") + xlab("TE Count") +
  theme_gray(base_size = 14) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+theme(axis.title.x = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(size = 15)) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))  +
  facet_wrap(~Subpopulation)

quantile(gpi_gene_region_summary[which(gpi_gene_region_summary$Subpopulation == "mixed"),]$te_count_assigned)
quantile(gpi_gene_region_summary[which(gpi_gene_region_summary$Subpopulation == "nss"),]$te_count_assigned)
quantile(gpi_gene_region_summary[which(gpi_gene_region_summary$Subpopulation == "popcorn"),]$te_count_assigned)
quantile(gpi_gene_region_summary[which(gpi_gene_region_summary$Subpopulation == "ss"),]$te_count_assigned)
quantile(gpi_gene_region_summary[which(gpi_gene_region_summary$Subpopulation == "sweet"),]$te_count_assigned)
quantile(gpi_gene_region_summary[which(gpi_gene_region_summary$Subpopulation == "ts"),]$te_count_assigned)

########################################################################################
# take top 10 of each category of ins_prop and ins_per_1kb for each gene region
load("gpi_gene_region_summary.RData")

order_by_count = gpi_gene_region_summary[order(-gpi_gene_region_summary$te_count_assigned),]

sink("top10_bottom10.txt")
print("Inbreds with highest insertion count:")
order_by_count$Inbred[1:10]
order_by_count[1:10,]

print("Inbreds with lowest insertion count:")
order_by_count$Inbred[261:271]
order_by_count[261:271,]
sink()


########################################################################################
# linear model

lm1 = glm(te_count_assigned~ Subpopulation, data = gpi_gene_region_summary)
summary(lm1) # population is not significant factor

# te_count_assigned is only significant for proportion of promoter insertions, no other gene regions
#sink(file = "ins_prop_linear_models.txt")
lm2 = glm(promoter_ins_prop~te_count_assigned + Subpopulation, data = gpi_gene_region_summary)
summary(lm2)

lm22 = glm(promoter_ins_prop~te_count_assigned + Subpopulation + te_count_assigned*Subpopulation, data = gpi_gene_region_summary)
summary(lm22) # no interaction effect between population and te count

lm3 = glm(utr5_ins_prop~te_count_assigned + Subpopulation, data = gpi_gene_region_summary)
summary(lm3)

lm4 = glm(exon_ins_prop~te_count_assigned + Subpopulation, data = gpi_gene_region_summary)
summary(lm4)

lm5 = glm(intron_ins_prop~te_count_assigned + Subpopulation, data = gpi_gene_region_summary)
summary(lm5)

lm6 = glm(intergenic_ins_prop~te_count_assigned + Subpopulation, data = gpi_gene_region_summary)
summary(lm6)
#sink()

########################################################################################
# EXAMINE DATA BY GENE REGION, BOXPLOTS

name_vector = rep(gpi_gene_region_summary$Inbred, 6)
gene_region_vector = rep(c("Promoter", "UTR5", "Exon", "Intron", "UTR3", "Intergenic"), each = 271)
ins_prop_vector = c(gpi_gene_region_summary$promoter_ins_prop, gpi_gene_region_summary$utr5_ins_prop, gpi_gene_region_summary$exon_ins_prop,
                    gpi_gene_region_summary$intron_ins_prop, gpi_gene_region_summary$utr3_ins_prop, gpi_gene_region_summary$intergenic_ins_prop)
ins_per_kb_vector = c(gpi_gene_region_summary$promoter_ins_per_1kb, gpi_gene_region_summary$utr5_ins_per_1kb, gpi_gene_region_summary$exon_ins_per_1kb,
                      gpi_gene_region_summary$intron_ins_per_1kb, gpi_gene_region_summary$utr3_ins_per_1kb, gpi_gene_region_summary$intergenic_ins_per_1kb)

boxplot_df = data.frame(inbred = name_vector, gene_region = gene_region_vector, ins_prop = ins_prop_vector,ins_per_kb = ins_per_kb_vector)
boxplot_df$gene_region = factor(boxplot_df$gene_region, levels = c("Promoter", "UTR5", "Exon", "Intron", "UTR3", "Intergenic"))

boxplot_df = boxplot_df[which(boxplot_df$ins_per_kb<0.001),]

ggplot(data = boxplot_df, aes(x=gene_region, y=ins_per_kb, fill=gene_region)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Mutator Insertions per 1Kb by Gene Region and Population") + ylab("Insertions per 1000 bp") + xlab("") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw()

########################################################################################
#     EXAMINE DATA BY GENE REGION AND POPULATION TYPE GROUPS, BARCHARTS

group_size = as.numeric(table(gpi_gene_region_summary$Subpopulation))
group_names = names(table(gpi_gene_region_summary$Subpopulation))
population_avgs = aggregate(.~gpi_gene_region_summary$Subpopulation, data = gpi_gene_region_summary[,10:ncol(gpi_gene_region_summary)], mean)
colnames(population_avgs)[1] = "population"
save(population_avgs,file= "population_avgs.RData")

population_sd = aggregate(.~gpi_gene_region_summary$Subpopulation, data = gpi_gene_region_summary[,10:ncol(gpi_gene_region_summary)],sd)

population_se = data.frame()
for(x in 1:length(group_size)) {
  temp_vector = population_sd[x,2:ncol(population_sd)]/rep(sqrt(group_size[x]), ncol(population_sd)-1)
  population_se = rbind(population_se, temp_vector)
}

population_se$population = group_names

save(population_se,file= "population_se.RData")

name_vector = rep(group_names, 6)
gene_region_vector = rep(c("Promoter", "UTR5", "Exon", "Intron", "UTR3", "Intergenic"), each = 6)
ins_prop_vector = c(population_avgs$promoter_ins_prop, population_avgs$utr5_ins_prop, population_avgs$exon_ins_prop,
                    population_avgs$intron_ins_prop, population_avgs$utr3_ins_prop, population_avgs$intergenic_ins_prop)

ins_prop_se_vector = c(population_se$promoter_ins_prop,population_se$utr5_ins_prop,population_se$exon_ins_prop,
                       population_se$intron_ins_prop,population_se$utr3_ins_prop,population_se$intergenic_ins_prop)

ins_per_kb_vector = c(population_avgs$promoter_ins_per_1kb, population_avgs$utr5_ins_per_1kb,population_avgs$exon_ins_per_1kb,
                      population_avgs$intron_ins_per_1kb,population_avgs$utr3_ins_per_1kb,population_avgs$intergenic_ins_per_1kb)

ins_per_kb_se_vector = c(population_se$promoter_ins_per_1kb, population_se$utr5_ins_per_1kb,population_se$exon_ins_per_1kb,
                         population_se$intron_ins_per_1kb,population_se$utr3_ins_per_1kb,population_se$intergenic_ins_per_1kb)

plot_gr_summary_df = data.frame(population = name_vector, gene_region = gene_region_vector, ins_prop = ins_prop_vector,
                                ins_prop_se = ins_prop_se_vector, ins_per_kb = ins_per_kb_vector, ins_per_kb_se = ins_per_kb_se_vector)

plot_gr_summary_df$population = factor(plot_gr_summary_df$population, levels = c("mixed", "nss", "popcorn","ss", "sweet", "ts"))
plot_gr_summary_df$gene_region = factor(plot_gr_summary_df$gene_region, levels = c("Promoter", "UTR5", "Exon", "Intron", "UTR3", "Intergenic"))

rm(name_vector, gene_region_vector, ins_prop_se_vector, ins_per_kb_vector, ins_prop_vector, ins_per_kb_se_vector)
save(plot_gr_summary_df, file = "plot_gr_summary_df.RData")

# barplots

ins_prop_plot = ggplot(data = plot_gr_summary_df, aes(fill=population, y=ins_prop, x=gene_region)) +
  geom_bar(position="dodge", stat="identity") +
  labs(title = "Proportion of Mutator Insertions by Gene Region and Population") + ylab("Average proportion of insertions") + xlab("") +
  theme_gray(base_size = 14) + geom_text(aes(label=round(ins_prop, digits = 3)), position=position_dodge(width=0.9), vjust=-0.5, size = 2) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 15, face = "bold", color = "black")) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))     

kb_plot = ggplot(data = plot_gr_summary_df, aes(fill=population, y=ins_per_kb, x=gene_region)) +
  geom_bar(position="dodge", stat="identity") +
  labs(title = "Mutator Insertions per 1Kb by Gene Region and Population") + ylab("Insertions per 1000 bp") + xlab("") +
  theme_gray(base_size = 14) + geom_text(aes(label=round(ins_per_kb, digits = 8)), position=position_dodge(width=0.9), vjust=-0.5, size = 1) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 15, face = "bold", color = "black")) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))


########################################################################################
# PCA AND ALLELIC FREQUENCY, INS_MAT ROUNDED TO NEAREST 1000
########################################################################################

wd = "/Users/gcampidilli/Documents/research/mu_tir/step_6"
setwd(wd)

load("gpi_gene_region_summary.RData")
library(ggfortify)

ins_mat = read.csv("../step_4/jan21_insertion_matrix_thousand.csv", head = T, row.names = 1)
ins_mat = ins_mat[,-c(grep("scaf_", colnames(ins_mat)))]
colnames(ins_mat) = substring(colnames(ins_mat), first=4)

# assign population to inbreds (?)
row.names(ins_mat) == gpi_gene_region_summary$Inbred # rows are not in same order, use merge() func
ins_mat$Inbred = row.names(ins_mat)
ins_mat_merged = merge(gpi_gene_region_summary[,c(1,9)], ins_mat, by = "Inbred")
ins_mat_merged$Subpopulation = factor(ins_mat_merged$Subpopulation, levels = c("mixed", "nss", "popcorn", "ss", "sweet", "ts"))
ins_mat = ins_mat[,-c(ncol(ins_mat))] # remove inbred column from ins_mat
pca_dat = prcomp(ins_mat_merged[,c(3:ncol(ins_mat_merged))])

# if want to see specific inbred names
row.names(ins_mat_merged) = ins_mat_merged$Inbred
autoplot(pca_dat, data = ins_mat_merged, colour = 'Subpopulation', shape = F,label.size = 3)
autoplot(pca_dat,x=1,y=2, data = ins_mat_merged, colour = 'Subpopulation')
autoplot(pca_dat,x=1,y=2, data = ins_mat_merged, colour = 'Subpopulation', shape = F,label.size = 3)


autoplot(pca_dat,x=2,y=3, data = ins_mat_merged, colour = 'Subpopulation', shape = F,label.size = 3)
autoplot(pca_dat,x=2,y=3, data = ins_mat_merged, colour = 'Subpopulation')

# no inbred names, just dots
autoplot(pca_dat,x=1,y=3, data = ins_mat_merged, colour = 'Subpopulation', shape = F,label.size = 3)
autoplot(pca_dat,x=1,y=3, data = ins_mat_merged, colour = 'Subpopulation')



########################################################################################
# insertion matrix - look at locations that are most frequently inserted into across inbreds

# Allele frequency spectrum (x axis allele frequency colSum/#individuals, y axis count at each bin)
# - are these evolving under neutral genetic drift, or is selection pushing them to lower frequencies
# identify genes at high and low allelic frequencies


# ins_mat = read.csv("../step_4/jan21_insertion_matrix_exact.csv", head = T, row.names = 1)
# colnames(ins_mat) = substring(colnames(ins_mat), first=4)

allele_freq = unlist(lapply(1:ncol(ins_mat), function(x){
  print(x)
  sum(ins_mat[,x])/nrow(ins_mat)
}))

allele_freq = round(allele_freq, digits = 5)
allele_location = colnames(ins_mat)

allele_df = data.frame(location = allele_location, freq = allele_freq)
write.csv(allele_df, file="allele_df_thousand.csv")

ggplot(data = allele_df, aes(x=freq)) + geom_histogram(bins = 40)+
  labs(title = "Mutator TE Allele Frequency Spectrum - thousand location") + ylab("Count") + xlab("Allele Frequency") +
  theme_gray(base_size = 14) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+theme(axis.title.x = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(size = 15)) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))  


########################################################################################
# allelic frequency plotted against chromosomal location

library(dplyr)
library(tidyr)
library(stringr)

# allele_df_ordered = read.csv("allele_df_ordered.csv", head = T)
# allele_df_ordered = allele_df_ordered[,-c(1)]

location_cols = str_split_fixed(allele_df$location, "_", 2)
allele_df_loc = data.frame(chrom = as.numeric(location_cols[,1]), bp = as.numeric(location_cols[,2]), freq = allele_df$freq)
allele_df_loc = allele_df_loc[order(allele_df_loc$chrom, allele_df_loc$bp),]
allele_df_loc = allele_df_loc[!is.na(allele_df_loc$chrom),]

allele_df_loc$chrom = factor(allele_df_loc$chrom, levels = c(1,2,3,4,5,6,7,8,9,10))
ggplot(data = allele_df_loc, aes(x=bp, y = freq)) + geom_point()+
  labs(title = "Mutator TE Allele Frequency by chromosome and bp location - thousand location") + ylab("Frequency") + xlab("Basepair") +
  theme_gray(base_size = 14) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+theme(axis.title.x = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))  +
  facet_wrap(~chrom)

