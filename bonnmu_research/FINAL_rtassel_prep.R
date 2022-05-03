####### Intro rTASSEL with BonnMu Genotype and Phenotype Data - prepping data for GLM with rTASSEL software#######

wd = "/Users/gcampidilli/Documents/research"
setwd(wd)

library(rTASSEL)
library(data.table)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(stringr)
library(rtracklayer)

# this value will be used to determine the number of different F2 families that must have an insertion in the same gene
# for that gene to be included in downstream analyses
FREQ = 4


### PREPARE AND FORMAT GENOTYPE DATA TO BE LOADED INTO RTASSEL SOFTWARE ####

load("cleaned_bonnMu.RData")
genes = as.data.frame(table(cleaned_bonnMu$GM))
genes = genes[order(-genes$Freq),]

# 4207 GENES ARE INSERTED INTO BY 5+ DIFFERENT F2 FAMILIES
# +1908 4 families
# +2728 3 families
# +3741 2 families
num_genes = length(genes$Freq[genes$Freq > FREQ])
num_genes

# GENES OF INTEREST = GOI
GOI = genes[1:num_genes,]$Var1

# GIVEN 4207 GENES WE ARE GOING TO:
# FIGURE OUT THE CHROMOSOMES TO WHICH THEY BELONG
# FIND AVERAGE POSITION OF INSERSTION INTO THOSE GENES (ROUNDED TO NEAREST INT)
# CREATE HAPMAP DATAFRAME WHICH WE WILL USE IN THE PHENO-GENO ASSOCIATION ANALYSES
# HAPMAP HAS GENES AS ROWS AND F2 FAMILIES AS COLUMNS
# EACH CELL WILL HAVE EITHER "TT" IF INSERTION HAS OCCURED OR "AA" IF NO INSERTION HAS OCCURED

# 1st - for each gene in list, calculate average insertion position 
## FIND OUT HOW TO DO THIS W BRACKET OPERATIONS
position_range = vector()
i = 1
for(i in 1:length(GOI)) {
  position_range[i] = mean(start(ranges(cleaned_bonnMu[grepl(GOI[i], cleaned_bonnMu$GM)])))
} 
position_range = as.integer(trunc(position_range))

## FIND OUT HOW TO DO THIS W BRACKET OPERATIONS
chromosomes = vector()
for(i in 1:length(GOI)) {
  chromosomes[i] = seqnames(cleaned_bonnMu[grepl(GOI[i], cleaned_bonnMu$GM)])
} 

chromosomes = str_replace_all(chromosomes, "Chr", "")
chromosomes = as.integer(chromosomes)

# 2nd - create dataframe of stock names 
fam_names = unique(sort(cleaned_bonnMu$Stock))
#fam_names = str_replace_all(fam_names, ".", "-")
fam_names_df = data.frame(matrix(ncol = length(fam_names), nrow = length(GOI)))
colnames(fam_names_df) = fam_names

# 3rd - create overarching dataframe that will be loaded into rTASSEL as a hapmap
hm_total = data.frame(rs = GOI, alleles = "A/T", chrom = chromosomes, pos = position_range, strand = '+',
                      assembly = NA, center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA)

colnames(hm_total)[1] = "rs#"
colnames(hm_total)[6] = "assembly#"

hm_total = hm_total[order(hm_total$chrom, hm_total$pos),]
hm_total = cbind(hm_total, fam_names_df)
safe = hm_total

rows_hm = nrow(hm_total)

x = 1
for(x in 1:rows_hm) {hm_total[x, unique(cleaned_bonnMu[grepl(hm_total$`rs#`[x], cleaned_bonnMu$GM)]$Stock)] ="TT"}
families = hm_total[,12:ncol(hm_total)]
families[is.na(families)] = "AA"
hm_total = cbind(hm_total[,1:11],families)

colnames(hm_total)[12:ncol(hm_total)] = str_replace_all(colnames(hm_total)[12:ncol(hm_total)], "-",".")

#check that insertions were labeled correctly in hm
summary(grepl("TT", hm_total[20,]))
summary((grepl("Zm00001d027351", cleaned_bonnMu$GM))) # equal to eachother = 9
colnames(hm_total[(grepl("TT", hm_total[20,]))]) # provides list of families from hm
cleaned_bonnMu[grepl("Zm00001d027351", cleaned_bonnMu$GM)]$Stock # list of families from cleaned_bonnMu - lists have same contents!

summary(grepl("TT", hm_total[1,]))
summary((grepl("Zm00001d027242", cleaned_bonnMu$GM))) # equal to eachother = 14

#save hapmap!
save(hm_total, file = "bonnmu_hapmap_FREQ_2PLUS.RData")
write.table(hm_total, "bonnmu_hapmap_FREQ_2PLUS.hmp.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


### CREATE AND LOAD PHENOTYPE DATA INTO RTASSEL  ####

phenoDF <- read.csv("/Users/gcampidilli/Documents/research/bonnmu_phenotypes.csv", header = TRUE)
phenoDF = phenoDF[1:1671,2:4] ## 1680 families

## use cleaned_bonnMu to create easy list of stock  names
bm = data.frame(cleaned_bonnMu) %>% group_by(Stock, img, Genetic_background) %>% dplyr::summarize(count=n())
dim(bm) ## 1671 families

# all same families/images are present
setdiff(phenoDF$Image_name, bm$img)

# assign stock/taxon names to phenoDF
phenoDF_updated = cbind(Taxon = bm$Stock, phenoDF[2:3])

## make taxon IDs consistent across geno and pheno data
phenoDF_updated$Taxon = str_replace_all(phenoDF_updated$Taxon, "-", ".")
#save(phenoDF_updated, file = "RAW_bonnmu_pheno_data.RData")

### transform pheno data with logit link function since our phenotypes are proportions ##

load("bonnmu_pheno_data.RData")
phenoDF_updated$Taxon = str_replace_all(phenoDF_updated$Taxon, "-", ".")
colnames(phenoDF_updated)[1] = "Taxa"
logit_germ = unlist(lapply(phenoDF_updated$prop_ungerminated, function(x)(log(x))/(1-x)))
logit_chlor = unlist(lapply(phenoDF_updated$prop_chlorotic, function(x)(log(x))/(1-x)))
logit_germ = str_replace_all(logit_germ, "-Inf", '0')
logit_chlor = str_replace_all(logit_chlor, "-Inf", '0')
phenoDF_updated$logit_germ = logit_germ
phenoDF_updated$logit_chlorotic = logit_chlor
logit_phenoDF = phenoDF_updated[,c(1,4,5)]
colnames(logit_phenoDF)[3] = "logit_chlor"
logit_phenoDF$logit_germ = as.numeric(logit_phenoDF$logit_germ)
logit_phenoDF$logit_chlor = as.numeric(logit_phenoDF$logit_chlor)


save(logit_phenoDF, file = "logit_phenoDF.RData")
write.table(logit_phenoDF, file = "logit_phenoDF.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# Add following header to txt file otherwise rTASSEL analyses will not work:
# <Phenotype>
# taxa  data  data

# move to rTASSEL software
# 1) Union join logit_phenoDF and geno hapmap
# 2) Select logit_phenoDF + geno_hapmap and perform GLM analysis
# 3) pvalues are provided in GLM_stat table result. Additive model is used

# The table shows the F-statistics and p-values for the requested F-tests for the main and additive models, 
# It also contains marker_Rsq, mean squares (MS) and degrees of freedom (DF) for the marker effect, for the model (corrected for the mean), 
# and for error. If taxa are replicated (across reps or environments), then the markers are tested using the taxa within marker mean square. 
# If taxa are unreplicated, then the residual mean square is used. 
# Marker_Rsq is the marginal R-squared for the marker calculated as SS Marker (after fitting all other model terms) / SS Total, where SS stands for sum of squares. 












