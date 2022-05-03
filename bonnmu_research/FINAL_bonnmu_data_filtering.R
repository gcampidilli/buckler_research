# bonnMu data filtering

## COME BACK TO THIS
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

bonnMu =import.gff3("BonnMu.fromMaizeGDBbrowser.gff3.gz", format = "gff3")

####### FILTER BONNMU GENO DATASET  ########

## 260 rows with "B73V4" in seqnames col, no images to go along with these -> 5-B through 6-C
## 1 row with "Mt" in seqnames col

## remove rows with B73V4 and Mt as seqname ID, update levels #####
intermediate = bonnMu[-c(which(grepl("B73V4", seqnames(bonnMu)))),] 
cleaned_bonnMu = intermediate[-c(which(grepl("Mt", seqnames(intermediate)))),]
seqlevels(cleaned_bonnMu) = c(paste("Chr", 1:10, sep = ""))

cleaned_bonnMu = cleaned_bonnMu[-c(which(grepl("3-E",cleaned_bonnMu$Stock))),]
cleaned_bonnMu = cleaned_bonnMu[-c(which(grepl("None",cleaned_bonnMu$img))),]

cleaned_bonnMu_rm_double = cleaned_bonnMu[!duplicated(start(ranges(cleaned_bonnMu))), ]


save(cleaned_bonnMu_rm_double, file = "cleaned_bonnMu_rm_double.RData")

### SUMMARY OF CURRENT GENO DATA ####
length(unique(cleaned_bonnMu$Stock)) #1671 families
length(cleaned_bonnMu$GM) #60,000 total insertions into 
length(unique(cleaned_bonnMu$GM)) # 18,510 genes
length(unique(ranges(cleaned_bonnMu))) #57,240 different insertion positions

###### crossreference images  between bonnMu geno dataset and bonnMu pheno dataset ######
###### half of 5-B, all of 6-B and 6-C listed in BonnMu geno data dont have images ########

images = read.csv("/Users/gcampidilli/Documents/research/bonnmu_phenotypes.csv", header = TRUE)
images = images[,1]

########### EXAMINE ERROR/INCONSISTENCIES BETWEEN GENO AND PHENO DATA #################

# 204 images that are unique to the image database & NOT in BonnMu geno dataset. Primarily 3-E, a few 4-A, 4-B, 5-B
###### out of 204 images that are NOT in geno dataset, 195 of these images of 3-E stocks ####
unique_images = setdiff(images$Image_name, cleaned_bonnMu$img)
length(unique_images[grepl("3-E", unique_images)]) # = 195

###### pheno data provides 520 images for 3-E, where each image is representative of a 3-E stock #####
length(images[grepl("3-E", images$Image_name),])

###### take closer look at 3-E, all 3-E is geno Co125, whereas the rest of the dataset is B73 ######
co125 = bonnMu[grepl("Co125", bonnMu$Genetic_background)]
length(unique(co125$Stock)) # = 358 

e_3 = bonnMu[grepl("3-E", bonnMu$Stock)]
length(unique(e_3$Stock)) # = 358 

###### geno data includes 358 3-E stocks, with only 325 of those stocks having pheno data ####
co125_no_img = co125[grepl("None", co125$img),]
length(unique(co125_no_img$Stock)) # = 33 (358-33=325)

###### 325 3-E images in geno dataset + 195 3-E images NOT in geno dataset = 520 images total, which is consistent with pheno data ####


