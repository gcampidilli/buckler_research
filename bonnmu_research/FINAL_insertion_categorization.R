###### CONSTRUCT CBM_CAT_UPDATED.RDATA from B73 reference genome transcript cleaned_bonnMu_rm_double ####
###### cleaned_bonnMu_rm_double created at the end of this script #####
###### this script resolves 1) transcripts from reference genome that overlap in 2+ categorical regions and 2) insertions that span 2+ categorical regions


library(data.table)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(rtracklayer)
library(xml2)
library(GenomicFeatures)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")

wd = "/Users/gcampidilli/Documents/research"
setwd(wd)


# go to updated_insertion_count.R for construction of cleaned_bonnMu_rm_double.RData
load("cleaned_bonnMu_rm_double.RData")
head(cleaned_bonnMu_rm_double)

# cbm_cat created in insertion_categorization.R
load("cbm_cat_updated.RData")
head(cbm_cat_updated)

#####################
# original code sourced from gene_region_geno_pheno_assoc.R and updated_insertion_count.R for breakdown of intergenic category
## Load maize GFF data into txdb object
# Transcript database
txdb <- makeTxDbFromGFF("Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz", format="gff3")
gff_db =import.gff3("Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz", format = "gff3")

# refine txdb object so it only includes chromosomes 1:10
seqlevels(txdb) = as.character(1:10)

# transcripts = set of exons, a gene is composed of alternatively spliced transcripts
transcripts <- transcripts(txdb)
first_transcripts = transcripts[grepl("T001", transcripts$tx_name)] # 40,923

# use.names = T, group names are transcript names "tx_name" (i.e. Zm00001d027235_T001), 
# if use.names is false, group names are "tx_id" (i.e. 67598)
introns_bt = unlist(intronsByTranscript(txdb, use.names = T)) #1mil+ introns
introns_bt01 = introns_bt[grepl("T001", names(introns_bt))] #152472 introns

utr5_bt = unlist(fiveUTRsByTranscript(txdb, use.names = T)) #248164 utr5
utr5_bt01 = utr5_bt[grepl("T001", names(utr5_bt))] # 35228 utr5

utr3_bt = unlist(threeUTRsByTranscript(txdb, use.names = T)) # 251982 utr3
utr3_bt01 = utr3_bt[grepl("T001", names(utr3_bt))] # 33721 utr3

cds_bt = unlist(cdsBy(txdb, use.names = T)) # 918634 cds
cds_bt01 = cds_bt[grepl("T001", names(cds_bt))] # 170269 cds

# can make analysis more specific by using the promoter region -500 to 0 bp and the region upstream of promoter -2000 to -501 bp
# as a complement to intergenic category

# promoter is at correct location at + and - strand
promoters_bt = promoters(txdb, upstream=500, downstream = 0, use.names=TRUE) #133940 promoter regions
promoters_bt01 = promoters_bt[grepl("T001", names(promoters_bt))] #41277 promoter regions
length(reduce(promoters_bt01)) # 40979

upstream_of_promoterbt = promoters(txdb, upstream = 2000, downstream = 0, use.names = TRUE)
upstream_of_promoterbt_01 = upstream_of_promoterbt[grepl("T001", names(upstream_of_promoterbt))]
# if on positive: 
end(ranges(upstream_of_promoterbt_01[strand(upstream_of_promoterbt_01) == "+"])) = 
  end(ranges(upstream_of_promoterbt_01[strand(upstream_of_promoterbt_01) == "+"])) - 500
#end(ranges(upstream_of_promoterbt_01)) = end(ranges(upstream_of_promoterbt_01))-500
# if on negative: 
start(ranges(upstream_of_promoterbt_01[strand(upstream_of_promoterbt_01) == "-"])) = 
  start(ranges(upstream_of_promoterbt_01[strand(upstream_of_promoterbt_01) == "-"])) + 500
width(ranges(upstream_of_promoterbt_01)[40000]) # 1500 bp

table(width(ranges(upstream_of_promoterbt_01)))
intersect(promoters_bt01, upstream_of_promoterbt_01) # 859 intersections
intersect(upstream_of_promoterbt_01, promoters_bt01) # 0 intersections

chromInfo <- seqinfo(txdb)
chrom_ranges_gff = gff_db[grepl("chromosome", gff_db$type)]
chrom_ranges_gff = chrom_ranges_gff[order(seqnames(chrom_ranges_gff))]
chromRanges <- GRanges(1:10, ranges(chrom_ranges_gff[1:10,]))

tx <- reduce(transcripts)
strand(tx) <- "*"
intergenic <- setdiff(chromRanges, tx) # 39784
# how many total basepairs are intergenic? = 1944445180! ~85% of maize genome (2.3e9 bp)
sum(width(intergenic))
strand(first_transcripts) = "*"
intergenic_bt = setdiff(chromRanges, first_transcripts)


# remove overlap among transcript info for each genomic category
remove_overlap = function(a,b,c,d,e,f,g) {
  a = setdiff(a,b)
  a = setdiff(a,c)
  a = setdiff(a,d)
  a = setdiff(a,e)
  a = setdiff(a,f)
  a = setdiff(a,g)
  return(a)
}

# add param for intergenic
remove_overlap_intergenic = function(a,b,c,d,e,f,g) {
  a = setdiff(a,b, ignore.strand = T)
  a = setdiff(a,c, ignore.strand = T)
  a = setdiff(a,d, ignore.strand = T)
  a = setdiff(a,e, ignore.strand = T)
  a = setdiff(a,f, ignore.strand = T)
  a = setdiff(a,g, ignore.strand = T)
  return(a)
}

temp_upstream_bt01 = remove_overlap(upstream_of_promoterbt_01, promoters_bt01, utr5_bt01, cds_bt01, 
                                    introns_bt01, utr3_bt01, intergenic_bt)

temp_promoters_bt01 = remove_overlap(promoters_bt01, utr5_bt01, cds_bt01, introns_bt01,
                                     utr3_bt01, intergenic_bt, upstream_of_promoterbt_01)

temp_utr5_bt01 = remove_overlap(utr5_bt01, cds_bt01, introns_bt01, utr3_bt01, 
                                intergenic_bt, upstream_of_promoterbt_01, promoters_bt01)

temp_cds_bt01 = remove_overlap(cds_bt01, introns_bt01, utr3_bt01, intergenic_bt, 
                               upstream_of_promoterbt_01, promoters_bt01, utr5_bt01)

temp_introns_bt01 = remove_overlap(introns_bt01, utr3_bt01, intergenic_bt, upstream_of_promoterbt_01, 
                                   promoters_bt01, utr5_bt01, cds_bt01)

temp_utr3_bt01 = remove_overlap(utr3_bt01, intergenic_bt, upstream_of_promoterbt_01, promoters_bt01, 
                                utr5_bt01, cds_bt01, introns_bt01)

temp_intergenic_bt = remove_overlap_intergenic(intergenic_bt, upstream_of_promoterbt_01, promoters_bt01, utr5_bt01, 
                                               cds_bt01, introns_bt01, utr3_bt01)

rm(upstream_of_promoterbt_01, promoters_bt01, utr5_bt01, cds_bt01, introns_bt01, utr3_bt01, intergenic_bt)
# reassign to variables to maintain consistency
upstream_of_promoterbt_01 = temp_upstream_bt01
promoters_bt01 = temp_promoters_bt01
utr5_bt01 = temp_utr5_bt01
cds_bt01 = temp_cds_bt01
introns_bt01 = temp_introns_bt01
utr3_bt01 = temp_utr3_bt01
intergenic_bt = temp_intergenic_bt

############################################
# efficacy measure
# NOW:
# examine intersections of transcript assignments...
# no overlap should exist
intersect(utr5_bt01, cds_bt01)
intersect(utr5_bt01, introns_bt01)
intersect(utr5_bt01, utr3_bt01)
intersect(utr5_bt01, upstream_of_promoterbt_01)
intersect(utr5_bt01, intergenic_bt)

intersect(cds_bt01, introns_bt01)
intersect(cds_bt01, utr3_bt01)
intersect(cds_bt01, promoters_bt01)
intersect(cds_bt01, upstream_of_promoterbt_01)
intersect(cds_bt01, intergenic_bt)

intersect(introns_bt01, utr3_bt01)
intersect(introns_bt01, promoters_bt01)
intersect(introns_bt01, upstream_of_promoterbt_01)
intersect(introns_bt01, intergenic_bt)

intersect(utr3_bt01, promoters_bt01)
intersect(utr3_bt01, upstream_of_promoterbt_01)
intersect(utr3_bt01, intergenic_bt)


# load bonnMu genotype data 
load("cleaned_bonnMu_rm_double.RData")
seqlevels(cleaned_bonnMu_rm_double) = str_replace_all(seqlevels(cleaned_bonnMu_rm_double), "Chr", "")
cleaned_bonnMu = cleaned_bonnMu[order(seqnames(cleaned_bonnMu_rm_double), ranges(cleaned_bonnMu_rm_double))]

# utr5 overlap with insertions in cleaned_bonnMu
utr5_overlap = findOverlaps(cleaned_bonnMu_rm_double, utr5_bt01, type = "within")
count_utr5 = sum(countOverlaps(cleaned_bonnMu_rm_double, utr5_bt01, type = "within")) # 21918
bp_utr5 = sum(width(utr5_bt01)) # 7367575 bp

insertion_vector = queryHits(findOverlaps(cleaned_bonnMu_rm_double, utr5_bt01, type="within"))
table(duplicated(ins_vector)) # 130 items are in the vector more than once!
ins_vector[duplicated(ins_vector)] # 21918-130 = 21788, which correlates with the number we get when assigning to cbm_cat


# cds overlap with insertions in cleaned_bonnM
cds_overlap = findOverlaps(cleaned_bonnMu_rm_double, cds_bt01, type = "within")
count_cds = sum(countOverlaps(cleaned_bonnMu_rm_double, cds_bt01, type = "within")) # 13683
bp_cds = sum(width(cds_bt01)) # 41374458 bp

#intron overlap with insertions in cleaned_bonnM
intron_overlap = findOverlaps(cleaned_bonnMu_rm_double, introns_bt01, type = "within")
count_introns = sum(countOverlaps(cleaned_bonnMu_rm_double, introns_bt01, type = "within")) # 11099
bp_introns = sum(width(introns_bt01)) # 87657969 bp

# utr3 overlap with insertions in cleaned_bonnM
utr3_overlap = findOverlaps(cleaned_bonnMu_rm_double, utr3_bt01, type = "within")
count_utr3 = sum(countOverlaps(cleaned_bonnMu_rm_double, utr3_bt01, type = "within")) # 1567
bp_utr3 = sum(width(utr3_bt01)) # 11038316 bp


intergenic_overlap = findOverlaps(cleaned_bonnMu_rm_double, intergenic_bt, type = "within")
count_intergenic = sum(countOverlaps(cleaned_bonnMu_rm_double, intergenic_bt, type = "within")) # 3561
bp_intergenic = sum(width(intergenic_bt)) # 1959932612 bp
#count_remaining_intergenic = count_intergenic

promoter_overlap = findOverlaps(cleaned_bonnMu_rm_double, promoters_bt01, type = "within")
count_promoter = sum(countOverlaps(cleaned_bonnMu_rm_double, promoters_bt01, type = "within")) # 3449
bp_promoter = sum(width(promoters_bt01))

upstream_overlap = findOverlaps(cleaned_bonnMu_rm_double, upstream_of_promoterbt_01, type = "within")
count_upstream = sum(countOverlaps(cleaned_bonnMu_rm_double, upstream_of_promoterbt_01, type = "within")) # 2847
bp_upstream = sum(width(upstream_of_promoterbt_01))

#insertioncount = t(data.frame(count_utr5,count_cds, count_introns, count_utr3, count_promoter, count_upstream, count_intergenic))
insertioncount = t(data.frame(count_promoter = count_promoter + count_upstream, count_utr5,count_cds, count_introns, count_utr3, count_intergenic))
save(insertioncount, file = "insertion_count.RData")

sum(insertioncount) # 58124 - this number is not accurate, as it includes the double counts explained in line 157
# the cleaned_bonnMu_rm_double only has 57247 insertions - this is reflected in the code below



## quick bp maths
# calculate proportions of where insertions are occurring
bp_total = sum(width(chromRanges))
prop_utr5 = bp_utr5/bp_total
prop_cds = bp_cds/bp_total
prop_intron = bp_introns/bp_total
prop_utr3 = bp_utr3/bp_total
prop_intergenic = bp_intergenic/bp_total
prop_promoter = (bp_promoter+bp_upstream)/bp_total

# go down to ibp for comment abt including prop_upstream20000

bp_proportions = t(data.frame(prop_promoter, prop_utr5, prop_cds, prop_intron, prop_utr3, prop_intergenic))
# bp_count = t(data.frame(bp_upsream2000bp = (bp_promoter+bp_upstream), bp_utr5, bp_cds, bp_introns, bp_utr3, bp_intergenic))
bp_count = t(data.frame(bp_promoter = bp_promoter+bp_upstream,bp_utr5, bp_cds, bp_introns, bp_utr3, bp_intergenic))
save(bp_count, file = "bp_count.RData")
save(bp_proportions, file = "bp_proportions.RData")

sum(bp_proportions)
sum(bp_utr5,bp_cds,bp_introns,bp_utr3,bp_intergenic) - bp_total #1032813 excess basepairs = `0.05%`> than 100%

upstream2000bp_ibp = (count_promoter+count_upstream)/(bp_promoter + bp_upstream)
promoter_ibp = upstream2000bp_ibp
utr5_ibp = count_utr5/bp_utr5
cds_ibp = count_cds/bp_cds
introns_ibp = count_introns/bp_introns
utr3_ibp = count_utr3/bp_utr3
intergenic_ibp = count_intergenic/bp_intergenic

# how to transform these small values
# don't include upstream2000bp_ibp bc not significant proportion of insertions in that area compared to utr5, makes point unclear to reviewers
insertionbp_prop = t(data.frame(promoter_ibp,utr5_ibp, cds_ibp, introns_ibp, utr3_ibp, intergenic_ibp))
ibp_1000 = 1000*insertionbp_prop

save(ibp_1000, file = "insertions_per_kb.RData")

bp_insertion_summary = cbind(bp_proportions, insertioncount, ibp_1000)
colnames(bp_insertion_summary) = c("proportion_of_genome", "insertion_count", "insertions_per_1kb")
rownames(bp_insertion_summary) = c("promoter", "utr5", "cds", "intron", "utr3", "intergenic")
save(bp_insertion_summary, file = "bp_insertion_summary.RData")

##################
load("cleaned_bonnMu_rm_double.RData")
# following code was sourced from insertion_categorization.R
# categorize each insertion/row in bonnMu dataset by the genomic region it inserts into
# cbm_cat = cleaned bonnMu categorized

insertion_category = matrix(nrow=length(cleaned_bonnMu_rm_double), ncol = 1)
cbm_cat_updated = cbind(data.frame(cleaned_bonnMu_rm_double), insertion_category)
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, utr5_bt01, type="within")),]$insertion_category = "utr5" # 21788
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, cds_bt01, type="within")),]$insertion_category = "cds"
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, introns_bt01, type="within")),]$insertion_category = "intron"
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, utr3_bt01, type="within")),]$insertion_category = "utr3"
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, promoters_bt01, type="within")),]$insertion_category = "promoter"
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, upstream_of_promoterbt_01, type="within")),]$insertion_category = "upstream_of_promoter"
cbm_cat_updated[queryHits(findOverlaps(cleaned_bonnMu_rm_double, intergenic_bt, type="within")),]$insertion_category = "intergenic"
cbm_cat_updated[is.na(cbm_cat_updated$insertion_category),] = "not_assigned"
cbm_cat_updated = cbm_cat_updated[order(cbm_cat_updated$Stock),]

table(cbm_cat_updated$insertion_category)
sum(table(cbm_cat_updated$insertion_category)) #57240 - correct! 57240 - 2185 (not assigned) = 55055

save(cbm_cat_updated, file="cbm_cat_updated.RData")




# ### Construction of cleaned_bonnMu_rm_double - remove double count insertions from bonnMu data so cbm can be constructed ##
# 
# utr5_utr3_insertionoverlap=intersect(queryHits(utr5_overlap),queryHits(utr3_overlap)) 
# utr5_cds_insertionoverlap=intersect(queryHits(utr5_overlap),queryHits(cds_overlap)) 
# utr5_introns_insertionoverlap=intersect(queryHits(utr5_overlap),queryHits(intron_overlap)) 
# utr5_upstream_insertionoverlap=intersect(queryHits(utr5_overlap),queryHits(upstream_overlap)) 
# utr5_promoters_insertionoverlap=intersect(queryHits(utr5_overlap),queryHits(promoter_overlap)) 
# utr5_doublecount_insertions = temp #2079 insertions are counted 2+ times
# 
# cds_intron_insertionoverlap=intersect(queryHits(cds_overlap),queryHits(intron_overlap)) 
# cds_utr3_insertionoverlap=intersect(queryHits(cds_overlap),queryHits(utr3_overlap)) 
# cds_upstream_insertionoverlap=intersect(queryHits(cds_overlap),queryHits(upstream_overlap)) 
# cds_promoters_insertionoverlap=intersect(queryHits(cds_overlap),queryHits(promoter_overlap)) 
# 
# intron_utr3_insertionoverlap=intersect(queryHits(intron_overlap),queryHits(utr3_overlap)) 
# intron_upstream_insertionoverlap=intersect(queryHits(intron_overlap),queryHits(upstream_overlap)) 
# intron_promoters_insertionoverlap=intersect(queryHits(intron_overlap),queryHits(promoter_overlap)) 
# 
# utr3_upstream_insertionoverlap=intersect(queryHits(utr3_overlap),queryHits(upstream_overlap)) 
# utr3_promoter_insertionoverlap=intersect(queryHits(utr3_overlap),queryHits(promoter_overlap)) 
# 
# promoter_upstream_overlap = intersect(queryHits(promoter_overlap), queryHits(upstream_overlap))
# 
# cds_doublecount_insertions = union(cds_intron_insertionoverlap, cds_utr3_insertionoverlap)
# cds_doublecount_insertions = union(cds_doublecount_insertions, cds_upstream_insertionoverlap)
# cds_doublecount_insertions = union(cds_doublecount_insertions, cds_promoters_insertionoverlap)
# 
# intron_doublecount_insertions = union(intron_utr3_insertionoverlap, intron_upstream_insertionoverlap)
# intron_doublecount_insertions = union(intron_doublecount_insertions, intron_promoters_insertionoverlap)
# 
# utr3_doublecount_insertions = union(utr3_promoter_insertionoverlap, utr3_upstream_insertionoverlap)
# 
# intergenic_promoter_insertionoverlap = intersect(queryHits(intergenic_overlap), queryHits(promoter_overlap))
# intergenic_upstream_insertionoverlap = intersect(queryHits(intergenic_overlap), queryHits(upstream_overlap))
# 
# total_doublecount_insertions = union(utr5_doublecount_insertions, cds_doublecount_insertions)
# total_doublecount_insertions = union(total_doublecount_insertions, intron_doublecount_insertions)
# total_doublecount_insertions = union(total_doublecount_insertions, utr3_doublecount_insertions)
# total_doublecount_insertions = union(total_doublecount_insertions, promoter_upstream_overlap)
# total_doublecount_insertions = union(total_doublecount_insertions, intergenic_promoter_insertionoverlap)
# total_doublecount_insertions = union(total_doublecount_insertions, intergenic_upstream_insertionoverlap)
# total_doublecount_insertions_vector = runValue(total_doublecount_insertions)
# save(total_doublecount_insertions, file="total_doublecount_insertions.RData")
# 
# # 2753 insertions were in 2+ categories, so they were removed
# cleaned_bonnMu_rm_double = cleaned_bonnMu[-total_doublecount_insertions_vector,]
# save(cleaned_bonnMu_rm_double, file="cleaned_bonnMu_rm_double.RData")



