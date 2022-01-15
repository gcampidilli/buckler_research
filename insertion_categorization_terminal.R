###### INPUT: 1) blast output for given inbred against Mutator TE
######        2) B73 reference genome by gene region 
######           Possible gene regions: "Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic", "Not assigned"

##### OUTPUT: 1) data frame with each insertion categorized by the gene region it inserted into

##### Script outline:
#     PRE-PART 1 (BOTTOM OF SCRIPT): DIVIDE REFERENCE GENOME INTO REFERENCE GENOME BY GENIC CATEGORY (INTRON, UTR5, CDS, PROMOTER, ETC)
#     PART 1: LOAD PACKAGES + DATA, FILTER, AND FORMAT
#     PART 2: IDENTIFY AND REMOVE INSERTIONS THAT EITHER 1) MAP TO 2+ TRANSCRIPTS OF THE SAME GENE REGION CATEGORIES 2) MAP TO 2+ TRANSCRIPTS OF DIFFERENT GENE REGION CATEGORIES
#     PART 3: CONSTRUCT INSERTION SUMMARY BY GENE REGION
#     PART 4: ASSIGN EACH INSERTION A GENE REGION CATEGORY (INTRON, UTR5, CDS, PROMOTER, ETC)
#     PART 5: VISUALIZE GENE REGION INSERTION SUMMARY

########################################################################################
#     PART 1: LOAD PACKAGES + DATA, FILTER, AND FORMAT
args = commandArgs(trailingOnly=TRUE)
inbred_name = args[1]

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

wd = "/Users/gcampidilli/Documents/research/mu_tir/step_5"
setwd(wd)

#     LOAD INBRED INSERTION DATA AND TRANSFORM INTO GRANGE OBJECT
inbred_dat = read.csv(paste("../step_4/blast_out_flanking_filtered/",inbred_name, "_flanking_filtered.csv", sep = ""), row.names = 1)
colnames(inbred_dat)[c(2,9,10)] = c("chr", "start", "end")
#inbred_dat$chr= str_replace_all(inbred_dat$chr, "chr", "")

inbred_dat$strand = NA
#     ASSIGN STRAND
inbred_dat$strand = unlist(lapply(1:nrow(inbred_dat), function(x){
  if(inbred_dat$start[x]<inbred_dat$end[x]){inbred_dat$strand = "+"} else{inbred_dat$strand = "-"}}))
head(inbred_dat)
#     SWITCH START AND END FOR - STRAND TO GET IN CORRECT FORMAT FOR GRANGE OBJECT
neg_strand_index = which(inbred_dat$strand=="-")
neg_strand_temp = inbred_dat[neg_strand_index,]
inbred_dat[neg_strand_index,]$start = neg_strand_temp$end
inbred_dat[neg_strand_index,]$end = neg_strand_temp$start

inbred_dat_grange = makeGRangesFromDataFrame(inbred_dat, keep.extra.columns = T)
head(inbred_dat_grange)

#     LOAD B73 REFERENCE GENOME BY GENE REGION (PROMOTER, UTR5, CDS, INTRON, ETC)
#     CONTRUCTION OF THIS DATA IS AT THE BOTTOM OF THE SCRIPT (PRE-PART 1)

load("/Users/gcampidilli/Documents/research/mu_tir/step_6/upstream_of_promoter_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/promoter_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/utr5_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/cds_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/intron_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/utr3_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/intergenic_ref.RData")
load("/Users/gcampidilli/Documents/research/mu_tir/step_6/chromRanges.RData")

print("PART 1 DONE")
########################################################################################
#     PART 2: IDENTIFY AND REMOVE INSERTIONS THAT EITHER 
#             1) MAP TO 2+ TRANSCRIPTS OF THE SAME GENE REGION CATEGORIES
#             2) MAP TO 2+ TRANSCRIPTS OF DIFFERENT GENE REGION CATEGORIES

# utr5 overlap with insertions in cleaned_bonnMu
utr5_overlap = findOverlaps(inbred_dat_grange, utr5_ref, type = "within")
cds_overlap = findOverlaps(inbred_dat_grange, cds_ref, type = "within")
intron_overlap = findOverlaps(inbred_dat_grange, intron_ref, type = "within")
utr3_overlap = findOverlaps(inbred_dat_grange, utr3_ref, type = "within")
promoter_overlap = findOverlaps(inbred_dat_grange, promoter_ref, type = "within")
upstream_overlap = findOverlaps(inbred_dat_grange, upstream_of_promoter_ref, type = "within")
intergenic_overlap = findOverlaps(inbred_dat_grange, intergenic_ref, type = "within")

remove_same_region_overlap = function(region_overlap){
  remove_vector = vector()
  region_overlap_index = findOverlaps(inbred_dat_grange[queryHits(region_overlap),], type = "any", select = "first")
  numeric_vector = as.numeric(table(region_overlap_index))
  overlap_vector = which(numeric_vector>1)
  for(x in 1:length(overlap_vector)){
    duplicates = which(region_overlap_index == overlap_vector[x])
    duplicates_row_index = queryHits(region_overlap)[duplicates]
    if(length(unique(inbred_dat_grange$percent_identity[duplicates_row_index])) > 1) { 
      temp = which(inbred_dat_grange$percent_identity[duplicates_row_index] != max(inbred_dat_grange$percent_identity[duplicates_row_index]))
      remove = queryHits(region_overlap)[duplicates[temp]]
    } else{
      temp = which(inbred_dat_grange$alignment_length[duplicates_row_index] != max(inbred_dat_grange$alignment_length[duplicates_row_index]))
      remove = queryHits(region_overlap)[duplicates[temp]]
    }
    remove_vector = append(remove_vector, remove)
  }
  remove_vector
}

same_region_remove_vector = vector()
same_region_remove_vector = append(same_region_remove_vector, if(length(utr5_overlap) != 0){remove_same_region_overlap(utr5_overlap)})
same_region_remove_vector = append(same_region_remove_vector,if(length(cds_overlap) != 0){remove_same_region_overlap(cds_overlap)})
same_region_remove_vector = append(same_region_remove_vector,if(length(intron_overlap) != 0){remove_same_region_overlap(intron_overlap)})
same_region_remove_vector = append(same_region_remove_vector,if(length(utr3_overlap) != 0){remove_same_region_overlap(utr3_overlap)})
same_region_remove_vector = append(same_region_remove_vector,if(length(promoter_overlap) != 0){remove_same_region_overlap(promoter_overlap)})
same_region_remove_vector = append(same_region_remove_vector,if(length(upstream_overlap) != 0){remove_same_region_overlap(upstream_overlap)})
same_region_remove_vector = append(same_region_remove_vector,if(length(intergenic_overlap) != 0){remove_same_region_overlap(intergenic_overlap)})

all_region_query = vector()
all_region_query = append(all_region_query, queryHits(utr5_overlap))
all_region_query = append(all_region_query, queryHits(cds_overlap))
all_region_query = append(all_region_query, queryHits(intron_overlap))
all_region_query = append(all_region_query, queryHits(promoter_overlap))
all_region_query = append(all_region_query, queryHits(upstream_overlap))
all_region_query = append(all_region_query, queryHits(intergenic_overlap))

different_region_remove_vector = unique(all_region_query[duplicated(all_region_query)])

total_remove_vector = unique(union(same_region_remove_vector, different_region_remove_vector))
# save(total_doublecount_insertions, file="total_doublecount_insertions.RData")

if(length(total_remove_vector) != 0){ inbred_dat_grange_nd = inbred_dat_grange[-(total_remove_vector),]
} else { inbred_dat_grange_nd = inbred_dat_grange}
# save(inbred_dat_grange_nd, file="inbred_dat_grange_nd.RData")

print("PART 2 DONE")
########################################################################################
#     PART 3: CONSTRUCT INSERTION SUMMARY BY GENE REGION

# load("/Users/gcampidilli/Documents/research/inbred_dat_grange_nd.RData") # blast output of TE insertions for each inbred

# utr5 insertion count + total ref bp
count_utr5 = sum(countOverlaps(inbred_dat_grange_nd, utr5_ref, type = "within")) # 21918
bp_utr5 = sum(width(utr5_ref)) # 7367575 bp

# utr5 insertion count + total ref bp
count_cds = sum(countOverlaps(inbred_dat_grange_nd, cds_ref, type = "within")) # 13683
bp_cds = sum(width(cds_ref)) # 41374458 bp

# utr5 insertion count + total ref bp
count_introns = sum(countOverlaps(inbred_dat_grange_nd, intron_ref, type = "within")) # 11099
bp_introns = sum(width(intron_ref)) # 87657969 bp

# utr5 insertion count + total ref bp
count_utr3 = sum(countOverlaps(inbred_dat_grange_nd, utr3_ref, type = "within")) # 1567
bp_utr3 = sum(width(utr3_ref)) # 11038316 bp

# utr5 insertion count + total ref bp
count_promoter = sum(countOverlaps(inbred_dat_grange_nd, promoter_ref, type = "within")) # 3449
bp_promoter = sum(width(promoter_ref))

# utr5 insertion count + total ref bp
count_upstream = sum(countOverlaps(inbred_dat_grange_nd, upstream_of_promoter_ref, type = "within")) # 2847
bp_upstream = sum(width(upstream_of_promoter_ref))

# utr5 insertion count + total ref bp
count_intergenic = sum(countOverlaps(inbred_dat_grange_nd, intergenic_ref, type = "within")) # 3561
bp_intergenic = sum(width(intergenic_ref)) # 1959932612 bp
#count_remaining_intergenic = count_intergenic

total_count = count_promoter + count_upstream + count_utr5+count_cds+ count_introns+ count_utr3+ count_intergenic
count_notassigned = length(inbred_dat_grange_nd) - total_count

insertion_count = t(data.frame(count_promoter = count_promoter + count_upstream, count_utr5,count_cds, count_introns, count_utr3, count_intergenic, count_notassigned))
# save(insertioncount, file = "insertion_count.RData")

# calculate basepair proportions by gene region category
bp_total = sum(width(chromRanges))
prop_utr5 = bp_utr5/bp_total
prop_cds = bp_cds/bp_total
prop_intron = bp_introns/bp_total
prop_utr3 = bp_utr3/bp_total
prop_intergenic = bp_intergenic/bp_total
prop_promoter = (bp_promoter+bp_upstream)/bp_total

bp_proportions = t(data.frame(prop_promoter, prop_utr5, prop_cds, prop_intron, prop_utr3, prop_intergenic, NA))
bp_count = t(data.frame(bp_promoter = bp_promoter+bp_upstream,bp_utr5, bp_cds, bp_introns, bp_utr3, bp_intergenic, NA))

# sum(bp_proportions)
# sum(bp_utr5,bp_cds,bp_introns,bp_utr3,bp_intergenic) - bp_total #1032813 excess basepairs = `0.05%`> than 100%

# calculate ratio of insertions per basepair by gene region category
promoter_ibp = (count_promoter+count_upstream)/(bp_promoter + bp_upstream)
utr5_ibp = count_utr5/bp_utr5
cds_ibp = count_cds/bp_cds
introns_ibp = count_introns/bp_introns
utr3_ibp = count_utr3/bp_utr3
intergenic_ibp = count_intergenic/bp_intergenic

insertionbp_prop = t(data.frame(promoter_ibp, utr5_ibp, cds_ibp, introns_ibp, utr3_ibp, intergenic_ibp, NA))

# transform to insertions per 1000 bp 
ibp_1000 = 1000*insertionbp_prop

inbred_name_vector = data.frame(rep(inbred_name, 7))
region_vector = c("Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic", "Not assigned")
insertion_proportion =  round(insertion_count/total_count, digits = 3)
insertion_proportion

insertion_summary_by_generegion = cbind(inbred_name_vector, region_vector, bp_count, bp_proportions, insertion_count, insertion_proportion, ibp_1000)
colnames(insertion_summary_by_generegion) = c("inbred_name_vector", "region", "bp_count","bp_proportion", "insertion_count", "insertion_proportion", "insertions_per_1kb")
insertion_summary_by_generegion$bp_proportion = round(insertion_summary_by_generegion$bp_proportion, digits = 3)
rownames(insertion_summary_by_generegion) = NULL

head(insertion_summary_by_generegion)
write.csv(insertion_summary_by_generegion, file = args[2])

print("PART 3 DONE")
########################################################################################
#     PART 4: ASSIGN EACH INSERTION A GENE REGION CATEGORY (INTRON, UTR5, CDS, PROMOTER, ETC)

# load("/Users/gcampidilli/Documents/research/inbred_dat_grange_nd.RData") # blast output of TE insertions for each inbred

insertion_category = matrix(nrow=length(inbred_dat_grange_nd), ncol = 1)
inbred_dat_cat = cbind(data.frame(inbred_dat_grange_nd), insertion_category)
print("test")
if(count_utr5 != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, utr5_ref, type="within")),]$insertion_category = "utr5"}
if(count_cds != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, cds_ref, type="within")),]$insertion_category = "cds"}
if(count_introns != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, intron_ref, type="within")),]$insertion_category = "intron"}
if(count_utr3 != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, utr3_ref, type="within")),]$insertion_category = "utr3"}
if(count_promoter != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, promoter_ref, type="within")),]$insertion_category = "promoter"}
if(count_upstream != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, upstream_of_promoter_ref, type="within")),]$insertion_category = "upstream_of_promoter"}
if(count_intergenic != 0){inbred_dat_cat[queryHits(findOverlaps(inbred_dat_grange_nd, intergenic_ref, type="within")),]$insertion_category = "intergenic"}
inbred_dat_cat[is.na(inbred_dat_cat$insertion_category),]$insertion_category = "not_assigned"
inbred_dat_cat = inbred_dat_cat[order(inbred_dat_cat$seqnames, inbred_dat_cat$start),]
print("test2")

# table(inbred_dat_cat$insertion_category)
# sum(table(inbred_dat_cat$insertion_category)) 
tail(inbred_dat_cat) # here we see that the last 2 insertions are likely the same insertion, so devise additional filtering

write.csv(inbred_dat_cat, file = args[3])

print("PART 4 DONE")
##########################################################################################
#     PART 5: VISUALIZE GENE REGION INSERTION SUMMARY

# load("insertion_summary_by_generegion.RData")
insertion_summary_by_generegion$region = factor(insertion_summary_by_generegion$region, levels = c("Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic", "Not Assigned"))

per1kb_barplot =
  ggplot(data = insertion_summary_by_generegion[1:6,], aes(x = region, y = insertions_per_1kb)) +
  geom_bar(stat = "identity", fill = "black") +
  labs(title = paste(inbred_name,": Frequency of Mutator insertions in each gene region", sep = "")) + ylab("# insertions per 1000 base pairs") + xlab("") +
  theme_gray(base_size = 14) + geom_text(aes(label=round(insertions_per_1kb, digits = 8)), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() + theme(axis.title.y = element_text(size =15, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 15, face = "bold", color = "black")) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

png(args[4], width=800, height=600)
per1kb_barplot
dev.off()

bp_proportion_barplot =
  ggplot(data = insertion_summary_by_generegion[1:6,], aes(x = region, y = bp_proportion)) +
  geom_bar(stat = "identity", fill = "black") +
  labs(title = "Proportion of genome in each gene region") + ylab("% of genome") + xlab("") +
  theme_gray(base_size = 14) + geom_text(aes(label=bp_proportion), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() + theme(axis.title.y = element_text(size =12, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 12, face = "bold", color = "black")) + theme(axis.text.y = element_text(size = 15)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))

print("PART 5 DONE")







# ########################################################################################
# #     PRE: PART 1: DIVIDE REFERENCE GENOME INTO REFERENCE GENOME BY GENIC CATEGORY (INTRON, UTR5, CDS, PROMOTER, ETC)
# Load maize GFF data into txdb object
# Transcript database
# txdb <- makeTxDbFromGFF("../Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", format="gff3")
# gff_db =import.gff3("../Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz", format = "gff3")
# 
# # refine txdb object so it only includes chromosomes 1:10
# seqlevels(txdb) = c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")
# 
# # transcripts = set of exons, a gene is composed of alternatively spliced transcripts
# transcripts <- transcripts(txdb)
# first_transcripts = transcripts[grepl("T001", transcripts$tx_name)] # 40,923
# 
# # use.names = T, group names are transcript names "tx_name" (i.e. Zm00001d027235_T001),
# # if use.names is false, group names are "tx_id" (i.e. 67598)
# introns_bt = unlist(intronsByTranscript(txdb, use.names = T)) #1mil+ introns
# introns_ref = introns_bt[grepl("T001", names(introns_bt))] #152472 introns
# 
# utr5_bt = unlist(fiveUTRsByTranscript(txdb, use.names = T)) #248164 utr5
# utr5_ref = utr5_bt[grepl("T001", names(utr5_bt))] # 35228 utr5
# 
# utr3_bt = unlist(threeUTRsByTranscript(txdb, use.names = T)) # 251982 utr3
# utr3_ref = utr3_bt[grepl("T001", names(utr3_bt))] # 33721 utr3
# 
# cds_bt = unlist(cdsBy(txdb, use.names = T)) # 918634 cds
# cds_ref = cds_bt[grepl("T001", names(cds_bt))] # 170269 cds
# 
# # can make analysis more specific by using the promoter region -500 to 0 bp and the region upstream of promoter -2000 to -501 bp
# # as a complement to intergenic category
# 
# # promoter is at correct location at + and - strand
# promoters_bt = promoters(txdb, upstream=500, downstream = 0, use.names=TRUE) #133940 promoter regions
# promoters_ref = promoters_bt[grepl("T001", names(promoters_bt))] #41277 promoter regions
# length(reduce(promoters_ref)) # 40979
# 
# upstream_of_promoterbt = promoters(txdb, upstream = 2000, downstream = 0, use.names = TRUE)
# upstream_of_promoter_ref = upstream_of_promoterbt[grepl("T001", names(upstream_of_promoterbt))]
# # if on positive:
# end(ranges(upstream_of_promoter_ref[strand(upstream_of_promoter_ref) == "+"])) =
#   end(ranges(upstream_of_promoter_ref[strand(upstream_of_promoter_ref) == "+"])) - 500
# #end(ranges(upstream_of_promoter_ref)) = end(ranges(upstream_of_promoter_ref))-500
# # if on negative:
# start(ranges(upstream_of_promoter_ref[strand(upstream_of_promoter_ref) == "-"])) =
#   start(ranges(upstream_of_promoter_ref[strand(upstream_of_promoter_ref) == "-"])) + 500
# 
# table(width(ranges(upstream_of_promoter_ref)))
# GenomicRanges::intersect(promoters_ref, upstream_of_promoter_ref) # 859 GenomicRanges::intersections
# 
# chromInfo <- seqinfo(txdb)
# chrom_ranges_gff = gff_db[grepl("chromosome", gff_db$type)]
# chrom_ranges_gff = chrom_ranges_gff[order(seqnames(chrom_ranges_gff))]
# chromRanges <- GRanges(c("chr1", "chr2", "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"),
#                        ranges(chrom_ranges_gff[1:10,]))
# 
# tx <- reduce(transcripts)
# strand(tx) <- "*"
# 
# intergenic <- GenomicRanges::setdiff(chromRanges, tx) # 39784
# # how many total basepairs are intergenic? = 1944445180! ~85% of maize genome (2.3e9 bp)
# sum(width(intergenic))
# strand(first_transcripts) = "*"
# intergenic_bt = GenomicRanges::setdiff(chromRanges, first_transcripts)
# 
# 
# # remove overlap among transcript info for each genomic category
# remove_overlap = function(a,b,c,d,e,f,g) {
#   a = GenomicRanges::setdiff(a,b)
#   a = GenomicRanges::setdiff(a,c)
#   a = GenomicRanges::setdiff(a,d)
#   a = GenomicRanges::setdiff(a,e)
#   a = GenomicRanges::setdiff(a,f)
#   a = GenomicRanges::setdiff(a,g)
#   return(a)
# }
# 
# # add param for intergenic
# remove_overlap_intergenic = function(a,b,c,d,e,f,g) {
#   a = GenomicRanges::setdiff(a,b, ignore.strand = T)
#   a = GenomicRanges::setdiff(a,c, ignore.strand = T)
#   a = GenomicRanges::setdiff(a,d, ignore.strand = T)
#   a = GenomicRanges::setdiff(a,e, ignore.strand = T)
#   a = GenomicRanges::setdiff(a,f, ignore.strand = T)
#   a = GenomicRanges::setdiff(a,g, ignore.strand = T)
#   return(a)
# }
# 
# temp_upstream_ref = remove_overlap(upstream_of_promoter_ref, promoters_ref, utr5_ref, cds_ref,
#                                     introns_ref, utr3_ref, intergenic_bt)
# 
# temp_promoters_ref = remove_overlap(promoters_ref, utr5_ref, cds_ref, introns_ref,
#                                      utr3_ref, intergenic_bt, upstream_of_promoter_ref)
# 
# temp_utr5_ref = remove_overlap(utr5_ref, cds_ref, introns_ref, utr3_ref,
#                                 intergenic_bt, upstream_of_promoter_ref, promoters_ref)
# 
# temp_cds_ref = remove_overlap(cds_ref, introns_ref, utr3_ref, intergenic_bt,
#                                upstream_of_promoter_ref, promoters_ref, utr5_ref)
# 
# temp_introns_ref = remove_overlap(introns_ref, utr3_ref, intergenic_bt, upstream_of_promoter_ref,
#                                    promoters_ref, utr5_ref, cds_ref)
# 
# temp_utr3_ref = remove_overlap(utr3_ref, intergenic_bt, upstream_of_promoter_ref, promoters_ref,
#                                 utr5_ref, cds_ref, introns_ref)
# 
# temp_intergenic_bt = remove_overlap_intergenic(intergenic_bt, upstream_of_promoter_ref, promoters_ref, utr5_ref,
#                                                cds_ref, introns_ref, utr3_ref)
# 
# rm(upstream_of_promoter_ref, promoters_ref, utr5_ref, cds_ref, introns_ref, utr3_ref, intergenic_bt)
# # reassign to variables to maintain consistency
# upstream_of_promoter_ref = temp_upstream_ref
# promoter_ref = temp_promoters_ref
# utr5_ref = temp_utr5_ref
# cds_ref = temp_cds_ref
# intron_ref = temp_introns_ref
# utr3_ref = temp_utr3_ref
# intergenic_ref = temp_intergenic_bt
# 
# # FINAL: efficacy measure
# # examine GenomicRanges::intersections of transcript assignments...
# # no overlap should exist
# GenomicRanges::intersect(utr5_ref, cds_ref)
# GenomicRanges::intersect(utr5_ref, intron_ref)
# GenomicRanges::intersect(utr5_ref, utr3_ref)
# GenomicRanges::intersect(utr5_ref, upstream_of_promoter_ref)
# GenomicRanges::intersect(utr5_ref, intergenic_ref)
# 
# GenomicRanges::intersect(cds_ref, intron_ref)
# GenomicRanges::intersect(cds_ref, utr3_ref)
# GenomicRanges::intersect(cds_ref, promoters_ref)
# GenomicRanges::intersect(cds_ref, upstream_of_promoter_ref)
# GenomicRanges::intersect(cds_ref, intergenic_ref)
# 
# GenomicRanges::intersect(intron_ref, utr3_ref)
# GenomicRanges::intersect(intron_ref, promoter_ref)
# GenomicRanges::intersect(intron_ref, upstream_of_promoter_ref)
# GenomicRanges::intersect(intron_ref, intergenic_ref)
# 
# GenomicRanges::intersect(utr3_ref, promoters_ref)
# GenomicRanges::intersect(utr3_ref, upstream_of_promoter_ref)
# GenomicRanges::intersect(utr3_ref, intergenic_ref)
# 
# 
# save(upstream_of_promoter_ref, file="upstream_of_promoter_ref.RData")
# save(promoter_ref, file="promoter_ref.RData")
# save(utr5_ref, file="utr5_ref.RData")
# save(cds_ref, file="cds_ref.RData")
# save(intron_ref, file="intron_ref.RData")
# save(utr3_ref, file="utr3_ref.RData")
# save(intergenic_ref, file="intergenic_ref.RData")
# save(chromRanges, file="chromRanges.RData")


