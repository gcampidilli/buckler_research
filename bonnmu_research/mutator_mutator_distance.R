# For each family, what is the minimum and average Mutator-Mutator distance in that individual? 
# Are insertion positions overdispersed in individuals, preventing deleterious chromosome breaks?

library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)

wd = "/Users/gcampidilli/Documents/research"
setwd(wd)

# go to FINAL_insertion_count.R for construction of cleaned_bonnMu_rm_double.RData
load("cleaned_bonnMu_rm_double.RData")
head(cleaned_bonnMu_rm_double)

# obtain minimum and average Mutator-Mutator distance by 
# subsetting GRanges and using distanceToNearest

# do this for x in 1:num_insertions times
# return vector with list of nearest_neighbor_distances
# from vector, calculate mean and min nearest_neighbor_distance for given family

nearest_distance = function(x, fam_ranges){
  nearest_range = distanceToNearest(fam_ranges,fam_ranges[x])
  distances = mcols(nearest_range)$distance
  nearest_neighbor_index = order(distances)[2]
  nearest_neighbor_distance = distances[nearest_neighbor_index]
  
  return(nearest_neighbor_distance)
}

fam_summary = data.frame(cleaned_bonnMu_rm_double) %>% group_by(Stock) %>% dplyr::summarize(count=n())

fam_list = fam_summary$Stock
num_fams = length(fam_list)
mu_distance_list = list()
mu_distance_mean = list()
mu_distance_min = list()

for(x in 1:num_fams) {
  fam = cleaned_bonnMu_rm_double[cleaned_bonnMu_rm_double$Stock == fam_list[x]]
  num_insertions = length(fam)
  current_fam_ranges = ranges(fam)
  
  distance_vector = unlist(lapply(1:num_insertions, nearest_distance, fam_ranges = current_fam_ranges))
  mu_distance_list[[x]] = distance_vector
  mu_distance_mean[[x]] = mean(distance_vector)
  mu_distance_min[[x]] = min(distance_vector)
}


mu_distance_summary = data.frame(family = fam_list, mean_mu_distance = unlist(mu_distance_mean), min_mu_distance = unlist(mu_distance_min))
save(mu_distance_summary, file = "mu_distance_summary.RData")








