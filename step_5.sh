#!/bin/bash
echo "Categorizing TE insertions by the gene region it inserted into"

cd /Users/gcampidilli/Documents/research/mu_tir/step_5

while read LINE
do

INBRED=$LINE
echo "Inbred is $INBRED"

echo "Running insertion categorization Rscript on $INBRED"
Rscript --vanilla insertion_categorization_terminal.R ${INBRED} generegion_summary/${INBRED}_gr_summary.csv insertion_dat_categorized/${INBRED}_dat_cat.csv generegion_summary/graphs/${INBRED}_per1kb.png

done < ./inbred_list.txt

echo "Combining all inbred generegion_summary dataframes "
Rscript --vanilla combine_summaries.R total_generegion_summary.csv
