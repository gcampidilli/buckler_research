#!/bin/bash
echo "Filtering and combining TE flanking sequence data for all inbreds."
echo "Constructing prematrix dataframe from which insertion matrix will be constructed"
echo "All flanking blast outputs are in blast_out_flanking/"

cd /Users/gcampidilli/Documents/research/mu_tir/step_4

CPU=4 ## adjust as needed

while read LINE
do

INBRED=$LINE
echo "Inbred is ${INBRED}"

echo "Running filtering Rscript on ${INBRED}"
Rscript --vanilla prematrix_df_setup.R ${INBRED} blast_out_flanking_filtered/${INBRED}_flanking_filtered.csv prematrix_df/exact/${INBRED}_prematrix_df_exact.csv prematrix_df/hundred/${INBRED}_prematrix_df_hundred.csv prematrix_df/thousand/${INBRED}_prematrix_df_thousand.csv

done < ./inbred_list.txt

echo "Combining all inbred prematrix dataframes "
Rscript --vanilla combine_prematrix_df.R total_prematrix_df.csv
