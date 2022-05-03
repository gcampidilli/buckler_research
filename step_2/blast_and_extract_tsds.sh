#!/bin/bash
echo "All fastas are in ~/goodman_fastas"

CPU=24 ## adjust as needed

while read LINE
do

INBRED=$LINE
FASTA="${INBRED}_out.fa"

echo "Inbred is ${INBRED}"
echo "Fasta name is ${FASTA}"

echo "Blasting ${INBRED} fasta against mu group query, outputting blast results in ~/goodman_blast_out"
blastn -db mu_groups_norc.fa -query ~/goodman_fastas/${FASTA} -out ~/goodman_blast_out/${INBRED}_mu_blast_out.txt -outfmt 6

echo "Extracting TSDs from ${INBRED} blast results and reference bam df"
Rscript --vanilla ./extract_tsds.R ${INBRED} ~/goodman_insertion_df/${INBRED}_insertion_df.csv

echo "Filtering insertions from TSD output file such that we exclude all TSDs that aren't 9bp, keep TSD duplicates bc multiple times coverage"
Rscript --vanilla ./filter_tsds.R ${INBRED} ~/filtered_goodman_insertion_df/${INBRED}_insertion_df_9TSD.csv

done < ./inbred_list.txt
