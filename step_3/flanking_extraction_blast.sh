#!/bin/bash
echo "All bams and bais are in mu_aln/bam"

CPU=24 ## adjust as needed

while read LINE
do

INBRED=$LINE
echo "Inbred is ${INBRED}"

echo "Running Rscript, which extracts flank + TSD sequence for each Mu hit in ${INBRED}"
Rscript --vanilla extract_tsd_flank.R ${INBRED} ~/goodman_flanking_fa/${INBRED}_flanking.fa

echo "Blasting ${INBRED}_flanking_df.csv against B73 database"
blastn -db Zm-B73-REFERENCE-NAM-5.0.fa -query ~/goodman_flanking_fa/${INBRED}_flanking.fa -out ~/blast_out_flanking/${INBRED}_blast_out_flanking.txt -outfmt 6

done < ./inbred_list.txt
