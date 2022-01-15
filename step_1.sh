#!/bin/bash

#
CPU=24 ## adjust as needed

echo "Downloading bam with iRods from the cloud"
while read LINE
do
echo "Directory is ${LINE}"
echo "Downloading bam with iRods"

iget $LINE

ORIGBAM=$(basename $LINE)
ORIGSAMPLE=$(basename $ORIGBAM .bam) ## get just the sample name without bam suffix

echo "File is ${ORIGBAM}"
echo "Sample is ${ORIGSAMPLE}"

## convert bam to fastq, map with bwa mem
echo "Converting to fastq, then mapping, filtering, and sorting all in one step"

samtools bamshuf -uO $ORIGBAM | \
samtools bam2fq - | \
bwa mem -t $CPU -p mu_groups_norc.fa - | \
samtools view -Sbu -F 260 - | \
samtools sort -@ $CPU - > ${ORIGSAMPLE}.sortedfiltered.MuTIR.bam

samtools index ${ORIGSAMPLE}.sortedfiltered.MuTIR.bam

echo "Copying files to homedir and removing original bam"
scp ${ORIGSAMPLE}.sortedfiltered.MuTIR.bam* ~/mu_aln/bams
rm $ORIGBAM

done < ~/inbred_list.txt

# all inbreds have more than 1 bam files due to multiple coverage short-reads
# combine bam files of the same inbred, and make format more usable
while read LINE
do

INBRED=$LINE
NUMBAMS=$(find . -name "*${INBRED}*.bam" | wc -l)
BAMARRAY=$(find . -name "*${INBRED}*.bam")

echo "Inbred is ${INBRED}"
echo "There are ${NUMBAMS} bam files for ${INBRED}, they are:"
echo "${BAMARRAY}"

echo "Combining bams to df"
Rscript --vanilla ~/combine_bams_to_fasta.R ${INBRED}  ~/goodman_df/${INBRED}_bam_df.csv ${BAMARRAY}

done < ~/inbred_list.txt
