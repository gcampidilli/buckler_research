# buckler_research

Buckler Lab for Maize Genetics and Diversity: https://www.maizegenetics.net/ <br>

My research questions: <br>
How are Mutator transposable elements (TEs) conserved across maize genotypes and across the andropogoneae clade? What is the impact of Mutator insertions on plant fitness? <br>

For each of ~270 maize inbreds, as well as active transposon line BonnMu:

| Step | Description|
| ---- | -----------|
| 1 | Map short-read bam files to Mutator TE sequence |
| 2 | Blast bam files against Mutator TE sequence and filter outputs|
| 3 | Construct TE insertion matrix |
| 4 | Identify what gene region each TE inserted into (ex. "Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic")|
| 5 | Expression Association Analysis |
| 6 | Phenotype Association Analysis |


#### Map short-read bam files to Mutator TE sequence
whileloop_map_to_mu_from_iRods.sh
  - Map the SR bams to the TIR sequence of Mutator
  - Sort and filter the output to produce \*.sortedfiltered.MuTIR.bam for each of the SR bams on iRods
combine_bams_to_fasta.R
  - Combine \*.sortedfiltered bams of the same genotype, as each inbred has multiple bams

