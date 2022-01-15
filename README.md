# buckler_research

Buckler Lab for Maize Genetics and Diversity: https://www.maizegenetics.net/ <br>

My research questions: <br>
How are Mutator transposable elements (TEs) conserved across maize genotypes and across the andropogoneae clade? And how does climate impact accumulation and distribution of Mutator TEs across both taxa? <br>

For each of ~270 maize inbreds, as well as active transposon lines BonnMu and UniformMu:

| Step | Description|
| ---- | -----------|
| 1 | Map short-read bam files to Mutator TE sequence |
| 2 | Blast bam files against Mutator TE sequence |
| 3 | Filter blast outputs|
| 4 | Construct insertion matrix |
| 5 | Identify what gene region each TE inserted into (ex. "Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic")|
| 6 | Combine with phenotype data and visualize results | <br>

** NOTE: Only a small sample of code used as part of this project is included in this repository. All scripts will be uploaded after findings are published. <br>

step_1.sh includes combine_bams_to_fasta.R <br>
step_5.sh includes insertion_categorization_terminal.R and combine_summaries.R <br>
