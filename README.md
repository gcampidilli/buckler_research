# buckler_research

Buckler Lab for Maize Genetics and Diversity: https://www.maizegenetics.net/ <br>

My research questions: <br>
How are Mutator transposable elements (TEs) conserved across maize genotypes? What is the impact of Mutator insertions on plant fitness? <br>

As a primer, Mutator TE behavior was investigated in the active mutator line BonnMu. These files are in bonnmu_research.<br/><br/>
Then, for all ~270 maize inbreds in the Goodman Maize Diversity Panel the following steps were carried out:

| Step | Description|
| ---- | -----------|
| 1 | Map short-read bam files to Mutator TE sequence |
| 2 | Blast bam files against Mutator TE sequence and filter outputs|
| 3 | Construct TE insertion matrix |
| 4 | Identify gene region at location of each TE insertion (ex. "Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic")|
| 5 | Expression Association Analysis |
| 6 | Phenotype Association Analysis |


#### Map short-read bam files to Mutator TE sequence
1. whileloop_map_to_mu_from_iRods.sh
    - Map the SR bams to the TIR sequence of Mutator
    - Sort and filter the output to produce \*.sortedfiltered.MuTIR.bam for each of the SR bams on iRods
2. combine_bams_to_fasta.R
    - Combine \*.sortedfiltered bams of the same genotype, as each inbred has multiple bams
#### Blast bam files against Mutator TE sequence and filter outputs
1. blast_and_extract_tsds.sh
    - Shell script that applies R scripts to all ~270 inbred fastas
2. goodman_fastas.zip
3. extract_tsds.R
    - Retrieve target site duplication sequences for each insertion, this acts as the 'ID' of the insertion
4. filter_tsds.R
    - Filters insertions from TSD output file such that we exclude all TSDs that aren't 9bp, keep TSD duplicates bc multiple times coverage
5. filtered_goodman_insertion_df.zip folder
6. flanking_extraction_blast.sh
    - utilizes extract_tsd_blank.R
7. extract_tsd_blank.R
    - extract bp flanking TSDs
8. blast_out_flanking_filtered.zip folder
#### Construct TE insertion matrix
1. insertion_matrix_pt1.sh
    - Utilizes prematrix_df_setup.R and combine_prematrix_df.R
2. prematrix_df_setup.R
    - Creates precurser to insertion matrix for each inbred
3. combine_prematrix_df.R
    - Combines  all inbred prematrices into one total prematix, from which the insertion matrix is constructed
4. pre_insertion_matricies folder
    - Has preinsertion matricies for the exact insertion locations (within 18bp of eachother), insertion locations rounded to the nearest hundred bp, and insertion locations rounded to nearest thousand bp
5. insertion_matrix_pt2.R
    - Construct insertion matricies given the pre-insertion matricies
6. insertion_matricies folder
    - Has insertion matricies (V5)for the exact insertion locations (within 18bp of eachother), insertion locations rounded to the nearest hundred bp, and insertion locations rounded to nearest thousand bp
    - V4 exact insertion matrix
7. convert_to_v4.R
    -  convert insertion matrix to v4 gene list from v5
8. ins_mat_to_hapmap.R
    - convert insertion matrix from csv to hapmap format
#### Identify gene region at location of each TE insertion
1. insertion_categorization.sh
    - Utilizes assign_insertion_categories.R and combine_summaries.R
2. assign_insertion_categories.R
    - Use B73 reference genome to categorize which gene region each insertion inserted into (ex. "Promoter","UTR5", "Exon", "Intron", "UTR3", "Intergenic")
3. combine_summaries.R
    - Create summary total_generegion_summary.csv
4. total_generegion_summary.csv
    - For each inbred, has statistics about # of insertion and proportion of insertions in each gene region category
5. insertion_category_graphs_by_inbred folder
    - Barplot of proportion of insertions in each gene region for each inbred

