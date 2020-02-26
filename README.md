# Gene-Expression-Dataset

In this project, I analyzed  a real gene expression dataset. This was downloaded from the Gene Expression Omnibus (GEO, NCBI). In this link you’ll see the description of the experiment and the data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77466

There are a couple input files: 

1) The dataset GSE77466_matrix_table_raw_counts.tsv.gz. is renamed as gene_counts.tsv. The first column is the gene (represented by its ENSEMBL ID). The other columns have the reads mapped to each gene, for each corresponding sample. 

2) samples_data.tsv : A tab delimited file, which contains some details about the samples, extracted from the description of the experiment.
This contains the following columns:
- id: the name of the sample
- accession: the accession code of the sample
- type: the type of sample, being either L (lesional derived keratinocytes) or NL (non lesional derived
keratinocytes).

3) results_L-vs-NL_full.csv : A comma-delimited file, which contains the results of the differential expression analysis of the two conditions under study (L vs NL) using DEseq2. The following three columns are needed for this exercise:
- gene_ID: the gene ENSEMBL ID
- log2FoldChange: the log2 of the Fold Change between conditions
- padj: FDR corrected p-value

Each python script analyzes or visualizes the data 
1) normalize_read_counts.py: normalizes the read counts per gene by dividing the read counts to each library size
(library size = total counts in that sample).
2) heatmap.py:  using the normalized counts, make a Heatmap of the samples. The package seaborn has been installed for data visualization 
3) annotation.py: annotates the table of results with the HUGO gene name and the gene description. PANDAS is used here.
4) volcano_plot.py: makes a volcano plot. A line is dded at FDR = 0.05 and Log2(fold change) = 2. Good candidates genes are identifed and labeled with their gene names in the plot. Seaborn is used here again for data visualization. 
5) scatter_and_violin_plots.py: makes a scatter plots of the candidate genes (sample type vs –log2 of normalized gene counts)
and violin plots of the candidate genes (sample type vs –log2 of normalized gene counts). Seaborn is used here. 
6) test_client.py: simply runs each and every single function to get the outputs

Here were my observations: 

Based on the scatter plots and violin plots I can tell that for genes in the upper left hand side of the volcano plot (MCAM, KCNIP3, KRT75, CPA4, LMX1B), their normalized gene counts tend to be much higher in the NL sample (non lesional derived keratinocytes). But for genes in the upper right hand side (ADAM8, CDTI, RINL, TFPI2) their normalized gene counts tend to be much higher in the L sample (lesional derived keratinocytes). This may indicate that genes MCAM, KCNIP3, KRT75, CPA4, LMX1B are significant for non-lesional derived keratinocytes while genes, while genes ADAM8, CDTI, RINL, TFPI2 are more significant for lesional derived keratinocytes
 
I chose these candidate genes because they are the outliers on this graph (located on the upper left and right hand side) and thus they represent the most highly differentially expressed genes. In volcano plots, we can visualize genes with large fold changes that are also statistically significant. The genes in the upper left and right hand side show the largest and statistically most significant changes in expression.
