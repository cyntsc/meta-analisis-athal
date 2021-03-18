# Files are:
<br>

**0_alignments_statis_all.ipynb** <br>
Integrate and build the basic stats for the alignments. 

**1_log2_transformation_massive_ctrl.ipynb** <br>
First, A.thaliana healthy quantification files are taken to transform the measure values to Log2+1.<br>
Second, get the individual descriptive stats by file (log2) and calculate the genes in ceros. <br>
Third, descriptive stats by each sample is integrated in a unique, adding a new column with the number of genes in ceros by sample. <br>

**1_log2_transformation_massive_infect.ipynb**<br>
Does the same that **1_log2_transformation_massive_ctrl.ipynb**, but with the A.thaliana infected quantification files.

**3_merge_dataframes_ctrl.ipynb**<br>
1) Merge the x dataframes in a single matrix<br>
2) Some basic statistics are built with the data matrix.  <br>
3) Zeros across all samples are dropped.  <br>

**3_merge_dataframes_infect.ipynb**<br>
Does the same that **3_merge_dataframes_ctrl.ipynb**, but with the A.thaliana healthy quantification files.

**8_df_correlations.ipynb** <br>
1) Correlate the expression values of RNASeq samples provided in a df <br>
2) Plot the resulted matrix in a heatmap <br>
3) Provided resources to extract data and plot specific samples of interest<br>

**Gene_length_extraction_from_GTF** 
Create a file with ID-Genes and lenghts from a GTF file <br>
For those concerned with converting raw counts from an expression array to TMP counts for any purpose. These scripts are an adaptation from the https://www.reneshbedre.com/blog/expression_units.html who uses https://github.com/reneshbedre/bioinfokit library. 
Note. Gene lengths are based on coding-gene features (CDS), if you wish a more complex approach like average exons or some-like, I encourage to do the necessary changes. 

**Transform_rawcounts_to_TPM** 
Take a expression matrix with raw counts and convert into TPM normalization.
Note. Need a table with the gene IDs and its lengths, if you do not have it please try with **Gene_length_extraction_from_GTF** first.
 