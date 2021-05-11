## Code to transform raw-count into TPM normalized and standarized values.
This code is complementary to my thesis report and is completely free. If you wish to use it as part of public paper material or online resources please give the corresponding credits as: Cynthia-Cardinault (2020). <br>

**0_alignments_statis_all.ipynb** <br>
Integration and construction of basic statistics for the alignments produced by STAR.<br>

**2_compouse_expr_matrix_raw_counts.ipynb** <br>
Compouse an expression matrix with raw counts (not normalized) produced by STAR <br>
Get some stats and distributions.<br>
Remove genes with ceros across all samples and plot again.

**3_Transform_rawcounts_to_TPM.ipynb**<br>
Take a expression matrix with raw counts and convert into TPM normalization.
Note. Need a table with the gene IDs and its lengths, if you do not have it please try with **Gene_length_extraction_from_GTF** first.

**4_TPM-Filter_Log2-TPM_transformations.ipynb** <br>
Standardization of a TPM-normalized expression matrix. <br>

**Gene_length_extraction_from_GTF** 
Create a file with ID-Genes and lenghts from a GTF file <br>
For those concerned with converting raw counts from an expression array to TMP counts for any purpose. These scripts are an adaptation from the https://www.reneshbedre.com/blog/expression_units.html who uses https://github.com/reneshbedre/bioinfokit library. 
Note. Gene lengths are based on coding-gene features (CDS), if you wish a more complex approach like average exons or some-like, I encourage to do the necessary changes. 

## Other resources

**8_df_correlations.ipynb** <br>
1) Correlate the expression values of RNASeq samples provided in a df <br>
2) Plot the resulted matrix in a heatmap <br>
3) Provided resources to extract data and plot specific samples of interest<br>

## DEPRECATED (Old versions)

**1_stats_raw_counts.ipynb** <br>
Getting basic statistics over raw count files. <br>

**1_log2_transformation_massive_ctrl.ipynb / 1_log2_transformation_massive_infect.ipynb** <br>
First, A.thaliana healthy or infected raw-counts are transformed to Log2+1 scale.<br>
Second, individual stats by file are generated, and the number of zeros by sample is calculated. <br>
Third, New stats are performed. <br>

**2_exploring_data_outliers_masive.ipynb / 2_exploring_data_outliers_masive_infect.ipynb** <br>
Getting basic statistics about the raw-count vs the (log2+1) distributions.

**3_merge_dataframes_ctrl.ipynb / 3_merge_dataframes_infect.ipynb**<br>
1) Merge a specific number of dataframes in a single matrix<br>
2) Some basic statistics are built with the data frame.  <br>
3) Zeros across all samples are dropped. <br>
