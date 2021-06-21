<h2> Scripts for transforming raw-count produced by HTSeq into a normalized and standardized TPM expression matrix </h2>
This code is complementary to my thesis work and will be fully free after this work be published. <br>
If you find useful all or some scripts published here, feel free to download them, just please give the corresponding credits as:<br> Cynthia-Cardinault (2020) <br><br>
Good luck and success<br>
Cynthia SC<br><br>

**STEP A <br> 
2_Matrix_A_integrate_raw_counts.ipynb** <br>
Create Matrix A: compouse an expression matrix with raw-counts (not normalized) produced by HTSeq <br>
Get some stats and distributions (text-files and plots).<br>
Remove genes with zeros across all samples and get the stats and plots again.  <br>

**STEP B <br> 
3_Matrix_B_TPM_normalization.ipynb**<br>
Create Matrix B: take the expression matrix A and transform into a TPM-expression matrix.<br>
Note. You need a file with the gene IDs and lengths of your spc, if you do not have it please try with **Gene_length_extraction_from_GTF** first, or perhaps you can be creative and integrate your data as convinient.<br>

**STEP C <br> 
4_Matrix_C_TPM_standardization.ipynb** <br>
Create Matrix C: take the expression matrix B and create a standardized TPM-expression matrix.<br>

**STEP D <br>
5_Matrix_D_E_Log2_Atypicals.ipynb** <br>
Create Matrix D: take the expression matrix C and transform into a Log2(TPM+1) expression matrix to reduce the effect of the extreme values in the data.<br>
Create Matrix E: take the expression matrix D and if exists remarkable atypical distributions, remove them from the expression matrix to reduce spurious association in a GCN. <br>

**CLUSTER DIFFERENTIATION <br>
8_modules_differentiation_mergedClusters.ipynb** <br>
Get the percentual difference between each module in the infected merged-ntw against the modules in the control merged-ntw. <br>

**Optionally<br>
5_Matrix_B_to_E.ipynb** <br>
Assumming that you do not want to filter the data, you can use this script to jump the step C, and just create a Log2(TPM+1) expression matrix.  <br>

**Other supplementary resources:** <br>
This files are not strictically part of the project, most of them were used to test some arrays (steps), thus these are not indispensable to build the matrix. But if you feel comfortable you can try some of them.

**1_percentual_alignments_statistics.ipynb** <br>
Extracts the coverage percentage of the alignments generated with STAR and plots them in a bar plot.<br>

**Gene_length_extraction_from_GTF** <br>
Create a file with ID-Genes and lenghts from a GTF file <br>
For those concerned with converting raw counts from an expression array to TMP counts for any purpose. These scripts are an adaptation from the https://www.reneshbedre.com/blog/expression_units.html who uses https://github.com/reneshbedre/bioinfokit library. <br>
Note. Gene lengths are based on coding-gene features (CDS), if you wish a more complex approach like average exons or some-like, I encourage to do the necessary changes. <br>

**6_Merge_two_matrix.ipynb** <br> 

**7_df_correlations.ipynb** <br> 

**Venn_diagram_genes_in_ceros.ipynb**<br>
