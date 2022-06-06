<h2> The below python scripts are the source code to buil a expression matrix for clustering since zero. Raw-count must be provided, such as those produced with HTSeq-count or its equivalent in Salmon. The expression matrix is composed by 4-Steps where a standardized log2(TPM+1) matrix is produced. </h2>
These scripts are complementary to my thesis work and will be fully free after be published. <br>
If you find some of them useful feel free to use them, just please give the corresponding credits as:<br> Cynthia-Cardinault (2022) <br><br>
Good luck and success<br>
Cynthia SC<br><br>

**Next files must be run sequentially:**
1_Step1_integrating_raw_counts.ipynb <br>
2_Step2_TPM_normalization.ipynb<br>
3_Step3_TPM_standardization.ipynb<br>
4_Step4_Log2_scale.ipynb<br>

*Description and instructions are provided in each step*

___________________________________________________________________________________________________________________


**Supplementary resource 1:** <br>
This file extracts the coverage percentage of the alignments generated with STAR and plots them in a bar plot.<br>
0_percentual_alignments_statistics.ipynb <br>

**Supplementary resource 2:** <br>
These files were used for testing additional assays or explore the additional information. Feel free to try some of them. <br>
5_preliminar_dataset_correlations.ipynb
6_modules_percentual_differentiation.ipynb

**Supplementary resource 3:** <br>
These files were used for deal with the clustering annotations recovered with the DAVID Tool. Feel free to try some of them. <br>
As we performed 4 runs, we have 1 file per run, RUN 0, 1, 3 y 6 <br>
Annotation_DAVID_Clusters_Decomposition_RUN_0.ipynb <br>
Annotation_DAVID_Clusters_Decomposition_RUN_1.ipynb <br>
Annotation_DAVID_Clusters_Decomposition_RUN_3.ipynb <br>
Annotation_DAVID_Clusters_Decomposition_RUN_6.ipynb <br>

**Supplementary resource 4:** <br>
There files were used for the visualization tasks (plots) to show the results in a compact way <br>
Annotation_DAVID_Clusters_Decomposition_Top15Genes.ipynb (this is a lm() model <br>
Annotation_Bubble_plot_for_DAVID_Clusters.ipynb (this is a bubble plot for the 4 runs <br>
Annotation_Venn_Diagram_for_DAVID_Clustering.ipynb (this is a Venn Diagram plot for the 4 runs <br>

**Supplementary resource 5:** <br>
These file is necessary for run the *Step 2* if you do not have a Gene lenght file.
Gene_length_extraction_from_GTF.ipynb <br>
Create a file with ID-Genes and lenghts from a GTF file <br>
For those concerned with converting raw counts from an expression array to TMP counts for any purpose. These scripts are an adaptation from the https://www.reneshbedre.com/blog/expression_units.html who uses https://github.com/reneshbedre/bioinfokit library. <br>
Note. Gene lengths are based on coding-gene features (CDS), if you wish a more complex approach like average exons or some-like, I encourage to do the necessary changes. <br>

