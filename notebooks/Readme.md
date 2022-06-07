<h3> Scripts for building an expression array for clustering from scratch.</h3>
1. The process to build the expression matrix is composed of 4 sequential steps where a standardized log2 matrix (TPM+1) is produced as a final product. <br>
2. The raw-count files to integrate the expression array must be ready and available before to try this pipeline.<br>
3. The raw-count used here were produced with the HTSeq-count tool from bulk RNASeq libraries. <br> <br>

These scripts are complementary to my thesis work and will be totally free after publication. <br>
If you find some of them useful, feel free to try them, just please give proper credits like:
<br> Cynthia-Cardinault (2022) <br><br>

<h3>Python scripts to build an expression array for clustering</h3>
<ol>
    <li>1_Step1_integrating_raw_counts.ipynb </li>
    <li>2_Step2_TPM_normalization.ipynb<br> </li>
    <li>3_Step3_TPM_standardization.ipynb<br> </li>
    <li>4_Step4_Log2_scale.ipynb<br> </li>
</ol>
*Descriptions and instructions are provided in each file*

___________________________________________________________________________________________________________________

### SUPLEMENTARY MATERIAL

**Supplementary resource 1: extract the coverage percentage** <br>
Script to extract and plot the coverage percentage stats from the alignment's files --gotten with STAR. <br>
A. 0_percentual_alignments_statistics.ipynb <br>

**Supplementary resource 2: build a gene length file** <br>
This file is required for the *Step 2*. You can create one yourself or try this script which extracts the gene lenghts from the corresponding GTF file. Some tweaking may be required.<br>
A. Gene_length_extraction_from_GTF.ipynb  <br>
The script is an adaptation from the https://www.reneshbedre.com/blog/expression_units.html who uses https://github.com/reneshbedre/bioinfokit library. <br>
Note. Gene lengths are based on coding-gene features (CDS), if you wish to try another genetic-feature, g.e. exons, sRNA, etc., I encourage to do try it adjusting the required parameters in the gene-feature target. <br>

**Supplementary resource 3: Exploratory-Data-Analysis (EDA)** <br>
Scripts to get information from metadata, clusters or simply explore preliminar results. <br>
A. Venn_diagram_genes_in_ceros.ipynb (build a Venn diagram to show the relationship between genes at zeros in 2 data sets) <br>
B. Venn_HighTPMs.ipynb (build a histogram & a venn diagram to show the relation between high TPM values in 2 datasets) <br>
C. Annotation_traking_gene_data.ipynb (get stats about gene-features in a meta-data file) <br>
D. Matrix_Healthy_Exploratory_HTMLrpt.ipynb (Use sweetviz library to get a dynamic HTML file for a fast EDA (df)) <br>
E. Matrix_Infected_Exploratory_HTMLrpt.ipynb (Use sweetviz library to get a dynamic HTML file for a fast EDA (df)) <br>
F. 5_preliminar_dataset_correlations.ipynb (build a coexpresion matrix for a fast overview) <br>
G. EDA_Genes_by_GO_Biological_Term.ipynb (plot a dataset (gene-cluster) based on GO-Terms) <br>
H. EDA_Genes_by_GO_Component_Term.ipynb (plot a dataset (gene-cluster) based on GO-Terms) <br>
G. Merge_two_datasets.ipynb (merge 2 datasets in one based on the gene-ID) <br>

**Supplementary resource 4: clustering annotation** <br>
Scripts to deal with the clusters & clustering annotation files queried in the DAVID platform website. <br>
*Makes logical comparissons to extract high differenciated clusters from the control grp* <br>
A. 6_modules_percentual_differentiation.ipynb <br>
*Separate in individual clusters the clustering output files from DAVID* <br>
As we performed 4 runs, we have 1 file per run, RUN 0, 1, 3 y 6 <br>
A. Annotation_DAVID_Clusters_Decomposition_RUN_0.ipynb <br>
B. Annotation_DAVID_Clusters_Decomposition_RUN_1.ipynb <br>
C. Annotation_DAVID_Clusters_Decomposition_RUN_3.ipynb <br>
D. Annotation_DAVID_Clusters_Decomposition_RUN_6.ipynb <br>

**Supplementary resource 5: categorical plots and lm for top genes** <br>
Scripts to plot & show the top clusters in a compact way. <br>
A. Annotation_DAVID_Clusters_Decomposition_Top15Genes.ipynb (this is a lm() model) <br>
B. Annotation_Bubble_plot_for_DAVID_Clusters.ipynb (this is a bubble plot for the 4 runs) <br>
C. Annotation_Venn_Diagram_for_DAVID_Clustering.ipynb (this is a Venn Diagram plot for the 4 runs) <br>

Good luck and success<br>
Cynthia SC<br><br>