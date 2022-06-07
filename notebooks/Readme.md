### Scripts for building an expression array for clustering from scratch
<ol>
    <li> The process to build the expression matrix is composed of 4 sequential steps where a standardized log2 matrix (TPM+1) is produced as a final product. </li>
    <li> The raw-count files to integrate the expression array must be ready and available before to try this pipeline. </li>
    <li> The raw-count used here were produced with the HTSeq-count tool from bulk RNASeq libraries. </li>
</ol>
    
These scripts are complementary to a thesis research and will be totally free after publication. <br>
If you find some of them useful, feel free to try them, just please give proper credits like: <br> <br> 
*Cynthia-Cardinault (2022)* <br>

<h3>Python scripts to build an expression array for clustering</h3>
<ul>
    <li>1_Step1_integrating_raw_counts.ipynb </li>
    <li>2_Step2_TPM_normalization.ipynb<br> </li>
    <li>3_Step3_TPM_standardization.ipynb<br> </li>
    <li>4_Step4_Log2_scale.ipynb<br> </li>
</ul>
*Descriptions and instructions are provided in each file*

___________________________________________________________________________________________________________________

### SUPLEMENTARY MATERIAL

**Supplementary resource 1: extract the coverage percentage** <br>
Script to extract and plot the coverage percentage stats from the alignment's files --gotten with STAR. <br>
<ul>
    <li> 0_percentual_alignments_statistics.ipynb </li>
</ul>

**Supplementary resource 2: build a gene length file** <br>
This file is required for the *Step 2*. You can create one yourself or try this script which extracts the gene lenghts from the corresponding GTF file. Some tweaking may be required.<br>
<ul>
    <li> Gene_length_extraction_from_GTF.ipynb  </li>
</ul>
The script is an adaptation from the https://www.reneshbedre.com/blog/expression_units.html who uses https://github.com/reneshbedre/bioinfokit library. <br>
Note. Gene lengths are based on coding-gene features (CDS), if you wish to try another genetic-feature, g.e. exons, sRNA, etc., I encourage to do try it adjusting the required parameters in the gene-feature target. <br><br>

**Supplementary resource 3: Exploratory-Data-Analysis (EDA)** <br>
Scripts to get information from metadata, clusters or simply explore preliminar results. <br>
<ul>
    <li> Venn_diagram_genes_in_ceros.ipynb (build a Venn diagram to show the relationship between genes at zeros in 2 datasets) </li>
    <li> Venn_HighTPMs.ipynb (build a histogram and a venn diagram to show the relation between high TPM values in 2 datasets) </li>
    <li> Annotation_traking_gene_data.ipynb (get stats about gene-features in a meta-data file) </li>
    <li> Matrix_Healthy_Exploratory_HTMLrpt.ipynb (Use sweetviz library to get a dynamic HTML file for a fast EDA (df)) </li>
    <li> Matrix_Infected_Exploratory_HTMLrpt.ipynb (Use sweetviz library to get a dynamic HTML file for a fast EDA (df)) </li>
    <li> 5_preliminar_dataset_correlations.ipynb (build a coexpresion matrix for a fast overview) </li>
    <li> EDA_Genes_by_GO_Biological_Term.ipynb (plot a dataset (gene-cluster) based on GO-Terms) </li>
    <li> EDA_Genes_by_GO_Component_Term.ipynb (plot a dataset (gene-cluster) based on GO-Terms) </li>
    <li> Merge_two_datasets.ipynb (merge 2 datasets in one based on the gene-ID) </li>
</ul>

**Supplementary resource 4: clustering annotation** <br>
Scripts to deal with the clusters & clustering annotation files queried in the DAVID platform website. <br><br>
*Makes logical comparissons to extract high differenciated clusters from the control grp* <br>
<ul><li> 6_modules_percentual_differentiation.ipynb </li></ul> <br>
*Separate into individual clusters the clustering output files from DAVID*<br><br>
As we performed 4 runs, we have 1 file per run, RUN 0, 1, 3 y 6 <br>
<ul> 
    <li> Annotation_DAVID_Clusters_Decomposition_RUN_0.ipynb  </li>
    <li> Annotation_DAVID_Clusters_Decomposition_RUN_1.ipynb  </li>
    <li> Annotation_DAVID_Clusters_Decomposition_RUN_3.ipynb  </li>
    <li> Annotation_DAVID_Clusters_Decomposition_RUN_6.ipynb  </li>
</ul>

**Supplementary resource 5: categorical plots and lm for top genes** <br>
Scripts to plot & show the top clusters in a compact way. <br>
<ul> 
    <li> Annotation_DAVID_Clusters_Decomposition_Top15Genes.ipynb (this is a lm() model)  </li>
    <li> Annotation_Bubble_plot_for_DAVID_Clusters.ipynb (this is a bubble plot for the 4 runs)  </li>
    <li> Annotation_Venn_Diagram_for_DAVID_Clustering.ipynb (this is a Venn Diagram plot for the 4 runs)  </li>
</ul> 
<
Good luck and success<br>
Cynthia SC