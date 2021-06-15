<h2>This repository containt the meta-analysis done to track the consensus response in Arabidopsis under several criterias of biotic stress caused by fungi</h2>
 
Data are bulk RNA-Seq libraries downloaded from the SRA-NCBI<br>
The ML approach was Gene-Coexpression-Ntw-Analysis ( Hierarchical Clustering ) <br>
All the results are genetic modules (hub-genes) found highly correlated in the included RNA samples <br>
Number of samples = 25, 17 infected and 8 muck treatments <br>
RNA was extracted from leaf among a range from 0 to 48 hpi<br>

***

**Folders containt**:<br>

**HTC_scripts_biotools:** scripts used to get the gene expression raw-counts.<br>
**R_scripts_WGCNA:** scripts written in R to build the genetic ntws.<br>
**grep_utilities:** grep utilities to get some stats.<br>
**meta-data:** mostly files to link external data to the genetic ntws.<br>
**notebooks:** phython 3+ scritps to build the expression matrices. <br>
**results-data:** statistical results, ntw results, intermedian results, all that can be considered a product comming from a process.  <br><br>
 
 ***

**WGCNA Coexpression Ntw(s) were built with the help of several online public resources. Here, I share you some links that hopefully can help you in your own project:**<br>
Please be aware you need to adjust code to your own needs <br>
https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/ <br>
http://pklab.med.harvard.edu/scw2014/WGCNA.html <br>
https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA <br>
https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/ <br><br>

 ***  

### If you just want to have a glance at the outputs of the main scripts of this project, you can access this information by following the links.

#### Script to build the expression matrices used in this project
**2_Matrix_A_integrate_raw_countsb** <br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/2_Matrix_A_integrate_raw_counts.ipynb <br>

**3_Matrix_B_TPM_normalization**<br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/3_Matrix_B_TPM_normalization.ipynb <br>

**4_Matrix_C_TPM_standardization**<br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/4_Matrix_C_TPM_standardization.ipynb<br>

**5_Matrix_D_E_Log2_Atypicals**<br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/5_Matrix_D_E_Log2_Atypicals.ipynb<br>

Some interactive HTML files about the expression matrices built can be downloaded from the **results-data/matrices_de_expresion/** folder: <br>
SWEETVIZ_RPT_ARABIDOPSIS_HEALTHY_B_C.html <br>
SWEETVIZ_RPT_ARABIDOPSIS_HEALTHY_D.html<br>
SWEETVIZ_RPT_ARABIDOPSIS_INFECTED_B_C.html <br>
SWEETVIZ_RPT_ARABIDOPSIS_INFECTED_D_E.html <br>
SWEETVIZ_RPT_ARABIDOPSIS_INFECTED_E.html <br>

#### **WGCNA scripts are in R_scripts_WGCNA folder**
#### **WGCNA genetic modules gotten are in results-data/wgcna folder**

Be happy and Enjoy!<br>
Cynthia SC





