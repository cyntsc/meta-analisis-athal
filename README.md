
<h2>Here is the metatranscriptomic analysis done to identify genetic patterns related with the consensus response in Arabidopsis under diverse  biotic stressors caused by fungi</h2>

**Experimental Design**
* The experiment was designed in 2 blocks. The first with 8 transcriptomes from a healthy plant arranged in 4 groups for control (healthy12, healthy18, healthy24 and healthy30), and the second with 17 transcriptomes from a plant infected with three types of Ascomycetes fungi (B=Botrytis cinerea, Ch=Colletotrichum higginsianum, and Ss=Sclerotinia sclerotiorum) arranged in 5 groups for the treatments (Bc12, Bc18, Bc24, Ch22, Ch40 and Ss30) --the digits correspond to the time of inoculation (hpi) with the fungus.<br>
* The treatments represent 68% of the included samples (32% for C higginsianum, 24% for B cinerea and 12% for S sclerotiorum), the remaining 32% corresponds to the controls.<br>
* RNA was extracted from leaf among a range from 0 to 48 hpi <br>
* Just transcriptomes sequenced in Illumina platforms with a lenght.seq>100 & reads>5G were considered. <br>

**Material & Methods**
* Data are bulk RNASeq libraries downloaded from the SRA-NCBI:<br>
ID Bioprojects: PRJNA148307 (SRR364389, SRR364390, SRR364391, SRR364392, SRR364400, SRR364401, SRR364398 and SRR364399) for arabidopsis infected with Colletotrichum higginsianum at 22 and 40 hpi, PRJNA315516 (SRR3383696, SRR3383697,  SRR3383779 and SRR3383780) arabidopsis infected with Botrytis cinerea at 12 and 18 hpi, PRJNA593073 (SRR10586397 and SRR10586399) with Botrytis cinerea at 24 hpi, and PRJNA418121 (SRR6283146, SRR6283147 and SRR6283148) arabidopsis infected with Sclerotinia sclerotiorum at 30 hpi.  The arabidopsis healthy RNASeq libraries were downloaded from the same repository under the ID Bioprojects: PRJNA315516 (SRR3383640 and SRR3383641) mock treatment at 12hr, (SRR3383782 and SRR3383783) mock treatment at 18 hr and (SRR3383821 and SRR3383822) mock treatment at 30hr, and PRJNA418121 (SRR6283144 and SRR6283145) mock treatment at 30hr. <br>
* For the RNASeq raw-counts, I followed the alignment to genetic reference approach. The GENOME TAIR10 (GenBank accessions CP002684 – CP002688) and the ARAPORTt11 annotation were used, with a target of 27655 Protein-Coding-Sequence (CDS).
* A Clustering (ML not-supervised) approach was used, made with WGCNA. Signed-networks were built with the *Pearson* method. A threshold 𝛃=0.80 for signed-ntw & merged at 0.1 euclidean distances was set. Genetic-Modules with corr ≷ 0.75 were extracted.<br>
* perks
* Genetic-Modules on infected plants of interest were identified through logical operations, extracting modules differentiated between 100 and 77% from the healthy plants.<br>

***

**Folder content**:<br>

**HTC_scripts_biotools:** scripts to get the raw-counts.<br>
**R_scripts_WGCNA:** scripts in R to build the genetic ntws with WGCNA.<br>
**grep_utilities:** grep build-in commands to extract data from files (g.e:alignment stats).<br>
**meta-data:** mostly files to link external data to the genetic ntws.<br>
**notebooks:** phython 3+ scritps to build the expression matrices. <br>
**results-data:** statistical results, ntw results, intermedian results, all that can be considered a product comming from a process.  <br><br>
**Supplementary material**
 
 ***

**WGCNA Coexpression Ntw(s) were built with the help of several online public resources. Here, I share you some links that hopefully can help you in your own project:**<br>
Please be aware you need to adjust code to your own needs <br>
https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/ <br>
http://pklab.med.harvard.edu/scw2014/WGCNA.html <br>
https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA <br>
https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/ <br><br>

 ***  

### If you just want to have a glance at the outputs of the main scripts of this project, you can access this information by following the links.

#### Script to build the expression matrices
**2_Matrix_A_integrate_raw_countsb** <br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/2_Matrix_A_integrate_raw_counts.ipynb <br>

**3_Matrix_B_TPM_normalization**<br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/3_Matrix_B_TPM_normalization.ipynb <br>

**4_Matrix_C_TPM_standardization**<br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/4_Matrix_C_TPM_standardization.ipynb<br>

**5_Matrix_D_E_Log2_Atypicals**<br>
https://nbviewer.jupyter.org/github/cyntsc/meta-analisis-athal/blob/master/notebooks/5_Matrix_D_E_Log2_Atypicals.ipynb<br>

Some interactive HTML files about the expression matrices built can be downloaded from the **results-data/matrices_de_expresion/** folder: (just download the file in your local PC and open it. <br>
[Interactive stats for ARABIDOPSIS HEALTHY MatrixD.html](https://github.com/cyntsc/meta-analisis-athal/blob/7c87b3532f85106127df3ff68ab47445221c5971/results-data/matrices_de_expresion/SWEETVIZ_RPT_ARABIDOPSIS_HEALTHY_D.html) 
<br>
[Interactive stats for ARABIDOPSIS INFECTED MatrixE.html](https://github.com/cyntsc/meta-analisis-athal/blob/7c87b3532f85106127df3ff68ab47445221c5971/results-data/matrices_de_expresion/SWEETVIZ_RPT_ARABIDOPSIS_INFECTED_E.html)
<br>
[Interactive stats for comparition in ARABIDOPSIS INFECTED MatrixD and E.html](https://github.com/cyntsc/meta-analisis-athal/blob/7c87b3532f85106127df3ff68ab47445221c5971/results-data/matrices_de_expresion/SWEETVIZ_RPT_ARABIDOPSIS_INFECTED_D_E.html)
<br>
[More stats...](https://github.com/cyntsc/meta-analisis-athal/tree/master/results-data/matrices_de_expresion)<br>

***

#### Script to get the genetic modules in the expression matrices

**WGCNA scripts:** [WGCNA scripts for Signed-Ntw (pearson)](https://github.com/cyntsc/meta-analisis-athal/tree/master/R_scripts_WGCNA)<br>

**WGCNA results:** [Genetic Modules gotten with WGCNA Signed-Ntw (pearson)](https://github.com/cyntsc/meta-analisis-athal/tree/master/results-data/wgcna)<br>
Preferred genetic modules results are in folders: <br>
Athal_healthy_mods_merged_MatrixD for A thaliana Healthy<br>
Athal_infected_mods_merged_MatrixE for A thaliana Infected<br>

***

Be happy and Enjoy!<br>
Cynthia SC



