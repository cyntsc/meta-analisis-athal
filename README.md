<h2>Meta-analysis to track the consensus genetic response in Arabidopsis under various criterias of biotic stress caused by fungi</h2>
 
Data are RNA-Seq raw libraries <br>
The methodology followed is Gene-Coexpression-Ntw-Analysis ( ML not supervised ) <br>
Results are genetic modules (hub-genes) found highly correlated in the selected samples <br>
N samples = 25, 17 infected samples and 8 control samples<br>
RNA extracted from leaf among 0 to 48 hpi<br>

  
**Folders containt**:<br>

**HTC_scripts_biotools:** scripts used to get the gene expression raw-counts.<br>
**R_scripts_WGCNA:** scripts written in R to build the genetic ntws.<br>
**grep_utilities:** grep utilities to get some stats.<br>
**meta-data:** mostly files to link external data to the genetic ntws.<br>
**notebooks:** phython 3+ scritps to build the expression matrices. <br>
**results-data:** statistical results, ntw results, intermedian results, all that can be considered a product comming from a process.  <br><br>

**WGCNA Coexpression Ntw(s) were built with the help of several online public resources. Here, I share you some links that hopefully can help you in your own project:**<br>
Please be aware you need to adjust code to your own needs <br>
https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/ <br>
http://pklab.med.harvard.edu/scw2014/WGCNA.html <br>
https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA <br>
https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/ <br><br>

Good luck!<br>
Cynthia SC





