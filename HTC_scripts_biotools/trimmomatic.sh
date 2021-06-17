#!/bin/bash

# Made on: CSC Sep-30-2020
# Last.md: Cynthia-Soto Feb-2-2021

######################################################################
#
#                          DATASET3 AThaliana Infectada
#			   Adding additional data to adjust correlation                          
#                         
######################################################################

# Load module
module load Trimmomatic/0.36

# Sintaxis


#for fn in Dataset2_At_infec_trimmo/SRR*;
#do
#samp=`basename ${fn}`
#echo "Processing sample ${samp}"
#	java -jar $TRIM/trimmomatic PE -phred33 -trimlog trimmo.log /data/run/cyntsc/tmp_trimmo/SRR10586399_1.fastq /data/run/cyntsc/tmp_trimmo/SRR10586399_2.fastq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_1_paired.fq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_1_unpaired.fq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_2_paired.fq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_2_unpaired.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40

#echo "done"
#done 


#java -jar $TRIM/trimmomatic PE -phred33 -trimlog trimmo.log /data/run/cyntsc/tmp_trimmo/SRR10586399_1.fastq /data/run/cyntsc/tmp_trimmo/SRR10586399_2.fastq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_1_paired.fq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_1_unpaired.fq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_2_paired.fq /data/run/cyntsc/Dataset2_At_sana_trimmo/SRR10586399_2_unpaired.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40

#java -jar $TRIM/trimmomatic SE -phred33 tmp_trimmo/SRR3383780.fastq /data/run/cyntsc/Dataset2_At_infec_trimmo/SRR3383780_trimmed.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40


java -jar $TRIM/trimmomatic SE -phred33 sra/SRR364391.fastq /data/run/cyntsc/Dataset3_At_infec_trimmo/SRR364391_trimmed.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40
java -jar $TRIM/trimmomatic SE -phred33 sra/SRR364392.fastq /data/run/cyntsc/Dataset3_At_infec_trimmo/SRR364392_trimmed.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:4
java -jar $TRIM/trimmomatic SE -phred33 sra/SRR364400.fastq /data/run/cyntsc/Dataset3_At_infec_trimmo/SRR364400_trimmed.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40
java -jar $TRIM/trimmomatic SE -phred33 sra/SRR364401.fastq /data/run/cyntsc/Dataset3_At_infec_trimmo/SRR364401_trimmed.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40
java -jar $TRIM/trimmomatic SE -phred33 sra/SRR6283146.fastq /data/run/cyntsc/Dataset3_At_infec_trimmo/SRR6283146_trimmed.fq ILLUMINACLIP:../../software/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40

