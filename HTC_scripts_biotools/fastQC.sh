#!/bin/bash

# Made by: Cynthia-Soto Sep-2020

######################################################################
#
#                          DATASET2 AThaliana Infectada
#                          Datasets adicionales 1-Feb 2020 / csc                          
#                         
######################################################################

# Load module
module load FastQC/0.11.3

# Sintaxis:
# fastqc --extract --outdir=</path/to/my/output_directory> -f fastq <raw_sequences_file1.fastq>


for fn in Dataset3_At_infect_trimmo/*.fq
do
	samp=`basename ${fn}`
	echo "Processing sample ${samp}"
	#fastqc --extract --outdir=/data/run/cyntsc/Dataset3_At_infect_trimmo -f fastq tmp_trimmo/${samp}
	fastqc --outdir=/data/run/cyntsc/Dataset3_At_infect_trimmo -f fastq Dataset3_At_infect_trimmo/${samp}
	echo "done"
done 

