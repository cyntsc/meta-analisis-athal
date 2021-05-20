#/bin/sh  
#Cynthia-Soto September 2020

######################################################################
#
#                          DATASET2 AThaliana infected 
#                          Counting  
#                          First need a index genome directory: see STAR_Idx.sh
#                         
######################################################################

# Files to process:
# SRR10586397_trimmedPEAligned.sortedByCoord.out.bam  SRR10586399_trimmedPEAligned.sortedByCoord.out.bam

# NOTE. HTSeq reports the number of alignments in the SAM file which are of reads that mapped to more than one place, whereas STAR reports the number of reads which align to more than one place. 
# Good practice recommendations encourage to use splicing-aware aligner such as STAR.

# --nonunique none (read with more than one feature not counted for any features)
# --stranded=no (default is yes)
# If your RNA-Seq data has not been made with a strand-specific protocol, this causes half of the reads to be lost. Hence, make sure to set the option --stranded=no unless you have strand-specific data!
# For paired-end data, the alignment have to be sorted either by read name or by alignment position.

module load HTSeq/0.11.1.Py3
# Directory to store the quantification files 
#mkdir athal_HTSeq_counts2

#IFS='/'

echo > /data/run/cyntsc/mylog_htseq2_PE.txt 
echo "Counting  SAalignment reads from STAR with HTSeq-count" >> /data/run/cyntsc/mylog_htseq2_PE.txt
echo "Made by Cynthia.Soto" >> /data/run/cyntsc/mylog_htseq2_PE.txt
now=$(date)
echo "Date: $now" >> /data/run/cyntsc/mylog_htseq2_PE.txt
echo "Count files on athal_HTSeq_counts2 Directory"  >> mylog_htseq2_PE.txt
echo "My log bash file: mylog_htseq2_PE.txt " >> mylog_htseq2_PE.txt

# PE reads need specify required ordered sam file by name or position. Default is name.
htseq-count --format=bam --order=pos --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2_PE/SRR10586397_trimmedPEAligned.sortedByCoord.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR10586397pe

htseq-count --format=bam --order=pos --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2_PE/SRR10586399_trimmedPEAligned.sortedByCoord.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR10586399pe


