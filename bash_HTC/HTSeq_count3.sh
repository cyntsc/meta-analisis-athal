#/bin/sh  
#Cynthia-Soto 4 Febrero 2021

######################################################################
#
#                          DATASET3 AThaliana infected 
#                          Counting  
#                          First need a index genome directory: see STAR_Idx.sh
#                         
######################################################################

# Files to process:


# Good practice recommendations encourage to use splicing-aware aligner such as STAR.

# --nonunique none (read with more than one feature not counted for any features)
# --stranded=no (default is yes)
# If your RNA-Seq data has not been made with a strand-specific protocol, this causes half of the reads to be lost. Hence, make sure to set the option --stranded=no unless you have strand-specific data!
# For paired-end data, the alignment have to be sorted either by read name or by alignment position.

module load HTSeq/0.11.1.Py3
# Directory to store the quantification files 
#mkdir athal_HTSeq_counts3

#IFS='/'

echo > /data/run/cyntsc/mylog_htseq3.txt 
echo "Counting  SAalignment reads from STAR with HTSeq-count" >> /data/run/cyntsc/mylog_htseq3.txt
echo "Made by Cynthia.Soto" >> /data/run/cyntsc/mylog_htseq3.txt
now=$(date)
echo "Date: $now" >> /data/run/cyntsc/mylog_htseq3.txt
echo "Count files on  Directory"  >> mylog_htseq3.txt
echo "My log bash file: mylog_htseq3.txt " >> mylog_htseq3.txt

#htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments3/SRR364391_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts3/SRR364391

#htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments3/SRR364392_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts3/SRR364392

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments3/SRR364400_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts3/SRR364400

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments3/SRR364401_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts3/SRR364401

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments3/SRR6283146_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts3/SRR6283146



