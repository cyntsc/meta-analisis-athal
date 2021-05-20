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
#SRR3383696
#SRR3383697
#SRR3383779
#SRR3383780
#SRR10586397 PE
#SRR10586399 PE
#SRR364389
#SRR364390
#SRR364398
#SRR364399

# NOTE. HTSeq reports the number of alignments in the SAM file which are of reads that mapped to more than one place, whereas STAR reports the number of reads which align to more than one place. 
# Good practice recommendations encourage to use splicing-aware aligner such as STAR.

# --nonunique none (read with more than one feature not counted for any features)
# --stranded=no (default is yes)
# If your RNA-Seq data has not been made with a strand-specific protocol, this causes half of the reads to be lost. Hence, make sure to set the option --stranded=no unless you have strand-specific data!
# For paired-end data, the alignment have to be sorted either by read name or by alignment position.

module load HTSeq/0.11.1.Py3
# Directory to store the quantification files 
mkdir athal_HTSeq_counts2

#IFS='/'

echo > /data/run/cyntsc/mylog_htseq2.txt 
echo "Counting  SAalignment reads from STAR with HTSeq-count" >> /data/run/cyntsc/mylog_htseq2.txt
echo "Made by Cynthia.Soto" >> /data/run/cyntsc/mylog_htseq2.txt
now=$(date)
echo "Date: $now" >> /data/run/cyntsc/mylog_htseq2.txt
echo "Count files on athal_HTSeq_counts2 Directory"  >> mylog_htseq2.txt
echo "My log bash file: mylog_htseq2.txt " >> mylog_htseq2.txt

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR3383696.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR3383696

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR3383697.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR3383697

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR3383779_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR3383779

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR3383780_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR3383780

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR364389_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR364389

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR364390_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR364390

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR364398_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR364398

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR364399_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR364399

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR6283147_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR6283147

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments2/SRR6283148_trimmed.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts2/SRR6283148


