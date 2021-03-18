#/bin/sh  
#Cynthia-Soto August-2020

######################################################################
#
#                          DATASET1 AThaliana SANA
#                          Counting  
#                          First need a index genome directory: see STAR_Idx.sh
#                         
######################################################################


# NOTE. HTSeq reports the number of alignments in the SAM file which are of reads that mapped to more than one place, whereas STAR reports the number of reads which align to more than one place. 
# Good practice recommendations encourage to use splicing-aware aligner such as STAR.

# --nonunique none (read with more than one feature not counted for any features)
# --stranded=no (default is yes)
# If your RNA-Seq data has not been made with a strand-specific protocol, this causes half of the reads to be lost. Hence, make sure to set the option --stranded=no unless you have strand-specific data!
# For paired-end data, the alignment have to be sorted either by read name or by alignment position.

module load HTSeq/0.11.1.Py3
# Directory to store the quantification files 
#mkdir athal_HTSeq_counts

#IFS='/'

echo > /data/run/cyntsc/mylog_htseq.txt 
echo "Counting  SAalignments from STAR with HTSeq-count" >> /data/run/cyntsc/mylog_htseq.txt
echo "Made by Cynthia.Soto" >> /data/run/cyntsc/mylog_htseq.txt
echo "Date: Augost 17, 2020" >> /data/run/cyntsc/mylog_htseq.txt


#for file in /data/run/cyntsc/athal_STAR_alignments/SRR*.fqAligned.out.bam 
#do
#	sam_f=/data/run/cyntsc/athal_STAR_alignments/"${file##*/}"
#	gtf=/data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf
#	htseq_f=/data/run/cyntsc/athal_HTSeq_counts/"${file##*/}"

#	echo "${sam_f}" "${gtf}"\>"${htseq_f}"  
	
#	htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS "${sam_f}" "${gtf}" \> "${htseq_f}"  
#	#issue: read the last two files as one  :(	

#	echo "Executed job for file:" $file >> mylog_htseq.txt
#done

echo "Count files on athal_HTSeq_counts Directory"  >> mylog_htseq.txt
echo "My log bash file: mylog_htseq.txt " >> mylog_htseq.txt

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR3383640_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR3383640

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR3383641_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR3383641

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR3383782_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR3383782

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR3383783_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR3383783

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR3383821_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR3383821

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR3383822_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR3383822

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR6283144_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR6283144

htseq-count --mode=union --nonunique none --stranded=no --format=bam --idattr=gene_id --type=CDS /data/run/cyntsc/athal_STAR_alignments/SRR6283145_trimmo.fqAligned.out.bam /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf>/data/run/cyntsc/athal_HTSeq_counts/SRR6283145




