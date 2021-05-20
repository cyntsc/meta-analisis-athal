#made by  Cynthia-Soto August-2020

######################################################################
#
#                          DATASET1 AThaliana SANA
#                          Steps: 1
#                          Salmon ALignment Mode Quantification
######################################################################

# Load module
module purge
module load Salmon/1.3.0

# Sintaxis:
# you’ve prepared your alignments using your favorite aligner and the results are in the file aln.bam, and assume that the sequence of the transcriptome you want to quantify is in the file transcripts.fa. You would run Salmon as follows:
# salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant

# NOTES:
# -l A or --libType A for library layout autodetection

salmon quant -t /data/run/cyntsc/Dataset1_At_sana_trimmo/SRR3383640_trimmo.fq  -l A -a /data/run/cyntsc/athal_STAR_alignments/SRR3383640_sanaAligned.out.sam
 -o /data/run/cyntsc/athal_SALMONam_quants/SRR3383640

# Opc con BAM file
salmon quant -t /data/run/cyntsc/Dataset1_At_sana_trimmo/SRR3383640_trimmo.fq  -l A -a /data/run/cyntsc/athal_STAR_alignments/SRR3383640_BAMUnsortedAligned.sortedByCoord.out.bam -o /data/run/cyntsc/athal_SALMONam_quants/SRR3383640

