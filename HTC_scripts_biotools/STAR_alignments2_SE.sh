#Cynthia-Soto September-2020

######################################################################
#
#                          DATASET2 AThaliana infected with fungi
#                          Step 2/2 Creating Alignments for RNASeq SE layout  
#                          First need a index genome directory: see STAR_Idx.sh
#                         
######################################################################

# Load module
module load STAR/2.7.5a

# Sintaxis:
# The basic options for aligning reads to the genome using STAR are:
#    --runThreadN: number of threads / cores
#    --readFilesIn: /path/to/FASTQ_file
#    --genomeDir: /path/to/genome_indices_directory
#    --outFileNamePrefix: prefix for all output files

# Listed below are additional parameters that we will use in our command:
#    --outSAMtype: output filetype (SAM default)
#    --outSAMunmapped: what to do with unmapped reads

# NOTES:
# Default filtering is applied in which the maximum number of multiple alignments allowed for a read is set to 10. If a read exceeds this number there is no alignment output. To change the default you can use --outFilterMultimapNmax
# “STAR’s default parameters are optimized for mammalian genomes. Other species may require significant modifications of some alignment parameters; in particular, the maximum and minimum intron sizes have to be reduced for organisms with smaller introns”
# '--outSAMtype SAM' only, or with --outSAMtype BAM Unsorted|SortedByCoordinate

echo > mylog_STAR_alignments2_SE.txt 
echo "Making  SAalignments with STAR/2.7.5a" >> mylog_STAR_alignments2_SE.txt
echo "Made by Cynthia.Soto" >> mylog_STAR_alignments2_SE.txt
echo "Date:September 15, 2020" >> mylog_STAR_alignments2_SE.txt

mkdir athal_STAR_alignments2

for file in Dataset2_At_infec_trimmo_SE/*.fq
do
	echo $file
	STAR --genomeDir athal_STAR_idx_genome/ --runThreadN 12 --readFilesIn $file --outFileNamePrefix athal_STAR_alignments2/$file -alignIntronMin 8 --alignIntronMax 25000 --outSAMtype BAM Unsorted #SortedByCoordinate
        echo "Executed job for file:" $file >> mylog_STAR_alignments2_SE.txt
done

echo "SAM/BAM Outputs on athal_STAR_alignments2 Directory"  >> mylog_STAR_alignments2_SE.txt
echo "My log bash file: mylog_STAR_alignments2_SE.txt " >> mylog_STAR_alignments2_SE.txt
  

# Command line example
# STAR --genomeDir /data/run/cyntsc/athal_STAR_idx_genome/ --runThreadN 12 
#	--readFilesIn /data/run/cyntsc/Dataset2_At_infec_trimmo/SRR3383640_trimmo.fq 
#	--outFileNamePrefix /data/run/cyntsc/athal_STAR_alignments2/SRR3383640 
#	-alignIntronMin 8 --alignIntronMax 25000  --outSAMtype BAM Unsorted SortedByCoordinate 













