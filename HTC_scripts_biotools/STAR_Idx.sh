#made by  Cynthia-Soto August-2020

######################################################################
#
#                          DATASET1 AThaliana SANA
#                          Step 1/2 Creating a genome index 
#                          
#                         
######################################################################

# Load module
module load STAR/2.7.5a

# Sintaxis:
# The basic options to generate genome indices using STAR are as follows:
#    --runThreadN: number of threads
#    --runMode: genomeGenerate mode
#    --genomeDir: /path/to/store/genome_indices
#    --genomeFastaFiles: /path/to/FASTA_file
#    --sjdbGTFfile: /path/to/GTF_file
#    --sjdbOverhang: readlength -1  ( in the control dataset readlenght range from 93-125

# NOTES:
# For small genomes, you may need to add the following to create a proper index: "--genomeSAindexNbases 6"
# Star is memory intensive and requires at least 30 Gb to align to the human or mouse genomes
#  In case of reads of varying length, the ideal value for --sjdbOverhang is max(ReadLength)-1. In most cases, the default value of 100 will work similarly to the ideal value.

#Be sure the references are in the correct path


#Option 1:
STAR --runThreadN 12 --runMode genomeGenerate --genomeSAsparseD 3 --genomeSAindexNbases 7 --genomeDir /data/run/cyntsc/athal_STAR_idx_genome --genomeFastaFiles /data/run/cyntsc/athal_STAR_idx_genome/Athaliana_447_TAIR10.fa --sjdbGTFfile /data/run/cyntsc/athal_STAR_idx_genome/Araport11_GFF3_genes_transposons.201606.gtf --sjdbOverhang 124

