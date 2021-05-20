
# Cyntsc Ago.2020

# Validaciones para los alineamientos. 
# The Python script htseq-qa takes a file with sequencing reads (either raw or aligned reads) and produces a PDF file with useful plots to assess the technical quality of a run.

load module SAMTools/1.9 
load module matplotlib/3.1.3.Py3

# Primero convertimos el archivo bam a sam con SAMTools v1.9
samtools view -h SRR3383640Aligned.out.bam > SRR3383640Aligned.out.sam


# Crea un plot(PDF)  de secuencias alineadas y no alineadas
htseq-qa -t sam SRR3383640Aligned.out.sam

