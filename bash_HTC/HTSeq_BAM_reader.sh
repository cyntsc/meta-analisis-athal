#/bin/sh  
#Cynthia-Soto September 2020

######################################################################
#
#                          UTILITY: bam reader 
#                         
######################################################################


module load HTSeq/0.11.1.Py3
bam_reader = HTSeq.BAM_Reader( "/data/run/cyntsc/athal_STAR_alignments2/SRR3383696.fqAligned.out.bam" )
for a in itertools.islice( bam_reader, 5 ):
	print(a)

bam_writer = HTSeq.BAM_Writer.from_BAM_Reader( "region.bam", bam_reader )   #set-up˓→BAM_Writer with same header as reader
for a in bam_reader.fetch( region = "1:249000000-249200000" ):
	print("Writing Alignment", a, "to file", bam_writer.filename)
	bam_writer.write( a )
