

> my_statistics_STAR_alignment2_PE.txt
grep 'Uniquely mapped reads %' SRR*Log.final.out | sort   >> my_statistics_STAR_alignment2_PE.txt 
grep 'Average mapped length' SRR*Log.final.out | sort  >> my_statistics_STAR_alignment2_PE.txt 
grep 'Number of splices: Total' SRR*Log.final.out | sort  >> my_statistics_STAR_alignment2_PE.txt 
grep 'Mismatch rate per base, %' SRR*Log.final.out | sort  >> my_statistics_STAR_alignment2_PE.txt 
grep 'Deletion average length' SRR*Log.final.out | sort  >> my_statistics_STAR_alignment2_PE.txt 
grep '% of reads unmapped: too short' SRR*Log.final.out | sort >> my_statistics_STAR_alignment2_PE.txt
cat my_statistics_STAR_alignment2_PE.txt
