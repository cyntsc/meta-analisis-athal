

> my_statistics_STAR_alignment2_SE.txt
grep 'Uniquely mapped reads %' SRR*fqLog.final.out | sort   >> my_statistics_STAR_alignment2_SE.txt 
grep 'Average mapped length' SRR*fqLog.final.out | sort  >> my_statistics_STAR_alignment2_SE.txt 
grep 'Number of splices: Total' SRR*fqLog.final.out | sort  >> my_statistics_STAR_alignment2_SE.txt 
grep 'Mismatch rate per base, %' SRR*fqLog.final.out | sort  >> my_statistics_STAR_alignment2_SE.txt 
grep 'Deletion average length' SRR*fqLog.final.out | sort  >> my_statistics_STAR_alignment2_SE.txt 
grep '% of reads unmapped: too short' SRR*fqLog.final.out | sort >> my_statistics_STAR_alignment2_SE.txt
cat my_statistics_STAR_alignment2_SE.txt
