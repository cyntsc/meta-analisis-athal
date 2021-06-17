

> my_statistics_STAR_alignment.txt
grep 'Uniquely mapped reads %' SRR*_trimmo.fqLog.final.out | sort   >> my_statistics_STAR_alignment.txt 
grep 'Average mapped length' SRR*_trimmo.fqLog.final.out | sort  >> my_statistics_STAR_alignment.txt 
grep 'Number of splices: Total' SRR*_trimmo.fqLog.final.out | sort  >> my_statistics_STAR_alignment.txt 
grep 'Mismatch rate per base, %' SRR*_trimmo.fqLog.final.out | sort  >> my_statistics_STAR_alignment.txt 
grep 'Deletion average length' SRR*_trimmo.fqLog.final.out | sort  >> my_statistics_STAR_alignment.txt 
grep '% of reads unmapped: too short' SRR*_trimmo.fqLog.final.out | sort >> my_statistics_STAR_alignment.txt
cat my_statistics_STAR_alignment.txt
