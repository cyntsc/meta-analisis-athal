# Made by. Cynthia SC
# Get some statistics from all HTSeq output reports set in a specific directory. 

# IN: Just go to the <dir> where are the HTSeq output reports and run the code.
# OUT: A file with all the lines including the specified pattern.
 
> my_statistics_HTSeq_counts.txt
grep 'no_feature' SRR* | sort   >> my_statistics_HTSeq_counts.txt 
grep 'ambiguous' SRR* | sort  >> my_statistics_HTSeq_counts.txt
grep 'alignment_not_unique' SRR* | sort  >> my_statistics_HTSeq_counts.txt 
grep 'not_aligned' SRR* | sort >> my_statistics_HTSeq_counts.txt

