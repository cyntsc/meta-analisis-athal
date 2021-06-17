
# Get the statistics of the counts 

> my_statistics_HTSeq_counts.txt
grep 'no_feature' SRR* | sort   >> my_statistics_HTSeq_counts.txt 
grep 'ambiguous' SRR* | sort  >> my_statistics_HTSeq_counts.txt
grep 'alignment_not_unique' SRR* | sort  >> my_statistics_HTSeq_counts.txt 
grep 'not_aligned' SRR* | sort >> my_statistics_HTSeq_counts.txt

