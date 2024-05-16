# Alignments

+ *Tiliquini_AHE_alns.zip*: individual locus and concatenated alignments of all AHE loci, including partition file. To split the concatenated alignment into locus specific alignment files we recommend using the [*AMAS*](https://github.com/marekborowiec/AMAS) tool.  
You can quickly divide the alignment using the provided partition file with AMAS like this:
```
AMAS.py split -i Tiliquini_AHE_concat.fasta -f fasta -d dna -l Tiliquini_AHE_partitions.txt -u fasta
```
+ *Tiliquini_AHE_CDS_RBH.csv*: results of the reciprocal best blast hit search of AHE targets against *Anolis carolinensis* CDS.

+ *Tiliquini_PerLocus_Summary.csv*: the locus specific summary as output by *AMAS* (function: *summary*).

+ *Tiliquini_PerSample_Summary.csv*: the sample specific summary, including number of sequenced AHE targets, registration numbers, subfamily status, and original data source. 
