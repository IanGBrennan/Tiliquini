# Alignments

*Tiliquini_Concatenated.phylip.zip*: the concatenated alignment of all AHE loci. To split this alignment into locus specific alignment files we recommend using the [*AMAS*](https://github.com/marekborowiec/AMAS) tool.  
You can quickly divide the alignment using the provided partition file with AMAS like this:
```
AMAS.py split -i Tiliquini_Concatenated.phylip -f phylip -d dna -l Tiliquini_Concatenated_partitions.txt -u fasta
```

*Tiliquini_Concatenated_partitions.txt*: the locus partitions for the concatenated alignment file above.

*Tiliquini_PerLocus_Summary.csv*: the locus specific summary as output by *AMAS* (function: *summary*).