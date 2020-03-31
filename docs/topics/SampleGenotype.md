**SampleGenotype** - *Example genotype inferrence results*

Description
--------------------

A `data.frame` of genotype inference results from [inferGenotype](inferGenotype.md)
after novel allele detection via [findNovelAlleles](findNovelAlleles.md).
Source data was a collection of V(D)J immunoglobulin sequences derived from a single
individual (PGP1), sequenced on the Roche 454 platform, and assigned by
IMGT/HighV-QUEST to IGHV1 family alleles.






Format
-------------------

A `data.frame` where rows correspond to genes carried by an
individual and columns lists the alleles of those genes and their counts.


References
-------------------


1.  Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
sequencing data reveals a high frequency of novel immunoglobulin V gene 
segment alleles. PNAS. 112(8):E862-70.





See also
-------------------

See [inferGenotype](inferGenotype.md) for detailed column descriptions.






