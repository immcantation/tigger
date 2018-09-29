**SampleDb** - *Example human immune repertoire data*

Description
--------------------

A `data.frame` of example V(D)J immunoglobulin sequences derived from a 
single individual (PGP1), sequenced on the Roche 454 platform, and assigned by
IMGT/HighV-QUEST to IGHV1 family alleles.




Format
-------------------
A `data.frame` where rows correspond to unique V(D)J sequences and
columns include:

+  `"SEQUENCE_IMGT"`: IMGT-gapped V(D)J nucleotide sequence.
+  `"V_CALL"`: IMGT/HighV-QUEST V segment allele calls.
+  `"D_CALL"`: IMGT/HighV-QUEST D segment allele calls.
+  `"J_CALL"`: IMGT/HighV-QUEST J segment allele calls.
+  `"JUNCTION_LENGTH"`: Junction region length.


References
-------------------


1.  Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
sequencing data reveals a high frequency of novel immunoglobulin V gene 
segment alleles. PNAS. 112(8):E862-70.






