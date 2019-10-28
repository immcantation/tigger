**airrDb** - *Example human immune repertoire data*

Description
--------------------

A `data.frame` of example V(D)J immunoglobulin sequences derived from a 
single individual (PGP1), sequenced on the Roche 454 platform, and assigned by
IMGT/HighV-QUEST to IGHV1 family alleles.






Format
-------------------
A `data.frame` where rows correspond to unique V(D)J sequences and
columns include:

+  `"sequence_alignment"`: IMGT-gapped V(D)J nucleotide sequence.
+  `"v_call"`: IMGT/HighV-QUEST V segment allele calls.
+  `"d_call"`: IMGT/HighV-QUEST D segment allele calls.
+  `"j_call"`: IMGT/HighV-QUEST J segment allele calls.
+  `"junction_length"`: Junction region length.


References
-------------------


1.  Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
sequencing data reveals a high frequency of novel immunoglobulin V gene 
segment alleles. PNAS. 112(8):E862-70.





See also
-------------------

See [SampleDb](SampleDb.md) for Change-O formatted version of `airrDb`.






