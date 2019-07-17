**reassignAlleles** - *Correct allele calls based on a personalized genotype*

Description
--------------------

`reassignAlleles` uses a subject-specific genotype to correct
correct preliminary allele assignments of a set of sequences derived
from a single subject.


Usage
--------------------
```
reassignAlleles(data, genotype_db, v_call = "V_CALL",
seq = "SEQUENCE_IMGT", method = "hamming", path = NA,
keep_gene = c("gene", "family", "repertoire"))
```

Arguments
-------------------

data
:   a `data.frame` containing V allele calls from a
single subject and the sample IMGT-gapped V(D)J sequences under
`"SEQUENCE_IMGT"`.

genotype_db
:   a vector of named nucleotide germline sequences
matching the calls detailed in `allele_calls`
and personalized to the subject

v_call
:   name of the column in `data` with V allele
calls. Default is `"V_CALL"`.

seq
:   name of the column in `data` with the 
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is SEQUENCE_IMGT

method
:   the method to be used when realigning sequences to
the genotype_db sequences. Currently, only `"hammming"`
(for Hamming distance) is implemented.

path
:   directory containing the tool used in the
realignment method, if needed. Hamming distance does
not require a path to a tool.

keep_gene
:   a string indicating if the gene (`"gene"`), 
family (`"family"`) or complete repertoire
(`"repertoire"`) assignments should be performed. 
Use of `"gene"` increases speed by minimizing required number of 
alignments, as gene level assignments will be maintained when possible.




Value
-------------------

A modifed input `data.frame` containing the best allele call from 
among the sequences listed in `genotype_db` in the 
`V_CALL_GENOTYPED` column.


Details
-------------------

In order to save time, initial gene assignments are preserved and
the allele calls are chosen from among those provided in `genotype_db`,
based on a simple alignment to the sample sequence.



Examples
-------------------

```R
# Extract the database sequences that correspond to the genotype
genotype_db <- genotypeFasta(SampleGenotype, SampleGermlineIGHV, novel=SampleNovel)

# Use the personlized genotype to determine corrected allele assignments
output_db <- reassignAlleles(SampleDb, genotype_db)
```




