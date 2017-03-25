





**reassignAlleles** - *Correct allele calls based on a personalized genotype*

Description
--------------------

`reassignAlleles` uses a subject-specific genotype to correct
correct preliminary allele assignments of a set of sequences derived
from a single subject.


Usage
--------------------
```
reassignAlleles(clip_db, genotype_db, method = "hamming", path = NA,
keep_gene = TRUE)
```

Arguments
-------------------

clip_db
:   a `data.frame` containing V allele calls from a
single subject under `"V_CALL"` and the sample
IMGT-gapped V(D)J sequences under
`"SEQUENCE_IMGT"`

genotype_db
:   a vector of named nucleotide germline sequences
matching the calls detailed in `allele_calls`
and personalized to the subject

method
:   the method to be used when realigning sequences to
the genotype_db sequences. Currently only "hammming"
(for Hamming distance) is implemented.

path
:   directory containing the tool used in the
realignment method, if needed. Hamming distance does
not require a path to a tool.

keep_gene
:   logical indicating if gene assignments should be
maintained when possible. Increases speed by
minimizing required number of alignments. Currently
only "TRUE" is implemented.




Value
-------------------

a single-column `data.frame` corresponding to `clip.db`
and containing the best allele call from among the sequences
listed in `genotype_db`


Details
-------------------

In order to save time, initial gene assignments are preserved and
the allele calls are chosen from among those provided in `genotype_db`,
based on a simple alignment to the sample sequence.



Examples
-------------------

```R
# Load example data
data(germline_ighv)
data(sample_db)
data(genotype)
data(novel_df)

# Extract the database sequences that correspond to the genotype
genotype_seqs = genotypeFasta(genotype, germline_ighv, novel_df)

# Use the personlized genotype to determine corrected allele assignments
V_CALL_GENOTYPED = reassignAlleles(sample_db, genotype_seqs)
sample_db = cbind(sample_db, V_CALL_GENOTYPED)
```




