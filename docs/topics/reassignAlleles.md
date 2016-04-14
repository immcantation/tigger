





**reassignAlleles** - *Correct allele calls based on a personalized genotype*

Description
--------------------

`reassignAlleles` uses a subject-specific genotype to correct
correct preliminary allele assignments of a set of sequences derived
from a single subject.


Usage
--------------------
```
reassignAlleles(clip_db, genotype_db)
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

# Infer genotype from the sample
novel_df = findNovelAlleles(sample_db, germline_ighv)
geno = inferGenotype(sample_db, find_unmutated = TRUE,
germline_db = germline_ighv, novel_df = novel_df)

# Find the sequences that correspond to the genotype
genotype_seqs = genotypeFasta(geno, germline_ighv, novel_df)

# Use the personlized genotype to determine corrected allele assignments
V_CALL_GENOTYPED = reassignAlleles(sample_db, genotype_seqs)
sample_db = bind_cols(sample_db, V_CALL_GENOTYPED)
```




