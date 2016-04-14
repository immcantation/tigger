





**findNovelAlleles** - *Find novel alleles from repertoire sequencing data*

Description
--------------------

`findNovelAlleles` analyzes mutation patterns in sequences thought to
align to each germline allele in order to determine which positions
might be polymorphic.


Usage
--------------------
```
findNovelAlleles(clip_db, germline_db, germline_min = 200, nproc = 4,
min_seqs = 50, auto_mutrange = TRUE, mut_range = 1:10,
pos_range = 1:312, y_intercept = 0.125, alpha = 0.05, j_max = 0.15,
min_frac = 0.75)
```

Arguments
-------------------

clip_db
:   a `data.frame` in Change-O format. See details.

germline_db
:   a vector of named nucleotide germline sequences
matching the V calls in `clip_db`

germline_min
:   the minimum number of sequences that must have a
particular germline allele call for the allele to
be analyzed

nproc
:   the number of processors to use

min_seqs
:   the minimum number of total sequences (within the
desired mutational range and nucleotide range)
required for the samples to be considered

auto_mutrange
:   if `TRUE`, the algorithm will attempt to
determine the appropriate mutation range
automatically using the mutation count of the most
common sequence assigned to each allele analyzed

mut_range
:   the range of mutations that sampled may carry and
be considered by the algorithm

pos_range
:   the range of IMGT-numbered positions that should be
considered by the algorithm

y_intercept
:   the y-intercept above which positions should be
considered potentially polymorphic

alpha
:   the alpha cutoff to be used when constructing the
confidence interval for the y-intercept

j_max
:   the maximum fraction of sequences perfectly aligning
to a potential novel allele that are allowed to
utilize to a particular combination of junction
length and J gene

min_frac
:   the minimum fraction of sequences that must have
usable nucleotides in a given position for that
position to considered



Value
-------------------

a `data.frame` with a row for each known allele analyzed.
Besides metadata on the the parameters used in the search, each row will have
either a note as to where the polymorphism-finding algorithm exited or a
nucleotide sequence for the predicted novel allele.

Details
-------------------

A `data.frame` in Change-O format contains the following
columns:

+  `"SEQUENCE_IMGT"` containing the IMGT-gapped nucleotide sequence
+  `"V_CALL"` containing the IMGT/V-QUEST V allele call(s)
+  `"J_CALL"` containing the IMGT/V-QUEST J allele call(s)
+  `"JUNCTION_LENGTH"` containing the junction length

The TIgGER allele-finding algorithm, briefly, works as follows:
Mutations are determined through comparison to the provided germline.
Mutation frequency at each *position* is determined as a function of
*sequence-wide* mutation counts. Polymorphic positions exhibit a high
mutation frequency despite sequence-wide mutation count. False positive of
potential novel alleles resulting from clonally-related sequences are guarded
against by ensuring that sequences perfectly matching the potential novel
allele utilize a wide range of combinations of J gene and junction length.



Examples
-------------------

```R
# Load example data and germlines
data(sample_db)
data(germline_ighv)

# Find novel alleles and return relevant data
novel_df = findNovelAlleles(sample_db, germline_ighv)
```



See also
-------------------

[plotNovel](plotNovel.md) to visualize the data supporting any
novel alleles hypothesized to be present in the data and
[inferGenotype](inferGenotype.md) to determine if the novel alleles are frequent
enought to be included in the subject's genotype



