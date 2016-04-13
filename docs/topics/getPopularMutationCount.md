





**getPopularMutationCount** - *Find Frequent Sequences' Mutation Counts*

Description
--------------------

`getPopularMutationCount` determines which sequences occur frequently
for each V gene and returns the mutation count of those sequences.


Usage
--------------------
```
getPopularMutationCount(sample_db, germline_db, gene_min = 0.001,
seq_min = 50, seq_p_of_max = 1/8, full_return = FALSE)
```

Arguments
-------------------

sample_db
:   A Change-O db data frame. See
[findNovelAlleles](findNovelAlleles.md) for a list of required
columns.

germline_db
:   A named list of IMGT-gapped germline sequences.

gene_min
:   The portion of all unique sequences a gene must
constitute to avoid exclusion.

seq_min
:   The number of copies of the V that must be present for
to avoid exclusion.

seq_p_of_max
:   For each gene, fraction of the most common V sequence's
count that a sequence must meet to avoid exclusion.

full_return
:   If true, will return all `sample_db` columns and
will include sequences with mutation count < 1.



Value
-------------------

A data frame of genes that have a frequent sequence mutation count
above 1.



Examples
-------------------

```R
data(sample_db, germline_ighv)
getPopularMutationCount(sample_db, germline_ighv)
```


```
Source: local data frame [1 x 2]

   V_GENE MUTATION_COUNT
    (chr)          (int)
1 IGHV1-8              1

```



See also
-------------------

[getMutatedPositions](getMutatedPositions.md) can be used to find which positions
of a set of sequences are mutated.



