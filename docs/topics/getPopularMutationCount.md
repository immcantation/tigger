**getPopularMutationCount** - *Find mutation counts for frequency sequences*

Description
--------------------

`getPopularMutationCount` determines which sequences occur frequently
for each V gene and returns the mutation count of those sequences.


Usage
--------------------
```
getPopularMutationCount(
data,
germline_db,
v_call = "v_call",
seq = "sequence_alignment",
gene_min = 0.001,
seq_min = 50,
seq_p_of_max = 1/8,
full_return = FALSE
)
```

Arguments
-------------------

data
:   `data.frame` in the Change-O format. See
[findNovelAlleles](findNovelAlleles.md) for a list of required
columns.

germline_db
:   named list of IMGT-gapped germline sequences.

v_call
:   name of the column in `data` with V allele calls. 
Default is `v_call`.

seq
:   name of the column in `data` with the 
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is `sequence_alignment`.

gene_min
:   portion of all unique sequences a gene must
constitute to avoid exclusion.

seq_min
:   number of copies of the V that must be present for
to avoid exclusion.

seq_p_of_max
:   ror each gene, the fraction of the most common V sequence
count that a sequence must meet to avoid exclusion.

full_return
:   if `TRUE`, will return all `data` columns and
will include sequences with mutation count < 1.




Value
-------------------

A data frame of genes that have a frequent sequence mutation count
above 1.



Examples
-------------------

```R
getPopularMutationCount(AIRRDb, SampleGermlineIGHV)
```


```
# A tibble: 1 x 2
  v_gene  mutation_count
  <chr>            <int>
1 IGHV1-8              1

```



See also
-------------------

[getMutatedPositions](getMutatedPositions.md) can be used to find which positions
of a set of sequences are mutated.






