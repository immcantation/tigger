**subsampleDb** - *Subsample repertoire*

Description
--------------------

`subsampleDb` will sample the same number of sequences for each gene, family
or allele (specified with `mode`) in `data`. Samples or subjects can
be subsampled independent by setting `group`.


Usage
--------------------
```
subsampleDb(
data,
gene = "v_call",
mode = c("gene", "allele", "family"),
min_n = 1,
max_n = NULL,
group = NULL
)
```

Arguments
-------------------

data
:   `data.frame` containing repertoire data.

gene
:   name of the column in `data` with allele calls. Default
is `v_call`.

mode
:   one of `c("gene", "family", "allele")` defining the degree of
specificity regarding allele calls when subsetting sequences.
Determines how `data` will be split into subsets from
which the same number of sequences will be subsampled. See
also `group`.

min_n
:   minimum number of observations to sample from each group. A group with
less observations than the minimum is excluded.

max_n
:   maximum number of observations to sample for all `mode` groups.
If `NULL`, it will be set automatically to the size of
the smallest group. If `max_n` is larger than the available
number of sequences for any `mode` group, it will be
automatically adjusted and the effective `max_n` used
will be the size of the smallest `mode` group.

group
:   columns containing additional grouping variables, e.g. sample_id.
These groups will be subsampled independently. If
`max_n` is `NULL`, a `max_n` will be
automatically set for each `group`.




Value
-------------------

Subsampled version of the input `data`.


Details
-------------------

`data` will be split into gene, allele or family subsets (`mode`) from
which the same number of sequences will be subsampled. If `mode=gene`,
for each gene in the field `gene` from `data`, a maximum of
`max_n` sequences will be subsampled. Input sequences
that have multiple gene calls (ties), can be subsampled from any of their calls,
but these duplicated samplings will be removed, and the final
subsampled `data` will contain unique rows.



Examples
-------------------

```R
subsampleDb(AIRRDb)

```


```
# A tibble: 140 × 26
   sequence_id                 sequence rev_comp productive v_call d_call j_call
   <chr>                       <chr>    <lgl>    <lgl>      <chr>  <chr>  <chr> 
 1 044192_2498_0295_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD1… IGHJ4…
 2 212810_1480_1640_length=45… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6…
 3 008603_0178_2137_length=45… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ5…
 4 031387_2356_0163_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD1… IGHJ6…
 5 194287_1422_2294_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6…
 6 159420_1086_3690_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD4… IGHJ4…
 7 031876_0245_1885_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ4…
 8 185805_1084_2905_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ6…
 9 184495_3403_0466_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ5…
10 067831_0282_1020_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ5…
# ℹ 130 more rows
# ℹ 19 more variables: sequence_alignment <chr>, germline_alignment <lgl>,
#   junction <chr>, junction_aa <lgl>, v_cigar <lgl>, d_cigar <lgl>,
#   j_cigar <lgl>, vj_in_frame <lgl>, stop_codon <lgl>, v_germline_end <lgl>,
#   np1_length <dbl>, np2_length <dbl>, j_germline_end <lgl>,
#   junction_length <dbl>, mutated_invariant <lgl>, indels <lgl>,
#   d_5_trim <dbl>, d_3_trim <dbl>, j_5_trim <dbl>

```



See also
-------------------

[selectNovel](selectNovel.md)






