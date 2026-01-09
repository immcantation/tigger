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
 1 204043_1358_2442_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ3…
 2 171443_1072_2674_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ4…
 3 247875_1598_0850_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6…
 4 123799_0962_2523_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ3…
 5 116794_1023_0376_length=46… NNNNNNN… NA       TRUE       IGHV1… IGHD1… IGHJ6…
 6 052525_0411_0611_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD4… IGHJ3…
 7 108481_3020_1175_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ1…
 8 190412_1510_1034_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ5…
 9 198447_1406_3616_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD4… IGHJ6…
10 188582_1299_0211_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ6…
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






