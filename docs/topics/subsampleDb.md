**subsampleDb** - *Subsample repertoire*

Description
--------------------

`subsampleDb` will sample the same number of sequences for each gene, family
or allele (specified with `mode`) in `data`. Samples or subjects can
be subsampled indepently by setting `group`.


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
:   minimum number of observations to sample from each groupe. A group with 
less observations than the minimum is excluded.

max_n
:   maximum number of observations to sample for all `mode` groups.
If `NULL`, it will be set automatically to the size of 
the smallest group. If `max_n` is larger than the availabe 
number of sequences for any `mode` group, if will be 
automatically adjusted and the efective `max_n` used 
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
# A tibble: 139 × 26
   sequence_id                 sequence rev_comp productive v_call d_call j_call
   <chr>                       <chr>    <lgl>    <lgl>      <chr>  <chr>  <chr> 
 1 077294_2726_1483_length=42… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ5…
 2 240679_3698_3066_length=42… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ3…
 3 014026_2186_3677_length=45… NNNNNNN… NA       TRUE       IGHV1… IGHD4… IGHJ6…
 4 087866_2612_1813_length=47… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6…
 5 186255_1506_0138_length=45… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ6…
 6 094811_0577_1785_length=45… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6…
 7 200273_3493_1744_length=45… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ6…
 8 008329_2221_1510_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ6…
 9 098375_3215_3198_length=43… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ5…
10 011171_2123_0412_length=44… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ4…
# ℹ 129 more rows
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






