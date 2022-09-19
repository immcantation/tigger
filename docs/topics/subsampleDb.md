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
# A tibble: 140 × 26
   sequen…¹ seque…² rev_c…³ produ…⁴ v_call d_call j_call seque…⁵ germl…⁶ junct…⁷
   <chr>    <chr>   <lgl>   <lgl>   <chr>  <chr>  <chr>  <chr>   <lgl>   <chr>  
 1 145605_… NNNNNN… NA      TRUE    IGHV1… IGHD3… IGHJ3… CAGGTG… NA      TGTGCG…
 2 183698_… NNNNNN… NA      TRUE    IGHV1… IGHD3… IGHJ3… CAGGTG… NA      TGTGCG…
 3 227123_… NNNNNN… NA      TRUE    IGHV1… IGHD3… IGHJ5… CAGGTC… NA      TGTGCA…
 4 096429_… NNNNNN… NA      TRUE    IGHV1… IGHD3… IGHJ6… CAGGTG… NA      TGTGCG…
 5 127175_… NNNNNN… NA      TRUE    IGHV1… IGHD3… IGHJ6… CAGGTG… NA      TGTGCG…
 6 189175_… NNNNNN… NA      TRUE    IGHV1… IGHD1… IGHJ6… CAGGTC… NA      TGTGCA…
 7 179957_… NNNNNN… NA      TRUE    IGHV1… IGHD4… IGHJ4… CAGGTG… NA      TGTGCG…
 8 196538_… NNNNNN… NA      TRUE    IGHV1… IGHD1… IGHJ5… CAGGTG… NA      TGTGCG…
 9 071543_… NNNNNN… NA      TRUE    IGHV1… IGHD5… IGHJ6… CAGGTG… NA      TGTGCG…
10 075839_… NNNNNN… NA      TRUE    IGHV1… IGHD2… IGHJ4… CAGGTG… NA      TGTGCG…
# … with 130 more rows, 16 more variables: junction_aa <lgl>, v_cigar <lgl>,
#   d_cigar <lgl>, j_cigar <lgl>, vj_in_frame <lgl>, stop_codon <lgl>,
#   v_germline_end <lgl>, np1_length <dbl>, np2_length <dbl>,
#   j_germline_end <lgl>, junction_length <dbl>, mutated_invariant <lgl>,
#   indels <lgl>, d_5_trim <dbl>, d_3_trim <dbl>, j_5_trim <dbl>, and
#   abbreviated variable names ¹​sequence_id, ²​sequence, ³​rev_comp, ⁴​productive,
#   ⁵​sequence_alignment, ⁶​germline_alignment, ⁷​junction

```



See also
-------------------

[selectNovel](selectNovel.md)






