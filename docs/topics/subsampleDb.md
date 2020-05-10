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
:   a `data.frame` in Change-O format.

gene
:   name of the column in `data` with allele calls. Default
is `v_call`.

mode
:   one of c("gene", "family", "allele") defining the degree of
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

A `data.frame`, subsampled from `data`.


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
# A tibble: 140 x 26
   sequence_id sequence rev_comp productive v_call d_call j_call sequence_alignm… germline_alignm… junction junction_aa v_cigar d_cigar j_cigar
   <chr>       <chr>    <lgl>    <lgl>      <chr>  <chr>  <chr>  <chr>            <lgl>            <chr>    <lgl>       <lgl>   <lgl>   <lgl>  
 1 189598_135… NNNNNNN… NA       TRUE       IGHV1… IGHD4… IGHJ6… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
 2 084252_280… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ6… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
 3 011265_017… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ4… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
 4 180277_351… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ4… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
 5 029365_020… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ6… CAGGTCCAGCTGGTA… NA               TGTGCAA… NA          NA      NA      NA     
 6 122489_102… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ4… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
 7 187939_128… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6… CAGGTCCAGCTGGTA… NA               TGTGCAA… NA          NA      NA      NA     
 8 175428_336… NNNNNNN… NA       TRUE       IGHV1… IGHD2… IGHJ3… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
 9 121034_089… NNNNNNN… NA       TRUE       IGHV1… IGHD6… IGHJ3… CAGGTCCGGCTGGAA… NA               TGTGCAA… NA          NA      NA      NA     
10 096429_076… NNNNNNN… NA       TRUE       IGHV1… IGHD3… IGHJ6… CAGGTGCAGCTGGTG… NA               TGTGCGA… NA          NA      NA      NA     
# … with 130 more rows, and 12 more variables: vj_in_frame <lgl>, stop_codon <lgl>, v_germline_end <lgl>, np1_length <dbl>, np2_length <dbl>,
#   j_germline_end <lgl>, junction_length <dbl>, mutated_invariant <lgl>, indels <lgl>, d_5_trim <dbl>, d_3_trim <dbl>, j_5_trim <dbl>

```



See also
-------------------

[selectNovel](selectNovel.md)






