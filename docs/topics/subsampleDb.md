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
[38;5;246m# A tibble: 140 x 26[39m
   sequence_id       sequence           rev_comp productive v_call d_call j_call
   [3m[38;5;246m<chr>[39m[23m             [3m[38;5;246m<chr>[39m[23m              [3m[38;5;246m<lgl>[39m[23m    [3m[38;5;246m<lgl>[39m[23m      [3m[38;5;246m<chr>[39m[23m  [3m[38;5;246m<chr>[39m[23m  [3m[38;5;246m<chr>[39m[23m 
[38;5;250m 1[39m 088849_0682_2530… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD6… IGHJ4…
[38;5;250m 2[39m 101756_3039_2086… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD6… IGHJ6…
[38;5;250m 3[39m 020245_2263_3065… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ4…
[38;5;250m 4[39m 063394_2816_2054… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ5…
[38;5;250m 5[39m 110553_0911_1215… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD2… IGHJ4…
[38;5;250m 6[39m 213015_1355_0914… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ5…
[38;5;250m 7[39m 164458_1266_2351… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD5… IGHJ4…
[38;5;250m 8[39m 044962_0498_3654… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ1…
[38;5;250m 9[39m 072351_0586_3203… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD1… IGHJ6…
[38;5;250m10[39m 211468_3820_3687… NNNNNNNNNNNNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD1… IGHJ5…
[38;5;246m# … with 130 more rows, and 19 more variables: sequence_alignment <chr>,[39m
[38;5;246m#   germline_alignment <lgl>, junction <chr>, junction_aa <lgl>,[39m
[38;5;246m#   v_cigar <lgl>, d_cigar <lgl>, j_cigar <lgl>, vj_in_frame <lgl>,[39m
[38;5;246m#   stop_codon <lgl>, v_germline_end <lgl>, np1_length <dbl>,[39m
[38;5;246m#   np2_length <dbl>, j_germline_end <lgl>, junction_length <dbl>,[39m
[38;5;246m#   mutated_invariant <lgl>, indels <lgl>, d_5_trim <dbl>,[39m
[38;5;246m#   d_3_trim <dbl>, j_5_trim <dbl>[39m

```



See also
-------------------

[selectNovel](selectNovel.md)






