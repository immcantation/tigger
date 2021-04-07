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
[90m# A tibble: 140 x 26[39m
   sequence_id sequence rev_comp productive v_call d_call j_call
   [3m[90m<chr>[39m[23m       [3m[90m<chr>[39m[23m    [3m[90m<lgl>[39m[23m    [3m[90m<lgl>[39m[23m      [3m[90m<chr>[39m[23m  [3m[90m<chr>[39m[23m  [3m[90m<chr>[39m[23m 
[90m 1[39m 251112_188… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD6… IGHJ4…
[90m 2[39m 062828_048… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD1… IGHJ3…
[90m 3[39m 133739_332… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD6… IGHJ6…
[90m 4[39m 066458_268… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD2… IGHJ6…
[90m 5[39m 129724_305… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD6… IGHJ6…
[90m 6[39m 199692_345… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ3…
[90m 7[39m 182790_352… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ5…
[90m 8[39m 204947_141… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD2… IGHJ4…
[90m 9[39m 019164_015… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD3… IGHJ4…
[90m10[39m 256485_170… NNNNNNN… [31mNA[39m       TRUE       IGHV1… IGHD2… IGHJ1…
[90m# … with 130 more rows, and 19 more variables: sequence_alignment [3m[90m<chr>[90m[23m,[39m
[90m#   germline_alignment [3m[90m<lgl>[90m[23m, junction [3m[90m<chr>[90m[23m, junction_aa [3m[90m<lgl>[90m[23m, v_cigar [3m[90m<lgl>[90m[23m,[39m
[90m#   d_cigar [3m[90m<lgl>[90m[23m, j_cigar [3m[90m<lgl>[90m[23m, vj_in_frame [3m[90m<lgl>[90m[23m, stop_codon [3m[90m<lgl>[90m[23m,[39m
[90m#   v_germline_end [3m[90m<lgl>[90m[23m, np1_length [3m[90m<dbl>[90m[23m, np2_length [3m[90m<dbl>[90m[23m,[39m
[90m#   j_germline_end [3m[90m<lgl>[90m[23m, junction_length [3m[90m<dbl>[90m[23m, mutated_invariant [3m[90m<lgl>[90m[23m,[39m
[90m#   indels [3m[90m<lgl>[90m[23m, d_5_trim [3m[90m<dbl>[90m[23m, d_3_trim [3m[90m<dbl>[90m[23m, j_5_trim [3m[90m<dbl>[90m[23m[39m

```



See also
-------------------

[selectNovel](selectNovel.md)






