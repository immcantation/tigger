**inferGenotype** - *Infer a subject-specific genotype using a frequency method*

Description
--------------------

`inferGenotype` infers an subject's genotype using a frequency method.
The genotype is inferred by finding the minimum number set of alleles that 
can explain the majority of each gene's calls. The most common allele of 
each gene is included in the genotype first, and the next most common allele 
is added until the desired fraction of alleles can be explained. In this 
way, mistaken allele calls (resulting from sequences which
by chance have been mutated to look like another allele) can be removed.


Usage
--------------------
```
inferGenotype(data_db, v_call = "V_CALL", fraction_to_explain = 0.875,
gene_cutoff = 1e-04, find_unmutated = TRUE, germline_db = NA,
novel_df = NA)
```

Arguments
-------------------

data_db
:   a `data.frame` containing V allele
calls from a single subject. If
`find_unmutated` is `TRUE`, then
the sample IMGT-gapped V(D)J sequence should

v_call
:   column in `data_db` with V allele calls.
Default is `"V_CALL"`.                            
be provided in a column `"SEQUENCE_IMGT"`

fraction_to_explain
:   the portion of each gene that must be
explained by the alleles that will be included
in the genotype.

gene_cutoff
:   either a number of sequences or a fraction of
the length of `allele_calls` denoting the
minimum number of times a gene must be
observed in `allele_calls` to be included
in the genotype.

find_unmutated
:   if `TRUE`, use `germline_db` to
find which samples are unmutated. Not needed
if `allele_calls` only represent
unmutated samples.

germline_db
:   named vector of sequences containing the
germline sequences named in
`allele_calls`. Only required if
`find_unmutated` is `TRUE`.

novel_df
:   an optional `data.frame` of the type
novel returned by
[findNovelAlleles](findNovelAlleles.md) containing
germline sequences that will be utilized if
`find_unmutated` is `TRUE`. See
Details.




Value
-------------------

A `data.frame` of alleles denoting the genotype of the subject containing 
the following columns:


+  `GENE`: The gene name without allele.
+  `ALLELES`: Comma separated list of alleles for the given `GENE`.
+  `COUNTS`: Comma separated list of observed sequences for each 
corresponding allele in the `ALLELES` list.
+  `TOTAL`: The total count of observed sequences for the given `GENE`.
+  `NOTE`: Any comments on the inferrence.



Details
-------------------

Allele calls representing cases where multiple alleles have been
assigned to a single sample sequence are rare among unmutated
sequences but may result if nucleotides for certain positions are
not available. Calls containing multiple alleles are treated as
belonging to all groups. If `novel_df` is provided, all
sequences that are assigned to the same starting allele as any
novel germline allele will have the novel germline allele appended
to their assignent prior to searching for unmutated sequences.


Note
-------------------

This method works best with data derived from blood, where a large
portion of sequences are expected to be unmutated. Ideally, there
should be hundreds of allele calls per gene in the input.



Examples
-------------------

```R
# Infer IGHV genotype, using only unmutated sequences, including novel alleles
inferGenotype(SampleDb, find_unmutated=TRUE, germline_db=GermlineIGHV,
novel_df=SampleNovel)
```


```
        GENE     ALLELES      COUNTS TOTAL NOTE
1    IGHV1-2       02,04     664,302   966     
2    IGHV1-3          01         226   226     
3    IGHV1-8 01,02_G234T     467,370   837     
4   IGHV1-18          01        1005  1005     
5   IGHV1-24          01         105   105     
6   IGHV1-46          01         624   624     
7   IGHV1-58       01,02       23,18    41     
8   IGHV1-69    01,04,06 515,469,280  1279     
9 IGHV1-69-2          01          31    31     

```



See also
-------------------

[plotGenotype](plotGenotype.md) for a colorful visualization and
[genotypeFasta](genotypeFasta.md) to convert the genotype to nucleotide sequences.
See [inferGenotypeBayesian](inferGenotypeBayesian.md) to infer a subject-specific genotype 
using a Bayesian approach.



