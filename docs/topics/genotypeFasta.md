**genotypeFasta** - *Return the nucleotide sequences of a genotype*

Description
--------------------

`genotypeFasta` converts a genotype table into a vector of nucleotide
sequences.


Usage
--------------------
```
genotypeFasta(
genotype,
germline_db,
novel = NA,
include_unseen = FALSE,
strip_d = TRUE
)
```

Arguments
-------------------

genotype
:   `data.frame` of alleles denoting a genotype,
as returned by [inferGenotype](inferGenotype.md).

germline_db
:   vector of named nucleotide germline sequences
matching the alleles detailed in `genotype`.

novel
:   an optional `data.frame` containing putative
novel alleles of the type returned by
[findNovelAlleles](findNovelAlleles.md).

include_unseen
:   if `TRUE`, include germline database alleles for
genes that are not present in `genotype`. For
genes present in `genotype`, include only the
genotyped alleles.

strip_d
:   if `TRUE` (default) remove the "D" from the end of
gene annotations (denoting a duplicate gene in the locus)
when matching genotype alleles to `germline_db`. If
`FALSE`, alleles are matched exactly.




Value
-------------------

A named vector of strings containing the germline nucleotide
sequences of the alleles in the provided genotype.



Examples
-------------------

```R
# Find the sequences that correspond to the genotype
genotype_db <- genotypeFasta(SampleGenotype, SampleGermlineIGHV, SampleNovel)

```



See also
-------------------

[inferGenotype](inferGenotype.md)






