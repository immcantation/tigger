**genotypeFasta** - *Return the nucleotide sequences of a genotype*

Description
--------------------

`genotypeFasta` converts a genotype table into a vector of nucleotide
sequences.


Usage
--------------------
```
genotypeFasta(genotype, germline_db, novel = NA)
```

Arguments
-------------------

genotype
:   a `data.frame` of alleles denoting a genotype, 
as returned by [inferGenotype](inferGenotype.md).

germline_db
:   a vector of named nucleotide germline sequences
matching the alleles detailed in `genotype`.

novel
:   an optional `data.frame` containing putative
novel alleeles of the type returned by
[findNovelAlleles](findNovelAlleles.md).




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



