





**genotypeFasta** - *Return the nucleotide sequences of a genotype*

Description
--------------------

`genotypeFasta` converts a genotype table into a vector of nucleotide
sequences.


Usage
--------------------
```
genotypeFasta(genotype, germline_db, novel_df = NA)
```

Arguments
-------------------

genotype
:   a table of alleles denoting a genotype, as returned by
[inferGenotype](inferGenotype.md)

germline_db
:   a vector of named nucleotide germline sequences
matching the alleles detailed in `genotype`

novel_df
:   an optional `data.frame` containing putative
novel alleeles of the type returned by
[findNovelAlleles](findNovelAlleles.md)




Value
-------------------

A named vector of strings containing the germline nucleotide
sequences of the alleles in the provided genotype



Examples
-------------------

```R
# Load example data
data(germline_ighv)
data(novel_df)
data(genotype)

# Find the sequences that correspond to the genotype
genotype_seqs = genotypeFasta(genotype, germline_ighv, novel_df)
```



See also
-------------------

[inferGenotype](inferGenotype.md)



