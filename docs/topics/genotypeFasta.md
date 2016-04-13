





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
data(sample_db)

# Infer and view a genotype from the sample
novel_df = findNovelAlleles(sample_db, germline_ighv)
geno = inferGenotype(sample_db, find_unmutated = TRUE,
germline_db = germline_ighv, novel_df = novel_df)
print(geno)

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


```R

# Find the sequences that correspond to the genotype
genotype_seqs = genotypeFasta(geno, germline_ighv, novel_df)
```



See also
-------------------

[inferGenotype](inferGenotype.md)



