





**selectNovel** - *Select rows containing novel alleles*

Description
--------------------

`selectNovel` takes the result from [findNovelAlleles](findNovelAlleles.md) and
selects only the rows containing unique, novel alleles.


Usage
--------------------
```
selectNovel(novel_df, keep_alleles = FALSE)
```

Arguments
-------------------

novel_df
:   A `data.frame` of the type returned by
[findNovelAlleles](findNovelAlleles.md)

keep_alleles
:   A `logical` indicating if different alleles
leading to the same novel sequence should be kept.
See details.




Value
-------------------

A `data.frame` containing only unique, novel alleles (if any)
that were in the input.


Details
-------------------

If, for instance, subject has in his genome IGHV1-2*02 and a novel 
allele equally close to IGHV1-2*02 and IGHV1-2*05, the novel allele may be
detected by analyzing sequences that best align to either of these alleles.
If `keep_alleles` is `TRUE`, both polymorphic allele calls will
be retained. In the case that multiple mutation ranges are checked for the
same allele, only one mutation range will be kept in the output.



Examples
-------------------

```R
data(novel_df)
novel = selectNovel(novel_df)
```




