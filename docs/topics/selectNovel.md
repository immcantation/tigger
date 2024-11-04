**selectNovel** - *Select rows containing novel alleles*

Description
--------------------

`selectNovel` takes the result from [findNovelAlleles](findNovelAlleles.md) and
selects only the rows containing unique, novel alleles.


Usage
--------------------
```
selectNovel(novel, keep_alleles = FALSE)
```

Arguments
-------------------

novel
:   `data.frame` of the type returned by
[findNovelAlleles](findNovelAlleles.md).

keep_alleles
:   `logical` indicating if different alleles
leading to the same novel sequence should be kept.
See Details.




Value
-------------------

A `data.frame` containing only unique, novel alleles (if any)
that were in the input.


Details
-------------------

If, for instance, subject has in his genome `IGHV1-2*02` and a novel 
allele equally close to `IGHV1-2*02` and `IGHV1-2*05`, the novel allele may be
detected by analyzing sequences that best align to either of these alleles.
If `keep_alleles` is `TRUE`, both polymorphic allele calls will
be retained. In the case that multiple mutation ranges are checked for the
same allele, only one mutation range will be kept in the output.



Examples
-------------------

```R
novel <- selectNovel(SampleNovel)

```








