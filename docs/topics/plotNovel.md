





**plotNovel** - *Visualize evidence of novel V alleles*

Description
--------------------

`plotNovel` is be used to visualize the evidence of any novel V
alleles found using [findNovelAlleles](findNovelAlleles.md).


Usage
--------------------
```
plotNovel(clip_db, novel_df_row, ncol = 1)
```

Arguments
-------------------

clip_db
:   a `data.frame` in Change-O format. See
[findNovelAlleles](findNovelAlleles.md) for details.

novel_df_row
:   a single row from a data frame as output by
[findNovelAlleles](findNovelAlleles.md) that contains a
polymorphism-containing germline allele

ncol
:   number of columns to use when laying out the plots





Examples
-------------------

```R
# Load example data and germlines
data(sample_db)
data(germline_ighv)

# Find novel alleles and return relevant data
novel_df = findNovelAlleles(sample_db, germline_ighv)
# Plot the evidence for the first (and only) novel allele in the example data
novel = selectNovel(novel_df)
plotNovel(sample_db, novel[1,])
```

![2](plotNovel-2.png)



