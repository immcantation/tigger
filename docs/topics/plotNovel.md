**plotNovel** - *Visualize evidence of novel V alleles*

Description
--------------------

`plotNovel` is be used to visualize the evidence of any novel V
alleles found using [findNovelAlleles](findNovelAlleles.md). It can also be used to
visualize the results for alleles that did


Usage
--------------------
```
plotNovel(clip_db, novel_df_row, ncol = 1, v_call = "V_CALL")
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

v_call
:   name of the column in `clip_db` with V allele
calls. Default is "V_CALL"




Details
-------------------

The first panel in the plot shows, for all sequences which align to a particular 
germline allele, the mutation frequency at each postion along the aligned 
sequece as a function of the sequence-wide mutation. Sequences that pass 
the novel allele test are colored red, while sequences that don't pass
the test are colored yellow.

The second panel shows the nucleotide usage at the positions 
as a function of sequence-wide mutation count. 

To avoid cases where a clonal expansion might lead to a false positive, tigger examines
the combinations of J gene and junction length among sequences which perfectly match the proposed
germline allele.



Examples
-------------------

```R
# Plot the evidence for the first (and only) novel allele in the example data
novel <- selectNovel(SampleNovel)
plotNovel(SampleDb, novel[1,])
```

![2](plotNovel-2.png)



