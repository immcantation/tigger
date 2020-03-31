**plotNovel** - *Visualize evidence of novel V alleles*

Description
--------------------

`plotNovel` is be used to visualize the evidence of any novel V
alleles found using [findNovelAlleles](findNovelAlleles.md). It can also be used to
visualize the results for alleles that did


Usage
--------------------
```
plotNovel(
data,
novel_row,
v_call = "v_call",
j_call = "j_call",
seq = "sequence_alignment",
junction = "junction",
junction_length = "junction_length",
ncol = 1,
multiplot = TRUE
)
```

Arguments
-------------------

data
:   a `data.frame` in AIRR or Change-O format. See
[findNovelAlleles](findNovelAlleles.md) for details.

novel_row
:   a single row from a data frame as output by
[findNovelAlleles](findNovelAlleles.md) that contains a
polymorphism-containing germline allele

v_call
:   name of the column in `data` with V allele
calls. Default is `v_call`..

j_call
:   name of the column in `data` with J allele calls. 
Default is `j_call`.

seq
:   name of the column in `data` with the 
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is `sequence_alignment`.

junction
:   Junction region nucleotide sequence, which includes
the CDR3 and the two flanking conserved codons. Default
is `junction`.

junction_length
:   Number of junction nucleotides in the junction sequence.
Default is `junction_length`.

ncol
:   number of columns to use when laying out the plots

multiplot
:   Whether to return one single plot (`TRUE`) or a list 
with the three individual plots (`FALSE`).




Details
-------------------

The first panel in the plot shows, for all sequences which align to a particular 
germline allele, the mutation frequency at each postion along the aligned 
sequece as a function of the sequence-wide mutation. Sequences that pass 
the novel allele test are colored red, while sequences that don't pass
the test are colored yellow. The second panel shows the nucleotide usage at the 
positions as a function of sequence-wide mutation count.

To avoid cases where a clonal expansion might lead to a false positive, tigger examines
the combinations of J gene and junction length among sequences which perfectly 
match the proposed germline allele.



Examples
-------------------

```R
# Plot the evidence for the first (and only) novel allele in the example data
novel <- selectNovel(SampleNovel)
plotNovel(airrDb, novel[1, ], 
v_call="v_call", j_call="j_call", 
seq="sequence_alignment", 
junction="junction", junction_length="junction_length", multiplot=TRUE)
```

![2](plotNovel-2.png)







