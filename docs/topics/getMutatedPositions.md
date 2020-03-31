**getMutatedPositions** - *Find the location of mutations in a sequence*

Description
--------------------

`getMutatedPositions` takes two vectors of aligned sequences and
compares pairs of sequences. It returns a list of the nucleotide positions of
any differences.


Usage
--------------------
```
getMutatedPositions(
samples,
germlines,
ignored_regex = "[\\.N-]",
match_instead = FALSE
)
```

Arguments
-------------------

samples
:   a vector of strings respresenting aligned sequences

germlines
:   a vector of strings respresenting aligned sequences
to which `samples` will be compared. If only
one string is submitted, it will be used for all
`samples`.

ignored_regex
:   a regular expression indicating what characters
should be ignored (such as gaps and N nucleotides).

match_instead
:   if `TRUE`, the function returns the positions
that are the same instead of those that are
different.




Value
-------------------

A list of the nucleotide positions of any differences between the
input vectors.



Examples
-------------------

```R
# Create strings to act as a sample sequences and a reference sequence
seqs <- c("----GATA", "GAGAGAGA", "TANA")
ref <- "GATAGATA"

# Find the differences between the two
getMutatedPositions(seqs, ref)
```


```
[[1]]
integer(0)

[[2]]
[1] 3 7

[[3]]
[1] 1


```








