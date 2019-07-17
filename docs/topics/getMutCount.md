**getMutCount** - *Determine the mutation counts from allele calls*

Description
--------------------

`getMutCount` takes a set of nucleotide sequences and their allele calls
and determines the distance between that seqeunce and any germline alleles
contained within the call


Usage
--------------------
```
getMutCount(samples, allele_calls, germline_db)
```

Arguments
-------------------

samples
:   a vector of IMGT-gapped sample V sequences

allele_calls
:   a vector of strings respresenting Ig allele calls for
the sequences in `samples`, where multiple
calls are separated by a comma

germline_db
:   a vector of named nucleotide germline sequences
matching the calls detailed in `allele_calls`




Value
-------------------

A list equal in length to `samples`, containing the Hamming
distance to each germline allele contained within each call within
each element of `samples`



Examples
-------------------

```R
# Insert a mutation into a germline sequence
s2 <- s3 <- SampleGermlineIGHV[1]
stringi::stri_sub(s2, 103, 103) <- "G"
stringi::stri_sub(s3, 107, 107) <- "C"

sample_seqs <- c(SampleGermlineIGHV[2], s2, s3)

# Pretend that one sample sequence has received an ambiguous allele call
sample_alleles <- c(paste(names(SampleGermlineIGHV[1:2]), collapse=","),
names(SampleGermlineIGHV[2]),
names(SampleGermlineIGHV[1]))

# Compare each sequence to its assigned germline(s) to determine the distance
getMutCount(sample_seqs, sample_alleles, SampleGermlineIGHV)
```


```
[[1]]
[[1]][[1]]
[1] 1

[[1]][[2]]
[1] 0


[[2]]
[1] 2

[[3]]
[1] 1


```




