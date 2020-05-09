**findUnmutatedCalls** - *Determine which calls represent an unmutated allele*

Description
--------------------

`findUnmutatedCalls` determines which allele calls would represent a 
perfect match with the germline sequence, given a vector of allele calls and
mutation counts. In the case of multiple alleles being assigned to a
sequence, only the subset that would represent a perfect match is returned.


Usage
--------------------
```
findUnmutatedCalls(allele_calls, sample_seqs, germline_db)
```

Arguments
-------------------

allele_calls
:   a vector of strings respresenting Ig allele calls,
where multiple calls are separated by a comma.

sample_seqs
:   V(D)J-rearranged sample sequences matching the order
of the given `allele_calls`.

germline_db
:   a vector of named nucleotide germline sequences




Value
-------------------

A vector of strings containing the members of `allele_calls`
that represent unmutated sequences.



Examples
-------------------

```R
# Find which of the sample alleles are unmutated
calls <- findUnmutatedCalls(AIRRDb$v_call, AIRRDb$sequence_alignment, 
germline_db=SampleGermlineIGHV)
```

**Error in gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), **: object 'AIRRDb' not found






