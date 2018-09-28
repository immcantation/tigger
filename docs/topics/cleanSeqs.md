**cleanSeqs** - *Clean up nucleotide sequences*

Description
--------------------

`cleanSeqs` capitalizes nucleotides, replaces "." with "-", and then
replaces all characters besides ACGT- with "N".


Usage
--------------------
```
cleanSeqs(seqs)
```

Arguments
-------------------

seqs
:   a vector of nucleotide sequences




Value
-------------------

A vector of nucleotide sequences



Examples
-------------------

```R
# Create messy nucleotide sequences
seqs <- c("AGAT.taa-GAG...ATA",
"GATACAGTXXXXXAGNNNPPPACA")
# Clean them up
cleanSeqs(seqs)
```


```
[1] "AGAT-TAA-GAG---ATA"       "GATACAGTNNNNNAGNNNNNNACA"

```



See also
-------------------

[sortAlleles](sortAlleles.md) and [updateAlleleNames](updateAlleleNames.md) can
help format a list of allele names.



