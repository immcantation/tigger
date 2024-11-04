**cleanSeqs** - *Clean up nucleotide sequences*

Description
--------------------

`cleanSeqs` capitalizes nucleotides and replaces all characters 
besides `c("A", "C", "G", "T", "-", ".")` with `"N"`.


Usage
--------------------
```
cleanSeqs(seqs)
```

Arguments
-------------------

seqs
:   vector of nucleotide sequences.




Value
-------------------

A modified vector of nucleotide sequences.



Examples
-------------------

```R
# Clean messy nucleotide sequences
seqs <- c("AGAT.taa-GAG...ATA", "GATACAGTXXZZAGNNPPACA")
cleanSeqs(seqs)

```


```
[1] "AGAT.TAA-GAG...ATA"    "GATACAGTNNNNAGNNNNACA"

```



See also
-------------------

[sortAlleles](sortAlleles.md) and [updateAlleleNames](updateAlleleNames.md) can
help format a list of allele names.






