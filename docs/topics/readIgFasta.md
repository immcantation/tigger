**readIgFasta** - *Read immunoglobulin sequences*

Description
--------------------

`readIgFasta` reads a fasta-formatted file of immunoglobulin (Ig)
sequences and returns a named vector of those sequences.


Usage
--------------------
```
readIgFasta(fasta_file, strip_down_name = TRUE, force_caps = TRUE)
```

Arguments
-------------------

fasta_file
:   fasta-formatted file of immunoglobulin sequences.

strip_down_name
:   if `TRUE`, will extract only the allele name
from the strings fasta file's sequence names.

force_caps
:   if `TRUE`, will force nucleotides to
uppercase.




Value
-------------------

Named vector of strings representing Ig alleles.



Examples
-------------------

```R
### Not run:
# germlines <- readIgFasta("ighv.fasta")

```



See also
-------------------

[writeFasta](writeFasta.md) to do the inverse.






