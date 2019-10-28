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
:   fasta-formatted file of immunoglobuling sequences.

strip_down_name
:   if `TRUE`, will extract only the allele name
from the strings fasta file's sequence names.

force_caps
:   if `TRUE`, will force nucleotides to
uppercase.




Value
-------------------

Named vector of strings respresenting Ig alleles.




See also
-------------------

[writeFasta](writeFasta.md) to do the inverse.






