**writeFasta** - *Write to a fasta file*

Description
--------------------

`writeFasta` writes a named vector of sequences to a file in fasta
format.


Usage
--------------------
```
writeFasta(named_sequences, file, width = 60, append = FALSE)
```

Arguments
-------------------

named_sequences
:   vector of named string representing sequences

file
:   the name of the output file.

width
:   the number of characters to be printed per line.
if not between 1 and 255, width with be infinite.

append
:   `logical` indicating if the output should be
appended to `file` instead of overwriting it




Value
-------------------

A named vector of strings respresenting Ig alleles.



Examples
-------------------

```R
### Not run:
# writeFasta(germlines, "ighv.fasta")
```



See also
-------------------

[readIgFasta](readIgFasta.md) to do the inverse.






