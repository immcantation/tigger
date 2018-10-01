**insertPolymorphisms** - *Insert polymorphisms into a nucleotide sequence*

Description
--------------------

`insertPolymorphisms` replaces nucleotides in the desired locations of a
provided sequence.


Usage
--------------------
```
insertPolymorphisms(sequence, positions, nucleotides)
```

Arguments
-------------------

sequence
:   starting nucletide sequence.

positions
:   numeric vector of positions which to be changed.

nucleotides
:   character vector of nucletides to which to change the
positions.




Value
-------------------

A sequence with the desired nucleotides in the provided locations.



Examples
-------------------

```R
insertPolymorphisms("HUGGED", c(1, 6, 2), c("T", "R", "I"))
```


```
[1] "TIGGER"

```




