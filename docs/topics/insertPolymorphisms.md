





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
:   the starting nucletide sequence

positions
:   a vector of positions which to be changed

nucleotides
:   a vector of nucletides to which to change the
positions




Value
-------------------

a sequence with the desired nucleotides in provided locations



Examples
-------------------

```R
insertPolymorphisms("hugged", c(1,6,2), c("t","r","i"))
```


```
[1] "tigger"

```




