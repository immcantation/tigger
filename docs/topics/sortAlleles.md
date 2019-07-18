**sortAlleles** - *Sort allele names*

Description
--------------------

`sortAlleles` returns a sorted vector of strings respresenting Ig allele
names. Names are first sorted by gene family, then by gene, then by allele.
Duplicated genes have their alleles are sorted as if they were part of their
non-duplicated counterparts (e.g. `IGHV1-69D*01` comes after `IGHV1-69*01` 
but before `IGHV1-69*02`), and non-localized genes (e.g. `IGHV1-NL1*01`) 
come last within their gene family.


Usage
--------------------
```
sortAlleles(allele_calls, method = c("name", "position"))
```

Arguments
-------------------

allele_calls
:   a vector of strings respresenting Ig allele names.

method
:   a string defining the method to use when sorting alleles.
If `"name"` then sort in lexicographic order. If
`"position"` then sort by position in the locus, as
determined by the final two numbers in the gene name.




Value
-------------------

A sorted vector of strings respresenting Ig allele names.



Examples
-------------------

```R
# Create a list of allele names
alleles <- c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
"IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
"IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")

# Sort the alleles by name
sortAlleles(alleles)

```


```
 [1] "IGHV1-2"               "IGHV1-2*01"            "IGHV1-2*01,IGHV1-2*05" "IGHV1-2*02"           
 [5] "IGHV1-69*01"           "IGHV1-69D*01"          "IGHV1-69*02"           "IGHV1-69-2*01"        
 [9] "IGHV1-NL1*01"          "IGHV2-5*01"           

```


```R

# Sort the alleles by position in the locus
sortAlleles(alleles, method="pos")
```


```
 [1] "IGHV1-NL1*01"          "IGHV1-69-2*01"         "IGHV1-69*02"           "IGHV1-69*01"          
 [5] "IGHV1-69D*01"          "IGHV2-5*01"            "IGHV1-2*02"            "IGHV1-2*01"           
 [9] "IGHV1-2*01,IGHV1-2*05" "IGHV1-2"              

```



See also
-------------------

Like `sortAlleles`, [updateAlleleNames](updateAlleleNames.md) can help
format a list of allele names.



