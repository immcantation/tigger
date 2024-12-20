**updateAlleleNames** - *Update IGHV allele names*

Description
--------------------

`updateAlleleNames` takes a set of IGHV allele calls and replaces any
outdated names (e.g. IGHV1-f) with the new IMGT names.


Usage
--------------------
```
updateAlleleNames(allele_calls)
```

Arguments
-------------------

allele_calls
:   vector of strings representing IGHV allele names.




Value
-------------------

Vector of strings representing updated IGHV allele names.


Note
-------------------

IMGT has removed `IGHV2-5*10` and `IGHV2-5*07` as it has determined they
are actually alleles `02` and `04`, respectively. The updated allele
names are based on IMGT release 2014-08-4.


References
-------------------


1.  Xochelli et al. (2014) Immunoglobulin heavy variable (IGHV) genes
and alleles: new entities, new names and implications for research and
prognostication in chronic lymphocytic leukaemia. Immunogenetics. 67(1):61-6




Examples
-------------------

```R
# Create a vector that uses old gene/allele names.
alleles <- c("IGHV1-c*01", "IGHV1-f*02", "IGHV2-5*07")

# Update the alleles to the new names
updateAlleleNames(alleles)

```


```
[1] "IGHV1-38-4*01" "IGHV1-69-2*02" "IGHV2-5*04"   

```



See also
-------------------

Like `updateAlleleNames`, [sortAlleles](sortAlleles.md) can help
format a list of allele names.






