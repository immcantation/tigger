# TIgGER #

High-throughput sequencing of B cell immunoglobulin receptors is providing unprecedented insight into adaptive immunity. A key step in analyzing these data involves assignment of the germline V, D and J gene segment alleles that comprise each immunoglobulin sequence by matching them against a database of known V(D)J alleles. However, this process will fail for sequences that utilize previously undetected alleles, whose frequency in the population is unclear.

**TIgGER is a computational method that significantly improves V(D)J allele assignments by first determining the complete set of gene segments carried by an individual (including novel alleles) from V(D)J-rearrange sequences. TIgGER can then infer a subject's genotype from these sequences, and use this genotype to correct the initial V(D)J allele assignments.**

The application of TIgGER continues to identify a surprisingly high frequency of novel alleles in humans, highlighting the critical need for this approach. (TIgGER, however, can and has been used with data from other species.)

## Core Abilities ##

* Detecting novel alleles
* Inferring a subject's genotype
* Correcting preliminary allele calls

## Required Input ##

* A table of sequences from a single individual, with columns containing the following:
    * V(D)J-rearranged nucleotide sequence (in IMGT-gapped format)
    * Preliminary V allele calls
    * Preliminary J allele calls
    * Length of the junction region
* Germline Ig sequences in IMGT-gapped fasta format (e.g., as those downloaded from [IMGT/GENE-DB](http://www.imgt.org/genedb))

The former can be created through the use of [IMGT/HighV-QUEST](http://www.imgt.org) and [Change-O](http://changeo.readthedocs.io).

## Contact ##

For help, questions, or suggestions, please contact [Daniel Gadala-Maria](mailto:daniel.gadala-maria@yale.edu) or use the [issue tracker](http://bitbucket.org/kleinstein/tigger/issues).

## Related Projects ##

* [pRESTO](http://presto.readthedocs.io) - 
  Raw read assembly, quality control and UMI processing 
* [Change-O](http://changeo.readthedocs.io) - 
  V(D)J alignment standardization and clonal clustering
* [Alakazam](http://alakazam.readthedocs.io) - 
  Lineage reconstruction, V(D)J gene usage, repertoire diversity and 
  physicochemical property analysis
* [SHazaM](http://shazam.readthedocs.io) - 
  Mutation profiling and selection strength quantification

Dependencies
---------------

**Depends:** alakazam, shazam, ggplot2  
**Imports:** tidyr, dplyr, doParallel, foreach, graphics, grid, iterators, parallel, stats  
**Suggests:** knitr

Authors
---------------

[Daniel Gadala-Maria](mailto:daniel.gadala-maria@yale.edu) (aut, cre)  
[Jason Vander Heiden](mailto:jason.vanderheiden@yale.edu) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)

Citing
---------------


Gadala-Maria D, Yaari G, Uduman M and Kleinstein Ss (2015). “Automated analysis of high-throughput B cell
sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles.” _Proceedings of the
National Academy of Sciency of the United States of America_, pp. E862-70.

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Automated analysis of high-throughput B cell sequencing data
         reveals a high frequency of novel immunoglobulin V gene segment
         alleles.},
    author = {Daniel Gadala-Maria and Gur Yaari and Mohamed Uduman and Steven H. steven.kleinstein@yale.edu Kleinstein},
    year = {2015},
    journal = {Proceedings of the National Academy of Sciency of the United
         States of America},
    number = {122},
    pages = {E862-70},
  }

