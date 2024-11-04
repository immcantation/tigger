[![](http://cranlogs.r-pkg.org/badges/grand-total/tigger)](https://www.r-pkg.org/pkg/tigger)
[![](https://cranlogs.r-pkg.org/badges/tigger)](https://www.r-pkg.org/pkg/tigger)
[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

!!! important "2025 Immcantation Users Group Meeting"
    *Are you an Immcantation user and/or interested in adaptive immune receptor repertoire analysis?*
    
    Register now for the upcoming Immcantation Users Group Meeting!
    It will be held virtually on **January 30th, 2025, from 10 to 1:30pm (ET)**.
    All talks will be from user-submitted abstracts.

    Full information here: [https://immcantation.github.io/users-meeting/](https://immcantation.github.io/users-meeting/)

TIgGER
-------------------------------------------------------------------------------

High-throughput sequencing of B cell immunoglobulin receptors is providing unprecedented insight into adaptive immunity. A key step in analyzing these data involves assignment of the germline V, D and J gene segment alleles that comprise each immunoglobulin sequence by matching them against a database of known V(D)J alleles. However, this process will fail for sequences that utilize previously undetected alleles, whose frequency in the population is unclear.

TIgGER is a computational method that significantly improves V(D)J allele assignments by first determining the complete set of gene segments carried by an individual (including novel alleles) from V(D)J-rearrange sequences. TIgGER can then infer a subject's genotype from these sequences, and use this genotype to correct the initial V(D)J allele assignments.

The application of TIgGER continues to identify a surprisingly high frequency of novel alleles in humans, highlighting the critical need for this approach. (TIgGER, however, can and has been used with data from other species.)

Core Abilities
-------------------------------------------------------------------------------

* Detecting novel alleles
* Inferring a subject's genotype
* Correcting preliminary allele calls

Required Input
-------------------------------------------------------------------------------

* A table of sequences from a single individual, with columns containing the following:
    * V(D)J-rearranged nucleotide sequence (in IMGT-gapped format)
    * Preliminary V allele calls
    * Preliminary J allele calls
    * Length of the junction region
* Germline Ig sequences in IMGT-gapped fasta format (e.g., as those downloaded from [IMGT/GENE-DB](https://www.imgt.org/genedb/)

The former can be created through the use of [IMGT/HighV-QUEST](https://www.imgt.org) and [Change-O](http://changeo.readthedocs.io).

Contact
-------------------------------------------------------------------------------

For help, questions, or suggestions, please contact the [Immcantation Group](mailto:immcantation@googlegroups.com) or use the [issue tracker](https://bitbucket.org/kleinstein/tigger/issues?status=new&status=open).
