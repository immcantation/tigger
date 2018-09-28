**tigger** - *tigger*

Description
--------------------

Here we provide a **T**ool for **I**mmuno**g**lobulin
**G**enotype **E**lucidation via
**R**ep-Seq (TIgGER). TIgGER inferrs the set of Ig alleles carried by an
individual (including any novel alleles) and then uses this set of alleles to
correct the initial assignments given to sample sequences by existing tools.




Details
-------------------

Immunoglobulin Repertoire-Sequencing (Rep-Seq) data is currently the
subject of much study. A key step in analyzing these data involves assigning
the closest known V(D)J germline alleles to the (often somatically mutated)
sample sequences using a tool such as IMGT/HighV-QUEST. However, if the
sample utilizes alleles not in the germline database used for alignment, this
step will fail. Additionally, this alignment has an associated error rate of
~5 percent, notably among sequences carrying a large number of somatic
mutations. The purpose of TIgGER is to address these issues.


Core tigger functions
-------------------



+ [findNovelAlleles](findNovelAlleles.md):   Detect novel alleles
+ [plotNovel](plotNovel.md):          Plot evidence of novel alleles
+ [inferGenotype](inferGenotype.md):      Infer an Ig genotype
+ [plotGenotype](plotGenotype.md):       A colorful genotype visualization
+ [genotypeFasta](genotypeFasta.md):      Convert a genotype to sequences
+ [reassignAlleles](reassignAlleles.md):    Correct allele calls



Mutation-related functions
-------------------



+ [getMutatedPositions](getMutatedPositions.md):      Find mutation locations
+ [getMutCount](getMutCount.md):              Find distance from germline
+ [findUnmutatedCalls](findUnmutatedCalls.md):       Subset unmutated sequences
+ [getPopularMutationCount](getPopularMutationCount.md):  Find most common sequence's
mutation count
+ [insertPolymorphisms](insertPolymorphisms.md):      Insert SNPs into a sequence



Input and formatting
-------------------



+ [readIgFasta](readIgFasta.md):        Read a fasta file of Ig sequences
+ [updateAlleleNames](updateAlleleNames.md):  Correct outdated allele names
+ [sortAlleles](sortAlleles.md):        Sort allele names intelligently
+ [cleanSeqs](cleanSeqs.md):          Standardize sequence format



References
-------------------

Gadala-Maria *et al*. (2015) Automated analysis of
high-throughput B cell sequencing data reveals a high frequency of novel
immunoglobulin V gene segment alleles. *PNAS*. 112(8):E862-70.





