**findNovelAlleles** - *Find novel alleles from repertoire sequencing data*

Description
--------------------

`findNovelAlleles` analyzes mutation patterns in sequences thought to
align to each germline allele in order to determine which positions
might be polymorphic.


Usage
--------------------
```
findNovelAlleles(data, germline_db, v_call = "V_CALL",
j_call = "J_CALL", seq = "SEQUENCE_IMGT", junction = "JUNCTION",
junction_length = "JUNCTION_LENGTH", germline_min = 200,
min_seqs = 50, auto_mutrange = TRUE, mut_range = 1:10,
pos_range = 1:312, y_intercept = 0.125, alpha = 0.05,
j_max = 0.15, min_frac = 0.75, nproc = 1)
```

Arguments
-------------------

data
:   a `data.frame` in Change-O format. See details.

germline_db
:   a vector of named nucleotide germline sequences
matching the V calls in `data`.

v_call
:   name of the column in `data` with V allele calls. 
Default is V_CALL.

j_call
:   name of the column in `data` with J allele calls. 
Default is J_CALL.

seq
:   name of the column in `data` with the 
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is SEQUENCE_IMGT.

junction
:   Junction region nucleotide sequence, which includes
the CDR3 and the two flanking conserved codons. Default
is JUNCTION.

junction_length
:   Number of junction nucleotides in the junction sequence.
Default is JUNCTION_LENGTH.

germline_min
:   the minimum number of sequences that must have a
particular germline allele call for the allele to
be analyzed

min_seqs
:   the minimum number of total sequences (within the
desired mutational range and nucleotide range)
required for the samples to be considered

auto_mutrange
:   if `TRUE`, the algorithm will attempt to
determine the appropriate mutation range
automatically using the mutation count of the most
common sequence assigned to each allele analyzed

mut_range
:   the range of mutations that samples may carry and
be considered by the algorithm

pos_range
:   the range of IMGT-numbered positions that should be
considered by the algorithm

y_intercept
:   the y-intercept threshold above which positions should be
considered potentially polymorphic

alpha
:   the alpha value used for determining whether the 
fit y-intercept is greater than the `y_intercept`
threshold

j_max
:   the maximum fraction of sequences perfectly aligning
to a potential novel allele that are allowed to
utilize to a particular combination of junction
length and J gene

min_frac
:   the minimum fraction of sequences that must have
usable nucleotides in a given position for that
position to considered

nproc
:   the number of processors to use




Value
-------------------

A `data.frame` with a row for each known allele analyzed.
Besides metadata on the the parameters used in the search, each row will have
either a note as to where the polymorphism-finding algorithm exited or a
nucleotide sequence for the predicted novel allele, along with columns providing
additional evidence.

The output contains the following columns:

+  `GERMLINE_CALL`: The input (uncorrected) V call.
+  `NOTE`: Comments regarding the inferrence.
+  `POLYMORPHISM_CALL`: The novel allele call.
+  `NT_SUBSTITUTIONS`: Mutations identified in the novel allele, relative
to the reference germline (`GERMLINE_CALL`)
+  `NOVEL_IMGT`: The novel allele sequence.
+  `NOVEL_IMGT_COUNT`:  The number of times the sequence `NOVEL_IMGT` 
is found in the input data. Considers the subsequence of `NOVEL_IMGT` 
in the `pos_range`.
+  `NOVEL_IMGT_UNIQUE_J`: Number of distinct J calls associated to `NOVEL_IMGT` 
in the input data. Considers the subsequence of `NOVEL_IMGT` in the `pos_range`.       
+  `NOVEL_IMGT_UNIQUE_CDR3`: Number of distinct CDR3 sequences associated
with `NOVEL_IMGT` in the input data. Considers the subsequence of `NOVEL_IMGT` 
in the `pos_range`.                                              
+  `PERFECT_MATCH_COUNT`: Final number of sequences retained to call the new 
allele. These are unique sequences that have V segments that perfectly match 
the predicted germline in the `pos_range`.
+  `PERFECT_MATCH_FREQ`: `PERFECT_MATCH_COUNT / GERMLINE_CALL_COUNT`
+  `GERMLINE_CALL_COUNT`: The number of sequences with the `GERMLINE_CALL` 
in the input data that were initially considered for the analysis.
+  `GERMLINE_CALL_FREQ`: The fraction of sequences with the `GERMLINE_CALL` 
in the input data initially considered for the analysis.              
+  `GERMLINE_IMGT`: Germline sequence for `GERMLINE_CALL`.
+  `GERMLINE_IMGT_COUNT`: The number of times the `GERMLINE_IMGT` 
sequence is found in the input data.
+  `MUT_MIN`: Minimum mutation considered by the algorithm.
+  `MUT_MAX`: Maximum mutation considered by the algorithm.
+  `MUT_PASS_COUNT`: Number of sequences in the mutation range.
+  `POS_MIN`: First position of the sequence considered by the algorithm (IMGT numbering).
+  `POS_MAX`: Last position of the sequence considered by the algorithm (IMGT numbering).
+  `Y_INTERCEPT`: The y-intercept above which positions were considered 
potentially polymorphic.
+  `Y_INTERCEPT_PASS`: Number of positions that pass the `Y_INTERCEPT` threshold.
+  `SNP_PASS`: Number of sequences that pass the `Y_INTERCEPT` threshold and are
within the desired nucleotide range (`min_seqs`).
+  `UNMUTATED_COUNT`: Number of unmutated sequences.
+  `UNMUTATED_FREQ`: Number of unmutated sequences over `GERMLINE_IMGT_COUNT`.
+  `UNMUTATED_SNP_J_GENE_LENGTH_COUNT`: Number of distinct combinations
of SNP, J gene, and junction length.     
+  `SNP_MIN_SEQS_J_MAX_PASS`: Number of SNPs that pass both the `min_seqs` 
and `j_max` thresholds.
+  `ALPHA`: Significance threshold to be used when constructing the 
confidence interval for the y-intercept.
+  `MIN_SEQS`: Input `min_seqs`. The minimum number of total sequences 
(within the desired mutational range and nucleotide range) required 
for the samples to be considered.
+  `J_MAX`: Input `j_max`. The maximum fraction of sequences perfectly 
aligning to a potential novel allele that are allowed to utilize to a particular 
combination of junction length and J gene.
+  `MIN_FRAC`: Input `min_frac`. The minimum fraction of sequences that must
have usable nucleotides in a given position for that position to be considered.


The following comments can appear in the `NOTE` column:


+  *Novel allele found*: A novel allele was detected.
+  *Plurality sequence too rare*: No sequence is frequent enough to pass 
the J test (`j_max`).
+  *A J-junction combination is too prevalent*: Not enough J diversity (`j_max`).
+  *No positions pass y-intercept test*: No positions above `y_intercept`.
+  *Insufficient sequences in desired mutational range*: 
`mut_range` and `pos_range`.
+  *Not enough sequences*: Not enough sequences in the desired mutational 
range and nucleotide range (`min_seqs`).
+  *No unmutated versions of novel allele found*: All observed variants of the 
allele are mutated.



Details
-------------------

The TIgGER allele-finding algorithm, briefly, works as follows:
Mutations are determined through comparison to the provided germline.
Mutation frequency at each *position* is determined as a function of
*sequence-wide* mutation counts. Polymorphic positions exhibit a high
mutation frequency despite sequence-wide mutation count. False positive of
potential novel alleles resulting from clonally-related sequences are guarded
against by ensuring that sequences perfectly matching the potential novel
allele utilize a wide range of combinations of J gene and junction length.



Examples
-------------------

```R
# Find novel alleles and return relevant data
novel <- findNovelAlleles(SampleDb, SampleGermlineIGHV)
selectNovel(novel)

```


```
  GERMLINE_CALL                NOTE POLYMORPHISM_CALL NT_SUBSTITUTIONS
1    IGHV1-8*02 Novel allele found!  IGHV1-8*02_G234T           234G>T
                                                                                                                                                                                                                                                                                                                        NOVEL_IMGT
1 CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATTACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG
  NOVEL_IMGT_COUNT NOVEL_IMGT_UNIQUE_J NOVEL_IMGT_UNIQUE_CDR3 PERFECT_MATCH_COUNT
1              657                   6                    626                 661
  PERFECT_MATCH_FREQ GERMLINE_CALL_COUNT GERMLINE_CALL_FREQ MUT_MIN MUT_MAX MUT_PASS_COUNT
1          0.7295806                 906              0.052       1      10            760
                                                                                                                                                                                                                                                                                                                     GERMLINE_IMGT
1 CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG
  GERMLINE_IMGT_COUNT POS_MIN POS_MAX Y_INTERCEPT Y_INTERCEPT_PASS SNP_PASS UNMUTATED_COUNT
1                   0       1     312       0.125                1      754             661
  UNMUTATED_FREQ UNMUTATED_SNP_J_GENE_LENGTH_COUNT SNP_MIN_SEQS_J_MAX_PASS ALPHA MIN_SEQS
1      0.7295806                                83                       1  0.05       50
  J_MAX MIN_FRAC
1  0.15     0.75

```


```R
# Note: In this example, with SampleGermlineIGHV,
# which contains reference germlines retrieved on August 2018,
# TIgGER finds the allele IGHV1-8*02_G234T. This allele
# was added to IMGT as IGHV1-8*03 on March 28, 2018.
```



See also
-------------------

[selectNovel](selectNovel.md) to filter the results to show only novel alleles.
[plotNovel](plotNovel.md) to visualize the data supporting any
novel alleles hypothesized to be present in the data and
[inferGenotype](inferGenotype.md) to determine if the novel alleles are frequent
enought to be included in the subject's genotype.



