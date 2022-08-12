**findNovelAlleles** - *Find novel alleles from repertoire sequencing data*

Description
--------------------

`findNovelAlleles` analyzes mutation patterns in sequences thought to
align to each germline allele in order to determine which positions
might be polymorphic.


Usage
--------------------
```
findNovelAlleles(
data,
germline_db,
v_call = "v_call",
j_call = "j_call",
seq = "sequence_alignment",
junction = "junction",
junction_length = "junction_length",
germline_min = 200,
min_seqs = 50,
auto_mutrange = TRUE,
mut_range = 1:10,
pos_range = 1:312,
pos_range_max = NULL,
y_intercept = 0.125,
alpha = 0.05,
j_max = 0.15,
min_frac = 0.75,
nproc = 1
)
```

Arguments
-------------------

data
:   `data.frame` containing repertoire data. See details.

germline_db
:   vector of named nucleotide germline sequences
matching the V calls in `data`. These should be 
the gapped reference germlines used to make the V calls.

v_call
:   name of the column in `data` with V allele calls. 
Default is `v_call`.

j_call
:   name of the column in `data` with J allele calls. 
Default is `j_call`.

seq
:   name of the column in `data` with the 
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is `sequence_alignment`.

junction
:   Junction region nucleotide sequence, which includes
the CDR3 and the two flanking conserved codons. Default
is `junction`.

junction_length
:   Number of junction nucleotides in the junction sequence.
Default is `junction_length`.

germline_min
:   the minimum number of sequences that must have a
particular germline allele call for the allele to
be analyzed.

min_seqs
:   minimum number of total sequences (within the
desired mutational range and nucleotide range)
required for the samples to be considered.

auto_mutrange
:   if `TRUE`, the algorithm will attempt to
determine the appropriate mutation range
automatically using the mutation count of the most
common sequence assigned to each allele analyzed.

mut_range
:   range of mutations that samples may carry and
be considered by the algorithm.

pos_range
:   range of IMGT-numbered positions that should be
considered by the algorithm.

pos_range_max
:   Name of the column in `data` with the ending
positions of the V alignment in the germline 
(usually `v_germline_end`). The end of the alignment will
be used to limit the range of positions to be 
considered to count mutations. With `NULL` all 
positions in the IMGT V region will be considered. In
this case, in sequences where the V was trimmed 
on the 3', mutated nucleotides could include nucleotides
from the CDR3.

y_intercept
:   y-intercept threshold above which positions should be
considered potentially polymorphic.

alpha
:   alpha value used for determining whether the 
fit y-intercept is greater than the `y_intercept`
threshold.

j_max
:   maximum fraction of sequences perfectly aligning
to a potential novel allele that are allowed to
utilize to a particular combination of junction
length and J gene. The closer to 1, the less strict 
the filter for junction length and J gene diversity
will be.

min_frac
:   minimum fraction of sequences that must have
usable nucleotides in a given position for that
position to considered.

nproc
:   number of processors to use.




Value
-------------------

A `data.frame` with a row for each known allele analyzed.
Besides metadata on the the parameters used in the search, each row will have
either a note as to where the polymorphism-finding algorithm exited or a
nucleotide sequence for the predicted novel allele, along with columns providing
additional evidence.

The output contains the following columns:

+  `germline_call`: The input (uncorrected) V call.
+  `note`: Comments regarding the inferrence.
+  `polymorphism_call`: The novel allele call.
+  `nt_substitutions`: Mutations identified in the novel allele, relative
to the reference germline (`germline_call`)
+  `novel_imgt`: The novel allele sequence.
+  `novel_imgt_count`:  The number of times the sequence `novel_imgt` 
is found in the input data. Considers the subsequence of `novel_imgt` 
in the `pos_range`.
+  `novel_imgt_unique_j`: Number of distinct J calls associated to `novel_imgt` 
in the input data. Considers the subsequence of `novel_imgt` in the `pos_range`.       
+  `novel_imgt_unique_cdr3`: Number of distinct CDR3 sequences associated
with `novel_imgt` in the input data. Considers the subsequence of `novel_imgt` 
in the `pos_range`.                                              
+  `perfect_match_count`: Final number of sequences retained to call the new 
allele. These are unique sequences that have V segments that perfectly match 
the predicted germline in the `pos_range`.
+  `perfect_match_freq`: `perfect_match_count / germline_call_count`
+  `germline_call_count`: The number of sequences with the `germline_call` 
in the input data that were initially considered for the analysis.
+  `germline_call_freq`: The fraction of sequences with the `germline_call` 
in the input data initially considered for the analysis.              
+  `germline_imgt`: Germline sequence for `germline_call`.
+  `germline_imgt_count`: The number of times the `germline_imgt` 
sequence is found in the input data.
+  `mut_min`: Minimum mutation considered by the algorithm.
+  `mut_max`: Maximum mutation considered by the algorithm.
+  `mut_pass_count`: Number of sequences in the mutation range.
+  `pos_min`: First position of the sequence considered by the algorithm (IMGT numbering).
+  `pos_max`: Last position of the sequence considered by the algorithm (IMGT numbering).
+  `y_intercept`: The y-intercept above which positions were considered 
potentially polymorphic.
+  `y_intercept_pass`: Number of positions that pass the `y_intercept` threshold.
+  `snp_pass`: Number of sequences that pass the `y_intercept` threshold and are
within the desired nucleotide range (`min_seqs`).
+  `unmutated_count`: Number of unmutated sequences.
+  `unmutated_freq`: Number of unmutated sequences over `germline_imgt_count`.
+  `unmutated_snp_j_gene_length_count`: Number of distinct combinations
of SNP, J gene, and junction length.     
+  `snp_min_seqs_j_max_pass`: Number of SNPs that pass both the `min_seqs` 
and `j_max` thresholds.
+  `alpha`: Significance threshold to be used when constructing the 
confidence interval for the y-intercept.
+  `min_seqs`: Input `min_seqs`. The minimum number of total sequences 
(within the desired mutational range and nucleotide range) required 
for the samples to be considered.
+  `j_max`: Input `j_max`. The maximum fraction of sequences perfectly 
aligning to a potential novel allele that are allowed to utilize to a particular 
combination of junction length and J gene.
+  `min_frac`: Input `min_frac`. The minimum fraction of sequences that must
have usable nucleotides in a given position for that position to be considered.


The following comments can appear in the `note` column:


+  *Novel allele found*: A novel allele was detected.
+  *Same as:*: The same novel allele sequence
has been identified multiple times. If this happens, the function
will also throw the message 'Duplicated polymorphism(s) found'.
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
# Note: In this example, with SampleGermlineIGHV,
# which contains reference germlines retrieved on August 2014,
# TIgGER finds the allele IGHV1-8*02_G234T. This allele
# was added to IMGT as IGHV1-8*03 on March 28, 2018.

# Find novel alleles and return relevant data
novel <- findNovelAlleles(AIRRDb, SampleGermlineIGHV)
selectNovel(novel)
```


```
  germline_call                note polymorphism_call nt_substitutions
1    IGHV1-8*02 Novel allele found!  IGHV1-8*02_G234T           234G>T
                                                                                                                                                                                                                                                                                                                        novel_imgt
1 CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATTACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG
  novel_imgt_count novel_imgt_unique_j novel_imgt_unique_cdr3
1              657                   6                    626
  perfect_match_count perfect_match_freq germline_call_count germline_call_freq
1                 661          0.7295806                 906              0.052
  mut_min mut_max mut_pass_count
1       1      10            760
                                                                                                                                                                                                                                                                                                                     germline_imgt
1 CAGGTGCAGCTGGTGCAGTCTGGGGCT...GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTC............ACCAGCTATGATATCAACTGGGTGCGACAGGCCACTGGACAAGGGCTTGAGTGGATGGGATGGATGAACCCTAAC......AGTGGTAACACAGGCTATGCACAGAAGTTCCAG...GGCAGAGTCACCATGACCAGGAACACCTCCATAAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGG
  germline_imgt_count pos_min pos_max y_intercept y_intercept_pass snp_pass
1                   0       1     312       0.125                1      754
  unmutated_count unmutated_freq unmutated_snp_j_gene_length_count
1             661      0.7295806                                83
  snp_min_seqs_j_max_pass alpha min_seqs j_max min_frac
1                       1  0.05       50  0.15     0.75

```



See also
-------------------

[selectNovel](selectNovel.md) to filter the results to show only novel alleles.
[plotNovel](plotNovel.md) to visualize the data supporting any
novel alleles hypothesized to be present in the data and
[inferGenotype](inferGenotype.md) and [inferGenotypeBayesian](inferGenotypeBayesian.md) to determine if the novel alleles are frequent
enought to be included in the subject's genotype.






