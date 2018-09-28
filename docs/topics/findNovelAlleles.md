**findNovelAlleles** - *Find novel alleles from repertoire sequencing data*

Description
--------------------

`findNovelAlleles` analyzes mutation patterns in sequences thought to
align to each germline allele in order to determine which positions
might be polymorphic.


Usage
--------------------
```
findNovelAlleles(clip_db, germline_db, v_call = "V_CALL",
germline_min = 200, min_seqs = 50, auto_mutrange = TRUE,
mut_range = 1:10, pos_range = 1:312, y_intercept = 0.125,
alpha = 0.05, j_max = 0.15, min_frac = 0.75, nproc = 1)
```

Arguments
-------------------

clip_db
:   a `data.frame` in Change-O format. See details.

germline_db
:   a vector of named nucleotide germline sequences
matching the V calls in `clip_db`

v_call
:   name of the column in clip_db with V allele calls. 
Default is V_CALL.

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

a `data.frame` with a row for each known allele analyzed.
Besides metadata on the the parameters used in the search, each row will have
either a note as to where the polymorphism-finding algorithm exited or a
nucleotide sequence for the predicted novel allele, along with columns providing
additional evidence.

Messages in the field `NOTE`:


+  *Novel allele found!* 
+  *Plurality sequence too rare*: No sequence is frequent 
enough to pass the J test (`j_max`).
+  *a J-junction combination is too prevalent*: Not enough
J diversity (`j_max`).
+  *No positions pass y-intercept test*. No positions above
`y_intercept`.
+  *Insufficient sequences in desired mutational range*. 
`mut_range` and `pos_range`.
+  *not enough sequences*:  not enough sequences in the 
desired mutational range and nucleotide range.
`min_seqs`
+  *no unmutated versions of novel allele found*


Other fields:

+  *GERMLINE_CALL*: The input V call
+  *POLYMORPHISM_CALL*: The new allele call
+  *NT_SUBSTITUTIONS*: Mutations identified in the new allele, relative
to the reference germline (`GERMLINE_CALL`)
+  *NOVEL_IMGT*: New allele
+  *NOVEL_IMGT_COUNT*:  the number of times the sequence 
`NOVEL_IMGT` is found in the input data 
`clip_db`. Considers the subsequence of 
`NOVEL_IMGT` in the range `pos_range`.  
+  *NOVEL_IMGT_UNIQUE_J*: Number of disctinct J calls associated
to `NOVEL_IMGT` in `clip_db`. Considers
the subsequence of `NOVEL_IMGT` in the range 
`pos_range`.       
+  *NOVEL_IMGT_UNIQUE_CDR3*: Number of disctinct CDR3 associated
to `NOVEL_IMGT` in `clip_db`. Considers
the subsequence of `NOVEL_IMGT` in the range 
`pos_range`.                                              
+  *PERFECT_MATCH_COUNT*: Final number of sequences retained to call 
the new allele. These are unique sequences that have 
Vs that perfect match the predicted germline in the 
range `pos_range`.
+  *PERFECT_MATCH_FREQ*: `PERFECT_MATCH_COUNT`/`GERMLINE_CALL_COUNT`
+  *GERMLINE_CALL_COUNT*: the number of sequences with the particular
`GERMLINE_CALL` in `clip_db` initially
considered for the analysis.
+  *GERMLINE_CALL_PERC*: the percent of sequences with the particular
`GERMLINE_CALL` in `clip_db` initially      
considered for the analysis.              
+  *GERMLINE_IMGT*: Germline sequence for `GERMLINE_CALL`
+  *GERMLINE_IMGT_COUNT*:  the number of times the sequence sequence
`GERMLINE_IMGT` is found in the input data 
`clip_db`.     
+  *MUT_MIN*: Minimum mutation considered by the algorithm
+  *MUT_MAX*: Maximum mutation considered by the algorithm
+  *MUT_PASS_COUNT*: Number of sequences in the mutation range
+  *POS_MIN*: First position of the sequence considered by the
algorithm (IMGT numbering)
+  *POS_MAX*: Last position of the sequence considered by the 
algorithm (IMGT numbering)
+  *Y_INTERCEPT*: The y-intercept above which positions were 
considered potentially polymorphic
+  *Y_INTERCEPT_PASS*: Number of positions that pass Y_INTERCEPT  
+  *SNP_PASS*: Number of sequences that pass Y_INTERCEPT and are
within the desired nucleotide range (`min_seqs`)   
+  *UNMUTATED_COUNT*: Number of unmutated sequences
+  *UNMUTATED_FREQ*: Number of unmutated sequences over 
`GERMLINE_IMGT_COUNT`
+  *UNMUTATED_SNP_J_GENE_LENGTH_COUNT*: Number of distinct combinations
of SNP, J gene and junction length.     
+  *SNP_MIN_SEQS_J_MAX_PASS*: Number of SNP strings that pass both the 
`min_seqs` and `j_max` thresholds                                                  
+  *ALPHA*: Significance cutoff to be used when constructing the 
confidence interval for the y-intercept
+  *MIN_SEQS*: Input `min_seqs`. The minimum number of total sequences (within the 
desired mutational range and nucleotide range) required 
for the samples to be considered
+  *J_MAX*: Input `j_max`. The maximum fraction of sequences perfectly aligning to 
a potential novel allele that are allowed to utilize to a 
particular combination of junction length and J gene
+  *MIN_FRAC*: Input `min_frac`. The minimum fraction of sequences that must have 
usable nucleotides in a given position for that position to 
be considered



Details
-------------------

A `data.frame` in Change-O format contains the following
columns:

+  `"SEQUENCE_IMGT"` containing the IMGT-gapped nucleotide sequence
+  `"V_CALL"` containing the IMGT/V-QUEST V allele call(s)
+  `"J_CALL"` containing the IMGT/V-QUEST J allele call(s)
+  `"JUNCTION_LENGTH"` containing the junction length

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
### Not run:
# Find novel alleles and return relevant data
# novel_df <- findNovelAlleles(SampleDb, GermlineIGHV)
```



See also
-------------------

[plotNovel](plotNovel.md) to visualize the data supporting any
novel alleles hypothesized to be present in the data and
[inferGenotype](inferGenotype.md) to determine if the novel alleles are frequent
enought to be included in the subject's genotype



