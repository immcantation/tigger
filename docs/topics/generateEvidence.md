**generateEvidence** - *Generate evidence*

Description
--------------------

`generateEvidence` builds a table of evidence metrics for the final novel V 
allele detection and genotyping inferrences.


Usage
--------------------
```
generateEvidence(data, novel, genotype, genotype_db, germline_db,
j_call = "j_call", junction = "junction", fields = NULL)
```

Arguments
-------------------

data
:   a `data.frame` containing sequence data that has been
passed through [reassignAlleles](reassignAlleles.md) to correct the allele 
assignments.

novel
:   the `data.frame` returned by [findNovelAlleles](findNovelAlleles.md).

genotype
:   the `data.frame` of alleles generated with [inferGenotype](inferGenotype.md) 
denoting the genotype of the subject.

genotype_db
:   a vector of named nucleotide germline sequences in the genotype.
Returned by [genotypeFasta](genotypeFasta.md).

germline_db
:   the original uncorrected germline database used to by
[findNovelAlleles](findNovelAlleles.md) to identify novel alleles.

j_call
:   name of the column in `data` with J allele calls. 
Default is `j_call`.

junction
:   Junction region nucleotide sequence, which includes
the CDR3 and the two flanking conserved codons. Default
is JUNCTION

fields
:   character vector of column names used to split the data to 
identify novel alleles, if any. If `NULL` then the data is 
not divided by grouping variables.




Value
-------------------

Returns the `genotype` input `data.frame` with the following additional columns 
providing supporting evidence for each inferred allele:


+  `FIELD_ID`: Data subset identifier, defined with the input paramter `fields`.
+  A variable number of columns, specified with the input parameter `fields`.
+  `POLYMORPHISM_CALL`: The novel allele call.
+  `NOVEL_IMGT`: The novel allele sequence.
+  `CLOSEST_REFERENCE`: The closest reference gene and allele in 
the `germline_db` database.
+  `CLOSEST_REFERENCE_IMGT`: Sequence of the closest reference gene and 
allele in the `germline_db` database.
+  `GERMLINE_CALL`: The input (uncorrected) V call.
+  `GERMLINE_IMGT`: Germline sequence for `GERMLINE_CALL`.
+  `NT_DIFF`: Number of nucleotides that differ between the new allele and
the closest reference (`CLOSEST_REFERENCE`) in the `germline_db` database.
+  `NT_SUBSTITUTIONS`: A comma separated list of specific nucleotide 
differences (e.g. `112G>A`) in the novel allele.
+  `AA_DIFF`: Number of amino acids that differ between the new allele and the closest 
reference (`CLOSEST_REFERENCE`) in the `germline_db` database.
+  `AA_SUBSTITUTIONS`: A comma separated list with specific amino acid 
differences (e.g. `96A>N`) in the novel allele.
+  `SEQUENCES`: Number of sequences unambiguosly assigned to this allele.
+  `UNMUTATED_SEQUENCES`: Number of records with the unmutated novel allele sequence.
+  `UNMUTATED_FREQUENCY`: Proportion of records with the unmutated novel allele 
sequence (`UNMUTATED_SEQUENCES / SEQUENCE`).
+  `ALLELIC_PERCENTAGE`: Percentage at which the (unmutated) allele is observed 
in the sequence dataset compared  to other (unmutated) alleles.
+  `UNIQUE_JS`: Number of unique J sequences found associated with the 
novel allele. The sequences are those who have been unambiguously assigned 
to the novel allelle (`POLYMORPHISM_CALL`).
+  `UNIQUE_CDR3S`: Number of unique CDR3s associated with the inferred allele.
The sequences are those who have been unambiguously assigned to the 
novel allelle (POLYMORPHISM_CALL).
+  `MUT_MIN`: Minimum mutation considered by the algorithm.
+  `MUT_MAX`: Maximum mutation considered by the algorithm.
+  `POS_MIN`: First position of the sequence considered by the algorithm (IMGT numbering).
+  `POS_MAX`: Last position of the sequence considered by the algorithm (IMGT numbering).
+  `Y_INTERCEPT`: The y-intercept above which positions were considered 
potentially polymorphic.
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
+  `NOTE`: Comments regarding the novel allele inferrence.




Examples
-------------------

```R
# Generate input data
novel <- findNovelAlleles(airrDb, SampleGermlineIGHV,
v_call="v_call", j_call="j_call", junction="junction", 
junction_length="junction_length", seq="sequence_alignment")
genotype <- inferGenotype(airrDb, find_unmutated=TRUE, 
germline_db=SampleGermlineIGHV,
novel=novel,
v_call="v_call", seq="sequence_alignment")
genotype_db <- genotypeFasta(genotype, SampleGermlineIGHV, novel)
data_db <- reassignAlleles(airrDb, genotype_db, 
v_call="v_call", seq="sequence_alignment")

# Assemble evidence table
evidence <- generateEvidence(data_db, novel, genotype, 
genotype_db, SampleGermlineIGHV,
j_call = "j_call", 
junction = "junction")
```



See also
-------------------

See [findNovelAlleles](findNovelAlleles.md), [inferGenotype](inferGenotype.md) and [genotypeFasta](genotypeFasta.md) 
for generating the required input.






