**generateEvidence** - *Generate evidence*

Description
--------------------

`generateEvidence` builds a table of evidence metrics for the final novel V 
allele detection and genotyping inferrences.


Usage
--------------------
```
generateEvidence(
data,
novel,
genotype,
genotype_db,
germline_db,
j_call = "j_call",
junction = "junction",
fields = NULL
)
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
is `junction`.

fields
:   character vector of column names used to split the data to 
identify novel alleles, if any. If `NULL` then the data is 
not divided by grouping variables.




Value
-------------------

Returns the `genotype` input `data.frame` with the following additional columns 
providing supporting evidence for each inferred allele:


+  `field_id`: Data subset identifier, defined with the input paramter `fields`.
+  A variable number of columns, specified with the input parameter `fields`.
+  `polymorphism_call`: The novel allele call.
+  `novel_imgt`: The novel allele sequence.
+  `closest_reference`: The closest reference gene and allele in 
the `germline_db` database.
+  `closest_reference_imgt`: Sequence of the closest reference gene and 
allele in the `germline_db` database.
+  `germline_call`: The input (uncorrected) V call.
+  `germline_imgt`: Germline sequence for `germline_call`.
+  `nt_diff`: Number of nucleotides that differ between the new allele and
the closest reference (`closest_reference`) in the `germline_db` database.
+  `nt_substitutions`: A comma separated list of specific nucleotide 
differences (e.g. `112G>A`) in the novel allele.
+  `aa_diff`: Number of amino acids that differ between the new allele and the closest 
reference (`closest_reference`) in the `germline_db` database.
+  `aa_substitutions`: A comma separated list with specific amino acid 
differences (e.g. `96A>N`) in the novel allele.
+  `sequences`: Number of sequences unambiguosly assigned to this allele.
+  `unmutated_sequences`: Number of records with the unmutated novel allele sequence.
+  `unmutated_frequency`: Proportion of records with the unmutated novel allele 
sequence (`unmutated_sequences / sequences`).
+  `allelic_percentage`: Percentage at which the (unmutated) allele is observed 
in the sequence dataset compared  to other (unmutated) alleles.
+  `unique_js`: Number of unique J sequences found associated with the 
novel allele. The sequences are those who have been unambiguously assigned 
to the novel allelle (`polymorphism_call`).
+  `unique_cdr3s`: Number of unique CDR3s associated with the inferred allele.
The sequences are those who have been unambiguously assigned to the 
novel allelle (polymorphism_call).
+  `mut_min`: Minimum mutation considered by the algorithm.
+  `mut_max`: Maximum mutation considered by the algorithm.
+  `pos_min`: First position of the sequence considered by the algorithm (IMGT numbering).
+  `pos_max`: Last position of the sequence considered by the algorithm (IMGT numbering).
+  `y_intercept`: The y-intercept above which positions were considered 
potentially polymorphic.
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
+  `note`: Comments regarding the novel allele inferrence.




Examples
-------------------

```R
# Generate input data
novel <- findNovelAlleles(AIRRDb, SampleGermlineIGHV,
v_call="v_call", j_call="j_call", junction="junction", 
junction_length="junction_length", seq="sequence_alignment")

```

**Error in eval(lhs, parent, parent)**: object 'AIRRDb' not found
```R
genotype <- inferGenotype(AIRRDb, find_unmutated=TRUE, 
germline_db=SampleGermlineIGHV,
novel=novel,
v_call="v_call", seq="sequence_alignment")

```

**Error in gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), **: object 'AIRRDb' not found
```R
genotype_db <- genotypeFasta(genotype, SampleGermlineIGHV, novel)

```

**Error in gsub("[Dd]\\*", "*", genotype$gene)**: object 'genotype' not found
```R
data_db <- reassignAlleles(AIRRDb, genotype_db, 
v_call="v_call", seq="sequence_alignment")

```

**Error in reassignAlleles(AIRRDb, genotype_db, v_call = "v_call", seq = "sequence_alignment")**: object 'AIRRDb' not found
```R

# Assemble evidence table
evidence <- generateEvidence(data_db, novel, genotype, 
genotype_db, SampleGermlineIGHV,
j_call = "j_call", 
junction = "junction")
```

**Error in eval(lhs, parent, parent)**: object 'genotype' not found

See also
-------------------

See [findNovelAlleles](findNovelAlleles.md), [inferGenotype](inferGenotype.md) and [genotypeFasta](genotypeFasta.md) 
for generating the required input.






