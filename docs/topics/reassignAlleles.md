**reassignAlleles** - *Correct allele calls based on a personalized genotype*

Description
--------------------

`reassignAlleles` uses a subject-specific genotype to correct
correct preliminary allele assignments of a set of sequences derived
from a single subject.


Usage
--------------------
```
reassignAlleles(
data,
genotype_db,
v_call = "v_call",
seq = "sequence_alignment",
method = "hamming",
path = NA,
keep_gene = c("gene", "family", "repertoire"),
trim_seq = FALSE,
overwrite = FALSE,
ignored_regex = "[\\.N-]",
treat_multigene_as_uncalled = FALSE,
top_k = NULL,
top_by = c("alphabetical", "mutation_count"),
strip_d = TRUE,
reassign_uncalled = TRUE
)
```

Arguments
-------------------

data
:   `data.frame` containing V allele calls from a
single subject and the sample IMGT-gapped V(D)J sequences under
`seq`.

genotype_db
:   vector of named nucleotide germline sequences
matching the calls detailed in `allele_calls`
and personalized to the subject

v_call
:   name of the column in `data` with V allele
calls. Default is `v_call`.

seq
:   name of the column in `data` with the
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is SEQUENCE_IMGT

method
:   method to use when realigning sequences to
the genotype_db sequences. Currently, only `"hamming"`
(for Hamming distance) is implemented.

path
:   directory containing the tool used in the
realignment method, if needed. Hamming distance does
not require a path to a tool.

keep_gene
:   string indicating if the gene (`"gene"`),
family (`"family"`) or complete repertoire
(`"repertoire"`) assignments should be performed.
Use of `"gene"` increases speed by minimizing required number of
alignments, as gene level assignments will be maintained when possible.

trim_seq
:   if `TRUE`, trim sample and germline sequences
to the segment boundaries before calculating Hamming
distance. Boundaries are determined from the segment
prefix of `v_call`, such as `v_*`,
`d_*`, or `j_*` columns.

overwrite
:   if `TRUE`, replace `v_call` with reassigned
calls instead of writing a `*_call_genotyped`
column.

ignored_regex
:   regular expression indicating characters to ignore
when comparing sequences. May also be `TRUE` to
ignore nothing (every position counts), as used for
D and J segments.

treat_multigene_as_uncalled
:   if `TRUE`, sequences whose call
spans more than one gene are treated as uncalled and
realigned against the whole genotype rather than kept
at their first gene. Only applies when `keep_gene`
is `"gene"` or `"repertoire"`.

top_k
:   maximum number of equally-best alleles to report per
sequence. `NULL` (default) keeps all ties.

top_by
:   how to break ties when more than `top_k` alleles
are equally close. `"alphabetical"` keeps the
first `top_k` by name; `"mutation_count"`
keeps all ties.

strip_d
:   if `TRUE` (default) remove the "D" from the end of
gene annotations (denoting a duplicate gene in the locus)
when grouping the `genotype_db` alleles by gene. If
`FALSE`, the "D" is kept, so the genotype grouping
and the sequence calls are matched consistently (use
together with `genotypeFasta(strip_d=FALSE)`).

reassign_uncalled
:   if `TRUE` (default), sequences whose call is
empty or `NA` are also realigned against the whole
genotype. If `FALSE`, they are left unassigned.




Value
-------------------

A modified input `data.frame` containing the best allele call from
among the sequences listed in `genotype_db` in the
`*_call_genotyped` column, or in `v_call` when
`overwrite=TRUE`.


Details
-------------------

In order to save time, initial gene assignments are preserved and
the allele calls are chosen from among those provided in `genotype_db`,
based on a simple alignment to the sample sequence.



Examples
-------------------

```R
# Extract the database sequences that correspond to the genotype
genotype_db <- genotypeFasta(SampleGenotype, SampleGermlineIGHV, novel=SampleNovel)

# Use the personalized genotype to determine corrected allele assignments
output_db <- reassignAlleles(AIRRDb, genotype_db, v_call="v_call",
seq="sequence_alignment")

```








