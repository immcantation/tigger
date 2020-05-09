**inferGenotype** - *Infer a subject-specific genotype using a frequency method*

Description
--------------------

`inferGenotype` infers an subject's genotype using a frequency method.
The genotype is inferred by finding the minimum number set of alleles that 
can explain the majority of each gene's calls. The most common allele of 
each gene is included in the genotype first, and the next most common allele 
is added until the desired fraction of alleles can be explained. In this 
way, mistaken allele calls (resulting from sequences which
by chance have been mutated to look like another allele) can be removed.


Usage
--------------------
```
inferGenotype(
data,
germline_db = NA,
novel = NA,
v_call = "v_call",
seq = "sequence_alignment",
fraction_to_explain = 0.875,
gene_cutoff = 1e-04,
find_unmutated = TRUE
)
```

Arguments
-------------------

data
:   a `data.frame` containing V allele
calls from a single subject. If
`find_unmutated` is `TRUE`, then
the sample IMGT-gapped V(D)J sequence should

germline_db
:   named vector of sequences containing the
germline sequences named in
`allele_calls`. Only required if
`find_unmutated` is `TRUE`.

novel
:   an optional `data.frame` of the type
novel returned by
[findNovelAlleles](findNovelAlleles.md) containing
germline sequences that will be utilized if
`find_unmutated` is `TRUE`. See
Details.

v_call
:   column in `data` with V allele calls.
Default is `"v_call"`.

seq
:   name of the column in `data` with the 
aligned, IMGT-numbered, V(D)J nucleotide sequence.
Default is `sequence_alignment`.

fraction_to_explain
:   the portion of each gene that must be
explained by the alleles that will be included
in the genotype.

gene_cutoff
:   either a number of sequences or a fraction of
the length of `allele_calls` denoting the
minimum number of times a gene must be
observed in `allele_calls` to be included
in the genotype.

find_unmutated
:   if `TRUE`, use `germline_db` to
find which samples are unmutated. Not needed
if `allele_calls` only represent
unmutated samples.




Value
-------------------

A `data.frame` of alleles denoting the genotype of the subject containing 
the following columns:


+  `gene`: The gene name without allele.
+  `alleles`: Comma separated list of alleles for the given `gene`.
+  `counts`: Comma separated list of observed sequences for each 
corresponding allele in the `alleles` list.
+  `total`: The total count of observed sequences for the given `gene`.
+  `note`: Any comments on the inferrence.



Details
-------------------

Allele calls representing cases where multiple alleles have been
assigned to a single sample sequence are rare among unmutated
sequences but may result if nucleotides for certain positions are
not available. Calls containing multiple alleles are treated as
belonging to all groups. If `novel` is provided, all
sequences that are assigned to the same starting allele as any
novel germline allele will have the novel germline allele appended
to their assignent prior to searching for unmutated sequences.


Note
-------------------

This method works best with data derived from blood, where a large
portion of sequences are expected to be unmutated. Ideally, there
should be hundreds of allele calls per gene in the input.



Examples
-------------------

```R
# Infer IGHV genotype, using only unmutated sequences, including novel alleles
inferGenotype(AIRRDb, germline_db=SampleGermlineIGHV, novel=SampleNovel,
find_unmutated=TRUE)
```

**Error in gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), **: object 'AIRRDb' not found

See also
-------------------

[plotGenotype](plotGenotype.md) for a colorful visualization and
[genotypeFasta](genotypeFasta.md) to convert the genotype to nucleotide sequences.
See [inferGenotypeBayesian](inferGenotypeBayesian.md) to infer a subject-specific genotype 
using a Bayesian approach.






