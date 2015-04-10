# sortAlleles -------------------------------------------------------------
#' Sort allele names
#'
#' \code{sortAlleles} returns a sorted vector of strings respresenting Ig allele
#' names. Names are first sorted by gene family, then by gene, then by allele.
#' Duplicated genes have their alleles are sorted as if they were part of their
#' non-duplicated counterparts (e.g. IGHV1-69D*01 comes after IGHV1-69*01 but
#' before IGHV1-69*02), and non-localized genes (e.g. IGHV1-NL1*01) come last
#' within their gene family.
#' 
#' @param    allele_calls  a vector of strings respresenting Ig allele names
#' @return   A sorted vector of strings respresenting Ig allele names
#' 
#' @examples
#' # Create a list of allele names
#' alleles = c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
#' "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*02,IGHV1-2*05", "IGHV1-2",
#' "IGHV1-2*02", "IGHV1-69*02")
#' 
#' # Sort the alleles
#' sortAlleles(alleles)
#' 
#' @export
sortAlleles <- function(allele_calls) {  
  
  # Standardize format of submitted alleles, first
  SUBMITTED_CALLS = getAllele(allele_calls, first = FALSE)
  
  allele_df = data.frame(SUBMITTED_CALLS, stringsAsFactors = FALSE) %>%
    # Sort to help with the Ds later
    arrange(SUBMITTED_CALLS) %>%
    # Determine the family
    mutate(FAMILY = getFamily(SUBMITTED_CALLS)) %>%
    # Determine the gene (exclude family); convert letters to numbers for sort
    mutate(GENE = getGene(SUBMITTED_CALLS)) %>%
    mutate(GENE1 = gsub("[^-]+-([^-\\*D]+).*","\\1",SUBMITTED_CALLS)) %>%
    mutate(GENE1 = as.numeric(gsub("NL|a|b|f", "99", GENE1))) %>%
    # If there is a second gene number, determine that, too
    mutate(GENE2 = gsub("[^-]+-[^-]+-?","",GENE)) %>%
    mutate(GENE2 = as.numeric(gsub("NL|a|b|f", "99", GENE2))) %>%
    mutate(ALLELE = as.numeric(gsub(".+\\*?","",getAllele(SUBMITTED_CALLS))))
  
  # Convert missing values to 0, sort data frame
  allele_df[is.na(allele_df)] = 0
  sorted_df = arrange(allele_df, FAMILY, GENE1, GENE2, ALLELE)
  
  return(sorted_df$SUBMITTED_CALLS)

}


# assignAlleleGroups ------------------------------------------------------
#' Find indicies of allele calls
#'
#' \code{assignAlleleGroups} determines the locations of unique alleles within a
#' mixed group. 
#' 
#' @param    allele_calls     \code{character} vector respresenting Ig allele
#'                            calls. Calls may consist of multiple
#'                            comma-separated alleles.
#' @param    allele_min       \code{numeric} indicating the minimum fraction of
#'                            \code{allele_calls} that must contain an allele
#'                            for it to be retained. Integers of 1 or greater
#'                            are interprted as a minumum sequence count.
#' @param    binomial_cutoff  \code{logical} indicating if an \code{allele_min}
#'                            cutoff < 1 should be applied in a binomial manner.
#' @param    alpha            \code{numeric} indicating the alpha cutoff used
#'                            when applying a binomial cutoff of
#'                            \code{allele_min}.
#' @return   \code{list} of indicies in \code{allele_calls} where each unique
#' input allele can be found.
#' 
#' @examples
#' # Create a sample vector of allele calls
#' allele_names = c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
#' "IGHV2-5*01","IGHV1-NL1*01","IGHV1-2*02", "IGHV1-69*02")
#' allele_counts = c(24, 15, 26, 36, 15, 43, 2, 42)
#' alleles = rep(allele_names, allele_counts)
#' 
#' # Find how many of each allele there are
#' assignAlleleGroups(alleles)
#' 
#' @export
assignAlleleGroups <- function(allele_calls, allele_min=1e-4,
                               binomial_cutoff=TRUE, alpha=0.05){
  
  # Find which calls have which alleles
  allele_calls = getAllele(allele_calls, first=F)
  unique_alleles =  sortAlleles(unique(unlist(strsplit(allele_calls, ","))))
  alleles_i = sapply(unique_alleles, grep, allele_calls, fixed=T)
  
  # Determine which alleles are too rare
  if(binomial_cutoff){
    if(allele_min > 1){ error("Allele_min must be < 1 for binomial cutoff") }
    cutoff = qbinom(1-alpha, length(allele_calls), allele_min, lower.tail=FALSE)
  } else if (allele_min < 1){
    cutoff = round(length(allele_calls)*allele_min)
  } else {
    cutoff = allele_min
  }
  allele_counts = sapply(alleles_i, length)
  return(alleles_i[allele_counts >= cutoff])
}


# getMutatedPositions -----------------------------------------------------
#' Find the location of mutations in a sequence
#'
#' \code{getMutatedPositions} takes two vectors of aligned sequences and
#' compares pairs of sequences. It returns a list of the nucleotide positions of
#' any differences.
#' 
#' @param    samples        a vector of strings respresenting aligned sequences
#' @param    germlines      a vector of strings respresenting aligned sequences
#'                          to which \code{samples} will be compared. If only
#'                          one string is submitted, it will be used for all
#'                          \code{samples}.
#' @param    ignored_regex  a regular expression indicating what characters
#'                          should be ignored (such as gaps and N nucleotides).
#' @param    match_instead  if \code{TRUE}, the function returns the positions
#'                          that are the same instead of those that are
#'                          different.
#' @return   A list of the nucleotide positions of any differences between the
#'           input vectors.
#' 
#' @examples
#' # Create strings to act as a sample sequences and a reference sequence
#' seqs = c("----GATA","GAGAGAGA","TANA")
#' ref = "GATAGATA"
#' 
#' # Find the differences between the two
#' getMutatedPositions(seqs, ref)
#' 
#' @export
getMutatedPositions <- function(samples, germlines, ignored_regex="[\\.N-]",
                                match_instead=FALSE) {
  
  # If only one germline sequence is given, use it for all the sample seqs
  if(length(germlines) == 1){ germlines = rep(germlines, length(samples)) }
  if(length(samples) != length(germlines)) {
    stop("Number of input sequences does not match number of germlines.")
  }
  
  # Truncate each pair of sequences to the length of the shorter
  germ_mins = lapply(germlines, nchar)
  samp_mins = lapply(samples, nchar)
  min_lens = mapply(min, germ_mins, samp_mins)
  germ = toupper(mapply(substr, germlines, 1, min_lens, SIMPLIFY=FALSE))
  samp = toupper(mapply(substr, samples, 1, min_lens, SIMPLIFY=FALSE))
  
  # Calculate poisitions of mutations (or matches), ignoring gaps, Ns, and CDR3
  samp_char = strsplit(samp,"")
  germ_char = strsplit(germ,"")
  if(!match_instead){
    muts = lapply(mapply("!=", samp_char, germ_char, SIMPLIFY=FALSE), which)
  } else {
    muts = lapply(mapply("==", samp_char, germ_char, SIMPLIFY=FALSE), which)
  }
  ignore_germ = gregexpr(ignored_regex, germ)
  ignore_samp = gregexpr(ignored_regex, samp)
  ignore = mapply(c, ignore_germ, ignore_samp, SIMPLIFY=FALSE)
  
  muts = mapply(function(x, y) x[!x%in%y], muts, ignore, SIMPLIFY=F)
  return(muts)
}


# summarizeMutations ------------------------------------------------------
#' Find positional mutation counts vs sequence-wide mutation count
#'
#' \code{summarizeMutations} takes the positions of that are similar and that
#' are different between a set of sequences and a germline, and returns a pair
#' of tables summarizing the positional mutation counts.
#' 
#' @param    mut_list   a list of the nucleotide positions of any differences
#'                      between the two vectors of sequences, as generated by
#'                      \code{getMutatedPositions}.
#' @param    match_list a list of the nucleotide positions of any
#'                      similarities between the two vectors of sequences,
#'                      as generated by \code{getMutatedPositions} where
#'                      \code{match_instead = TRUE}.
#' @return   A list containing two matricies. The first details counts of
#'           sequences mutated at given positions (rows) mutated as a
#'           function of sequence-wide mutation count (columns). The second is
#'           details how many usable nucleotides (i.e., not gaps or Ns) were
#'           found for each combination of position and sequence-wide mutation
#'           count.
#' @seealso \code{\link{getMutatedPositions}}
#' 
#' @examples
#' # Create strings to act as a sample sequences and a reference sequence
#' seqs = c("----GATA","GAGAGAGA","TANA")
#' ref = "GATAGATA"
#' 
#' # Find the differences/similarities between the two
#' muts = getMutatedPositions(seqs, ref)
#' matches = getMutatedPositions(seqs, ref, match_instead =TRUE)
#' 
#' # Find positional mutation and nucleotide counts
#' summarizeMutations(muts, matches)
#' 
#' @export
summarizeMutations <- function(mut_list, match_list) { 
  
  pos_list = mapply(c, mut_list, match_list)
  # Find what the longest (usable portion) of sequence is
  most_distant = max(unlist(pos_list))
  pos_range = as.character(1:most_distant)
  # Find the mutation counts of each sequence and their range
  mut_counts = sapply(mut_list, length)
  mut_range = min(unique(mut_counts)):max(unique(mut_counts))
  # Find which sequences have m mutations
  muts_i = sapply(mut_range, function(x) which(x == mut_counts))
  names(muts_i) = mut_range
  # For each sequence-wide mutation count, find the mutations at each position
  pos_muts = sapply(muts_i, function(x) table(unlist(mut_list[x]))[pos_range])
  rownames(pos_muts) = pos_range
  pos_muts[is.na(pos_muts)] = 0
  # For each sequence-wide mutation count, find the usable nts at each position
  pos_pos = sapply(muts_i, function(x) table(unlist(pos_list[x]))[pos_range])
  rownames(pos_pos) = pos_range
  pos_pos[is.na(pos_pos)] = 0
  
  return(list(pos_muts,pos_pos))
}


# trimMutMatrix -----------------------------------------------------------
#' Trim a mutation summary
#'
#' \code{trimMutMatrix} takes a pair of lists as returned by
#' \code{summarizeMutations} and returns a matrix of mutation frequencies at
#' given positions (rows) as a function of sequence-wide mutation count
#' (columns) for the desired ranges of each.
#' 
#' @param    mut_summary    a pair of lists as returned by
#'                          \code{summarizeMutations}
#' @param    mut_min        the minimum number of sequence-wide mutations for
#'                          sequences that will be included in the returned
#'                          matrix
#' @param    mut_max        the maximum number of sequence-wide mutations for
#'                          sequences that will be included in the returned
#'                          matrix
#' @param    nt_min         the first nucleotide position to be included in the
#'                          returned matrix
#' @param    nt_max         the last nucleotide position to be included in the
#'                          returned matrix
#' @param    min_seqs       the minimum number of total sequences (within the
#'                          desired mutational range and nucleotide range)
#'                          required for the function to return a value
#' @param    min_frac       the minimum fraction of sequences that must have
#'                          usable nucleotides in a given position for that
#'                          position to not be made \code{NA}
#' @param    verbose        if \code{TRUE}, a message will be printed when the
#'                          input causes a value of \code{NULL} to be returned
#' @return   A matrix of mutation frequencies at given positions (rows) as a
#'           function of sequence-wide mutation count (columns) for the desired
#'           ranges of each. If there is a problem with the number of sequences,
#'           etc., \code{NULL} will be returned.
#' @seealso \code{\link{summarizeMutations}}
#' 
#' @examples
#' # Create strings to act as a sample sequences and a reference sequence
#' seqs = c("----GATA","GAGAGAGA","GATAGGGA","TANA")
#' ref = "GATAGATA"
#' 
#' # Find the differences/similarities between the two
#' muts = getMutatedPositions(seqs, ref)
#' matches = getMutatedPositions(seqs, ref, match_instead =TRUE)
#' 
#' # Find positional mutation and nucleotide counts
#' mut_mat = summarizeMutations(muts, matches)
#' 
#' # Summarize the frequency for counts above one
#' trimMutMatrix(mut_mat, mut_max=2, nt_max=8, min_seqs=0, min_frac=0)
#' 
#' @export
trimMutMatrix <- function(mut_summary, mut_min=1, mut_max=10,
                          nt_min = 1, nt_max = 312, min_seqs=50, 
                          min_frac = 0.75, verbose=F){
  
  # Ensure the matrix covers the desired mutational range
  if(sum(!(mut_min:mut_max %in% colnames(mut_summary[[1]]))) > 0) {
    if(verbose){
      cat("Insufficient data within mutation range ",
          mut_min, " to ", mut_max,"; skipped.\n", sep="")
    }
    return(NULL)
  }
  
  # Ensure the matrix covers the desired nucleotide range
  if(sum(!(nt_min:nt_max %in% rownames(mut_summary[[1]]))) > 0) {
    if(verbose){
      cat("Insufficient data within nucleotide range ",
          nt_min, " to ", nt_max,"; skipped.\n", sep="")
    }
    return(NULL)
  }
  
  # Trim the matrix to the desired range of mutations
  mut_range = as.character(mut_min:mut_max)
  nt_range = as.character(nt_min:nt_max)
  mut_trim = list(mut_summary[[1]][nt_range,mut_range],
                  mut_summary[[2]][nt_range,mut_range])
  
  # Ensure there are actually sequences to work with in the range
  count_vs_muts = apply(mut_trim[[2]], 2, max, na.rm=TRUE)
  if(sum(count_vs_muts) < min_seqs){
    if(verbose){
      cat("Insufficient data within mutation range ",
          mut_min, " to ", mut_max,"; skipped.\n", sep="")
    }
    return(NULL)
  }
  
  # Find any outliers in the sequence counts.
  outs = boxplot(count_vs_muts,plot=F)$out
  big_outs = outs[which(outs > mean(count_vs_muts))]
  start_at = as.numeric(names(which(count_vs_muts == big_outs)))
  # If you find outliers, trim the matrix again
  if(length(start_at) > 0){
    mut_min = min(5, start_at[1])
    mut_range = as.character(mut_min:mut_max)
    mut_trim = list(mut_trim[[1]][,mut_range],
                    mut_trim[[2]][,mut_range])
  }
  
  mut_fracs = mut_trim[[1]]/mut_trim[[2]]
  
  # Ignore positions that are not found in enough sequences
  position_counts = apply(mut_trim[[2]], 1, sum)
  insufficient_counts = which(position_counts < (min_frac*max(position_counts)))
  mut_fracs[insufficient_counts,] = NaN
  
  return(mut_fracs)
  
}


# findIntercepts ----------------------------------------------------------
#' Find which y-intercepts are above a threshhold
#'
#' \code{findIntercepts} takes a matrix as returned by \code{trimMutMatrix},
#' where mutation frequencies at given positions (rows) are calculated as a
#' function of sequence-wide mutation count (columns), and determines which
#' y-intercepts (positionalal mutation frequency when sequence-wide mutation
#' counts equals zero) are predicted to be above a certain threshhold.
#' 
#' @details The y-intercept is tested by using a t-test to determine the
#'          range of values likely to contain the intercept. If the lower-bound
#'          of this confidence interval is greater than the y-intercept cutoff,
#'          the position will be returned. This is the key step in
#'          \code{\link{detectNovelV}}.
#' 
#' @param    mut_fracs    a matrix as returned by \code{trimMutMatrix}, where
#'                        rows are named with nucleotide position and columns
#'                        are named with sequence-wide mutation count
#' @param    y_intercept  the y-intercept above which positions should be
#'                        returned
#' @param    alpha        the alpha cutoff that should be used in calculating
#'                        the confidence interval for the y-intercept
#' @return   A matrix of mutation frequencies at given positions (rows) as a
#'           function of sequence-wide mutation count (columns) for the desired
#'           ranges of each. If there is a problem with the number of sequences,
#'           etc., \code{NULL} will be returned.
#' @seealso \code{\link{trimMutMatrix}}, \code{\link{detectNovelV}}
#' 
#' @examples
#' # Invent some mutation matrix
#' mut_fracs = matrix(c(sapply(1:10, function(x) max(rnorm(1, x),0)/20),
#'                      sapply(1:10, function(x) max(rnorm(1, x),0)/20),
#'                      sapply(1:10, function(x) max(rnorm(1, x),0)/20),
#'                      sapply(rep(5,10), function(x) max(rnorm(1, x),0)/10)),
#'                    nrow=4, byrow=T)
#' colnames(mut_fracs) = 1:10; rownames(mut_fracs) = 1:4
#' # Test to see if any have a y-intercept above 1/8 (position 4 should)
#' findIntercepts(mut_fracs)
#' 
#' @export
findIntercepts <- function(mut_fracs, y_intercept=1/8, alpha=0.05){
  
  # Clear out any NAs
  mut_fracs[which(is.na(mut_fracs))] = 0
  
  # Calculate y intercepts for each position
  mut_range = as.numeric(colnames(mut_fracs))
  intercepts = apply(mut_fracs, 1, function(x)
    confint(lm(x ~ mut_range), level=1-2*alpha)[1])
  
  # Determine which positions meet the cutoff
  return(intercepts[which(intercepts > y_intercept)])
}


# findNucletoideUsage -----------------------------------------------------
#' Determine nucleotide usage at a given position
#'
#' \code{findNucletoideUsage} determines the nucleotide distribution at a given
#' IMGT-numbered position as a function of sequence-wide mutation counts. 
#' 
#' @param    position   an integer representing the IMGT-numbered position of
#'                      interest
#' @param    samples    a vector of sample nucleotide sequences
#' @param    germline   a string with the germline nucleotide sequence
#' @param    mut_counts a vector containing the mutation count of each sample
#' @param    mut_min    the minimum number of sequence-wide mutations for
#'                      sequences that will be included in the returned matrix
#' @param    mut_max    the maximum number of sequence-wide mutations for
#'                      sequences that will be included in the returned matrix
#' @return   A table of nucleotide usage at a given position as a function of 
#'           sequence-wide mutation counts, with the germline base on top and
#'           the most frequent mutated-to base on the bottom
#' 
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Single out calls utilizing a particular germline sequence
#' germ = germline_ighv[1]
#' matches = grep(names(germ), sample_db$V_CALL, fixed=TRUE)
#' samples = sample_db$SEQUENCE_IMGT[matches]
#' 
#' # Find mutation counts in those sequences
#' mut_counts = getMutCount(samples, rep(names(germ), length(samples)), germ)
#' 
#' # Find nucleotide usage at position 2 as a function of mutation count
#' findNucletoideUsage(2, samples, germ, mut_counts)
#' 
#' @export
findNucletoideUsage <- function(position, samples, germline, mut_counts,
                                mut_min = 1, mut_max = 10){
  # Make everything uppercase
  samples = toupper(samples)
  germline = toupper(germline)
  nts = c("A","C","G","T")
  # Find nucleotide usage at requested position as a function of mutation count
  mut_range = mut_min:mut_max
  mut_counts_i = lapply(mut_range, function(x) which(mut_counts %in% x))
  nucs = lapply(mut_counts_i, function(x) substr(samples[x], position, position))
  nuc_mat = sapply(nucs, function(x) table(x)[nts])
  rownames(nuc_mat) = nts
  colnames(nuc_mat) = mut_range
  nuc_mat[is.na(nuc_mat)] = 0
  # Sort matrix so that the polymorphism is on the bottom and germline on top
  germ_nt = substr(germline, position, position)
  non_germ_nt = setdiff(nts, germ_nt)
  ordered_nt = names(sort(apply(nuc_mat[non_germ_nt,], 1, sum)))
  nuc_mat = nuc_mat[c(germ_nt,ordered_nt),]
  return(nuc_mat)
}


# insertPolymorphisms -----------------------------------------------------
#' Insert polymorphisms into a nucleotide sequence
#'
#' \code{insertPolymorphisms} replaces nucleotides in the desired locations of a
#' provided sequence.
#' 
#' 
#' @param    sequence     the starting nucletide sequence
#' @param    positions    a vector of positions which to be changed
#' @param    nucleotides  a vector of nucletides to which to change the
#'                        positions
#' @return   a sequence with the desired nucleotides in provided locations
#' 
#' @examples
#' insertPolymorphisms("hugged", c(1,2,6), c("t","i","r")) 
#' 
#' @export
insertPolymorphisms <- function(sequence, positions, nucleotides){
  
  if(length(positions) != length(nucleotides)){
    stop("Number of nucleotides and number of positions do not match.")
  }
  names(positions) = nucleotides
  for (i in 1:length(positions)){
    substr(sequence, positions[i], positions[i]) = names(positions[i])
  }
  
  return(sequence)
}


# createGermlines ---------------------------------------------------------
#' Create sequences with each combination of polymorphisms
#'
#' \code{createGermlines} inserts nucleotides in the desired locations of a
#' provided sequence, for each combination of possible insertions.
#' 
#' @details  The nucletoide sequence provided should be named, and will serve
#'           as the basis of the resulting names. For example, a sequence named
#'           IGHV1-2*02 with position 163 mutated to a C (from a T) will be
#'           named 1-2*01_T163C.
#' 
#' @param    sequence     the starting nucletide sequence
#' @param    positions    a vector of positions which to be changed
#' @param    nucleotides  a vector of nucletides to which to change the
#'                        positions
#' @return   Each combination of sequences with the desired nucleotides in
#'           provided locations, with names indicating the insertion(s) in the
#'           form _[germline_nucleotide][position][inserted_nucleotide].
#' 
#' @seealso  \code{\link{insertPolymorphisms}}
#' 
#' @examples
#' # Insert each combination of letters at the listed locations
#' # Note that 3! is 6
#' createGermlines("hugged", c(1,2,6), c("t","i","r")) 
#' 
#' @export
createGermlines <- function(germline, positions, nucleotides){
  allele = names(germline)
  n_pos = length(positions)
  germ_nts = sapply(positions, function(x) substr(germline, x,x))
  # There are 2^n possibilities, but one is no SNPs
  sequences = list()
  for (comb in 1:(2^n_pos-1)){
    # Find the indicies of each possibility of SNPs
    bits = as.numeric(intToBits(comb)[1:n_pos])
    inds = which(as.logical(bits))
    # Figure out what to call the new sequence
    snp_names = mapply(paste, germ_nts[inds], positions[inds],
                       nucleotides[inds], sep="")
    allele_name = paste(allele, paste(snp_names, collapse="_"), sep="_")
    # Introduce the polymorphisms
    nt_seq = insertPolymorphisms(germline, positions[inds], nucleotides[inds])
    sequences[[allele_name]] = nt_seq
  }
  return(unlist(sequences))
}


# findNovelAlleles --------------------------------------------------------
#' Find novel alleles in sequences thought to utilize one particular allele
#'
#' \code{findNovelAlleles} analyzes mutation patterns in sequences thought to
#' align to a particular germline allele in order to determine which positions
#' might be polymorphic.
#' 
#' @details  Mutations are determined through comparison to the provided
#' germline and the mutation frequency at each *position* is determined as a
#' function of *sequence-wide* mutation counts. Polymorphic positions are
#' expected to exhibit a high mutation frequency despite sequence-wide mutation
#' count. False positive of potential novel alleles resulting from clonally-
#' related sequences are guarded against by ensuring that sequences perfectly
#' matching the potential novel allele utilize a wide range of combinations of J
#'  gene and junction length.
#' 
#' @param    samples      a vector of IMGT-gappedsample V sequences thought to
#'                        be utilizing the same germline V allele
#' @param    germline     the germline V sequence utilized by the samples
#' @param    j_genes      a vector of J gene names utilized by the samples
#' @param    junc_lengths a vector of the junction lengths of the sample
#' @param    y_intercept  the y-intercept above which positions should be
#'                        considered potentially polymorphic, as utilized by
#'                        \code{\link{findIntercepts}}
#' @param    nt_min       the first nucleotide position to be considered, as
#'                        utilized by \code{\link{trimMutMatrix}}
#' @param    nt_max       the last nucleotide position to be considered, as
#'                        utilized by \code{\link{trimMutMatrix}}
#' @param    mut_min      the minimum number of sequence-wide mutations for
#'                        sequences that will be used in analysis, as utilized
#'                        by \code{\link{trimMutMatrix}}
#' @param    mut_max      the maximum number of sequence-wide mutations for
#'                        sequences that will be used in analysis, as utilized
#'                        by \code{\link{trimMutMatrix}}
#' @param    j_max        the maximum fraction of sequences perfectly aligning
#'                        to a potential novel allele that are allowed to
#'                        utilize to a particular combination of junction
#'                        length and J gene
#' @param    min_seqs     the minimum number of total sequences (within the
#'                        desired mutational range and nucleotide range)
#'                        required for the samples to be considered, as utilized
#'                        by \code{\link{trimMutMatrix}}
#' @param    min_frac     the minimum fraction of sequences that must have
#'                        usable nucleotides in a given position for that
#'                        position to considered, as utilized by
#'                        \code{\link{trimMutMatrix}}
#' @param    verbose      if \code{TRUE}, a message will be printed when the
#'                        samples do not meet the required parameters
#' @return   For each potential novel allele, a list of length five is returned 
#'           containing (1) the (named) germline sequence, (2) the y-intercepts
#'           of the position(s) which passed the y-intercept threshhold (with
#'           names indicating the positions themselves), (3) a matrix containing
#'           the fraction of sequences mutated at each nucleotide position
#'           (columns) as a function of sequence-wide mutation count (rows), (4)
#'           table(s) indicating the nucletoide usage at each polymorphic
#'           position as a function of mutation count, and (5) a table detailing
#'           the number of unmutated versions of the novel allele found to use
#'           each combination of J gene (columns) and junction length (rows).
#' 
#' @seealso  \code{\link{findIntercepts}}, \code{\link{summarizeMutations}},
#'           \code{\link{trimMutMatrix}}
#' 
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Single out calls utilizing a particular germline sequence and extract info
#' germ = germline_ighv[41]
#' matches = grep(names(germ), sample_db$V_CALL, fixed=TRUE)
#' samples = sample_db$SEQUENCE_IMGT[matches]
#' j_genes = getGene(sample_db$J_CALL[matches])
#' junc_lengths = sample_db$JUNCTION_LENGTH[matches]
#' 
#' # Find novel alleles and return relevant data
#' novel = findNovelAlleles(samples, germ, j_genes, junc_lengths)
#' # Print the nucleotide sequence of the first (if any) novel alleles found
#' novel[[1]][[1]]
#' 
#' @export
findNovelAlleles  <- function(samples, germlines, j_genes, junc_lengths,
                              y_intercept = 1/8, nt_min=1, nt_max=312,
                              mut_min=1, mut_max=10, j_max = 0.1,
                              min_seqs = 50, min_frac = 3/4, verbose=FALSE){
  
  #   args = as.list(match.call())
  #   mcd_options = names(which(nchar(formals(modifyChangeoDb))>0))
  #   
  #   # Functions to which arguments may be passed with ...
  #   gmp_options = names(which(nchar(formals(getMutatedPositions))>0))
  #   tmm_options = names(which(nchar(formals(trimMutMatrix))>0))
  #   fi_options = names(which(nchar(formals(findIntercepts))>0))
  #   fnu_options = names(which(nchar(formals(findNucletoideUsage))>0))
  
  # Find the positions of differences and similarities between sequences
  mut_list = getMutatedPositions(samples, germlines, match_instead = FALSE)
  mut_counts = sapply(mut_list, length)
  match_list = getMutatedPositions(samples, germlines, match_instead = TRUE)
  mut_summary = summarizeMutations(mut_list, match_list)
  mut_matrix = trimMutMatrix(mut_summary, mut_min, mut_max, nt_min, nt_max,
                             min_seqs, min_frac, verbose)
  if (is.null(mut_matrix)){ return(NULL) }
  intercepts = findIntercepts(mut_matrix, y_intercept)
  polymorphs = as.numeric(names(intercepts))
  if (length(polymorphs) > 0) {
    nuc_usage = lapply(polymorphs, findNucletoideUsage, samples, germlines,
                       mut_counts, mut_min, mut_max)
    names(nuc_usage) = polymorphs
    snp_nucs = sapply(nuc_usage, function(x) rownames(x)[4])
    putative = createGermlines(germlines, polymorphs, snp_nucs)
    put_mut_locs = lapply(putative, function(x) getMutatedPositions(samples, x))
    perfect_matches = lapply(put_mut_locs, function(x) which(sapply(x,length) == 0))
    j_junc_tables = lapply(perfect_matches, function(x) table(junc_lengths[x], j_genes[x]))
    pass_test = suppressWarnings(which(sapply(j_junc_tables, function(x) abs(max(x/sum(x),na.rm = T))) < j_max))
    if(length(pass_test) > 0){
      allele_summary = list()
      for(i in 1:length(pass_test)){
        germ_name = names(pass_test[i])
        positions = gsub("[A-Z]","",strsplit(germ_name, "_")[[1]][-1]) 
        allele_summary[[germ_name]] = list(putative[pass_test[i]],
                                           intercepts[positions],
                                           mut_matrix,  nuc_usage[positions],
                                           j_junc_tables[[pass_test[i]]])
      }
      return(allele_summary)
    }
  }
  return(NULL)
}  


# detectNovelV ------------------------------------------------------------
#' Find novel alleles from repertoire sequencing data
#'
#' \code{detectNovelV} analyzes mutation patterns in sequences thought to
#' align to each germline allele in order to determine which positions
#' might be polymorphic.
#' 
#' @details  \code{detectNovelV} applies \code{\link{findNovelAlleles}} to each
#' allele call found in the names of \code{allele_groups}. Mutations are
#' determined through comparison to the provided germline and the mutation
#' frequency at each *position* is determined as a function of *sequence-wide*
#' mutation counts. Polymorphic positions are expected to exhibit a high
#' mutation frequency despite sequence-wide mutation count. False positive of
#' potential novel alleles resulting from clonally-related sequences are guarded
#' against by ensuring that sequences perfectly matching the potential novel
#' allele utilize a wide range of combinations of J gene and junction length.
#' 
#' @param    v_sequences    a vector of IMGT-gapped sample V sequences
#' @param    j_genes        a vector of J gene names utilized by the samples
#' @param    junc_lengths   a vector of the junction lengths of the sample
#' @param    allele_groups  a list whose names match the alle names in
#'                          \code{germline_db} and the contents of which are the
#'                          indicies of \code{v_sequences} that are assigned to
#'                          those alleles. See \code{\link{assignAlleleGroups}}.
#' @param    germline_db    a vector of named nucleotide germline sequences
#'                          matching the calls detailed in \code{allele_groups}          
#' @param    y_intercept    the y-intercept above which positions should be
#'                          considered potentially polymorphic, as utilized by
#'                          \code{\link{findIntercepts}}
#' @param    nt_min         the first nucleotide position to be considered, as
#'                          utilized by \code{\link{trimMutMatrix}}
#' @param    nt_max         the last nucleotide position to be considered, as
#'                          utilized by \code{\link{trimMutMatrix}}
#' @param    mut_min        the minimum number of sequence-wide mutations for
#'                          sequences that will be used in analysis, as utilized
#'                          by \code{\link{trimMutMatrix}}
#' @param    mut_max        the maximum number of sequence-wide mutations for
#'                          sequences that will be used in analysis, as utilized
#'                          by \code{\link{trimMutMatrix}}
#' @param    j_max          the maximum fraction of sequences perfectly aligning
#'                          to a potential novel allele that are allowed to
#'                          utilize to a particular combination of junction
#'                          length and J gene
#' @param    min_seqs       the minimum number of total sequences (within the
#'                          desired mutational range and nucleotide range)
#'                          required for the samples to be considered, as
#'                          utilized by \code{\link{trimMutMatrix}}
#' @param    min_frac       the minimum fraction of sequences that must have
#'                          usable nucleotides in a given position for that
#'                          position to considered, as utilized by
#'                          \code{\link{trimMutMatrix}}
#' @param    verbose        if \code{TRUE}, a message will be printed when the
#'                          samples do not meet the required parameters
#' @return   for each potential novel allele, a list of length five is returned 
#'           containing (1) the (named) germline sequence, (2) the y-intercepts
#'           of the position(s) which passed the y-intercept threshhold (with
#'           names indicating the positions themselves), (3) a matrix containing
#'           the fraction of sequences mutated at each nucleotide position
#'           (columns) as a function of sequence-wide mutation count (rows), (4)
#'           table(s) indicating the nucletoide usage at each polymorphic
#'           position as a function of mutation count, and (5) a table detailing
#'           the number of unmutated versions of the novel allele found to use
#'           each combination of J gene (columns) and junction length (rows).
#' 
#' @seealso  \code{\link{findNovelAlleles}}, \code{\link{findIntercepts}},
#'           \code{\link{summarizeMutations}}, \code{\link{trimMutMatrix}}
#' 
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Single out calls utilizing particular germline sequences and extract info
#' germs = germline_ighv[c(1,41)]
#' matches = lapply(names(germs),grep, sample_db$V_CALL, fixed=TRUE)
#' names(matches) = names(germs)
#' samples = sample_db$SEQUENCE_IMGT
#' j_genes = getGene(sample_db$J_CALL)
#' junc_lengths = sample_db$JUNCTION_LENGTH
#' 
#' # Find novel alleles and return relevant data
#' novel = detectNovelV(samples, j_genes, junc_lengths, matches, germline_ighv)
#' # Print the nucleotide sequence of the first (if any) novel alleles found
#' novel[[1]][[1]]
#' 
#' @export
detectNovelV <- function(v_sequences, j_genes, junc_lengths, allele_groups,
                         germline_db,  y_intercept =1/8, nt_min=1, nt_max = 312,
                         mut_min=1, mut_max=10, j_max = 0.10, min_seqs = 50, min_frac= 1/8,
                         verbose=FALSE, quiet=T){
  
 # args = as.list(match.call())
  
 # args_findNovelAlleles = names(formals(findNovelAlleles))
  
  novel=list()
  for (allele_name in names(allele_groups)) {
    indicies = allele_groups[[allele_name]]
    samples = v_sequences[indicies]
    germline = germline_db[allele_name]
    if (!is.na(germline)){
      if(!quiet){ cat(".") }
      fna = findNovelAlleles(samples, germline,
                             j_genes[indicies],
                             junc_lengths[indicies],
                             y_intercept,
                             nt_min, nt_max,
                             mut_min, mut_max, j_max, min_seqs, min_frac,
                             verbose)
      novel = c(novel, fna)
    }
  }
  return(novel)
}


# getAllele ---------------------------------------------------------------
#' Get Ig segment allele, gene and family names
#' 
#' \code{getSegment} performs generic matching of delimited segment calls with a custom regular 
#' expression. \code{getAllele}, \code{getGene} and \code{getFamily} extract the allele, gene and 
#' family names, respectively, from a character vector of immunoglobulin (Ig) segment allele calls 
#' in IMGT format.
#'
#' @param     segment_call    character vector containing segment calls delimited by commas.
#' @param     segment_regex   string defining the segment match regular expression.
#' @param     first           if \code{TRUE} return only the first call in \code{segment_call};
#'                            if \code{FALSE} return all calls delimited by commas.
#' @param     collapse        if \code{TRUE} check for duplicates and return only unique segment
#'                            assignments; if \code{FALSE} return all assignments (faster). 
#'                            Has no effect if \code{first=TRUE}.
#' @param     sep             character defining both the input and output segment call delimiter.
#' @return    A character vector containing allele, gene or family names
#' 
#' @seealso   Uses \code{\link{str_extract}}.
#' @references
#'   \url{http://imgt.org}
#' @examples
#' kappa_call <- c("Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F", "Homsap IGKJ5*01 F")
#'
#' getAllele(kappa_call)
#' getAllele(kappa_call, first=FALSE)
#' 
#' getGene(kappa_call)
#' getGene(kappa_call, first=FALSE)
#' 
#' getFamily(kappa_call)
#' getFamily(kappa_call, first=FALSE)
#' getFamily(kappa_call, first=FALSE, collapse=TRUE)
#'
#' @export
getSegment <- function(segment_call, segment_regex, first=TRUE, collapse=TRUE, sep=",") {
  # Define boundaries of individual segment calls
  edge_regex <- if (first) { ".*" } else { paste0("[^", sep, "]*") }
  
  # Extract calls
  r <- gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), "\\1", 
            segment_call, perl=T)
  
  # Collapse to unique set if required
  if (!first & collapse) {
    r <- sapply(strsplit(r, sep), function(x) paste(unique(x), collapse=sep))
  }
  
  return(r)
}

#' @rdname getSegment
#' @export
getAllele <- function(segment_call, first=TRUE, collapse=TRUE, sep=",") {
  allele_regex <- '((IG[HLK]|TR[ABGD])[VDJ]\\d+[-/\\w]*[-\\*][\\.\\w]+)'
  r <- getSegment(segment_call, allele_regex, first=first, collapse=collapse, sep=sep)
  
  return(r)
}

#' @rdname getSegment
#' @export
getGene <- function(segment_call, first=TRUE, collapse=TRUE, sep=",") {
  gene_regex <- '((IG[HLK]|TR[ABGD])[VDJ]\\d+[-/\\w]*)'
  r <- getSegment(segment_call, gene_regex, first=first, collapse=collapse, sep=sep)
  
  return(r)
}


#' @rdname getSegment
#' @export
getFamily <- function(segment_call, first=TRUE, collapse=TRUE, sep=",") {
  family_regex <- '((IG[HLK]|TR[ABGD])[VDJ]\\d+)'
  r <- getSegment(segment_call, family_regex, first=first, collapse=collapse, sep=sep)
  
  return(r)
}

