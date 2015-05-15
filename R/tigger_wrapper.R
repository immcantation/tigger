#' Read a germline database
#'
#' \code{readGermlineDb} reads a fasta-formatted file of immunoglobulin (Ig)
#' sequences and returns a named vector of those sequences.
#' 
#' @param    fasta_file       fasta-formatted file of immunoglobuling sequences
#' @param    strip_down_name  if \code{TRUE}, will extract only the allele name
#'                            from the strings fasta file's sequence names
#' @param    force_caps       if \code{TRUE}, will force nucleotides to
#'                            uppercase
#' @return   a named vector of strings respresenting Ig alleles
#' 
#' @examples
#' \dontrun{ 
#' ## Read an imaginary file called "foo.fasta"
#' foo = readGermlineDb("foo.fasta")
#' }
#' 
#' @export
readGermlineDb <- function(fasta_file, 
                           strip_down_name = TRUE,
                           force_caps = TRUE){
  all_char = readChar(fasta_file, file.info(fasta_file)$size)
  split_by_sequence = strsplit(all_char, "[ \t\r\n\v\f]?>")
  add_name_break = sapply(split_by_sequence, function(x) sub("[\r\n]",">",x))
  cleaned_up = sapply(add_name_break, function(x) gsub("[ \t\r\n\v\f]", "", x))
  broken_names = sapply(cleaned_up, strsplit, ">")
  seqs = sapply(broken_names, "[", 2)
  seq_names = sapply(broken_names, "[", 1)
  if(force_caps){ seqs = toupper(seqs) }
  if(strip_down_name){ seq_names = getAllele(seq_names) }
  names(seqs) = seq_names
  return(seqs[which(!is.na(seqs))])
}


# Tool for Immunoglobulin Genotype Elucidation via Rep-Seq

#' Infer genotype (including novel alleles) and correct V calls
#'
#' \code{runTigger} takes a table of sample sequences from a single subject and a
#' vector of database germline sequences. It then performs the following:
#' (1) Infers the presence of novel IGHV alleles not in the germline database.
#' (2) Infers the individuals V genotype.
#' (3) Corrects the IGHV allele calls of the samples based on the IGHV genotype.
#' The sample sequences should be a data frame where each row is a sequence and
#' each column contains data about that sequence. The database germlines should
#' be a vector of sequences with names matching those in the table of sample
#' sequences.
#' 
#' @param    sample_db            a data frame with the colums described in
#'                                Details below.
#' @param    germline_db          a vector of named nucleotide germline
#'                                sequences matching the calls in
#'                                \code{sample_db}
#' @param    find_novel           logical. Should novel alleles be searched for?
#' @param    find_genotype        logical. Should the genotype be inferred?
#' @param    correct_calls        logical. Should the allele calls be corrected?
#' @param    allele_min           a number < 1 representing the minimum fraction
#'                                of sequences required for an allele to not be
#'                                excluded from analysis. or a number >= 1
#'                                representing the minimum count for sequences.
#'                                See \code{\link{assignAlleleGroups}}.
#' @param    y_intercept          the y-intercept above which positions should
#'                                be considered potentially polymorphic. See
#'                                \code{\link{detectNovelV}}.
#' @param    nt_min               the first nucleotide position to be considered
#'                                in intercept calculations. See
#'                                \code{\link{detectNovelV}}.
#' @param    nt_max               the last nucleotide position to be considered
#'                                in intercept calculations. See
#'                                \code{\link{detectNovelV}}.
#' @param    mut_min              the minimum number of mutations carried by
#'                                sequences used in in intercept calculations.
#'                                See \code{\link{detectNovelV}}.
#' @param    mut_max              the maximum number of mutations carried by
#'                                sequences used in in intercept calculations.
#'                                See \code{\link{detectNovelV}}.
#' @param    j_max                the maximum fraction of sequences perfectly
#'                                aligning to a potential novel allele that are
#'                                allowed to utilize to a particular combination
#'                                of junction length and J gene. See
#'                                \code{\link{detectNovelV}}.
#' @param    min_seqs             the minimum number of total sequences (within
#'                                the desired mutational range and nucleotide
#'                                range) required for the samples to be
#'                                analyzed for polymorphisms. See
#'                                \code{\link{detectNovelV}}.
#' @param    min_frac             the maxmium number of total sequences (within
#'                                the desired mutational range and nucleotide
#'                                range) required for the samples to be
#'                                analyzed for polymorphisms. See
#'                                \code{\link{detectNovelV}}.
#' @param    fraction_to_explain  the portion of each gene that must be
#'                                explained by the alleles that will be included
#'                                in the genotype. See
#'                                \code{\link{inferGenotype}}.
#' @param    gene_cutoff          the minimum fraction of the unmutated
#'                                sequences that must be attributed to a gene
#'                                in order for it to be included in the
#'                                genotype. See \code{\link{inferGenotype}}.                               
#' @param    seq_gap              the name of the column in \code{sample_db}
#'                                that includes the IMGT-gapped sequence
#' @param    v_call_col           the name of the column in \code{sample_db}
#'                                that includes the intial V call
#'                                in the column indicated by \code{seq_gap}
#' @param    v_length_col         the name of the column in \code{sample_db}
#'                                that includes the length of the V sequence
#'                                contained within \code{seq_gap}
#' @param    j_call_col           the name of the column in \code{sample_db}
#'                                that includes the initial J call
#' @param    junc_length_col      the name of the column in \code{sample_db}
#'                                that includes the junction length
#' @param    quiet                logical indicating if additional diagonostic
#'                                output will be suppressed
#' 
#' @return   a list containing data on new alleles, the inferred genotype, and
#' the corrected IGHV calls.
#' 
#' @details  The required columns that must be contained within \code{sample_db}
#' are detailed below:
#' \itemize{
#' \item{\code{SEQUENCE_IMGT}: }{V(D)J sequence in the IMGT gapped format}
#' \item{\code{V_CALL}: }{(Comma separated) name(s) of the nearest V allele(s)}
#' \item{\code{J_CALL}: }{(Comma separated) name(s) of the nearest J allele(s)}
#' \item{\code{JUNCTION_LENGTH}: }{Length of the junction region of the V(D)J sample}
#' }
#'
#' @seealso \code{\link{detectNovelV}}, \code{\link{inferGenotype}},
#'          \code{\link{reassignAlleles}}
#' 
#' @references Gadala-Maria D, Yaari G, Uduman M, Kleinstein SH (2015) Automated
#' analysis of high-throughput B cell sequencing data reveals a high frequency
#' of novel immunoglobulin V gene segment alleles. PNAS. 112(8):E862-70
#' 
#' @examples
#' \dontrun{
#' ## Load example data and run all aspects of TIgGER (takes a few minutes)
#' data(sample_db)
#' data(germline_ighv)
#' results = runTigger(sample_db, germline_ighv)
#' 
#' ## Summarize the detected novel alleles, add them to vector of all alleles
#' novel_sequences = novelSummary(results, seqs_to_return = "in genotype")
#' germline_ighv = c(germline_ighv, novel_sequences)
#' ## Plot positional mutation frequency versus sequence-wide mutation count
#' plotNovelLines(results$novel)
#' ## Plot nucleotide usage at polymorphic positions
#' plotNovelBars(results$novel)
#' ## Plot J and junction usage for sequences perfectly matching novel alleles
#' plotJunctionBars(results$novel)
#' 
#' ## View the inferred genotype
#' print(results$genotype)
#' ## Get the nucleotide sequences of all genotype alleles
#' genotype_sequences = genotypeFasta(results$genotype, germline_ighv)
#' 
#' ## Extract the corrected V allele calls and appened them to the data frame
#' V_CALL_GENOTYPED = results$new_calls
#' sample_db = cbind(sample_db, V_CALL_GENOTYPED)
#' }
#' 
#' @export
runTigger <- function(sample_db, germline_db,
                   find_novel = TRUE, find_genotype = TRUE, correct_calls = TRUE,
                   allele_min = 1e-4, y_intercept = 1/8, nt_min=1, nt_max = 312,
                   mut_min=1, mut_max=10, j_max = 0.1, min_seqs = 50, min_frac = 3/4,
                   fraction_to_explain = 7/8,  gene_cutoff = 0.001,
                   seq_gap = "SEQUENCE_IMGT",
                   v_call_col = "V_CALL", j_call_col = "J_CALL",
                   junc_length_col = "JUNCTION_LENGTH", quiet=FALSE){
  
  result = list(novel=NULL, genotype=NULL, new_calls=NULL)
  
  # EXTRACT USEFUL PORTIONS OF DB FILES

  v_calls = getAllele(sample_db[,v_call_col], first = FALSE)
  v_calls = updateAlleleNames(v_calls)
  # New version of sample db files are different, so check the columns names
  seqs = sample_db[,seq_gap]
  v_sequences = sapply(seqs, substr, 1, 312)

  # FIND NOVEL ALLELES

  if (find_novel){
    if(!quiet){ cat("Finding novel alleles...") }
    allele_groups = assignAlleleGroups(v_calls, allele_min)
    j_genes = getGene(sample_db[,j_call_col], first = FALSE)
    junc_lengths = sample_db[,junc_length_col]
    novel = detectNovelV(v_sequences, j_genes, junc_lengths, allele_groups,
                         germline_db,  y_intercept, nt_min, nt_max,
                         mut_min, mut_max, j_max, min_seqs, min_frac, 
                         verbose=FALSE, quiet=quiet)
    # Extract the nucleotide sequence portion
    fasta = unlist(unique(sapply(novel, "[", 1)))
    # In case we found the same allele multiple ways, ditch the duplicate
    fasta = fasta[which(!duplicated(fasta))]
    # Add the novel alleles to the germline_db
    germline_db = c(germline_db, fasta)
    result[["novel"]] = novel
    if(!quiet){ cat("done.\n") }
  }
  
  
  # DETERMINE GENOTYPE
  
  if (find_genotype){
    if(!quiet){ cat("Finding genotype...") }
    v_calls2 = v_calls
    # Paste novel alleles (if any) to all allele calls before determining dists
    if (find_novel & (length(fasta) > 0)){
      genes_novel = getGene(names(fasta))
      genes_groups = getGene(names(allele_groups))
      for(i in 1:length(fasta)){
        matching_groups = allele_groups[genes_groups %in% genes_novel[i]]
        indicies = unique(unlist(matching_groups))
        v_calls2[indicies] = sapply(v_calls2[indicies], paste, names(fasta)[i], sep=",")
      }
    }
    mut_counts = getMutCount(v_sequences, v_calls2, germline_db)
    unmutated_calls = findUnmutatedCalls(v_calls2, mut_counts)
    genotype = inferGenotype(unmutated_calls, fraction_to_explain,
                             gene_cutoff=gene_cutoff)
    genotype_db = genotypeFasta(genotype, germline_db)
    result[["genotype"]] = genotype
    if(!quiet){ cat("done.\n") }
  }
  
  
  # CORRECT ALLELE CALLS
  
  if (correct_calls) {
    if(!quiet){ cat("Correcting allele calls...") }
    if(!find_genotype){ genotype_db = germline_db }
    new_calls = reassignAlleles(v_calls, v_sequences, genotype_db)
    result[["new_calls"]] = new_calls
    if(!quiet){ cat("done.\n") }
  }
  
  return(result)
  
}


#' Return a summary of any novel alleles discovered
#'
#' \code{novelSummary} summaries the output of \code{\link{runTigger}}, stating
#' which novel alleles were included in the genotype. It returns the nucleotide
#' sequences of the novel alleles.
#' 
#' @param    tigger_result  the output of \code{\link{runTigger}}
#' @param    seqs_to_return either \code{"in genotype"} or \code{"all"},
#'                          indicating whether only those potential novel
#'                          alleles alleles in the genotype should be returned
#'                          or if all should be returned
#' @return   a named list of novel allele sequences, as well as text output
#'           indicating what number were detected versus included in the
#'           genotype
#' 
#' @examples
#' \dontrun{
#' ## Load example data and run all aspects of TIgGER (takes a few minutes)
#' data(sample_db)
#' data(germline_ighv)
#' results = runTigger(sample_db, germline_ighv)
#' 
#' ## Summarize the detected novel alleles, add them to vector of all alleles
#' novel_sequences = novelSummary(results, seqs_to_return = "in genotype")
#' germline_ighv = c(germline_ighv, novel_sequences)
#' ## Plot positional mutation frequency versus sequence-wide mutation count
#' plotNovelLines(results$novel)
#' ## Plot nucleotide usage at polymorphic positions
#' plotNovelBars(results$novel)
#' ## Plot J and junction usage for sequences perfectly matching novel alleles
#' plotJunctionBars(results$novel)
#' 
#' ## View the inferred genotype
#' print(results$genotype)
#' ## Get the nucleotide sequences of all genotype alleles
#' genotype_sequences = genotypeFasta(results$genotype, germline_ighv)
#' 
#' ## Extract the corrected V allele calls and appened them to the data frame
#' V_CALL_GENOTYPED = results$new_calls
#' sample_db = cbind(sample_db, V_CALL_GENOTYPED)
#' }
#' 
#' @export
novelSummary <- function(tigger_result,
                         seqs_to_return = c("in genotype", "all")[1]){
 
  r_names = names(which(sapply(tigger_result, length) > 0))
  
  if (!("novel" %in% r_names)){
    warning("Novel alleles not present in input.")
    return(NULL)
  }
  cat(length(tigger_result$novel), "potential novel alleles were detected.\n")
  
  if(seqs_to_return == "in genotype"){

    if(!("genotype" %in% r_names)){ stop("Genotype not present in input.") }
    splits = strsplit(tigger_result$genotype$alleles, ",")
    alleles = mapply(paste, tigger_result$genotype$gene, splits, sep="*")
    novel_in_geno = grep("_", unlist(alleles), value = TRUE)
    all_seqs = unlist(sapply(tigger_result$novel, "[", 1))
    all_names = names(tigger_result$novel)
    which_in_geno = which(gsub("D", "", all_names) %in% novel_in_geno)
    cat(length(which_in_geno),
        "novel alleles were common enough to be included in the genotype:\n",
        paste(all_names[which_in_geno], collapse=", "),"\n")
    seqs = all_seqs[which_in_geno]
    names(seqs) = all_names[which_in_geno]

  } else if (seqs_to_return == "all") {
  
    seqs = unlist(sapply(tigger_result$novel, "[", 1))
    names(seqs) = names(tigger_result$novel)

  }
  
  return(seqs)

}


#' Visualization of positional mutation frequencies
#'
#' \code{plotNovelLines} plots the mutation frequency of nucleotide positions as
#' a function of sequence-wide mutation count. Potentially polymorphic positions
#' are highlighted in red.
#' 
#' @param    novel  a list of the type returned by \code{\link{detectNovelV}}
#' 
#' @return   plot(s) the mutation frequency of nucleotide positions as
#' a function of sequence-wide mutation count.
#' 
#' @seealso  \code{\link{detectNovelV}}, \code{\link{runTigger}}
#' 
#' @examples
#' \dontrun{
#' ## Load example data and run all aspects of TIgGER (takes a few minutes)
#' data(sample_db)
#' data(germline_ighv)
#' results = runTigger(sample_db, germline_ighv)
#' 
#' ## Summarize the detected novel alleles, add them to vector of all alleles
#' novel_sequences = novelSummary(results, seqs_to_return = "in genotype")
#' germline_ighv = c(germline_ighv, novel_sequences)
#' ## Plot positional mutation frequency versus sequence-wide mutation count
#' plotNovelLines(results$novel)
#' ## Plot nucleotide usage at polymorphic positions
#' plotNovelBars(results$novel)
#' ## Plot J and junction usage for sequences perfectly matching novel alleles
#' plotJunctionBars(results$novel)
#' 
#' ## View the inferred genotype
#' print(results$genotype)
#' ## Get the nucleotide sequences of all genotype alleles
#' genotype_sequences = genotypeFasta(results$genotype, germline_ighv)
#' 
#' ## Extract the corrected V allele calls and appened them to the data frame
#' V_CALL_GENOTYPED = results$new_calls
#' sample_db = cbind(sample_db, V_CALL_GENOTYPED)
#' }
#' 
#' @export
plotNovelLines <- function(novel){
  for(n in novel){
    plot(NA, xlim = c(0,10), ylim = c(0,1), main = names(n[[1]]), las=1,
         xlab="Mutation Count (Sequence)",
         ylab="Mutation Frequency (Position)")
    apply(n[[3]], 1, function(x) lines(as.numeric(names(x)),x) )
    for (i in 1:length(n[[2]])){
      lines(colnames(n[[3]]),n[[3]][as.numeric(names(n[[2]]))[i],], col="red")
    }
    for (i in 1:length(n[[2]])){
      text(as.numeric(colnames(n[[3]])[1])+1-i,
           n[[3]][as.numeric(names(n[[2]]))[i],1],
           labels=names(n[[2]])[i],
           col="red",adj=c(1,1))
    }
  }
}


#' Visualization of nucleotide usage
#'
#' \code{plotNovelBars} shows the nucleotide usage at polymorphic positions as a
#' function of sequence-wide mutation count.
#' 
#' @param    novel  a list of the type returned by \code{\link{detectNovelV}}
#' 
#' @return   plot(s) of nucleotide usage at polymorphic positions as a
#' function of sequence-wide mutation count.
#' 
#' @seealso  \code{\link{detectNovelV}}, \code{\link{runTigger}}
#' 
#' @examples
#' \dontrun{
#' ## Load example data and run all aspects of TIgGER (takes a few minutes)
#' data(sample_db)
#' data(germline_ighv)
#' results = runTigger(sample_db, germline_ighv)
#' 
#' ## Summarize the detected novel alleles, add them to vector of all alleles
#' novel_sequences = novelSummary(results, seqs_to_return = "in genotype")
#' germline_ighv = c(germline_ighv, novel_sequences)
#' ## Plot positional mutation frequency versus sequence-wide mutation count
#' plotNovelLines(results$novel)
#' ## Plot nucleotide usage at polymorphic positions
#' plotNovelBars(results$novel)
#' ## Plot J and junction usage for sequences perfectly matching novel alleles
#' plotJunctionBars(results$novel)
#' 
#' ## View the inferred genotype
#' print(results$genotype)
#' ## Get the nucleotide sequences of all genotype alleles
#' genotype_sequences = genotypeFasta(results$genotype, germline_ighv)
#' 
#' ## Extract the corrected V allele calls and appened them to the data frame
#' V_CALL_GENOTYPED = results$new_calls
#' sample_db = cbind(sample_db, V_CALL_GENOTYPED)
#' }
#' 
#' @export
plotNovelBars <- function(novel){
  NUC_COLORS = c("#64F73F", "#FFB340", "#EB413C", "#3C88EE")
  names(NUC_COLORS) = c("A","C","G","T")
  
  for(n in novel){
    for (j in 1:length(n[[4]])){
      p = n[[4]][[j]]
      plot_name = paste(strsplit(names(n[[1]]),"_")[[1]][1])
      barplot(p[4:1,], col = NUC_COLORS[rownames(p)], las=1,
              xlab = "Mutation Count (Sequence)",
              ylab = "Sequence Count",
              main=plot_name)
      box()
      leg_names = paste(rownames(p), c("(germline)","","","(polymorphism)"))
      legend("topright",
             title = paste("Nucleotide at Position", names(n[[4]][j])),
             legend = leg_names,
             fill=NUC_COLORS[rownames(p)][4:1],
             bty = "n")
    }
  }
}


#' Visualization of J gene usage and junction length
#'
#' \code{plotJunctionBars} shows the frequency of each combination of J gene
#' junction length found among sequences representing unmutated versions of
#' potential novel alleles.
#' 
#' @param    novel  a list of the type returned by \code{\link{detectNovelV}}
#' 
#' @return   plot(s) of the frequency of each combination of J gene and
#' junction length among sequences using potential novel alleles
#' 
#' @seealso  \code{\link{detectNovelV}}, \code{\link{runTigger}}
#' 
#' @examples
#' \dontrun{
#' ## Load example data and run all aspects of TIgGER (takes a few minutes)
#' data(sample_db)
#' data(germline_ighv)
#' results = runTigger(sample_db, germline_ighv)
#' 
#' ## Summarize the detected novel alleles, add them to vector of all alleles
#' novel_sequences = novelSummary(results, seqs_to_return = "in genotype")
#' germline_ighv = c(germline_ighv, novel_sequences)
#' ## Plot positional mutation frequency versus sequence-wide mutation count
#' plotNovelLines(results$novel)
#' ## Plot nucleotide usage at polymorphic positions
#' plotNovelBars(results$novel)
#' ## Plot J and junction usage for sequences perfectly matching novel alleles
#' plotJunctionBars(results$novel)
#' 
#' ## View the inferred genotype
#' print(results$genotype)
#' ## Get the nucleotide sequences of all genotype alleles
#' genotype_sequences = genotypeFasta(results$genotype, germline_ighv)
#' 
#' ## Extract the corrected V allele calls and appened them to the data frame
#' V_CALL_GENOTYPED = results$new_calls
#' sample_db = cbind(sample_db, V_CALL_GENOTYPED)
#' }
#' 
#' @export
plotJunctionBars <- function(novel){
  for(n in novel){
    barplot(n[[5]], las=1, col = rainbow(nrow(n[[5]])),
            main = names(n[[1]]))
    legend("topleft", legend = rownames(n[[5]]), ncol = 3,
           title = "Junction Length",
           fill = rainbow(nrow(n[[5]])), bty="n")
    box()
  }
}


#' Standardize Sample Db data
#'
#' \code{modifyChangeoDb} takes a Change-O Sample Db and modifies it for use
#' with TIgGER.
#' 
#' @param  sample_db            A Change-O db data frame.
#' @param  seq_imgt_col         The name of the column in \code{sample_db}
#'                              containing the IMGT-gapped nucleotide sequence.
#' @param  v_call_col           The name of the column in \code{sample_db}
#'                              containingthe IMGT-assigned V allele call.
#' @param  j_call_col           The name of the column in \code{sample_db}
#'                              containing the IMGT-assigned J allele call.
#' @param  junc_len_col         The name of the column in \code{sample_db} 
#'                              containing the junction length.
#' @param  func_col             The name of the column in \code{sample_db}
#'                              indicating if the sequence is functional.
#' @param  cols_to_add          One or more of \code{"V_GENE"}, \code{"J_GENE"}, or
#'                              \code{"V_SEQUENCE_IMGT"}, indicating which columns
#'                              should be added. See details for more information.
#' @param  distinct_seq_only    Logical indicating if duplicate IMGT-gapped
#'                              sequences should be removed.
#' @param  functional_seq_only  Logical indicating if nonfunctional sequences
#'                              should be removed.
#' @param  nonempty_seq_only    Logical indicating if empty sequences should be
#'                              removed.
#' @param  standardize_calls    Logical indicating if IMGT-assigned allele calls
#'                              should be standardized.
#' @param  factor2char          Logical indicating if columns of type
#'                              \code{factor} should be converted to type
#'                              \code{character}.
#' @details  The supplied column names will be renamed to the current preferred
#' names (and utilized for the creation of new columns, if requested). Columns
#' \code{"V_GENE"}, \code{"J_GENE"}, and \code{"V_SEQUENCE_IMGT"}, if added,
#' will respectively contain the IMGT-assigned V gene call, the IMGT-assigned J
#' gene call, and the first 312 nucleotides (FWR1, CDR1, FWR2, CDR2, and FRW3)
#' of the IMGT-gapped nucleotide sequence.
#' 
#' @return  A corrected Change-O Db data frame
#' 
#' @examples
#' data(sample_db)
#' corrected_sample_db = modifySampleDb(sample_db, cols_to_add = c("V_GENE"))
#' 
#' @export
modifyChangeoDb <- function(sample_db,
                           seq_imgt_col = "SEQUENCE_IMGT",
                           v_call_col = "V_CALL",
                           j_call_col = "J_CALL",
                           junc_len_col = "JUNCTION_LENGTH", 
                           func_col = "FUNCTIONAL",
                           cols_to_add = c("V_GENE", "J_GENE", "V_SEQUENCE_IMGT"),
                           distinct_seq_only = TRUE,
                           nonempty_seq_only = TRUE,
                           functional_seq_only = FALSE,
                           standardize_calls = TRUE,
                           factor2char = TRUE
){
  
  # Convert factors to characters
  if(factor2char){
    i <- sapply(sample_db, is.factor)
    sample_db[,i] <- sapply(sample_db[,i], as.character)
  }
  
  # Rename columns, if needed
  sample_db = data.frame(sample_db, stringsAsFactors = FALSE)
  names(sample_db)[which(names(sample_db) == seq_imgt_col)] = "SEQUENCE_IMGT"
  names(sample_db)[which(names(sample_db) == v_call_col)] = "V_CALL"
  names(sample_db)[which(names(sample_db) == j_call_col)] = "J_CALL"
  names(sample_db)[which(names(sample_db) == junc_len_col)] = "JUNCTION_LENGTH"
  names(sample_db)[which(names(sample_db) == func_col)] = "FUNCTIONAL"
  
  # Filter out duplicates and non-functional seqs, if requested
  if(distinct_seq_only){ sample_db = distinct(sample_db, SEQUENCE_IMGT) }
  if(functional_seq_only & ("FUNCTIONAL" %in% names(sample_db))){
    sample_db = sample_db %>%
      mutate(FUNCTIONAL = FUNCTIONAL=="T" | FUNCTIONAL == "TRUE" | FUNCTIONAL == T) %>%
      filter(FUNCTIONAL == "TRUE")
  }
  # Standardize allele calls
  if(standardize_calls){
    sample_db = mutate(sample_db, V_CALL = updateAlleleNames(getAllele(V_CALL, first = F)))
    sample_db = mutate(sample_db, J_CALL = getAllele(J_CALL, first = F))
  }
  
  # Add columns, if requested
  if (("V_GENE" %in% cols_to_add) & !("V_GENE" %in% names(sample_db))){
    sample_db = mutate(sample_db, V_GENE = getGene(V_CALL))
  }
  if (("J_GENE" %in% cols_to_add) & !("J_GENE" %in% names(sample_db))){
    sample_db = mutate(sample_db, J_GENE = getGene(J_CALL))
  }
  if (("V_SEQUENCE_IMGT" %in% cols_to_add) & !("V_SEQUENCE_IMGT" %in% names(sample_db))){
    sample_db = mutate(sample_db, V_SEQUENCE_IMGT = substring(SEQUENCE_IMGT, 1, 312))
  }
  if (nonempty_seq_only){
    sample_db = filter(sample_db, nchar(SEQUENCE_IMGT) > 2)
  }
  
  return(sample_db)
  
}


#' Find Frequent Sequences' Mutation Counts
#'
#' \code{getPopularMutationCount} determines which sequences occur frequently
#' for each V gene and returns the mutation count of those sequences.
#' 
#' @param  sample_db     A Change-O db data frame. See \code{\link{runTigger}}
#'                       for a list of required columns.
#' @param  germline_db   A named list of IMGT-gapped germline sequences.
#' @param  gene_min      The portion of all unique sequences a gene must
#'                       constitute to avoid exclusion.
#' @param  seq_min       The number of copies of the V that must be present for
#'                       to avoid exclusion.
#' @param  seq_p_of_max  For each gene, fraction of the most common V sequence's
#'                       count that a sequence must meet to avoid exclusion.
#' @param  full_return   If true, will return all \code{sample_db} columns and
#'                       will include sequences with mutation count < 1.
#' @param  ...           Additional arguments to pass to
#'                       \code{\link{modifyChangeoDb}} prior to computation.
#' 
#' @return  A data frame of genes that have a frequent sequence mutation count
#' above 1.
#' 
#' @examples
#' data(sample_db, germline_ighv)
#' getPopularMutationCount(sample_db, germline_ighv)
#' 
#' @export
getPopularMutationCount <- function(sample_db,
                                    germline_db,
                                    gene_min = 1e-03,
                                    seq_min = 50,
                                    seq_p_of_max = 1/8,
                                    full_return = FALSE,
                                    ...){
  
  # Process dots for later
  args = as.list(match.call())
  mcd_options = names(which(nchar(formals(modifyChangeoDb))>0))
  
  modified_db = sample_db %>%
    # Correct common formatting problems, remove duplicates, add V_GENE column
    modifyChangeoDb(args[names(args) %in% mcd_options]) %>%
    # Count occurence of each unique IMGT-gapped V sequence
    group_by(V_GENE, V_SEQUENCE_IMGT) %>%
    mutate(V_SEQUENCE_IMGT_N = n()) %>%
    # Count occurence of each gene and determine count of most common sequence
    group_by(V_GENE) %>%
    mutate(V_GENE_N = n()) %>%
    mutate(V_SEQUENCE_IMGT_N_MAX = max(V_SEQUENCE_IMGT_N)) %>%
    # Remove rare V genes, rare sequences, and sequences not making up a
    # sufficient proportion of sequences as compared to the most common
    ungroup() %>%
    distinct(V_SEQUENCE_IMGT) %>%
    filter(V_GENE_N >= (nrow(sample_db)*gene_min)) %>%
    filter(V_SEQUENCE_IMGT_N >= seq_min) %>%
    mutate(V_SEQUENCE_IMGT_P_MAX = V_SEQUENCE_IMGT_N/V_SEQUENCE_IMGT_N_MAX) %>%
    filter(V_SEQUENCE_IMGT_P_MAX >= seq_p_of_max)
  
  # Determine the mutation counts of the V sequences and append them to the db
  MUTATION_COUNT = getMutCount(modified_db$V_SEQUENCE_IMGT,
                               modified_db$V_CALL,
                               germline_db) %>% 
    sapply(function(x) min(unlist(x)))
  merged_db = bind_cols(modified_db, data.frame(MUTATION_COUNT))
  
  # Strip down the data frame before returning it
  if (!full_return) {
    merged_db = merged_db %>%
      filter(MUTATION_COUNT > 0) %>%
      select(V_GENE, MUTATION_COUNT)
  }
  
  return(merged_db)
  
}














#          __  _-==-=_,-.
#         /--`' \_@-@.--<
#         `--'\ \   <___/.      The wonderful thing about Tiggers,
#             \ \\   " /        is Tiggers are wonderful things.
#               >=\\_/`<        Their tops are made out of rubber,
#   ____       /= |  \_/        their bottoms are made out of springs.
# _'    `\   _/=== \__/         They're bouncy, trouncy, flouncy, pouncy,
# `___/ //\./=/~\====\          Fun, fun, fun, fun, fun.
#     \   // /   | ===:         But the most wonderful thing about Tiggers is,
#      |  ._/_,__|_ ==:        __  I'm the only one.
#       \/    \\ \\`--|       / \\
#        |    _     \\:      /==:-\
#        `.__' `-____/       |--|==:
#           \    \ ===\      :==:`-'
#           _>    \ ===\    /==/
#          /==\   |  ===\__/--/
#         <=== \  /  ====\ \\/
#         _`--  \/  === \/--'
#        |       \ ==== |
#         -`------/`--' /
#                 \___-'
