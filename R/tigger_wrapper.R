# readGermlineDb ----------------------------------------------------------
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
#' @param    v_start_col          the name of the column in \code{sample_db}
#'                                that includes where the V nucleotides begin
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
#' \item{\code{SEQUENCE_GAP}: }{V(D)J sequence gapped in the IMGT gapped format}
#' \item{\code{V_CALL}: }{(Comma separated) name(s) of the nearest V allele(s)}
#' \item{\code{V_GERM_START}: }{Position in the germline sequence where sample V starts}
#' \item{\code{V_GAP_LENGTH}: }{Length (including gaps) of V sequence in \code{SEQUENCE_GAP}}
#' \item{\code{J_CALL}: }{(Comma separated) name(s) of the nearest J allele(s)}
#' \item{\code{JUNCTION_GAP_LENGTH}: }{Length of the junction region of the V(D)J sample}
#' }
#'
#' @seealso \code{\link{detectNovelV}}, \code{\link{inferGenotype}},
#'          \code{\link{reassignAlleles}}
#' 
#' @references http://clip.med.yale.edu/tigger/
#' 
#' @export
runTigger <- function(sample_db, germline_db,
                   find_novel = TRUE, find_genotype = TRUE, correct_calls = TRUE,
                   allele_min = 1e-4, y_intercept = 1/8, nt_min=1, nt_max = 312,
                   mut_min=1, mut_max=10, j_max = 0.1, min_seqs = 50, min_frac = 3/4,
                   fraction_to_explain = 7/8,  gene_cutoff = 0.001,
                   seq_gap = "SEQUENCE_GAP",
                   v_call_col = "V_CALL", v_start_col = "V_GERM_START",
                   v_length_col = "V_GAP_LENGTH", j_call_col = "J_CALL",
                   junc_length_col = "JUNCTION_GAP_LENGTH", quiet=FALSE){
  
  result = list(novel=NULL, genotype=NULL, new_calls=NULL)
  
  # EXTRACT USEFUL PORTIONS OF DB FILES

  v_calls = getAllele(sample_db[,v_call_col], first = FALSE)
  v_calls = updateAlleleNames(v_calls)
  # New version of sample db files are different, so check the columns names
  seqs = sample_db[,seq_gap]
  if ("V_GAP_LENGTH" %in% colnames(sample_db)){
    v_sequences = sapply(seqs, substr, 1, sample_db$V_GAP_LENGTH)
  } else {
    v_sequences = sapply(seqs, substr, 1, sample_db$V_GERM_START + sample_db$V_GERM_LENGTH)
  }
  
  
  # FIND NOVEL ALLELES

  if (find_novel){
    if(!quiet){ cat("Finding novel alleles...") }
    allele_groups = assignAlleleGroups(v_calls, allele_min)
    j_genes = getGene(sample_db[,j_call_col], first = FALSE)
    junc_lengths = sample_db[,junc_length_col]
    novel = detectNovelV(v_sequences, j_genes, junc_lengths, allele_groups,
                         germline_db,  y_intercept, nt_min, nt_max,
                         mut_min, mut_max, j_max, min_seqs, min_frac, 
                         verbose=FALSE)
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
    if (find_novel & (length(novel) > 0)){
      genes_novel = getGene(names(novel))
      genes_groups = getGene(names(allele_groups))
      for(i in 1:length(novel)){
        matching_groups = allele_groups[genes_groups %in% genes_novel[i]]
        indicies = unique(unlist(matching_groups))
        v_calls2[indicies] = sapply(v_calls2[indicies], paste, names(novel)[i], sep=",")
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



# novelSummary ------------------------------------------------------------
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
#' @export
novelSummary <- function(tigger_result,
                         seqs_to_return = c("in genotype", "all")[1]){
 
  r_names = names(which(sapply(tigger_result, length) > 0))
  
  if (!("novel" %in% r_names)){ stop("Novel alleles not present in input.") }
  cat(length(tigger_result$novel), "potential novel alleles were detected.\n")
  
  if(seqs_to_return == "in genotype"){

    if(!("genotype" %in% r_names)){ stop("Genotype not present in input.") }
    splits = strsplit(tigger_result$genotype$alleles, ",")
    alleles = mapply(paste, tigger_result$genotype$gene, splits, sep="*")
    novel_in_geno = grep("_", unlist(alleles), value = TRUE)
    cat(length(novel_in_geno),
        "novel alleles were common enough to be included in the genotype:\n",
        paste(novel_in_geno, collapse=", "),"\n")
    seqs = unlist(sapply(tigger_result$novel[novel_in_geno], "[", 1))
    names(seqs) = novel_in_geno

  } else if (seqs_to_return == "all") {
  
    seqs = unlist(sapply(tigger_result$novel, "[", 1))
    names(seqs) = names(tigger_result$novel)

  }
  
  return(seqs)

}

#




















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
