# updateAlleleNames -------------------------------------------------------
#' Update IGHV allele names
#'
#' \code{updateAlleleNames} takes a set of IGHV allele calls and replaces any
#' outdated names (e.g. IGHV1-f) with the new IMGT names.
#' @details The updated allele names are based on IMGT release 201408-4.
#' @note    IGMT has removed IGHV2-5*10 and IGHV2-5*07 as it has determined they
#'          are actually alleles *02 and *04, respectively.
#' 
#' @param    allele_calls  a vector of strings respresenting IGHV allele names
#' 
#' @return   A vector of strings respresenting updated IGHV allele names
#' 
#' @references Xochelli et al. (2014) Immunogenetics.
#' 
#' @examples
#' # Create a vector that uses old gene/allele names.
#' alleles = c("IGHV1-c*01", "IGHV1-f*02", "IGHV2-5*07")
#' 
#' # Update the alleles to the new names
#' updateAlleleNames(alleles)
#' 
#' @export
updateAlleleNames <- function(allele_calls){
  temporary_names = c("IGHV1-c*",
                      "IGHV1-f*",
                      "IGHV3-d*",
                      "IGHV3-h*",
                      "IGHV4-b*",
                      "IGHV5-a*",
                      "IGHV2-5*10",
                      "IGHV2-5*07")
  definitive_names = c("IGHV1-38-4*",
                       "IGHV1-69-2*",
                       "IGHV3-38-3*",
                       "IGHV3-69-1*",
                       "IGHV4-38-2*",
                       "IGHV5-10-1*",
                       "IGHV2-5*02",
                       "IGHV2-5*04")
  for (i in 1:length(temporary_names)){
    allele_calls = gsub(temporary_names[i], definitive_names[i], allele_calls,
                        fixed = TRUE)
  }
  return(allele_calls)
}


# getMutCount -------------------------------------------------------------
#' Determine the mutation counts from allele calls
#'
#' \code{getMutCount} takes a set of nucleotide sequences and their allele calls
#' and determines the distance between that seqeunce and any germline alleles
#' contained within the call
#' 
#' @param    samples       a vector of IMGT-gapped sample V sequences
#' @param    allele_calls  a vector of strings respresenting Ig allele calls for
#'                         the sequences in \code{samples}, where multiple
#'                         calls are separated by a comma
#' @param    germline_db   a vector of named nucleotide germline sequences
#'                         matching the calls detailed in \code{allele_calls}
#' 
#' @return   A list equal in length to \code{samples}, containing the Hamming
#'           distance to each germline allele contained within each call within
#'           each element of \code{samples}
#' 
#' @examples
#' # Load germline database
#' data(germline_ighv)
#' 
#' # Use createGermlines to insert a mutation into a germline sequence
#' sample_seqs = c(germline_ighv[2],
#'                 createGermlines(germline_ighv[1], 103, "G"),
#'                 createGermlines(germline_ighv[1], 107, "C"))
#' 
#' # Pretend that one sample sequence has received an ambiguous allele call
#' sample_alleles = c(paste(names(germline_ighv[1:2]), collapse=","),
#'                   names(germline_ighv[2]),
#'                   names(germline_ighv[1]))
#' 
#' # Compare each sequence to its assigned germline(s) to determine the distance
#' getMutCount(sample_seqs, sample_alleles, germline_ighv)
#' 
#' @export
getMutCount <- function(samples, allele_calls, germline_db){
  
  call_list = strsplit(allele_calls, ",")
  germline_list = lapply(call_list, function(x) germline_db[x])
  
  mut_pos_list = list()
  mut_count_list = list()
  # First, find mutations of all sequences with call count of 1
  call_count = sapply(germline_list, length)
  cc1 = which(call_count == 1)
  if (length(cc1) > 0) {
    mut_pos_list[cc1] = getMutatedPositions(samples[cc1],
                                            unlist(germline_list[cc1]))
    mut_count_list[cc1] = lapply(mut_pos_list[cc1], length)
  }
  # Then find mutations of all sequences with call count > 1
  ccm = which(call_count > 1)
  if (length(ccm) > 0){
    mut_pos_list[ccm] = mapply(getMutatedPositions,
                               germline_list[ccm], samples[ccm],
                               SIMPLIFY=FALSE)
    mut_count_list[ccm] = lapply(mut_pos_list[ccm],
                                 function(x) lapply(x,length))
  }
  
  return(mut_count_list)
}



# findUnmutatedCalls ------------------------------------------------------
#' Determine which calls represent an unmutated allele
#'
#' \code{findUnmutatedCalls} determines which allele calls would represent a 
#' perfect match with the germline sequence, given a vector of allele calls and
#' mutation counts. In the case of multiple alleles being assigned to a
#' sequence, only the subset that would represent a perfect match is returned.
#' 
#' @param    allele_calls   a vector of strings respresenting Ig allele calls,
#'                          where multiple calls are separated by a comma
#' @param    mut_counts     a list containing distance to each germline allele
#'                          call within \code{allele_calls}, as returned by
#'                          \code{\link{getMutCount}}
#' @param    only_unmutated if \code{TRUE}, calls where no allele that would
#'                          represent an unmutated sequence will be omitted from
#'                          the output
#' 
#' @return   A vector of strings containing the members of \code{allele_calls}
#'           that represent unmutated sequences
#' 
#' @examples
#' # Load germline database
#' data(germline_ighv)
#' 
#' # Use createGermlines to insert a mutation into a germline sequence
#' sample_seqs = c(germline_ighv[2],
#'                 createGermlines(germline_ighv[1], 103, "G"),
#'                 germline_ighv[1],
#'                 germline_ighv[2])
#' 
#' # Pretend that one sample sequence has received an ambiguous allele call
#' sample_alleles = c(paste(names(germline_ighv[1:2]), collapse=","),
#'                   names(germline_ighv[2]),
#'                   names(germline_ighv[1]),
#'                   names(germline_ighv[2]))
#' 
#' # Compare the sequence to a subset of the germlines
#' mut_counts = getMutCount(sample_seqs, sample_alleles, germline_ighv)
#' 
#' # Find which of the sample alleles are unmutated
#' findUnmutatedCalls(sample_alleles, mut_counts)
#' 
#' @export
findUnmutatedCalls <- function(allele_calls, mut_counts, only_unmutated = TRUE){
  
  # Find which seqs are unmutated and which of the allele calls that represents
  unmut_i = which(sapply(mut_counts, function(x) min(unlist(x))) == 0)
  which_no_muts = sapply(mut_counts, function(x) grep("^0$", unlist(x)) )
  unmut_alleles = rep("", length(allele_calls))
  
  # How many alleles represent perfect matches?
  n_gl_unmut = sapply(which_no_muts, length)
  
  one_unmut = which(n_gl_unmut == 1)
  split_names = strsplit(allele_calls, ",")
  if (length(one_unmut) > 0){
    inds = unlist(which_no_muts[one_unmut])
    unmut_alleles[one_unmut] = mapply("[", split_names[one_unmut], inds)
  }
  
  more_unmut = which(n_gl_unmut > 1)
  if (length(more_unmut) > 0){
    inds = which_no_muts[more_unmut]      
    unmut_multi = mapply(function(x,y) x[unlist(y)], split_names[more_unmut],
                         inds, SIMPLIFY = FALSE)
    unmut_alleles[more_unmut] = sapply(unmut_multi, paste, collapse=",")  
  }
  
  if(only_unmutated){unmut_alleles = unmut_alleles[unmut_i] }
  
  return(unmut_alleles)
  
}



# inferGenotype -----------------------------------------------------------
#' Infer a subject-specific genotype
#'
#' \code{inferGenotype} infers an subject's genotype by finding the minimum
#' number set of alleles that can explain the majority of each gene's calls. The
#' most common allele of each gene is included in the genotype first, and the
#' next most common allele is added until the desired fraction of alleles can be
#' explained. In this way, mistaken allele calls (resulting from sequences which
#' by chance have been mutated to look like another allele) can be removed.
#' 
#' @details  Allele calls representing cases where multiple alleles have been
#'           assigned to a single sample sequence are rare among unmutated
#'           sequences but may result if nucleotides for certain positions are
#'           not available. Calls containing multiple alleles are treated as
#'           belonging to all groups until one of those groups is included in
#'           the genotype.
#' 
#' @param    allele_calls         a vector of strings respresenting Ig allele
#'                                calls of unmutated sequences from a single
#'                                subject
#' @param    fraction_to_explain  the portion of each gene that must be
#'                                explained by the alleles that will be included
#'                                in the genotype
#' @param    gene_cutoff          either a number of sequences or a fraction of
#'                                the length of \code{allele_calls} denoting the
#'                                minimum number of times a gene must be
#'                                observed in \code{allele_calls} to be included
#'                                in the genotype
#' 
#' @return   A table of alleles denoting the genotype of the subject
#' 
#' @note     This method works best with data derived from blood, where a large
#'           portion of sequences are expected to be unmutated. Ideally, there
#'           should be hundreds of allele calls per gene in the input.
#' 
#' @examples
#' # Load example data; we'll pretend allele calls are unmutated
#' data(pgp1_example)
#' 
#' # Infer the V genotype
#' inferGenotype(pgp1_example[,"V_CALL"])
#' 
#' # Inger the J genotype
#' inferGenotype(pgp1_example[,"J_CALL"])
#' 
#' @export
inferGenotype <- function(allele_calls, # Calls of unique, unmutated sequences
                          fraction_to_explain = 7/8,
                          gene_cutoff = 1e-3 # Can be a no. of seqs or frac < 1
){
  # Standardize allele call names
  allele_calls = getAllele(allele_calls, first = FALSE)
  
  # Find the gene(s) of the allele calls; group duplicates (e.g. 1-69D*) as one
  gene_calls = getGene(allele_calls, first = FALSE, collapse = TRUE)
  gene_calls = gsub("([^H])D", "\\1", gene_calls)
  
  # If the sequences are assigned multiple genes, pick the more common gene
  # This should be very rare, since calls should be from unmutated sequences
  gctab = table(gene_calls)
  multigene = grep(",",names(gctab),value=TRUE)
  if (length(multigene) > 0){
    set_to = sapply(multigene, function(x) names(which.max(gctab[unlist(strsplit(x,","))])))
    for (i in 1:length(set_to)){
      if(!is.null(set_to[[i]])){
        gene_calls[which(gene_calls == names(set_to[i]))] = set_to[[i]]
      }
    }
  }
  
  # Remove genes that are too rare
  if(gene_cutoff < 1){ gene_cutoff = ceiling(length(gene_calls)*gene_cutoff) }
  rare_genes = names(which(table(gene_calls) < gene_cutoff))  
  df = data.frame(cbind(allele_calls,gene_calls),row.names=NULL,stringsAsFactors=F)
  exclude = which(gene_calls %in% rare_genes)
  if (length(exclude > 0)) { df = df[-exclude,] }
  
  #If after all that there still gene-ambiguous sequences, just keep the first gene
  stillmulti = grep(",",as.vector(df$gene_calls))
  if (length(stillmulti) > 0 ){
    df$gene_calls[stillmulti] = sapply(strsplit(as.vector(df$gene_calls[stillmulti]),","),"[",1)
  }
  # Make a table to store the resulting genotype
  gene = sortAlleles(as.character(unique(df$gene_calls)))
  gene = setdiff(gene, "")
  alleles = counts = rep("", length(gene))
  total = as.vector(table(as.character(df$gene_calls))[gene])
  genotype = cbind(gene, alleles,counts,total)
  
  # For each gene, find which alleles to include
  for (g in gene){
    
    ac = as.vector(df[df$gene_calls==g,"allele_calls"]) # find allele calls
    target = ceiling(fraction_to_explain*length(ac)) # how many we need to explain
    t_ac = table(ac) # table of allele calls
    potentials = unique(unlist(strsplit(names(t_ac),","))) # potential alleles
    # One allele? Easy!
    if (length(potentials) == 1 | length(t_ac) == 1){
      genotype[genotype[,"gene"]==g,"alleles"] = gsub("[^d\\*]*[d\\*]","",potentials )[1]
      genotype[genotype[,"gene"]==g,"counts"] = t_ac
    } else {
    # More alleles? Let's find the fewest that can explain the needed fraction
      
      # Make a table of whic alleles can explain which calls
      regexpotentials = paste(gsub("\\*","\\\\*", potentials),"$",sep="")
      regexpotentials = paste(regexpotentials,gsub("\\$",",",regexpotentials),sep="|")
      tmat = sapply(regexpotentials, function(x) grepl(x, names(t_ac),fixed=F))
      seqs_expl = as.data.frame(apply(tmat, 2, function(x) x*t_ac))
      colnames(seqs_expl) = potentials
      
      # Cycle through the table, including alleles to explain more sequences,
      # until we explain enough sequences
      included = counts = character(0)
      tot_expl = 0
      while(tot_expl < target){
        allele_tot = apply(seqs_expl, 2, sum)
        included = c(included, names(which.max(allele_tot)))
        counts = c(counts, max(allele_tot))
        tot_expl = max(allele_tot)  + tot_expl
        seqs_expl = seqs_expl[which(seqs_expl[,which.max(allele_tot)]==0),]
      }
      genotype[genotype[,"gene"]==g,"alleles"] = paste(gsub("[^d\\*]*[d\\*]","",included ),collapse=",")
      genotype[genotype[,"gene"]==g,"counts"] = paste(counts,collapse=",")
    }
    
  }
  
  return(as.data.frame(genotype, stringsAsFactors = FALSE))
}


# compareSepString --------------------------------------------------------
#' Compare two strings of separated values
#' 
#' \code{compareSepString} takes two strings, usually comma-separated, and
#' returns their intersection or difference in the form of a string using the
#' same separator.
#' 
#' @param    string1  a string of separated values, usually by a comma
#' @param    string2  a second string of separated values, usually by a comma
#' @param    value    what values to return. If "both" the intersection of the
#'                    values will be returned. "only1" and "only2" will return,
#'                    respectively, the values only in the first string or only
#'                    in the second string.
#' @param    sep      the separator that should be used to divide up the strings
#'                    before comparing the values they hold
#'                  
#' @return   A string of values representing the intersection or difference of
#'           of the input strings, separated in the same manner as the input
#'           
#' @examples
#' compareSepString("1,2,5,6", "1,5,7", value="both")
#' compareSepString("1,2,5,6", "1,5,7", value="only1")
#' compareSepString("1,2,5,6", "1,5,7", value="only2")
#' 
#' @export
compareSepString <- function(string1, string2, value="both", sep=",") {
  
  set1 = unlist(strsplit(string1, sep))
  set2 = unlist(strsplit(string2, sep))
  
  if (value == "both") {
    output = set1[set1 %in% set2]
  } else if (value == "only1") {
    output = set1[!(set1 %in% set2)]
  } else if (value == "only2") {
    output = set2[!(set2 %in% set1)]
  } else {
    output = set1[set1 %in% set2]
  }
  
  return(paste(output, collapse=sep))
  
}


# compareGenotypes --------------------------------------------------------
#' Compare two genotypes
#' 
#' \code{compareGenotypes} takes two genotypes, binds them together by gene,
#' and adds columns indicating the alleles only in the first, the alleles only
#' in the second, and the alleles shared between the two.
#' 
#' @param    genotype1  a genotype of the type returned by
#'                      \code{\link{inferGenotype}}
#' @param    genotype2  a genotype of the type returned by
#'                      \code{\link{inferGenotype}}
#'                  
#' @return   A data frame indicating which alleles are unique to each genotype
#'           or shared between then two
#'           
#' @seealso  \code{\link{inferGenotype}}          
#'           
#' @examples
#' # Load example data
#' data(pgp1_example)
#' 
#' # Determine a genotype
#' geno = geno2 = inferGenotype(pgp1_example[,"V_CALL"])
#' # Shuffle the gene names to make a different "genotype"
#' geno2$gene = sample(geno2$gene)
#' 
#' # Compare the two genotypes
#' compareGenotypes(geno, geno2)
#'                    
#' @export
compareGenotypes <- function(genotype1, genotype2){
  
  comparison = full_join(genotype1, genotype2, by  = "gene") %>%
    group_by(gene) %>%
    mutate(alleles.both = compareSepString(alleles.x, alleles.y)) %>%
    mutate(alleles.xonly = compareSepString(alleles.x, alleles.y, value="only1")) %>%
    mutate(alleles.yonly = compareSepString(alleles.x, alleles.y, value="only2"))
  
  return(comparison)
  
}


# genotypeFasta -----------------------------------------------------------
#' Return the nucleotide sequences of a genotype
#'
#' \code{genotypeFasta} converts a genotype table into a vector of nucleotide
#' sequences.
#' 
#' @param    genotype     a table of alleles denoting a genotype, as returned by
#'                        \code{\link{inferGenotype}}
#' @param    germline_db  a vector of named nucleotide germline sequences
#'                        matching the alleles detailed in \code{genotype} 
#' 
#' @return   A named vector of strings containing the germline nucleotide
#'           sequences of the alleles in the provided genotype
#' 
#' @seealso \code{\link{inferGenotype}}
#' 
#' @examples
#' # Load example data
#' data(germline_ighv)
#' data(pgp1_example)
#' 
#' # Infer and view a genotype from the sample
#' geno = inferGenotype(updateAlleleNames(pgp1_example[,"V_CALL"]))
#' geno
#' 
#' # Return the sequences that correspond to the genotype
#' genotypeFasta(geno, germline_ighv)
#' 
#' @export
genotypeFasta <- function(genotype, germline_db){
  g_names = names(germline_db)
  names(g_names) = gsub("D", "", names(germline_db))
  table_calls = mapply(paste, genotype$gene, strsplit(genotype$alleles, ","),
                       sep="*")
  seqs = germline_db[as.vector(g_names[unlist(table_calls)])]
  if(sum(is.na(seqs)) > 0){
    stop("The following genotype alleles were not found in germline_db: ",
         paste(unlist(table_calls)[which(is.na(seqs))], collapse = ", "))
  }
  return(seqs)
}


# writeFasta -----------------------------------------------------------
#' Write nucleotide sequences to a fasta file
#'
#' \code{writeFasta} write a vector of nucleotide sequences to a file
#' 
#' @param    named_sequences  a vector of nucleotide sequences
#' @param    file             a character string naming the file to write to
#' @param    char_per_line    how many characters should be printed per line
#' 
#' @return   Saves a fasta file containing the sequences of interest
#' 
#' @examples
#' \dontrun{
#' ## Load example IGHV germlines and write them to a fasta file
#' data(germline_ighv)
#' writeFasta(germline_ighv, file="germline_ighv.fasta")
#' }
#' 
#' @export
writeFasta <- function(named_sequences, file, char_per_line=60){
  # Ensure no sequences are too short or too long
  lens = sapply(named_sequences, nchar)
  if(min(lens) <= 0) { stop("All sequences must be of length > 0.") }
  if(max(lens) > 10^4) { stop("No sequences may have length > 10^4.") }
  # Break up sequences into lines of appropriate length
  split_seqs = lapply(named_sequences, substring,
                      seq(1, 10^4, char_per_line),
                      seq(char_per_line, 10^4, char_per_line))
  split_seqs = lapply(split_seqs, function(x) x[x!=""])
  seq_names = paste(">", names(named_sequences), sep="")
  to_print = as.vector(unlist(mapply(c, seq_names, split_seqs)))
  write(to_print, file=file)
}







