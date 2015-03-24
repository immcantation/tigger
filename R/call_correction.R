# reassignAlleles ---------------------------------------------------------
#' Correct allele calls based on a personalized genotype
#'
#' \code{reassignAlleles} uses a subject-specific genotype to correct
#' correct preliminary allele assignments of a set of sequences derived
#' from a single subject.
#' 
#' @param    v_calls       a vector of strings respresenting Ig allele calls for
#'                         the sequences in \code{v_sequences}, where multiple
#'                         calls are separated by a comma
#' @param    v_sequences   a vector of IMGT-gapped sample V sequences from a
#'                         single subject
#' @param    genotype_db   a vector of named nucleotide germline sequences
#'                         matching the calls detailed in \code{allele_calls}
#'                         and personalized to the subject
#' 
#' @return   a list equal in length to \code{v_calls}, best allele call from
#'           among the sequences listed in \code{genotype_db}
#' 
#' @examples
#' \dontrun{
#' ## Load example data and run all aspects of TIgGER (takes a few minutes)
#' data(pgp1_example)
#' data(germline_ighv)
#' results = runTigger(pgp1_example, germline_ighv)
#' 
#' ## Derive the subject-specific Ig sequences
#' novel_sequences = novelSummary(results, seqs_to_return = "in genotype")
#' germline_ighv = c(germline_ighv, novel_sequences)
#' genotype_db = genotypeFasta(sample_output$genotype, germline_ighv)
#' 
#' ## Extract the appropriate portions of example data
#' v_seqs = sapply(pgp1_example$SEQUENCE_GAP, substr, 1, 312)
#' 
#' ## Derive the vector of corrected calls
#' corrected_calls = reassignAlleles(pgp1_example$V_CALL, v_seqs, genotype_db)
#' }
#' 
#' @export
reassignAlleles <- function(v_calls, v_sequences, genotype_db){
  
  new_calls = rep("", length(v_calls))
  v_genes = getGene(v_calls, first = T)
  
  # Find which genotype genes are homozygous and assign those alleles first
  geno_genes = gsub("D", "", getGene(names(genotype_db)))
  names(geno_genes) = names(genotype_db)
  hetero_genes = unique(geno_genes[which(duplicated(geno_genes))])
  homo_genes = geno_genes[!(geno_genes %in% hetero_genes)]
  homo_alleles = names(homo_genes); names(homo_alleles) = homo_genes
  homo_calls_i = which(v_genes %in% homo_genes)
  new_calls[homo_calls_i] = homo_alleles[v_genes[homo_calls_i]]
  
  # Now realign the heterozygote sequences to each allele of that gene
  for (het_gene in hetero_genes){
    ind = which(v_genes %in% het_gene)
    if (length(ind) > 0){
      het_alleles = names(geno_genes[which(geno_genes == het_gene)])
      het_seqs = genotype_db[het_alleles]
      dists = lapply(het_seqs, function(x)
        sapply(getMutatedPositions(v_sequences[ind], x, match_instead=TRUE), length))
      dist_mat = matrix(unlist(dists), ncol = length(het_seqs))
      best_match = apply(dist_mat, 1, function(x) which(x == max(x)))
      best_alleles = sapply(best_match, function(x) het_alleles[x])   
      new_calls[ind] = sapply(best_alleles, paste, collapse=",")
    }
  }
  
  # Not realign the gene-not-in-genotype calls to every genotype allele
  hetero_calls_i = which(v_genes %in% hetero_genes)
  not_called = setdiff(1:length(v_genes), c(homo_calls_i, hetero_calls_i))
  dists = lapply(genotype_db, function(x)
    sapply(getMutatedPositions(v_sequences[not_called], x, match_instead=TRUE), length))
  dist_mat = matrix(unlist(dists), ncol = length(genotype_db))
  best_match = apply(dist_mat, 1, function(x) which(x == max(x)))
  best_alleles = sapply(best_match, function(x) names(genotype_db[x])) 
  new_calls[not_called] = sapply(best_alleles, paste, collapse=",")
  
  return(new_calls)
}


