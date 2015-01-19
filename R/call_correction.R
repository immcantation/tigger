

# reassignAlleles ---------------------------------------------------------

# genotype_db should be personalized to the individual
# the new v_calls will be limited to the names of the sequences in genotype_db
#' @export
reassignAlleles <- function(v_calls, v_sequences, genotype_db){
  
  new_calls = rep("", length(v_calls))
  v_genes = alakazam::getGene(v_calls, first = T)
  
  # Find which genotype genes are homozygous and assign those alleles first
  geno_genes = alakazam::getGene(names(genotype_db))
  names(geno_genes) = names(genotype_db)
  hetero_genes = geno_genes[which(duplicated(geno_genes))]
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
                               germline_list[ccm], samples[ccm])
    mut_count_list[ccm] = lapply(mut_pos_list[ccm],
                                 function(x) lapply(x,length))
  }
  
  return(mut_count_list)
}





