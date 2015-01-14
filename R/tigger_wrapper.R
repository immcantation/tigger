

tigger <- function(sample_db, germline_db,
                   find_novel = TRUE, find_genotype = TRUE, correct_calls = TRUE,
                   allele_min = 1e-4, y_intercept = 1/8, nt_min=1, nt_max = 312,
                   mut_min=1, mut_max=10, j_max = 0.1, min_seqs = 50, min_frac = 3/4,
                   fraction_to_explain = 7/8,
                   verbose=FALSE){
  
  result = list(novel="", genotype="", new_calls="")
  
  # EXTRACT USEFUL PORTIONS OF DB FILES

  v_calls = alakazam::getAllele(sample_db[,v_call_col], first = FALSE)
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
    if(verbose){ cat("Finding novel alleles...") }
    allele_groups = assignAlleleGroups(v_calls, allele_min)
    j_genes = alakazam::getGene(sample_db[,j_call_col], first = FALSE)
    junc_lengths = sample_db[,junc_length_col]
    novel = detectNovelV(v_sequences, j_genes, junc_lengths, allele_groups,
                         germline_db,  y_intercept, nt_min, nt_max,
                         mut_min, mut_max, j_max, min_seqs, min_frac, 
                         verbose)
    # Extract the nucleotide sequence portion
    fasta = unlist(unique(sapply(novel, "[", 1)))
    # In case we found the same allele multiple ways, ditch the duplicate
    fasta = fasta[which(!duplicated(fasta))]
    # Add the novel alleles to the germline_db
    germline_db = c(germline_db, fasta)
    result[["novel"]] = novel
    if(verbose){ cat("done.\n") }
  }
  
  
  # DETERMINE GENOTYPE
  
  if (find_genotype){
    if(verbose){ cat("Finding genotype...") }
    v_calls2 = v_calls
    # Paste novel alleles (if any) to all allele calls before determining dists
    if (find_novel){
      genes_novel = alakazam::getGene(names(novel))
      genes_groups = alakazam::getGene(names(allele_groups))
      for(i in 1:length(novel)){
        matching_groups = allele_groups[genes_groups %in% genes_novel[i]]
        indicies = unique(unlist(matching_groups))
        v_calls2[indicies] = sapply(v_calls2[indicies], paste, names(novel)[i], sep=",")
      }
    }
    
    mut_counts = getMutCount(v_sequences, v_calls2, germline_db)
    unmutated_calls = findUnmutatedCalls(v_calls2, mut_counts)
    genotype = inferGenotype(unmutated_calls, fraction_to_explain)
    genotype_db = genotypeFasta(genotype, germline_db)
    result[["genotype"]] = genotype
    if(verbose){ cat("done.\n") }
  }
  
  
  # CORRECT ALLELE CALLS
  
  if (correct_calls) {
    if(verbose){ cat("Correcting allele calls...") }
    if(!find_genotype){ genotype_db = germline_db }
    new_calls = reassignAlleles(v_calls, v_sequences, genotype_db)
    result[["new_calls"]] = new_calls
    if(verbose){ cat("done.\n") }
  }
  
  return(result)
  
}

















