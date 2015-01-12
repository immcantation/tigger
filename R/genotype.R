
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




# only_unmutated - if true, empty allele calls (meaining the sequence had no
# allele that would represent an unmutated sequence) will be omitted.

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


# Infer genotypefrom a list of allele calls
inferGenotype <- function(allele_calls, # Calls of unique, unmutated sequences
                          fraction_to_explain = 7/8,
                          gene_cutoff = 1e-3 # Can be a no. of seqs or frac < 1
){
  
  
  # Find the gene(s) of the allele calls; group duplicates (e.g. 1-69D*) as one
  gene_calls = alakazam::getGene(allele_calls, first = FALSE, collapse = TRUE)
  gene_calls = gsub("D", "", gene_calls)
  
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






#

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
















