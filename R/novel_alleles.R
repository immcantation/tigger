

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
#' @return   a sorted vector of strings respresenting Ig allele names
#' 
#' @examples
#' # Create a list of allele names
#' alleles = c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
#' "IGHV2-5*01","IGHV1-NL1*01","IGHV1-2*02", "IGHV1-69*02")
#' 
#' # Sort the alleles
#' sortAlleles(alleles)
#' 
#' @export
sortAlleles <- function(allele_calls) {  
  
  allele_calls = sort(allele_calls) # This will help align the Ds later.
  segment = gsub(".*(IG[HKL][VDJ]).*", "\\1", allele_calls)
  family = as.numeric(gsub(".*IG[HKL][VDJ]([0-9]*)[^0-9].*", "\\1", allele_calls))
  gene = gsub("[^-]*-([^-\\*D]*)[-\\*D].*", "\\1", allele_calls)
  gene = as.numeric(gsub("NL", "99", gene))
  gene2 = gsub("-", "", gsub("^[^-]*|[^-]*$", "", allele_calls))
  gene2 = as.numeric(gsub("D","",gene2))
  gene2[is.na(gene2)] = 0
  allele = as.numeric(gsub(".*\\*","",allele_calls))
  sorted_calls = allele_calls[order(segment, family, gene + gene2, allele)]
  
  return(sorted_calls)
}


# assignAlleleGroups ------------------------------------------------------

#' Find indicies of allele calls
#'
#' \code{assignAlleleGroups} takes a vector of allele calls (a call may consist
#' of multiple comma-separated alleles). It returns a list, named for the unique
#' alleles in the input vector, of indicies of the calls containing the alleles. 
#' 
#' @param    allele_calls  a vector of strings respresenting Ig allele calls
#' @param    allele_min  a number < 1 representing the minimum fraction of
#'           sequences--or a number >= 1 representing the minimum count for
#'           sequences--required for an allele to not be excluded.
#' @param    binomial_cutoff a logical indicating if an \code{allele_min} < 1
#'           should be applied in a binomial manner.
#' @param    alpha  the alpha cutoff used if \code{binomial_cutoff = TRUE}
#' @return   a list of indicies of calls that contain each unique input allele
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
  allele_counts = ldply(alleles_i, length)
  rare_i = allele_counts[,2] < cutoff
  
  # Print out some statistics
  if (sum(rare_i) > 0){
    cat(sum(rare_i), " alleles excluded due to rarity, representing ",
        sum(allele_counts[rare_i,2]), " sequences (",
        100*round(sum(allele_counts[rare_i,2])/sum(allele_counts[,2]),4),
        "% of input).\n", sep="")
  }
  cat(sum(!rare_i), "alleles detected and retained.\n")
  
  return(alleles_i[!rare_i])
}


# getMutatedPositions -----------------------------------------------------


#' Find the location of mutations in a sequence
#'
#' \code{getMutatedPositions} takes two vectors of aligned sequences and
#' compares pairs of sequences. It returns a list of the nucleotide positions of
#' any differences.
#' 
#' @param    samples  a vector of strings respresenting aligned sequences
#' @param    germlines  a vector of strings respresenting aligned sequences to
#'           which \code{samples} will be compared. If only one string is
#'           submitted, it will be used for all \code{samples}.
#' @param    ignored_regex a regular expression indicating what characters
#'           should be ignored (such as gaps and N nucleotides).
#' @param    match_instead  if \code{TRUE}, the function returns the positions
#'           that are the same instead of those that are different.
#' @return   a list of the nucleotide positions of any differences between the
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
#' \code{summarizeMutations} takes two vectors of aligned sequences and
#' compares pairs of sequences. It returns a matrix detailing which positions
#' are mutated as a function of sequence-wide mutation count.
#' 
#' @param    mut_list  a list of the nucleotide positions of any differences
#'           between the two vectors of sequences, as generated by
#'           \code{getMutatedPositions}.
#' @param    match_list  mut_list  a list of the nucleotide positions of any
#'           similarities between the two vectors of sequences, as generated by
#'           \code{getMutatedPositions} where \code{match_instead = TRUE}.
#' @return   a list containing two matricies. The first details counts of
#'           sequences mutated at given positions (rows) mutated as a
#'           function of sequence-wide mutation count (columns). The second is
#'           details how many usable nucleotides (i.e., not gaps or Ns) were
#'           found for each combination of position and sequence-wide mutation
#'           count.
#' 
#' @examples
#' # TO BE ADDED
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
#' If there is a problem with the number of sequences, etc., 
#'           a message will be printed and \code{NULL} will be returned.

trimMutMatrix <- function(mut_summary, mut_min=1, mut_max=10,
                          nt_min = 1, nt_max = 312, min_seqs=50){
  
  allele = names(germline)
  
  # Ensure the matrix covers the desired mutational range
  if(sum(!(mut_min:mut_max %in% colnames(mut_summary[[1]]))) > 0) {
    cat("Insufficient data for ", allele, " within mutation range ",
        mut_min, " to ", mut_max,"; skipped.\n", sep="")
    return(NULL)
  }
  
  # Ensure the matrix covers the desired nucleotide range
  if(sum(!(nt_min:nt_max %in% rownames(mut_summary[[1]]))) > 0) {
    cat("Insufficient data for ", allele, " within nucleotide range ",
        nt_min, " to ", nt_max,"; skipped.\n", sep="")
    return(NULL)
  }
  
  # Trim the matrix to the desired range of mutations
  mut_range = as.character(mut_min:mut_max)
  nt_range = as.character(nt_min:nt_max)
  mut_trim = list(mut_summary[[1]][nt_range,mut_range],
                  mut_summary[[2]][nt_range,mut_range])
  
  # Ensure there are actually sequences to work with in the range
  if(sum(count_vs_muts) < min_seqs){
    cat("Insufficient data for ", allele, " within mutation range ",
        mut_min, " to ", mut_max,"; skipped.\n", sep="")
    return(NULL)
  }
  
  # Find any outliers in the sequence counts.
  count_vs_muts = apply(mut_trim[[2]], 2, max)
  outs = boxplot(count_vs_muts,plot=F)$out
  big_outs = outs[which(outs > mean(count_vs_muts))]
  start_at = as.numeric(names(which(count_vs_muts == big_outs)))
  # If you find outliers, trim the matrix again
  if(length(start_at) > 0){
    mut_min = min(5, start_at[1])
    mut_range = as.character(mut_min:mut_max)
    mut_trim = list(mut_summary[[1]][,mut_range],
                    mut_summary[[2]][,mut_range])
  }
  
  mut_fracs = mut_trim[[1]]/mut_trim[[2]]
  
  return(mut_fracs)
  
}


# findIntercepts ----------------------------------------------------------

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
#' germline is a nt sequence and should have the allele name as its name

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

findNovelAlleles  <- function(samples, germline, j_genes, junc_lengths,
                              mut_min=1, mut_max=10, j_max = 0.1){
  # Find the positions of differences and similarities between sequences
  mut_list = getMutatedPositions(samples, germline)
  mut_counts = sapply(mut_list, length)
  match_list = getMutatedPositions(samples, germline, match_instead = TRUE)
  mut_summary = summarizeMutations(mut_list, match_list, )
  mut_matrix = trimMutMatrix(mut_summary, mut_min, mut_max)
  intercepts = findIntercepts(mut_matrix)
  polymorphs = as.numeric(names(intercepts))
  allele_summary = list()
  if (length(polymorphs) > 0) {
    nuc_usage = lapply(polymorphs, findNucletoideUsage, samples, germline,
                       mut_counts, mut_min, mut_max)
    names(nuc_usage) = polymorphs
    snp_nucs = sapply(nuc_usage, function(x) rownames(x)[4])
    putative = createGermlines(germline, polymorphs, snp_nucs)
    put_mut_locs = lapply(putative, function(x) getMutatedPositions(samples, x))
    perfect_matches = lapply(put_mut_locs, function(x) which(sapply(x,length) == 0))
    j_junc_tables = lapply(perfect_matches, function(x) table(junc_lengths[x], j_genes[x]))
    pass_test = which(sapply(j_junc_tables, function(x) max(x/sum(x))) < j_max)
    for(i in 1:length(pass_test)){
      germ_name = names(pass_test[i])
      positions = gsub("[A-Z]","",strsplit(germ_name, "_")[[1]][-1]) 
      allele_summary[[germ_name]] = list(putative[pass_test[i]],
                                         intercepts[positions],
                                         mut_matrix,  nuc_usage[positions],
                                         j_junc_tables[[pass_test[i]]])
    }
  }
  return(allele_summary)
}  


# old tigger stuff --------------------------------------------------------




# script ------------------------------------------------------------------


## START SCRIPT

# Constants

v_call_col = "V_CALL"
v_start_col = "V_GERM_START"
v_length_col = "V_GERM_LENGTH"
j_call_col = "J_CALL"
junc_length_col = "JUNCTION_GAP_LENGTH"

# Parameters you don't change
mut_min = 1
mut_max = 10


# Your part

# Data 

load("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/Bolen Twins/twinStudy_ClipTab_withFRdata_goodSeqs_deident.Rd")
Read_fasta <- function(fasta_file, quiet=F){
  # Reads a fasta file, returning named uppercase strings
  if(!quiet){ cat("Reading \"", fasta_file, "\"...", sep="") }
  seqs <- read.fasta(fasta_file, strip.desc=T, as.string=T)
  seq_names <- sapply(seqs, attr, "Annot")
  seqs = toupper(seqs)
  names(seqs) = seq_names
  if(!quiet){ cat("done.\n") }
  return(seqs)
}
germlines = Read_fasta("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/IMGT/IMGT Ig Sequences 2014-09-29.fasta")
dat = subset(combSeqData.noID, donor == "Donor_06" | donor == "Donor_09")

# Main

allele_groups = assignAlleleGroups(dat[,v_call_col])
v_seqs = sapply(dat$SEQUENCE_GAP, substr, 1,
                as.numeric(dat[,v_start_col]) + as.numeric(dat[,v_length_col]))
j_genes = alakazam::getGene(dat[,j_call_col], first = F)
junc_lengths = dat[,junc_length_col]


novel=list()

which(sapply(allele_groups, length) > 

for each allele with in size of something

  allele_gp = allele_groups[[2]]
  samples = v_seqs[]
  germline = 
  novel[] = findNovelAlleles(samples, germline, j_genes[allele_groups[[2]]],
                           junc_lengths[allele_groups[[2]]],
                           mut_min=1, mut_max=10, j_max = 0.1)
  



## END SCRIPT



  # Find hypothetical polymorphisms, on a per-allele basis
  

  muts_list = list(); pcts_list = list(); snps_list = list()
  for (gnam in names(igroups_list)){
    grp = igroups_list[[gnam]]
    v_samples = samples[grp]
    v_germline <- repseqs[gnam]
    muts_list[[gnam]] = Find_mutated_positions(v_samples, v_germline)
    pcts_list[[gnam]] = Find_mutation_by_position(muts_list[[gnam]], mmax, allele_cutoff)
    snps = Find_polymorphisms(v_samples, v_germline, muts_list[[gnam]],
                              pcts_list[[gnam]], y_intercept = y_i,
                              inter_alpha=i_a, mmin=mmin, mmax=mmax,
                              allele_cutoff=allele_cutoff)
    if(length(snps)>0){ snps_list[[gnam]] = snps }
  }
  # sum(sapply(muts_list, function(x) sum(table(unlist(x))/length(x) > 0.3)))3)))
  if (length(snps_list)>0){
    snps_table = as.table(Reduce(rbind, snps_list))
    if(!clean){
      write.table(snps_table, paste0(run, "/putative_SNPs_list.tab"),
                  sep="\t", row.names=F, quote=F)
    }
    # Create predicted germlines, adding multipy-SNP'd germlines if >1 SNP found in 1 allele
    seqs = c(unlist(Create_predicted_germlines(snps_table, repseqs)),
             unlist(Create_multiple_snp_germlines(snps_table, repseqs)))
    # Check the J/junction length distribtions
    jtab_list = list()
    for (s in names(seqs)){
      jtab_list[[s]] = Calc_jtab(seqs[s], igroups_list, samples,
                                 clip_tab[,j_call_col],
                                 clip_tab[,junc_col])
    }
    # Require 50+ perfect-matching sequences
    # over50 = which(sapply(jtab_list, function(x) sum(x)) > 50)
    # Require most common J-junction combinations to be < 10% among perfect sequences
    under15 = which(sapply(jtab_list, function(x) max(x/sum(x))) < jpct)
    #   pass = names(jtab_list[intersect(over50,under10)])
    pass = names(jtab_list[under15])
    seqs = pseqs = seqs[pass]
    # Check for/remove SNPd seqs that result in the same sequence as each other
    if (length(seqs) > 2){ seqs = Remove_duplicated_germlines(seqs, quiet=quiet) }
    # Merge new repseqs with old repseqs
    repseqs = c(repseqs, seqs)[ sort(unique(c(names(repseqs), names(seqs)))) ]
    # Write *unique* putative new alleles to fatsa
    if(!clean){
      write.fasta(lapply(repseqs, s2c), names(repseqs),
                  file.out=paste0(run, "/rep_+_putative.fasta"))
    }
    
  }
  if (!quiet){
    cat("done.\n", length(seqs), "novel allele(s) found:\n",
        paste(names(seqs),collapse="\n "),"\n")
  }
  



