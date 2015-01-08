#
#
#
#
#


# LOAD DEPENDENCIES

library(alakazam)
library(plyr)
source("C:/Users/Daniel Gadala-Maria/repos/tiggerpackage/R/novel_alleles.R")

# Parameters you don't change often
mut_min = 1
mut_max = 10
min_allele = 1e-4 # for allele detection and genotyping


# WHAT DO YOU WANT TO DO?

find_novel = TRUE
plot_novel = TRUE # will be effictively false if find_novel = FALSE
find_genotype = TRUE
correct_calls = TRUE
v_gap_length = c("V_GERM_LENGTH","V_GAP_LENGTH")[1]

# LOAD YOUR DATA

# clip_tab
load("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/HIVnew.RData")
ddply(hiv_annot, "SUBJECT", nrow)
#dat = subset(dat.ib, TIME %in% c("-1h", "-2d", "-8d") & FUNCTIONAL == "T")
subjects = unique(hiv_annot$SUBJECT); i = 1
dat = subset(hiv_annot, SUBJECT %in% subjects[i] & FUNCTIONAL == "TRUE")
dat = dat[which(!duplicated(dat$SEQUENCE_GAP)),] # REMOVE DUPLICATES
#rm(dat.ib)
germline_db_file = "C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/IMGT/IMGT Variable 2014-12-22.fasta"
germline_db = readGermlineDb(germline_db_file, strip_down_name = TRUE)


# EXTRACT SOME USEFUL PORTIONS

v_calls = alakazam::getAllele(dat[,v_call_col], first = FALSE)
v_calls = updateAlleleNames(v_calls)
# New version of Clip tabs are different, so this depends on the columns names
if (v_gap_length == "V_GAP_LENGTH"){
  v_sequences = sapply(dat$SEQUENCE_GAP, substr, 1, dat$V_GAP_LENGTH)
} else {
  v_sequences = sapply(dat$SEQUENCE_GAP, substr, 1, dat$V_GERM_START + dat$V_GERM_LENGTH)
}


# FIND NOVEL ALLELES

if (find_novel){
  allele_groups = assignAlleleGroups(v_calls, min_allele)
  j_genes = alakazam::getGene(dat[,j_call_col], first = F)
  junc_lengths = dat[,junc_length_col]
  novel = detectNovelV(v_sequences, j_genes, junc_lengths, allele_groups,
                       germline_db,  y_intercept =1/8, nt_min=1, nt_max = 312,
                       mut_min=1, mut_max=10, j_max = 0.15, verbose=FALSE)
  fasta = unlist(unique(sapply(novel, "[", 1)))
  # In case we found the same allele multiple ways, ditch the duplicate
  fasta = fasta[which(!duplicated(fasta))]
  # Add the novel alleles to the germline_db
  germline_db = c(germline_db, fasta)
  
  # Plot the supporting data for the novel alleles
  if (plot_novel){
    plotNovelLines(novel)
    plotNovelBars(novel)
    plotJunctionBars(novel) 
  }
}


# DETERMINE GENOTYPE

if (find_genotype){
  v_calls2 = v_calls
  if (find_novel){
    # Paste novel alleles to all allele calls before determining distance, if we
    # searched for novel alleles
    novel_names = paste(names(novel), collapse=",")
    v_calls2 = sapply(v_calls2, paste, novel_names, sep=",")
  }
  mut_counts = getMutCount(v_sequences, v_calls2, germline_db)
  unmutated_calls = findUnmutatedCalls(v_calls2, mut_counts)
  genotype = inferGenotype(unmutated_calls)
  genotype_db = genotypeFasta(genotype, germline_db)
  
}


# CORRECT ALLELE CALLS
if (correct_calls) {
  if(!find_genotype){ genotype_db = germline_db }
  new_calls = reassignAlleles(v_calls, v_sequences, genotype_db)
}














