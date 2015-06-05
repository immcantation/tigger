library(shm)
library(alakazam)
library(tigger)
  
verbose = TRUE

# DATA

db_file = "C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/fv_from_sonia/SG3M_Assembled_atleast-2_internalC-pass_HL_IGLV_multiFxDBs-F_jxngap-F.rdata"
load(db_file)
subj = toupper(gsub(".*/|\\.[^\\.]+","",db_file))
dat = dat.ib; rm(dat.ib); gc()
dat = dat.fv.l

# Clean up the dataset
dat = filter(dat, TIME %in% c("-1h","-8d","-2d"))
dat = dat %>%
  modifyChangeoDb(seq_imgt_col = "SEQUENCE_GAP", junc_len_col = "JUNCTION_GAP_LENGTH")


data(germline_ighv)
GERMLINE_CALL = names(germline_ighv)
TRUE_GERMLINE = as.character(germline_ighv)
germ_df = data.frame(GERMLINE_CALL, TRUE_GERMLINE, stringsAsFactors = FALSE)

#' @export
insertNVMutations <- function(sequence, mutation_count){
  mut_to = list("A" = c("C","G","T"),"G" = c("C","A","T"),
                "C" = c("A","G","T"),"T" = c("C","G","A"))
  sequence = toupper(sequence)
  if (mutation_count >0) {
    nt_positions = gregexpr("[ACGT]", sequence)[[1]]
    positions_to_mut = sample(nt_positions[nt_positions <= 312], mutation_count)
    current_nt = sapply(positions_to_mut, function(x) substring(sequence, x, x))
    new_nt = sapply(current_nt, function(x) sample(mut_to[[x]], 1))
    new_sequence = insertPolymorphisms(sequence, positions_to_mut, new_nt)
  } else {
    new_sequence = sequence
  }
  return(new_sequence)
}

tig = runTigger(dat, germline_ighv)
novel_sequences = novelSummary(tig, seqs_to_return = "in genotype")
genotype_sequences = genotypeFasta(tig$genotype, c(germline_ighv, novel_sequences))



for(nmuts in 0:10){
  cat("Starting on nmuts =", nmuts,"\n")
  germlines = sapply(genotype_sequences, insertNVMutations, nmuts)
  n <- dat %>%
    filter(ISOTYPE %in% c("IGHG", "IGHA")) %>%
    findNovelAlleles2(germlines, nproc = 4)
  if(nmuts == 0) {
    novel2 = n
  } else {
    novel2 = bind_rows(novel2,n)
  }
}

#load("ib distance 0 ro 50.rdata")
#novel = bind_rows(novel, novel2)

#save(novel, file="ib switched distance 0 to 10.rdata")


compateSets = function(s1, s2, is.in=TRUE){
  s1 = unlist(strsplit(as.character(s1), ","))
  s2 = unlist(strsplit(as.character(s2), ","))
  s1 = setdiff(s1, "integer(0)")
  if(is.in){
    return(paste(s1[s1 %in% s2], collapse=","))
  } else {
    return(paste(setdiff(s1, s2), collapse=","))
  }
}

sorter = function(a, b){
  if(a & b) { return("Perfect") }
  if(a) { return("Perfect and\nPolymorphic") }
  return("No Call")
}

result_df = novel %>%
  #  bind_rows(novel2) %>%
  # bind_rows(novel3) %>%
  left_join( germ_df, by = "GERMLINE_CALL") %>%
  group_by(1:n()) %>%
  mutate(GERMLINE_CALLED = nchar(NOVEL_IMGT) > 2 ) %>%
  mutate(N_DIFFERENCES = sapply(getMutatedPositions(GERMLINE_IMGT, TRUE_GERMLINE),length)) %>%
#  mutate(MISSED_POSITIONS = compateSets(getMutatedPositions(NOVEL_IMGT, TRUE_GERMLINE), positions) ) %>%
#  mutate(ADDED_POSITIONS = compateSets(getMutatedPositions(NOVEL_IMGT, TRUE_GERMLINE), positions, is.in=FALSE) ) %>%
  mutate(PERFECT = sapply(getMutatedPositions(NOVEL_IMGT, TRUE_GERMLINE),length) == 0 ) %>%
  mutate(PERFECT = PERFECT & GERMLINE_CALLED) %>%
  mutate(WRONG_POS = paste(getMutatedPositions(NOVEL_IMGT, TRUE_GERMLINE),paste=",")) %>%
  mutate(CALL_QUALITY = sorter(GERMLINE_CALLED, PERFECT)) %>%
  filter(grepl("_[0-9]", POLYMORPHISM_CALL) == FALSE)
#View(result_df)


rdf = result_df %>%
  filter(GERMLINE_CALL %in% names(genotype_sequences)) %>%
  filter(GERMLINE_CALL_COUNT >= 500) %>%
  filter(N_DIFFERENCES > 0)
rdf$WRONG_POS[grep("integer", rdf$WRONG_POS)] = ""
rdf$CALL_QUALITY[sapply(rdf$WRONG_POS, nchar)>0] = "False\nPositive"
rdf$CALL_QUALITY[grep("112,222,286", rdf$WRONG_POS)] = "Perfect"
rdf$CALL_QUALITY[grep("Perfect", rdf$CALL_QUALITY)] = "True\nPositive"
#rdf$N_DIFFERENCES = factor(rdf$N_DIFFERENCES, levels = 50:0)
ggplot(rdf, aes(factor(N_DIFFERENCES), fill = CALL_QUALITY)) +
  geom_bar(binwidth=1) +
  scale_x_discrete(breaks=c(1,seq(5, 50, 5)), labels = c(1,seq(5, 50, 5))) +
  theme_bw() + scale_fill_discrete(name="Call Quality") +
  xlab("Distance From IMGT Germline") +
  ylab("Count") + ylim(0,round(length(genotype_sequences),digits = -1))

ggplot(summarise(group_by(rdf,GERMLINE_CALL_COUNT),
                 Correct = sum(CALL_QUALITY=="Perfect")/n()),
       aes(GERMLINE_CALL_COUNT, Correct, color = Correct)) +
  geom_point() + ylim(0,1) +
  xlab("Gene Count") + ylab("Approximate Frequency of Detection")

ggplot(result_df , aes(CALL_QUALITY, GERMLINE_CALL_COUNT)) +
  geom_boxplot() + theme_bw() +
  xlab("Germline Call Quality") +
  ylab("Sequence Count")

ggplot(rdf, aes(GERMLINE_CALL_COUNT, MUT_MIN, color = CALL_QUALITY)) +
  geom_jitter()


View(novel[!is.na(novel$POLYMORPHISM_CALL),])

# Read in V and J germlines
v_germs = readGermlineDb("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/IMGT/IMGT Variable 2014-12-22.fasta")
j_germs = readGermlineDb("C:/Users/Daniel Gadala-Maria/Documents/Kleinstein/Datasets/IMGT/IMGT Joining 2015-02-06.fasta")


x = getPopularMutationCount(dat, v_germs, seq_imgt_col = "SEQUENCE_IMGT")


tigger_out = runTigger(dat, v_germs, j_max = .15)
v_germs = c(v_germs, novelSummary(tigger_out))

# MAKE SOME "CLONES"

dat2 = dat %>%
  mutate(CLONE = paste(V_GENE, J_GENE, JUNCTION, sep=",")) %>%
  mutate(CLONE = dense_rank(CLONE)) %>%
  group_by(CLONE) %>%
  mutate(CLONE_SIZE = n()) %>%
  group_by(V_SEQUENCE_IMGT) %>%
  mutate(V_COUNT = n()) %>%
  mutate(MEAN_CLONE_SIZE = mean(CLONE_SIZE)) %>%
  group_by(V_GENE) %>%
  mutate(MAX_V_COUNT = max(V_COUNT)) %>%
  mutate(V_FRAC_OF_MAX = V_COUNT/MAX_V_COUNT) %>%
  ungroup() %>%
  arrange(V_GENE, desc(V_COUNT))




test1 = function(x, y, ...){
  
  args = as.list(match.call())
  
  a_plot = names(formals(plot))
  a_axis = args[names(args) %in% names(formals(axis))]
  
  # print(args)
  
  do.call(plot, args[names(args) %in% a_plot])
  
  a_axis["side"] = 1
  do.call(axis, a_axis)
  a_axis["side"] = 3
  do.call(axis, a_axis)
  
  return(a_axis)
}


x = test1(1:10, 11:20, side = 3)

















suppressPackageStartupMessages(library(shm))
suppressPackageStartupMessages(library(tigger))


positionMutations <- function(db_subset, germline, pos_range){
  # Duplicate each sequence for all the positions to be analyzed
  pos_db = pos_range %>%
    length() %>%
    rep("db_subset", .) %>%
    paste(collapse=",") %>%
    paste("bind_rows(",., ")") %>%
    parse(text=.) %>%
    eval()
  pos_db$POSITION = c(sapply(pos_range, rep, nrow(db_subset)))
  # Find which positions are mutated
  pos_db = pos_db %>%
    mutate(NT = substring(SEQUENCE_IMGT, POSITION, POSITION)) %>%
    mutate(GERM_NT = substring(germline, POSITION, POSITION)) %>%
    mutate(MUTATED = (NT != GERM_NT & NT != "N" & NT != "-" & NT != "")) %>%
    mutate(OBSERVED = (NT != "-" & NT != ""))
  return(pos_db)
}

mutationRangeSubset <- function(db_subset, germline, mut_range, pos_range){
  pads = paste(rep("-", min(pos_range)-1), collapse="")
  db_subset$MUT_COUNT = db_subset$SEQUENCE_IMGT %>%
    substring(min(pos_range), max(pos_range)) %>%
    paste(pads, ., sep="") %>%
    getMutatedPositions(germline) %>%
    sapply(length)
  db_subset = db_subset %>%
    filter(MUT_COUNT %in% mut_range)
  return(db_subset)
}


findNovelAlleles2  <- function(clip_db, germlines,
                               germline_min = 200,
                               mut_min = 1, mut_range = 10,
                               auto_mutrange = TRUE,
                               nt_min = 1, nt_max = 312,
                               min_frac = 0.75,
                               y_intercept = 1/8, alpha = 0.05,
                               j_max = 0.10,
                               min_seqs = 50, nproc = 8){
  # Make seqs into caps, replace dots with dashes, make all non ACGT- into N
  germlines = germlines %>%
    toupper() %>%
    gsub(".", "-", . , fixed = TRUE) %>%
    gsub("[^ACGT-]", "N", .)
  db_strip = clip_db %>%
    select(SEQUENCE_IMGT, V_CALL, J_CALL, JUNCTION_LENGTH)
  db_strip$SEQUENCE_IMGT = db_strip$SEQUENCE_IMGT %>%
    toupper() %>%
    gsub(".", "-", . , fixed = TRUE) %>%
    gsub("[^ACGT-]", "N", .)
  
  # Find which rows' calls contain which germline alleles
  cutoff = ifelse(germline_min < 1, round(nrow(db_strip)*germline_min), germline_min)
  
  allele_groups = sapply(names(germlines), grep, db_strip$V_CALL, fixed=T)
  allele_groups = allele_groups[sapply(allele_groups, length) >= cutoff]
  # allele_groups = allele_groups[sortAlleles(names(allele_groups))]
  
  # Private function to find lower range of y intercept confidence interval
  findLowerY = function(x, y, mut_min, alpha){
    y = y+1-mut_min
    lowerY = suppressWarnings(confint(lm(x ~ y),level=1-2*alpha)[[1]])
    return(lowerY)
  }
  
  # Prepare for parallel processing!
  nproc <- min(nproc, getnproc()-1)
  # If user wants to paralellize this function and specifies nproc > 1, then
  # initialize and register slave R processes/clusters & 
  # export all nesseary environment variables, functions and packages.
  if( nproc==1 ) {
    # If needed to run on a single core/cpu then, regsiter DoSEQ 
    # (needed for 'foreach' in non-parallel mode)
    registerDoSEQ()
  } else {
    if(nproc != 0) { cluster <- makeCluster(nproc, type = "SOCK") }
    clusterExport( cluster, list( "germline_min",
                                  "mut_min", "mut_range",
                                  "auto_mutrange",
                                  "nt_min", "nt_max",
                                  "min_frac",
                                  "y_intercept", "alpha",
                                  "j_max",
                                  "min_seqs", "germlines", "db_strip",
                                  "allele_groups", "findLowerY",
                                  "getPopularMutationCount"), 
                   envir=environment() )
    clusterEvalQ(cluster, library(tigger))
    registerDoSNOW(cluster)
  }
  
  df_out <-
    foreach(a=icount(length(allele_groups)), .combine = rbind ) %dopar% {
      
      
      allele_name = names(allele_groups)[a]
      
      # Subset of data being analyzed
      germline = germlines[allele_name]
      indicies = allele_groups[[allele_name]]
      db_subset = slice(db_strip, indicies)
      
      # If mutrange is auto, find most popular mutation count and start from there
      gpm = db_subset %>%
        mutate(V_CALL = allele_name) %>% # For the use of the germline of interest
        getPopularMutationCount(germline,
                                gene_min=0, seq_min=min_seqs,
                                seq_p_of_max = 1/8, full_return = T)
      
      # Determine the mutation range(s) to scan
      if(auto_mutrange & sum(gpm$MUTATION_COUNT > 0) > 0 ){
        mut_mins = gpm$MUTATION_COUNT[gpm$MUTATION_COUNT > 0]
      } else {
        mut_mins = 1
      }

      for (mut_min in mut_mins) {
        mut_max = mut_min + (mut_range - 1)
        # Create the run's return object
        df_run = data.frame(GERMLINE_CALL = names(germline),
                            NOTE = "",
                            POLYMORPHISM_CALL = NA,
                            NOVEL_IMGT = NA,
                            PERFECT_MATCH_COUNT = NA,
                            GERMLINE_CALL_COUNT = length(indicies),
                            GERMLINE_IMGT = as.character(germline),
                            NT_MIN = nt_min,
                            NT_MAX = nt_max,
                            MUT_MIN = mut_min, 
                            MUT_MAX = mut_max,
                            Y_INTERCEPT = y_intercept,
                            ALPHA = alpha,
                            MIN_SEQS = min_seqs,
                            J_MAX = j_max,
                            stringsAsFactors = FALSE)
        
        # If gpm is zero no sequence is frequent enough to pass the J test so give up
        if(length(gpm) > 0) {
          
          # Add a mutation count column and filter out sequences not in our mut range
          db_subset = mutationRangeSubset(db_subset, germline, mut_range, pos_range)
          
          if(nrow(db_subset) > germline_min ){
            
            pos_range = nt_min:nt_max
            
            # Duplicate each sequence for all the positions to be analyzed
            # and find which positions are mutated
            pos_db = positionMutations(db_subset, germline, pos_range)
            
            
            # Find positional mut freq vs seq mut count
            pos_muts = pos_db %>%
              group_by(POSITION) %>%
              mutate(PASS = mean(OBSERVED) >= min_frac) %>%
              group_by(MUT_COUNT, POSITION) %>%
              summarise(POS_MUT_RATE = mean(MUTATED)*unique(PASS) ) %>% 
              ungroup()   
            
            # Calculate y intercepts
            pass_y = pos_muts %>%
              group_by(POSITION) %>%
              summarise(Y_INT_MIN = findLowerY(POS_MUT_RATE, MUT_COUNT, mut_min, alpha))
            
            # Find which positions have y-intercepts above the required cutoff
            pass_y = filter(pass_y, Y_INT_MIN > y_intercept)
            
            
            if(nrow(pass_y) > 0){
              
              # Make string to represent the nucleotides at positions of interest
              superSubstring = function(string, positions){
                chars = sapply(positions, function(x) substring(string, x, x))
                return(paste(chars, collapse=""))
              }
              
              gl_substring = superSubstring(germline, pass_y$POSITION)
              gl_minus_substring = insertPolymorphisms(germline, pass_y$POSITION,
                                                       rep("N", nrow(pass_y)))
              
              # Find the potential SNP positions and remove anything that matches
              # the germline at all those positions or any combo that is too rare
              db_y_subset = db_subset %>%
                group_by(1:n()) %>%
                mutate(SNP_STRING = superSubstring(SEQUENCE_IMGT, pass_y$POSITION)) %>%
                filter(SNP_STRING != gl_substring) %>%
                group_by(SNP_STRING) %>%
                mutate(STRING_COUNT = n()) %>%
                filter(STRING_COUNT >= min_seqs)
              if (nrow(db_y_subset) > 0 ){
                # Get mutation count at all positions that are not potential SNPs
                db_y_subset$MUT_COUNT_MINUS_SUBSTRING = db_y_subset$SEQUENCE_IMGT %>%
                  substring(nt_min, nt_max) %>%
                  paste(pads, ., sep="") %>% # padding is "-" repreated for nt_min-1
                  getMutatedPositions(gl_minus_substring) %>%
                  sapply(length)
                # Keep only unmutated seqences and then find the counts of J and
                # junction length for each of the SNP strings, and then check to
                # see which pass the j/junction and count requirements
                db_y_summary = db_y_subset %>%
                  filter(MUT_COUNT_MINUS_SUBSTRING == 0) %>%
                  mutate(J_GENE = getGene(J_CALL)) %>%
                  group_by(SNP_STRING, J_GENE, JUNCTION_LENGTH) %>%
                  summarise(COUNT = n()) %>%
                  group_by(SNP_STRING) %>%
                  mutate(FRACTION = COUNT/sum(COUNT)) %>%
                  summarise(TOTAL_COUNT = sum(COUNT), MAX_FRAC = max(FRACTION)) %>%
                  filter(TOTAL_COUNT >= min_seqs & MAX_FRAC <= j_max)
                
                if(nrow(db_y_summary) > 0){
                  germ_nts = unlist(strsplit(gl_substring,""))
                  for (r in 1:nrow(db_y_summary)) {
                    # Create the new germline
                    snp_nts = unlist(strsplit(db_y_summary$SNP_STRING[r],""))
                    germ = insertPolymorphisms(germline, pass_y$POSITION, snp_nts)
                    names(germ) = mapply(paste, germ_nts, pass_y$POSITION,
                                         snp_nts, sep="") %>%
                      paste(collapse="_") %>%
                      paste(names(germline), ., sep="_")
                    # Save the new germline to our data frame               
                    if (!is.na(df_run$NOVEL_IMGT[1])){
                      df_run = bind_rows(df_run[1,], df_run)
                    }
                    df_run$POLYMORPHISM_CALL[1] = names(germ)
                    df_run$NOVEL_IMGT[1] =  as.character(germ)
                    df_run$PERFECT_MATCH_COUNT[1] = db_y_summary$TOTAL_COUNT[r]
                  }
                } # end pass j/junction
              } # end nrow db_y_subset
            } # end pass_y
          } # end seq count after mutation filtering
        }
      } # end for each mutcount
      return(df_run)
    } # end foreach
  if(nproc > 1) { stopCluster(cluster) }
  return(df_out)
}



Rprof("tiggerProfile.out")
system.time( novel <- findNovelAlleles(dat, germline_ighv,y_intercept = 1/4) )
filter(novel, !is.na(POLYMORPHISM_CALL))
Rprof(NULL)
summaryRprof("tiggerProfile.out")
# 
# Plot positional mut freq vs seq mut count    
# Join y intercepts to positional mutation rate and proceed
pos_muts = left_join(pos_muts, pass_y, by = "POSITION")
#pos_muts$MUT_COUNT = factor(pos_muts$MUT_COUNT, levels = 1:10)
ggplot(pos_muts, aes(factor(MUT_COUNT), POS_MUT_RATE, group=POSITION,
                     color=POSITION %in% c(112, 222, 286))) +
  geom_line(show_guide=FALSE,
            size = 1, alpha = 0.25) +
  ylim(0,1) + 
 # scale_x_discrete(breaks = 1:10, labels = 1:10, limits = 1:10) +
  xlab("Mutation Count (Sequence)") + ylab("Mutation Frequency (Position)") +
  theme_bw()




plotTigger(dat, filter(novel_light, nchar(NOVEL_IMGT) > 3)[1,])

plotTigger <- function(clip_db, novel_df_row, plot_linear_model=TRUE){
  
  min_frac=0.75 # Need to integrate this into the tigger result
  
  # Use the data frame
  if(length(novel_df_row) > 0){
    if(is.data.frame(novel_df_row) & nrow(novel_df_row) == 1){
      pos_range = novel_df_row$NT_MIN:novel_df_row$NT_MAX
      germline = novel_df_row$GERMLINE_IMGT
      names(germline) = novel_df_row$GERMLINE_CALL
      mut_range = novel_df_row$MUT_MIN:novel_df_row$MUT_MAX
      y_intercept = novel_df_row$Y_INTERCEPT
      alpha = novel_df_row$ALPHA
      novel_imgt = novel_df_row$NOVEL_IMGT
    } else {
      stop("novel_df_row is not a data frame with only one row.")
    }
  }
  
  germline = germline %>%
    toupper() %>%
    gsub(".", "-", . , fixed = TRUE) %>%
    gsub("[^ACGT-]", "N", .)
  clip_db$SEQUENCE_IMGT = clip_db$SEQUENCE_IMGT %>%
    toupper() %>%
    gsub(".", "-", . , fixed = TRUE) %>%
    gsub("[^ACGT-]", "N", .)
  # Extract sequences assigned to the germline, determine which
  # have an appropriate range of mutations, and find the mutation
  # frequency of each position
  db_subset = clip_db %>%
    select(SEQUENCE_IMGT, V_CALL, J_CALL, JUNCTION_LENGTH) %>%
    filter(grepl(names(germline), V_CALL, fixed=T))
  pos_db = db_subset %>%  
    mutationRangeSubset(germline, mut_range, pos_range) %>%
    positionMutations(germline, pos_range)
  
  pos_muts = pos_db %>%
    group_by(POSITION) %>%
    mutate(PASS = mean(OBSERVED) >= min_frac) %>%
    group_by(MUT_COUNT, POSITION) %>%
    summarise(POS_MUT_RATE = mean(MUTATED)*unique(PASS) ) %>% 
    ungroup()
    
  # Calculate y intercepts
  # and find which positions have y-intercepts above the required cutoff
  pass_y = pos_muts %>%
    group_by(POSITION) %>%
    summarise(Y_INT_MIN = findLowerY(POS_MUT_RATE, MUT_COUNT, mut_min, alpha)) %>%
    filter(Y_INT_MIN > y_intercept)
  
  pos_muts = left_join(pos_muts, pass_y, by = "POSITION")
  
  if(nrow(pass_y) > 0){
    lms = lapply(pass_y$POSITION,
                 function(x) coef(lm(POS_MUT_RATE ~ MUT_COUNT,
                                     filter(pos_muts, POSITION == x)))) %>%
      unlist() %>%
      matrix(ncol = 2,byrow = TRUE) %>%
      as.data.frame()
  }

  db_subset$MUT_COUNT_NOVEL = db_subset$SEQUENCE_IMGT %>%
    substring(min(pos_range), max(pos_range)) %>%
    paste(pads, ., sep="") %>%
    getMutatedPositions(novel_imgt) %>%
    sapply(length)
  db_subset = db_subset %>%
    filter(MUT_COUNT_NOVEL == 0) %>%
    mutate(J_GENE = getGene(J_CALL))
  
  pos = pass_y$POSITION
  pos_muts = mutate(pos_muts, Polymorphic = POSITION %in% pos)
  # MAKE THE FIRST PLOT
  POLYCOLORS = setNames(DNA_COLORS[c(4,3)], c("FALSE", "TRUE"))
  p1 = ggplot(pos_muts, aes(factor(MUT_COUNT), POS_MUT_RATE, group=POSITION,
                           color=Polymorphic)) +
    geom_line(size = 0.75) +
    scale_color_manual(values = POLYCOLORS) +
    ylim(0,1) +
    xlab("Mutation Count (Sequence)") +
    ylab("Mutation Frequency (Position)") +
    theme_bw() +
    theme(legend.position=c(1,1), legend.justification=c(1,1)) +
    guides(color = guide_legend(ncol = 2))
  # MAKE THE SECOND PLOT
  p2 = ggplot(mutate(filter(pos_db, POSITION %in% pos),
                POSITION = paste("Position", POSITION)),
         aes(factor(MUT_COUNT), fill=NT)) +
    geom_bar(binwidth=1) +
    guides(fill = guide_legend("Nucleotide", ncol = 4)) +
    facet_grid(POSITION ~ .) +
    xlab("Mutation Count (Sequence)") + ylab("Sequence Count") +
    scale_fill_manual(values = DNA_COLORS, breaks=names(DNA_COLORS)) +
    theme_bw() +
    theme(legend.position=c(1,1), legend.justification=c(1,1))
  # MAKE THE THIRD PLOT
  p3 = ggplot(db_subset, aes(factor(J_GENE), fill=factor(JUNCTION_LENGTH))) +
    geom_bar() +
    guides(fill = guide_legend("Junction Length", ncol = 3)) +
    xlab("J Gene") + ylab("Unmutated Sequence Count") +
    theme_bw()
  
  multiplot(p1,p2,p3, cols = 3)

}

#   
# # Plot nucleoide usage vs seq mut count for given position
# DNA_COLORS <- c("A"="#64F73F", 
#                 "C"="#FFB340", 
#                 "G"="#EB413C", 
#                 "T"="#3C88EE",
#                 "N"="grey")
# p = pass_y$POSITION
p = c(112, 222, 286)
ggplot(mutate(filter(pos_db,POSITION %in% p),
              POSITION = paste("Position", POSITION)),
       aes(factor(MUT_COUNT), fill=NT)) +
  geom_bar(binwidth=1) +
  guides(fill = guide_legend("Nucleotide")) +
  facet_grid(POSITION ~ .) +
  xlab("Mutation Count (Sequence)") + ylab("Sequence Count") +
  scale_fill_manual(values = DNA_COLORS, breaks=names(DNA_COLORS)) +
  theme_bw()

# # Plot the J usage for sequences perfectly matching a given predicted allele
# ggplot(db_junc, aes(J_GENE, COUNT, fill=JUNCTION_LENGTH)) +
#   geom_bar(stat="identity") +
#   guides(fill = guide_legend("Junction Length")) +
#   xlab("J Gene") + ylab("Sequence Count") +
#   theme_bw()
# 


ops <- function(fun){ names(which(nchar(formals(fun))>0)) }

# Multiple plot function
# by  Winston Chang via Cookbook for R
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

