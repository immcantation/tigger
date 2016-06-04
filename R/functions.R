
# The TIgGER Trifecta -----------------------------------------------------

#' Find novel alleles from repertoire sequencing data
#'
#' \code{findNovelAlleles} analyzes mutation patterns in sequences thought to
#' align to each germline allele in order to determine which positions
#' might be polymorphic.
#' 
#' @details A \code{data.frame} in Change-O format contains the following
#' columns:
#' \itemize{
#'   \item \code{"SEQUENCE_IMGT"} containing the IMGT-gapped nucleotide sequence
#'   \item \code{"V_CALL"} containing the IMGT/V-QUEST V allele call(s)
#'   \item \code{"J_CALL"} containing the IMGT/V-QUEST J allele call(s)
#'   \item \code{"JUNCTION_LENGTH"} containing the junction length
#' }
#' The TIgGER allele-finding algorithm, briefly, works as follows:
#' Mutations are determined through comparison to the provided germline.
#' Mutation frequency at each *position* is determined as a function of
#' *sequence-wide* mutation counts. Polymorphic positions exhibit a high
#' mutation frequency despite sequence-wide mutation count. False positive of
#' potential novel alleles resulting from clonally-related sequences are guarded
#' against by ensuring that sequences perfectly matching the potential novel
#' allele utilize a wide range of combinations of J gene and junction length.
#' 
#' @param    clip_db        a \code{data.frame} in Change-O format. See details.
#' @param    germline_db    a vector of named nucleotide germline sequences
#'                          matching the V calls in \code{clip_db}
#' @param    germline_min   the minimum number of sequences that must have a
#'                          particular germline allele call for the allele to
#'                          be analyzed
#' @param    nproc          the number of processors to use
#' @param    auto_mutrange  if \code{TRUE}, the algorithm will attempt to
#'                          determine the appropriate mutation range
#'                          automatically using the mutation count of the most
#'                          common sequence assigned to each allele analyzed
#' @param    mut_range      the range of mutations that sampled may carry and
#'                          be considered by the algorithm
#' @param    pos_range      the range of IMGT-numbered positions that should be
#'                          considered by the algorithm
#' @param    alpha          the alpha cutoff to be used when constructing the
#'                          confidence interval for the y-intercept           
#' @param    y_intercept    the y-intercept above which positions should be
#'                          considered potentially polymorphic
#' @param    j_max          the maximum fraction of sequences perfectly aligning
#'                          to a potential novel allele that are allowed to
#'                          utilize to a particular combination of junction
#'                          length and J gene
#' @param    min_seqs       the minimum number of total sequences (within the
#'                          desired mutational range and nucleotide range)
#'                          required for the samples to be considered
#' @param    min_frac       the minimum fraction of sequences that must have
#'                          usable nucleotides in a given position for that
#'                          position to considered
#' @return   a \code{data.frame} with a row for each known allele analyzed.
#' Besides metadata on the the parameters used in the search, each row will have
#' either a note as to where the polymorphism-finding algorithm exited or a
#' nucleotide sequence for the predicted novel allele.
#' 
#' @seealso \link{plotNovel} to visualize the data supporting any
#' novel alleles hypothesized to be present in the data and
#' \link{inferGenotype} to determine if the novel alleles are frequent
#' enought to be included in the subject's genotype
#' 
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Find novel alleles and return relevant data
#' \dontrun{novel_df = findNovelAlleles(sample_db, germline_ighv)}
#' 
#' @export
findNovelAlleles  <- function(clip_db, germline_db,
                              germline_min = 200,
                              nproc = 4,
                              min_seqs = 50,
                              auto_mutrange = TRUE,
                              mut_range = 1:10,
                              pos_range = 1:312,
                              y_intercept = 0.125,
                              alpha = 0.05,
                              j_max = 0.15,
                              min_frac = 0.75){
  . = a = NULL
  # Keep only the columns we need and clean up the sequences
  missing = c("SEQUENCE_IMGT", "V_CALL", "J_CALL", "JUNCTION_LENGTH") %>%
    setdiff(colnames(clip_db))
  if (length(missing) != 0) {
    stop("Could not find required columns in clip_db:\n  ",
         paste(missing, collapse="\n  "))
  }
  empty_junctions = sum(clip_db$JUNCTION_LENGTH == 0, na.rm=TRUE)
  if (empty_junctions > 0) {
    stop(empty_junctions, " sequences have junction ", "length of zero. ",
         "Please remove these sequences.")
  }
  germlines = cleanSeqs(germline_db)
  names(germlines) = getAllele(names(germlines), first=FALSE, strip_d=FALSE)
  clip_db$SEQUENCE_IMGT = cleanSeqs(clip_db$SEQUENCE_IMGT)
  
  
  # Find which rows' calls contain which germline alleles
  cutoff =
    ifelse(germline_min < 1, round(nrow(clip_db)*germline_min), germline_min)
  allele_groups = sapply(names(germlines), grep, clip_db$V_CALL, fixed=TRUE,
                         simplify=FALSE)
  names(allele_groups) = names(germlines)
  allele_groups = allele_groups[sapply(allele_groups, length) >= cutoff]
  if(length(allele_groups) == 0){
    paste("Not enough sample sequences were assigned to any germline:\n",
          " (1) germline_min is too large or\n",
          " (2) sequences names don't match germlines.") %>%
      stop()
  }
  allele_groups = allele_groups[sortAlleles(names(allele_groups))]
  
  # Prepare for parallel processing
  nproc = ifelse(Sys.info()['sysname'] == "Windows",
                 Sys.getenv('NUMBER_OF_PROCESSORS'),
            ifelse(Sys.info()['sysname'] == "Darwin",  
                   system("sysctl -n hw.ncpu", intern=TRUE),
                   system("nproc", intern=TRUE))) %>%
    as.numeric() %>% 
    min(nproc, . - 1) %>%
    max(1, .)
  if(nproc == 1) {
    registerDoSEQ()
  } else {
    cluster <- makeCluster(nproc, type = "SOCK")
    clusterExport(cluster, list("allele_groups",
                                "germlines",
                                "clip_db",
                                "min_seqs",
                                "auto_mutrange",
                                "mut_range",
                                "pos_range",
                                "y_intercept",
                                "alpha",
                                "j_max",
                                "germline_min",
                                "min_frac"), 
                   envir=environment())
    clusterEvalQ(cluster, library(tigger))
    registerDoParallel(cluster)
  }
  
  df_out <- foreach(a=icount(length(allele_groups)), .combine=rbind) %dopar% {
    
    allele_name = names(allele_groups)[a]
    
    # Subset of data being analyzed
    germline = germlines[allele_name]
    indicies = allele_groups[[allele_name]]
    db_subset = slice_(clip_db, ~indicies)
    
    # If mutrange is auto, find most popular mutation count and start from there
    gpm = db_subset %>%
      mutate_(V_CALL = ~allele_name) %>%
      getPopularMutationCount(germline,
                              gene_min=0, seq_min=min_seqs,
                              seq_p_of_max=1/8, full_return=TRUE)
    
    # Determine the mutation range(s) to scan
    mut_mins = min(mut_range)
    if(auto_mutrange & sum(gpm$MUTATION_COUNT > 0) > 0 ){
      mut_mins = c(mut_mins, gpm$MUTATION_COUNT[gpm$MUTATION_COUNT > 0]) %>%
        unique() %>%
        sort()
    }
    
    # Create the run's return object
    df_run_empty = data.frame(GERMLINE_CALL = names(germline),
                              NOTE = "",
                              POLYMORPHISM_CALL = NA,
                              NOVEL_IMGT = NA,
                              PERFECT_MATCH_COUNT = NA,
                              GERMLINE_CALL_COUNT = length(indicies),
                              MUT_MIN = NA,
                              MUT_MAX = NA,
                              GERMLINE_IMGT = as.character(germline),
                              POS_MIN = min(pos_range),
                              POS_MAX = max(pos_range),
                              Y_INTERCEPT = y_intercept,
                              ALPHA = alpha,
                              MIN_SEQS = min_seqs,
                              J_MAX = j_max,
                              MIN_FRAC = min_frac,
                              stringsAsFactors = FALSE)
    
    for (mut_min in rev(mut_mins)) {
      
      if (mut_min == rev(mut_mins)[1]){
        df_run = df_run_empty
      } else {
        df_run = bind_rows(df_run_empty, df_run)
      }
      mut_max = mut_min + diff(range(mut_range))
      df_run$MUT_MIN[1] = mut_min
      df_run$MUT_MAX[1] = mut_max
      
      # If no sequence is frequent enough to pass the J test, give up now
      if(length(gpm) < 1) {
        df_run$NOTE[1] = "Plurality sequence too rare."
        if(mut_mins[1] == mut_min){
          return(df_run)
        } else {
          next
        }
      }
      
      # Add a mutation count column and filter out sequences not in our range
      db_subset_mm = mutationRangeSubset(db_subset, germline,
                                                  mut_min:mut_max, pos_range)
      
      if(nrow(db_subset_mm) < germline_min){
        df_run$NOTE[1] = "Insufficient sequences in desired mutational range."
        if(mut_mins[1] == mut_min){
          return(df_run)
        } else {
          next
        }
      }
      
      # Duplicate each sequence for all the positions to be analyzed
      # and find which positions are mutated
      pos_db = positionMutations(db_subset_mm, germline, pos_range)
      
      # Find positional mut freq vs seq mut count
      pos_muts = pos_db %>%
        group_by_(~POSITION) %>%
        mutate_(PASS = ~mean(OBSERVED) >= min_frac) %>%
        group_by_(~MUT_COUNT, ~POSITION) %>%
        summarise_(POS_MUT_RATE = ~ mean(MUTATED)*unique(PASS) ) %>% 
        ungroup()   
      
      # Calculate y intercepts, find which pass the test
      pass_y = pos_muts %>%
        group_by_(~POSITION) %>%
        summarise_(Y_INT_MIN = ~findLowerY(POS_MUT_RATE, MUT_COUNT,
                                                  mut_min, alpha)) %>%
        filter_(~Y_INT_MIN > y_intercept)
      
      if(nrow(pass_y) < 1){
        df_run$NOTE[1] = "No positions pass y-intercept test."
        if(mut_mins[1] == mut_min){
          return(df_run)
        } else {
          next
        }
      }
      
      gl_substring = superSubstring(germline, pass_y$POSITION)
      gl_minus_substring = insertPolymorphisms(germline, pass_y$POSITION,
                                               rep("N", nrow(pass_y)))
      
      # Find the potential SNP positions and remove anything that matches
      # the germline at all those positions or any combo that is too rare
      db_y_subset_mm = db_subset_mm %>%
        group_by(1:n()) %>%
        mutate_(SNP_STRING = ~superSubstring(SEQUENCE_IMGT,
                                                    pass_y$POSITION)) %>%
        filter_(~SNP_STRING != gl_substring) %>%
        group_by_(~SNP_STRING) %>%
        mutate_(STRING_COUNT = ~n()) %>%
        filter_(~STRING_COUNT >= min_seqs)
      
      if (nrow(db_y_subset_mm) < 1 ){
        df_run$NOTE[1] = paste("Position(s) passed y-intercept (",
                               paste(pass_y$POSITION, collapse = ","),
                               ") but the plurality sequence is too rare.",
                               sep="")
        if(mut_mins[1] == mut_min){
          return(df_run)
        } else {
          next
        }
      }
      
      # Get mutation count at all positions that are not potential SNPs
      pads = paste(rep("-", min(pos_range)-1), collapse="")
      db_y_subset_mm$MUT_COUNT_MINUS_SUBSTRING = db_y_subset_mm$SEQUENCE_IMGT %>%
        substring(min(pos_range), max(pos_range)) %>%
        paste(pads, ., sep="") %>% 
        getMutatedPositions(gl_minus_substring) %>%
        sapply(length)
      
      # Keep only unmutated seqences and then find the counts of J and
      # junction length for each of the SNP strings, and then check to
      # see which pass the j/junction and count requirements
      db_y_summary0 = db_y_subset_mm %>%
        filter_(~MUT_COUNT_MINUS_SUBSTRING == 0) %>%
        mutate_(J_GENE = ~getGene(J_CALL)) %>%
        group_by_(~SNP_STRING, ~J_GENE, ~JUNCTION_LENGTH) %>%
        summarise_(COUNT = ~n()) %>%
        group_by_(~SNP_STRING) %>%
        mutate_(FRACTION = ~COUNT/sum(COUNT)) %>%
        summarise_(TOTAL_COUNT = ~sum(COUNT), MAX_FRAC = ~max(FRACTION))
        
      if(nrow(db_y_summary0) < 1){
        df_run$NOTE[1] = paste("Position(s) passed y-intercept (",
                               paste(pass_y$POSITION, collapse = ","),
                               ") but no unmutated versions of novel allele",
                               " found.", sep="")
        if(mut_mins[1] == mut_min){
          return(df_run)
        } else {
          next
        }
      }
      
      db_y_summary = db_y_summary0 %>%
        filter_(~TOTAL_COUNT >= min_seqs & MAX_FRAC <= j_max)
      
      if(nrow(db_y_summary) < 1){
        df_run$NOTE[1] = paste("Position(s) passed y-intercept (",
                               paste(pass_y$POSITION, sep = ","),
                               ") but a ",
                               "J-junction combination is too prevalent (",
                               round(100*max(db_y_summary0$MAX_FRAC),1),
                               "% of sequences).",
                               sep="")
        df_run$PERFECT_MATCH_COUNT = max(db_y_summary0$TOTAL_COUNT)
        if(mut_mins[1] == mut_min){
          return(df_run)
        } else {
          next
        }
      }
      
      germ_nts = unlist(strsplit(gl_substring,""))
      for (r in 1:nrow(db_y_summary)) {
        if (r > 1){
          df_run = bind_rows(df_run[1,], df_run)
        }
        # Create the new germline
        snp_nts = unlist(strsplit(db_y_summary$SNP_STRING[r],""))
        remain_mut = db_y_summary$SNP_STRING[r] %>%
          getMutatedPositions(gl_substring) %>%
          unlist() %>%
          unique()
        germ = insertPolymorphisms(germline, pass_y$POSITION, snp_nts)
        names(germ) = mapply(paste, germ_nts[remain_mut],
                             pass_y$POSITION[remain_mut],
                             snp_nts[remain_mut], sep="") %>%
          paste(collapse="_") %>%
          paste(names(germline), ., sep="_")
        # Save the new germline to our data frame               
        df_run$POLYMORPHISM_CALL[1] = names(germ)
        df_run$NOVEL_IMGT[1] =  as.character(germ)
        df_run$PERFECT_MATCH_COUNT[1] = db_y_summary$TOTAL_COUNT[r]
        df_run$NOTE = "Novel allele found!"
      }
      
    } # end for each starting mutation counts
    return(df_run)
  } # end foreach allele
  if(nproc > 1) { stopCluster(cluster) }
  return(df_out)
}

#' Select rows containing novel alleles
#' 
#' \code{selectNovel} takes the result from \link{findNovelAlleles} and
#' selects only the rows containing unique, novel alleles.
#' 
#' @param   novel_df        A \code{data.frame} of the type returned by
#'                          \link{findNovelAlleles}
#' @param   keep_alleles    A \code{logical} indicating if different alleles
#'                          leading to the same novel sequence should be kept.
#'                          See details.
#'                          
#' @details  If, for instance, subject has in his genome IGHV1-2*02 and a novel 
#' allele equally close to IGHV1-2*02 and IGHV1-2*05, the novel allele may be
#' detected by analyzing sequences that best align to either of these alleles.
#' If \code{keep_alleles} is \code{TRUE}, both polymorphic allele calls will
#' be retained. In the case that multiple mutation ranges are checked for the
#' same allele, only one mutation range will be kept in the output.
#'                        
#' @return  A \code{data.frame} containing only unique, novel alleles (if any)
#' that were in the input.
#' 
#' @examples
#' data(novel_df)
#' novel = selectNovel(novel_df)
#' 
#' @export
selectNovel <- function(novel_df, keep_alleles=FALSE) {
  if (keep_alleles) {
    novel_df = novel_df %>% group_by_(~GERMLINE_CALL)
  }
  novel = novel_df %>%
    distinct_(~NOVEL_IMGT) %>%
    filter_(~nchar(NOVEL_IMGT) > 2)
  return(novel)
}

#' Visualize evidence of novel V alleles
#'
#' \code{plotNovel} is be used to visualize the evidence of any novel V
#' alleles found using \link{findNovelAlleles}.
#' 
#' @param    clip_db        a \code{data.frame} in Change-O format. See
#'                          \link{findNovelAlleles} for details.
#' @param    novel_df_row   a single row from a data frame as output by
#'                          \link{findNovelAlleles} that contains a
#'                          polymorphism-containing germline allele
#' @param    ncol           number of columns to use when laying out the plots            
#' @return   NULL
#' 
#' @examples
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Find novel alleles and return relevant data
#' \dontrun{novel_df = findNovelAlleles(sample_db, germline_ighv)}
#' data(novel_df)
#' # Plot the evidence for the first (and only) novel allele in the example data
#' novel = selectNovel(novel_df)
#' plotNovel(sample_db, novel[1,])
#' 
#' @export
plotNovel <- function(clip_db, novel_df_row, ncol = 1){
  . = NULL
    
  # Use the data frame
  if(length(novel_df_row) > 0){
    if(is.data.frame(novel_df_row) & nrow(novel_df_row) == 1){
      pos_range = novel_df_row$POS_MIN:novel_df_row$POS_MAX
      germline = novel_df_row$GERMLINE_IMGT
      names(germline) = novel_df_row$GERMLINE_CALL
      mut_range = novel_df_row$MUT_MIN[1]:novel_df_row$MUT_MAX[1]
      novel_imgt = novel_df_row$NOVEL_IMGT
      names(novel_imgt) = novel_df_row$POLYMORPHISM_CALL
      min_frac = novel_df_row$MIN_FRAC
    } else {
      stop("novel_df_row is not a data frame with only one row.")
    }
  }
  
  germline = cleanSeqs(germline)
  clip_db$SEQUENCE_IMGT = cleanSeqs(clip_db$SEQUENCE_IMGT)
  
  # Extract sequences assigned to the germline, determine which
  # have an appropriate range of mutations, and find the mutation
  # frequency of each position
  db_subset = clip_db %>%
    select_(~SEQUENCE_IMGT, ~V_CALL, ~J_CALL, ~JUNCTION_LENGTH) %>%
    filter_(~grepl(names(germline), V_CALL, fixed=TRUE))
  pos_db = db_subset %>%  
    mutationRangeSubset(germline, mut_range, pos_range) %>%
    positionMutations(germline, pos_range)
  pos_muts = pos_db %>%
    group_by_(~POSITION) %>%
    mutate_(PASS = ~mean(OBSERVED) >= min_frac) %>%
    group_by_(~MUT_COUNT, ~POSITION) %>%
    summarise_(POS_MUT_RATE = ~mean(MUTATED)*unique(PASS) ) %>% 
    ungroup()
  
  # Label the polymorphic positions as such
  pass_y = unlist(strsplit(names(novel_imgt), "_"))[-1] %>%
    gsub("[^0-9]", "", .) %>%
    as.numeric()
  p_y_f = unlist(strsplit(names(novel_imgt), "_"))[-1] %>%
    gsub("[0-9]+.", "", .)
  p_y_t = unlist(strsplit(names(novel_imgt), "_"))[-1] %>%
    gsub(".[0-9]+", "", .)
  to_from = paste(paste("Position", pass_y), paste(paste(p_y_f, "->"), p_y_t))
  names(to_from) = pass_y
  pos_muts = pos_muts %>%
    mutate_(Polymorphic = ~ifelse(POSITION %in% pass_y, "True", "False"))
  
  pads = paste(rep("-", min(pos_range)-1), collapse="")
  db_subset$MUT_COUNT_NOVEL = db_subset$SEQUENCE_IMGT %>%
    substring(min(pos_range), max(pos_range)) %>%
    paste(pads, ., sep="") %>%
    getMutatedPositions(novel_imgt) %>%
    sapply(length)
  db_subset = db_subset %>%
    filter_(~MUT_COUNT_NOVEL == 0) %>%
    mutate_(J_GENE = ~getGene(J_CALL))
  db_subset$JUNCTION_LENGTH = db_subset$JUNCTION_LENGTH %>%
    factor(levels=min(db_subset$JUNCTION_LENGTH):max(db_subset$JUNCTION_LENGTH))
  pos_muts$Polymorphic = pos_muts$Polymorphic %>%
    factor(levels = c("False", "True"))
  pos_db$NT = pos_db$NT %>%
    factor(levels = names(DNA_COLORS))
  pos_muts$GERMLINE = names(germline)
  
  # MAKE THE FIRST PLOT
  POLYCOLORS = setNames(DNA_COLORS[c(4,3)], c("False", "True"))
  p1 = ggplot(pos_muts, aes_(~factor(MUT_COUNT), ~POS_MUT_RATE, group=~POSITION,
                            color=~Polymorphic)) +
    geom_line(size = 0.75) +
    facet_grid(GERMLINE ~ .) +
    scale_color_manual(values = POLYCOLORS) +
    ylim(0,1) +
    xlab("Mutation Count (Sequence)") +
    ylab("Mutation Frequency (Position)") +
    theme_bw() +
    theme(legend.position=c(0.5,0.9), legend.justification=c(0.5,1),
          legend.background=element_rect(fill = "transparent")) +
    guides(color = guide_legend(ncol = 2, reverse = TRUE))
  # MAKE THE SECOND PLOT
  p2 = ggplot(mutate_(filter_(pos_db, ~POSITION %in% pass_y),
                     POSITION = ~to_from[as.character(POSITION)]),
              aes_(~factor(MUT_COUNT), fill=~NT)) +
    geom_bar(width=0.9) +
    guides(fill = guide_legend("Nucleotide", ncol = 4)) +
    facet_grid(POSITION ~ .) +
    xlab("Mutation Count (Sequence)") + ylab("Sequence Count") +
    scale_fill_manual(values = DNA_COLORS, breaks=names(DNA_COLORS),
                      drop=FALSE) +
    theme_bw() +
    theme(legend.position=c(1,1), legend.justification=c(1,1),
          legend.background=element_rect(fill = "transparent"))
  # MAKE THE THIRD PLOT
  p3 = ggplot(db_subset, aes_(~JUNCTION_LENGTH, fill=~factor(J_GENE))) +
    geom_bar(width=0.9) +
    guides(fill = guide_legend("J Gene", ncol = 2)) +
    xlab("Junction Length") + ylab("Unmutated Sequence Count") +
    theme_bw() +
    theme(legend.position=c(1,1), legend.justification=c(1,1),
          legend.background=element_rect(fill = "transparent"))
  
  multiplot(p1,p2,p3, cols = ncol)
}

#' Infer a subject-specific genotype
#'
#' \code{inferGenotype} infers an subject's genotype by finding the minimum
#' number set of alleles that can explain the majority of each gene's calls. The
#' most common allele of each gene is included in the genotype first, and the
#' next most common allele is added until the desired fraction of alleles can be
#' explained. In this way, mistaken allele calls (resulting from sequences which
#' by chance have been mutated to look like another allele) can be removed.
#' 
#' @param    clip_db              a \code{data.frame} containing V allele
#'                                calls from a single subject under
#'                                \code{"V_CALL"}. If
#'                                \code{find_unmutated} is \code{TRUE}, then
#'                                the sample IMGT-gapped V(D)J sequence should 
#'                                be provided in a column \code{"SEQUENCE_IMGT"}
#' @param    fraction_to_explain  the portion of each gene that must be
#'                                explained by the alleles that will be included
#'                                in the genotype
#' @param    gene_cutoff          either a number of sequences or a fraction of
#'                                the length of \code{allele_calls} denoting the
#'                                minimum number of times a gene must be
#'                                observed in \code{allele_calls} to be included
#'                                in the genotype
#' @param    find_unmutated       if \code{TRUE}, use \code{germline_db} to
#'                                find which samples are unmutated. Not needed
#'                                if \code{allele_calls} only represent
#'                                unmutated samples.
#' @param    germline_db          named vector of sequences containing the
#'                                germline sequences named in
#'                                \code{allele_calls}. Only required if
#'                                \code{find_unmutated} is \code{TRUE}.
#' @param    novel_df             an optional \code{data.frame} of the type
#'                                novel returned by
#'                                \link{findNovelAlleles} containing
#'                                germline sequences that will be utilized if
#'                                \code{find_unmutated} is \code{TRUE}. See
#'                                details.
#' @details  Allele calls representing cases where multiple alleles have been
#'           assigned to a single sample sequence are rare among unmutated
#'           sequences but may result if nucleotides for certain positions are
#'           not available. Calls containing multiple alleles are treated as
#'           belonging to all groups. If \code{novel_df} is provided, all
#'           sequences that are assigned to the same starting allele as any
#'           novel germline allele will have the novel germline allele appended
#'           to their assignent prior to searching for unmutated sequences.
#' 
#' @return   A table of alleles denoting the genotype of the subject
#' 
#' @note     This method works best with data derived from blood, where a large
#'           portion of sequences are expected to be unmutated. Ideally, there
#'           should be hundreds of allele calls per gene in the input.
#' 
#' @examples
#' # Infer the IGHV genotype, using only unmutated sequences, including any 
#' # novel alleles
#' data(sample_db)
#' data(germline_ighv)
#' data(novel_df)
#' inferGenotype(sample_db, find_unmutated = TRUE, germline_db = germline_ighv,
#'               novel_df = novel_df)
#' 
#' @seealso \link{plotGenotype} for a colorful visualization and
#'          \link{genotypeFasta} to convert the genotype to nucleotide sequences.
#' 
#' @export
inferGenotype <- function(clip_db, fraction_to_explain = 0.875,
                          gene_cutoff = 1e-4, find_unmutated = TRUE,
                          germline_db = NA, novel_df = NA){
  
  . = NULL
  allele_calls = getAllele(clip_db$V_CALL, first=FALSE, strip_d=FALSE)
  # Find the unmutated subset, if requested
  if(find_unmutated){
    if(is.na(germline_db[1])){
      stop("germline_db needed if find_unmutated is TRUE")
    }
    if(!is.null(nrow(novel_df))){
      novel_df = filter_(novel_df, ~!is.na(POLYMORPHISM_CALL)) %>%
        select_(~GERMLINE_CALL, ~POLYMORPHISM_CALL, ~NOVEL_IMGT)
      if(nrow(novel_df) > 0){
        # Extract novel alleles if any and add them to germline_db
        novel_gl = novel_df$NOVEL_IMGT
        names(novel_gl) = novel_df$POLYMORPHISM_CALL
        germline_db = c(germline_db, novel_gl)
        # Add the novel allele calls to allele calls of the same starting allele
        for(r in 1:nrow(novel_df)){
          ind = grep(novel_df$GERMLINE_CALL[r], allele_calls, fixed=TRUE)
          allele_calls[ind] = allele_calls[ind] %>%
            sapply(paste, novel_df$POLYMORPHISM_CALL[r], sep=",")
        }
      }
    }
    # Find unmutated sequences
    allele_calls = findUnmutatedCalls(allele_calls,
                                      as.character(clip_db$SEQUENCE_IMGT),
                                      germline_db)
    if(length(allele_calls) == 0){
      stop("No unmutated sequences found! Set 'find_unmutated' to 'FALSE'.")
    }
  }

  # Find which rows' calls contain which genes
  cutoff = ifelse(gene_cutoff < 1, length(allele_calls)*gene_cutoff, gene_cutoff)
  gene_regex = allele_calls %>% strsplit(",") %>% unlist() %>%
    getGene(strip_d=FALSE) %>%  unique() %>% paste("\\*", sep="")
  gene_groups = sapply(gene_regex, grep, allele_calls, simplify=FALSE)
  names(gene_groups) = gsub("\\*", "", gene_regex, fixed=TRUE)
  gene_groups = gene_groups[sapply(gene_groups, length) >= cutoff]
  gene_groups = gene_groups[sortAlleles(names(gene_groups))]
  
  # Make a table to store the resulting genotype
  GENE = names(gene_groups)
  ALLELES = COUNTS = NOTE = rep("", length(GENE))
  TOTAL = sapply(gene_groups, length)
  genotype = cbind(GENE, ALLELES, COUNTS, TOTAL, NOTE)
  
  # For each gene, find which alleles to include
  for (g in GENE){
    # Keep only the part of the allele calls that uses the gene being analyzed
    ac = allele_calls[gene_groups[[g]]] %>%
      strsplit(",") %>%
      lapply(function(x) x[grep(paste(g, "\\*", sep=""), x)]) %>%
      sapply(paste, collapse=",")
    target = ceiling(fraction_to_explain*length(ac)) # how many we need to explain
    t_ac = table(ac) # table of allele calls
    potentials = unique(unlist(strsplit(names(t_ac),","))) # potential alleles
    # One allele? Easy!
    if (length(potentials) == 1 | length(t_ac) == 1){
      genotype[genotype[,"GENE"]==g,"ALLELES"] =
        gsub("[^d\\*]*[d\\*]","",potentials )[1]
      genotype[genotype[,"GENE"]==g,"COUNTS"] = t_ac
    } else {
      # More alleles? Let's find the fewest that can explain the needed fraction
      # Make a table of which alleles can explain which calls
      regexpotentials = paste(gsub("\\*","\\\\*", potentials),"$",sep="")
      regexpotentials = 
        paste(regexpotentials,gsub("\\$",",",regexpotentials),sep="|")
      tmat = 
        sapply(regexpotentials, function(x) grepl(x, names(t_ac),fixed=FALSE))
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
      genotype[genotype[,"GENE"]==g,"ALLELES"] =
        paste(gsub("[^d\\*]*[d\\*]","",included ),collapse=",")
      genotype[genotype[,"GENE"]==g,"COUNTS"] =
        paste(counts,collapse=",")
    }
    
  }
  geno = as.data.frame(genotype, stringsAsFactors = FALSE)
  
  # Check for indistinguishable calls
  if(find_unmutated == TRUE){
    seqs = genotypeFasta(geno, germline_db)
    dist_mat = seqs %>%
      sapply(function(x) sapply((getMutatedPositions(seqs, x)), length))
    rownames(dist_mat) = colnames(dist_mat)
    for (i in 1:nrow(dist_mat)){ dist_mat[i,i] = NA }
    same = which(dist_mat == 0, arr.ind=TRUE)
    if (nrow(same) > 0 ) {
      for (r in 1:nrow(same)) {
        inds = as.vector(same[r,])
        geno[getGene(rownames(dist_mat)[inds][1]),]$NOTE =
          paste(rownames(dist_mat)[inds], collapse=" and ") %>%
          paste("Cannot distinguish", .)
      }
    }
  }
  rownames(geno) = NULL
  return(geno)
}

#' Show a colorful representation of a genotype
#'
#' \code{plotGenotype} plots a genotype table.
#' 
#' @param    genotype     a table of alleles denoting a genotype, as returned by
#'                        \link{inferGenotype}
#' @param    facet_by     a column name in \code{genotype} to facet the plot by. 
#'                        If \code{NULL}, then do not facet the plot. 
#' @param    gene_sort    a string defining the method to use when sorting alleles.
#'                        If \code{"name"} then sort in lexicographic order. If
#'                        \code{"position"} then sort by position in the locus, as
#'                        determined by the final two numbers in the gene name.
#' @param    text_size    the point size of the plotted text
#' @param    silent       if \code{TRUE} do not draw the plot and just return the ggplot
#'                        object; if \code{FALSE} draw the plot.
#' @param    ...          additional arguments to pass to ggplot2::theme.
#' 
#' @return  A ggplot object defining the plot.
#' 
#' @seealso \link{inferGenotype}
#' 
#' @examples
#' # Load example data
#' data(novel_df)
#' data(genotype)
#' 
#' # Plot genotype
#' plotGenotype(genotype)
#' 
#' # Facet by subject
#' genotypea = genotypeb = genotype
#' genotypea$SUBJECT = "A"
#' genotypeb$SUBJECT = "B"
#' geno_sub = rbind(genotypea, genotypeb)
#' plotGenotype(geno_sub, facet_by="SUBJECT", gene_sort="pos")
#' 
#' @export
plotGenotype = function(genotype, facet_by=NULL, gene_sort=c("name", "position"), 
                        text_size=12, silent=FALSE, ...) {
  # Check arguments
  gene_sort <- match.arg(gene_sort)
    
  # Split genes' alleles into their own rows
  alleles = strsplit(genotype$ALLELES, ",")
  geno2 = genotype
  r = 1
  for (g in 1:nrow(genotype)){
    for(a in 1:length(alleles[[g]])) {
      geno2[r, ] = genotype[g, ]
      geno2[r, ]$ALLELES = alleles[[g]][a]
      r = r + 1
    }
  }
  
  # Set the gene order
  geno2$GENE = factor(geno2$GENE, 
                      levels=rev(sortAlleles(unique(geno2$GENE), method=gene_sort)))
  
  # Create the base plot
  p = ggplot(geno2, aes_(x=~GENE, fill=~ALLELES)) +
    theme_bw() +
    theme(axis.ticks=element_blank(),
          axis.text.x=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          text=element_text(size=text_size),
          strip.background=element_blank(),
          strip.text=element_text(face="bold")) +
    geom_bar(position="fill") +
    coord_flip() + xlab("Gene") + ylab("") +
    scale_fill_hue(name="Allele", h=c(0, 270), h.start=10)
  
  # Plot, with facets by SUBJECT if that column is present
  if (!is.null(facet_by)) {
    p = p + facet_grid(paste0(".~", facet_by))
  }

  # Add additional theme elements
  p = p + do.call(theme, list(...))
  
  # Plot
  if (!silent) { plot(p) }
  
  invisible(p)
}

#' Return the nucleotide sequences of a genotype
#'
#' \code{genotypeFasta} converts a genotype table into a vector of nucleotide
#' sequences.
#' 
#' @param    genotype     a table of alleles denoting a genotype, as returned by
#'                        \link{inferGenotype}
#' @param    germline_db  a vector of named nucleotide germline sequences
#'                        matching the alleles detailed in \code{genotype}
#' @param    novel_df     an optional \code{data.frame} containing putative
#'                        novel alleeles of the type returned by
#'                        \link{findNovelAlleles}
#' 
#' @return   A named vector of strings containing the germline nucleotide
#'           sequences of the alleles in the provided genotype
#' 
#' @seealso \link{inferGenotype}
#' 
#' @examples
#' # Load example data
#' data(germline_ighv)
#' data(novel_df)
#' data(genotype)
#'                      
#' # Find the sequences that correspond to the genotype
#' genotype_seqs = genotypeFasta(genotype, germline_ighv, novel_df)
#' 
#' 
#' @export
genotypeFasta <- function(genotype, germline_db, novel_df=NA){
  if(!is.null(nrow(novel_df))){
    # Extract novel alleles if any and add them to germline_db
    novel_df = filter_(novel_df, ~!is.na(POLYMORPHISM_CALL)) %>%
      select_(~GERMLINE_CALL, ~POLYMORPHISM_CALL, ~NOVEL_IMGT)
    if(nrow(novel_df) > 0){
      novel_gl = novel_df$NOVEL_IMGT
      names(novel_gl) = novel_df$POLYMORPHISM_CALL
      germline_db = c(germline_db, novel_gl)
    }
  }
  
  genotype$GENE = gsub("D$|d$","",genotype$GENE)
  
  g_names = names(germline_db)
  names(g_names) = gsub("D", "", names(germline_db))
  table_calls = mapply(paste, genotype$GENE, strsplit(genotype$ALLELES, ","),
                       sep="*")
  seqs = germline_db[as.vector(g_names[unlist(table_calls)])]
  if(sum(is.na(seqs)) > 0){
    stop("The following genotype alleles were not found in germline_db: ",
         paste(unlist(table_calls)[which(is.na(seqs))], collapse = ", "))
  }
  return(seqs)
}

#' Correct allele calls based on a personalized genotype
#'
#' \code{reassignAlleles} uses a subject-specific genotype to correct
#' correct preliminary allele assignments of a set of sequences derived
#' from a single subject.
#' 
#' @details In order to save time, initial gene assignments are preserved and
#' the allele calls are chosen from among those provided in \code{genotype_db},
#' based on a simple alignment to the sample sequence.
#' 
#' @param    clip_db       a \code{data.frame} containing V allele calls from a
#'                         single subject under \code{"V_CALL"} and the sample
#'                         IMGT-gapped V(D)J sequences under
#'                         \code{"SEQUENCE_IMGT"}
#' @param    genotype_db   a vector of named nucleotide germline sequences
#'                         matching the calls detailed in \code{allele_calls}
#'                         and personalized to the subject
#' @param    method        the method to be used when realigning sequences to
#'                         the genotype_db sequences. Currently only "hammming"
#'                         (for Hamming distance) is implemented.
#' @param    path          directory containing the tool used in the
#'                         realignment method, if needed. Hamming distance does
#'                         not require a path to a tool.
#' @param    keep_gene     logical indicating if gene assignments should be
#'                         maintained when possible. Increases speed by
#'                         minimizing required number of alignments. Currently
#'                         only "TRUE" is implemented.
#' 
#' @return   a single-column \code{data.frame} corresponding to \code{clip.db}
#'           and containing the best allele call from among the sequences
#'           listed in \code{genotype_db}
#' 
#' @examples
#' # Load example data
#' data(germline_ighv)
#' data(sample_db)
#' data(genotype)
#' data(novel_df)
#'                      
#' # Extract the database sequences that correspond to the genotype
#' genotype_seqs = genotypeFasta(genotype, germline_ighv, novel_df)
#' 
#' # Use the personlized genotype to determine corrected allele assignments
#' V_CALL_GENOTYPED = reassignAlleles(sample_db, genotype_seqs)
#' sample_db = cbind(sample_db, V_CALL_GENOTYPED)
#' 
#' @export
reassignAlleles <- function(clip_db, genotype_db, method="hamming", path=NA,
                            keep_gene=TRUE){
  
  # Extract data subset and prepare output vector
  v_sequences = as.character(clip_db$SEQUENCE_IMGT)
  v_calls = getAllele(clip_db$V_CALL, first=FALSE, strip_d=FALSE)
  v_genes = getGene(v_calls, first = TRUE, strip_d=FALSE)
  V_CALL_GENOTYPED = rep("", length(v_calls))
  
  
  if(keep_gene){
    # Find which genotype genes are homozygous and assign those alleles first
    geno_genes = getGene(names(genotype_db),strip_d=TRUE)
    names(geno_genes) = names(genotype_db)
    hetero_genes = unique(geno_genes[which(duplicated(geno_genes))])
    homo_genes = geno_genes[!(geno_genes %in% hetero_genes)]
    homo_alleles = names(homo_genes); names(homo_alleles) = homo_genes
    homo_calls_i = which(v_genes %in% homo_genes)
    V_CALL_GENOTYPED[homo_calls_i] = homo_alleles[v_genes[homo_calls_i]]
    
    # Now realign the heterozygote sequences to each allele of that gene
    for (het_gene in hetero_genes){
      ind = which(v_genes %in% het_gene)
      if (length(ind) > 0){
        het_alleles = names(geno_genes[which(geno_genes == het_gene)])
        het_seqs = genotype_db[het_alleles]
        if(method == "hamming"){
        dists = lapply(het_seqs, function(x)
          sapply(getMutatedPositions(v_sequences[ind], x, match_instead=FALSE),
                 length))
        dist_mat = matrix(unlist(dists), ncol = length(het_seqs))
        } else {
          stop("Only Hamming distance is currently supported as a method.")
        }
        best_match = apply(dist_mat, 1, function(x) which(x == min(x)))
        best_alleles = sapply(best_match, function(x) het_alleles[x])   
        V_CALL_GENOTYPED[ind] = sapply(best_alleles, paste, collapse=",")
      }
    }
    
    # Now realign the gene-not-in-genotype calls to every genotype allele
    hetero_calls_i = which(v_genes %in% hetero_genes)
    not_called = setdiff(1:length(v_genes), c(homo_calls_i, hetero_calls_i))
    if(length(not_called)>1){
      if(method ==  "hamming"){
      dists = lapply(genotype_db, function(x)
        sapply(getMutatedPositions(v_sequences[not_called], x, match_instead=FALSE),
               length))
      dist_mat = matrix(unlist(dists), ncol = length(genotype_db))
      } else {
        stop("Only Hamming distance is currently supported as a method.")
      }
      best_match = apply(dist_mat, 1, function(x) which(x == min(x)))
      best_alleles = sapply(best_match, function(x) names(genotype_db[x])) 
      V_CALL_GENOTYPED[not_called] = sapply(best_alleles, paste, collapse=",")
    }
  } else {
    stop("Complete realignment is currently not supported.")
  }
  
  return(data.frame(V_CALL_GENOTYPED,stringsAsFactors=FALSE))
}


# Other Mutation-Related Functions ----------------------------------------

#' Find the location of mutations in a sequence
#'
#' \code{getMutatedPositions} takes two vectors of aligned sequences and
#' compares pairs of sequences. It returns a list of the nucleotide positions of
#' any differences.
#' 
#' @param    samples        a vector of strings respresenting aligned sequences
#' @param    germlines      a vector of strings respresenting aligned sequences
#'                          to which \code{samples} will be compared. If only
#'                          one string is submitted, it will be used for all
#'                          \code{samples}.
#' @param    ignored_regex  a regular expression indicating what characters
#'                          should be ignored (such as gaps and N nucleotides).
#' @param    match_instead  if \code{TRUE}, the function returns the positions
#'                          that are the same instead of those that are
#'                          different.
#' @return   A list of the nucleotide positions of any differences between the
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
  
  muts = mapply(function(x, y) x[!x%in%y], muts, ignore, SIMPLIFY=FALSE)
  return(muts)
}


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
#' #sample_seqs = c(germline_ighv[2],
#' #                createGermlines(germline_ighv[1], 103, "G"),
#' #                createGermlines(germline_ighv[1], 107, "C"))
#' 
#' # Pretend that one sample sequence has received an ambiguous allele call
#' #sample_alleles = c(paste(names(germline_ighv[1:2]), collapse=","),
#' #                  names(germline_ighv[2]),
#' #                  names(germline_ighv[1]))
#' 
#' # Compare each sequence to its assigned germline(s) to determine the distance
#' #getMutCount(sample_seqs, sample_alleles, germline_ighv)
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

#' Determine which calls represent an unmutated allele
#'
#' \code{findUnmutatedCalls} determines which allele calls would represent a 
#' perfect match with the germline sequence, given a vector of allele calls and
#' mutation counts. In the case of multiple alleles being assigned to a
#' sequence, only the subset that would represent a perfect match is returned.
#' 
#' @param    allele_calls   a vector of strings respresenting Ig allele calls,
#'                          where multiple calls are separated by a comma
#' @param    germline_db    a vector of named nucleotide germline sequences
#' @param    sample_seqs    V(D)J-rearranged sample sequences matching the order
#'                          of the given \code{allele_calls}
#' 
#' @return   A vector of strings containing the members of \code{allele_calls}
#'           that represent unmutated sequences
#' 
#' @examples
#' # Load data
#' data(germline_ighv)
#' data(sample_db)
#'
#' # Find which of the sample alleles are unmutated
#' findUnmutatedCalls(sample_db$V_CALL, sample_db$SEQUENCE_IMGT, germline_ighv)
#' 
#' @export
findUnmutatedCalls <- function(allele_calls, sample_seqs, germline_db){
  . = NULL
  allele_calls = getAllele(allele_calls, first = FALSE)
  sample_seqs = as.character(sample_seqs)
  
  # Remove calls not in germline_db
  not_in_db = allele_calls %>%
    strsplit(",") %>%
    unlist %>%
    setdiff(names(germline_db))
  no_call = which(allele_calls == "") 
  in_db = not_in_db %>%
    sapply(grep, allele_calls, fixed=TRUE) %>%
    unlist() %>%
    c(no_call) %>%
    unique() %>%
    setdiff(1:length(allele_calls), .)
  allele_calls = allele_calls[in_db]
  sample_seqs = sample_seqs[in_db]
  
  mut_counts = getMutCount(sample_seqs, allele_calls, germline_db)
  
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
  
  unmut_alleles = unmut_alleles[unmut_i]
  
  return(unmut_alleles)
  
}

#' Find Frequent Sequences' Mutation Counts
#'
#' \code{getPopularMutationCount} determines which sequences occur frequently
#' for each V gene and returns the mutation count of those sequences.
#' 
#' @param  sample_db     A Change-O db data frame. See
#'                       \link{findNovelAlleles} for a list of required
#'                       columns.
#' @param  germline_db   A named list of IMGT-gapped germline sequences.
#' @param  gene_min      The portion of all unique sequences a gene must
#'                       constitute to avoid exclusion.
#' @param  seq_min       The number of copies of the V that must be present for
#'                       to avoid exclusion.
#' @param  seq_p_of_max  For each gene, fraction of the most common V sequence's
#'                       count that a sequence must meet to avoid exclusion.
#' @param  full_return   If true, will return all \code{sample_db} columns and
#'                       will include sequences with mutation count < 1.
#' 
#' @return  A data frame of genes that have a frequent sequence mutation count
#'          above 1.
#' 
#' @seealso \link{getMutatedPositions} can be used to find which positions
#'          of a set of sequences are mutated.
#' 
#' @examples
#' data(sample_db, germline_ighv)
#' getPopularMutationCount(sample_db, germline_ighv)
#' 
#' @export
getPopularMutationCount <- function(sample_db, germline_db, gene_min = 1e-03,
                                    seq_min = 50, seq_p_of_max = 1/8,
                                    full_return = FALSE){
  modified_db = sample_db %>%
    mutate_(V_GENE = ~getGene(V_CALL)) %>%
    group_by_(~1:n()) %>%
    mutate_(V_SEQUENCE_IMGT = ~substring(SEQUENCE_IMGT, 1, 312)) %>%
    # Count occurence of each unique IMGT-gapped V sequence
    group_by_(~V_GENE, ~V_SEQUENCE_IMGT) %>%
    mutate_(V_SEQUENCE_IMGT_N = ~n()) %>%
    # Count occurence of each gene and determine count of most common sequence
    group_by_(~V_GENE) %>%
    mutate_(V_GENE_N = ~n()) %>%
    mutate_(V_SEQUENCE_IMGT_N_MAX = ~max(V_SEQUENCE_IMGT_N)) %>%
    # Remove rare V genes, rare sequences, and sequences not making up a
    # sufficient proportion of sequences as compared to the most common
    ungroup %>%
    distinct_(~V_SEQUENCE_IMGT) %>%
    filter_(~V_GENE_N >= (nrow(sample_db)*gene_min)) %>%
    filter_(~V_SEQUENCE_IMGT_N >= seq_min) %>%
    mutate_(V_SEQUENCE_IMGT_P_MAX = ~V_SEQUENCE_IMGT_N/V_SEQUENCE_IMGT_N_MAX) %>%
    filter_(~V_SEQUENCE_IMGT_P_MAX >= seq_p_of_max)
  # Determine the mutation counts of the V sequences and append them to the db
  MUTATION_COUNT = getMutCount(modified_db$V_SEQUENCE_IMGT,
                               modified_db$V_CALL,
                               germline_db) %>% 
    sapply(function(x) min(unlist(x)))
  merged_db = bind_cols(modified_db, data.frame(MUTATION_COUNT))
  # Strip down the data frame before returning it
  if (!full_return) {
    merged_db = merged_db %>%
      filter_(~MUTATION_COUNT > 0) %>%
      select_(~V_GENE, ~MUTATION_COUNT)
  }
  return(merged_db)
}

#' Insert polymorphisms into a nucleotide sequence
#'
#' \code{insertPolymorphisms} replaces nucleotides in the desired locations of a
#' provided sequence.
#' 
#' 
#' @param    sequence     the starting nucletide sequence
#' @param    positions    a vector of positions which to be changed
#' @param    nucleotides  a vector of nucletides to which to change the
#'                        positions
#' @return   a sequence with the desired nucleotides in provided locations
#' 
#' @examples
#' insertPolymorphisms("hugged", c(1,6,2), c("t","r","i")) 
#' 
#' @export
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

# Formatting and Cleanup --------------------------------------------------

#' Read immunoglobulin sequences
#'
#' \code{readIgFasta} reads a fasta-formatted file of immunoglobulin (Ig)
#' sequences and returns a named vector of those sequences.
#' 
#' @param    fasta_file       fasta-formatted file of immunoglobuling sequences
#' @param    strip_down_name  if \code{TRUE}, will extract only the allele name
#'                            from the strings fasta file's sequence names
#' @param    force_caps       if \code{TRUE}, will force nucleotides to
#'                            uppercase
#' @return   a named vector of strings respresenting Ig alleles
#' 
#' @seealso  \link{writeFasta} to do the inverse.
#' 
#' @export
readIgFasta <- function(fasta_file, 
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
  if(strip_down_name){ seq_names = getAllele(seq_names, strip_d=FALSE) }
  names(seqs) = seq_names
  return(seqs[which(!is.na(seqs))])
}

#' Write to a fasta file
#'
#' \code{writeFasta} writes a named vector of sequences to a file in fasta
#' format.
#' 
#' @param    named_sequences  a vector of named string representing sequences
#' @param    file             the name of the output file
#' @param    width            the number of characters to be printed per line.
#'                            If not between 1 and 255, width with be infinite.
#' @param    append           \code{logical} indicating if the output should be
#'                            appended to \code{file} instead of overwriting it
#' 
#' @return   a named vector of strings respresenting Ig alleles
#' 
#' @seealso  \link{readIgFasta} to do the inverse.
#' 
#' @export
writeFasta <- function(named_sequences, file, width=60, append=FALSE){
  . = NULL
  seq_names = names(named_sequences) %>%
    paste(">", ., "\n", sep="")
  seqs = as.character(named_sequences)
  if(is.numeric(width) & width > 0 & width < 256){
    width_regex = paste("(.{", width, ",", width, "})", sep="")
    seqs = gsub(width_regex, "\\1\n", seqs)
  }
  seqs = seqs %>%
    paste("\n", sep="") %>%
    gsub("\n\n", "\n", .)
  paste(seq_names, seqs, sep="", collapse="") %>%
    cat(file=file, append=append)
} 

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
#' @return   vector of strings respresenting updated IGHV allele names
#' 
#' @references Xochelli et al. (2014) Immunoglobulin heavy variable (IGHV) genes
#' and alleles: new entities, new names and implications for research and
#' prognostication in chronic lymphocytic leukaemia. Immunogenetics. 67(1):61-6
#' 
#' @seealso Like \code{updateAlleleNames}, \link{sortAlleles} can help
#'          format a list of allele names.
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
  . = NULL
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
    allele_calls = allele_calls %>%
      gsub(temporary_names[i], definitive_names[i], ., fixed = TRUE)
  }
  return(allele_calls)
}

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
#' @param    method        a string defining the method to use when sorting alleles.
#'                         If \code{"name"} then sort in lexicographic order. If
#'                         \code{"position"} then sort by position in the locus, as
#'                         determined by the final two numbers in the gene name.
#' @return   A sorted vector of strings respresenting Ig allele names
#' 
#' @seealso Like \code{sortAlleles}, \link{updateAlleleNames} can help
#'          format a list of allele names.
#' 
#' @examples
#' # Create a list of allele names
#' alleles = c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
#'             "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
#'             "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")
#' 
#' # Sort the alleles by name
#' sortAlleles(alleles)
#' 
#' # Sort the alleles by position in the locus
#' sortAlleles(alleles, method="pos")
#' 
#' @export
sortAlleles <- function(allele_calls, method=c("name", "position")) { 
  # Check arguments
  method <- match.arg(method)
    
  # Standardize format of submitted alleles, first
  SUBMITTED_CALLS = getAllele(allele_calls, first = FALSE, strip_d= FALSE) %>%
    sort()
  allele_df = data.frame(SUBMITTED_CALLS,stringsAsFactors = FALSE) %>%
    # Determine the family
    mutate_(FAMILY = ~getFamily(SUBMITTED_CALLS)) %>%
    # Determine the gene (exclude family); convert letters to numbers for sort
    mutate_(GENE = ~getGene(SUBMITTED_CALLS)) %>%
    mutate_(GENE1 = ~gsub("[^-]+-([^-\\*D]+).*","\\1",SUBMITTED_CALLS)) %>%
    mutate_(GENE1 = ~as.numeric(gsub("[^0-9]+", "99", GENE1))) %>%
    # If there is a second gene number, determine that, too
    mutate_(GENE2 = ~gsub("[^-]+-[^-]+-?","",GENE)) %>%
    mutate_(GENE2 = ~as.numeric(gsub("[^0-9]+", "99", GENE2))) %>%
    mutate_(ALLELE = ~as.numeric(sub("[^\\*]+\\*|[^\\*]+$","",
                                   getAllele(SUBMITTED_CALLS))))
  # Convert missing values to 0, sort data frame
  allele_df[is.na(allele_df)] = 0
  if (method == "name") {  
    sorted_df = arrange_(allele_df, ~FAMILY, ~GENE1, ~GENE2, ~ALLELE)
  } else if (method == "position") {
    sorted_df = arrange_(allele_df, ~desc(GENE1), ~desc(GENE2), ~FAMILY, ~ALLELE)
  }
  
  return(sorted_df$SUBMITTED_CALLS)
}

#' Clean up nucleotide sequences
#'
#' \code{cleanSeqs} capitalizes nucleotides, replaces "." with "-", and then
#' replaces all characters besides ACGT- with "N". 
#' 
#' @param    seqs  a vector of nucleotide sequences
#' @return   A vector of nucleotide sequences
#' 
#' @seealso \link{sortAlleles} and \link{updateAlleleNames} can
#'          help format a list of allele names.
#' 
#' @examples
#' # Create messy nucleotide sequences
#' seqs = c("AGAT.taa-GAG...ATA",
#'          "GATACAGTXXXXXAGNNNPPPACA")
#' # Clean them up
#' cleanSeqs(seqs)
#' 
#' @export
cleanSeqs <- function(seqs){
  . = NULL
  seqs %>%
    toupper %>%
    gsub(".", "-", . , fixed = TRUE) %>%
    gsub("[^ACGT-]", "N", .) %>%
    return
}


# Private Functions -------------------------------------------------------

# Find muations-by-position compared to a germline
#
# \code{positionMutations} duplicates the rows of a data frame for each
# position to be analyzed and determines if each sample is mutated at that
# position
# 
# @param  clip_db       A Change-O db data frame. See
#                       \link{findNovelAlleles} for a list of required
#                       columns.
# @param  germline      The germline to which all the sequences should be
#                       compared
# @param  pos_range     The range of positions within the sequence for which
#                       the rows should be duplicated and checked for mutation
# 
# @return  A data frame with rows duplicated for all the positions to be
# analyzed and a column indicating whether the position is mutated in
# comparison to the germline
#
positionMutations <- function(clip_db, germline, pos_range){
  . = NULL
  pos_db = pos_range %>%
    length() %>%
    rep("clip_db", .) %>%
    paste(collapse=",") %>%
    paste("bind_rows(",., ")") %>%
    parse(text=.) %>%
    eval()
  pos_db$POSITION = c(sapply(pos_range, rep, nrow(clip_db)))
  # Find which positions are mutated
  pos_db = pos_db %>%
    mutate_(NT = ~substring(SEQUENCE_IMGT, POSITION, POSITION)) %>%
    mutate_(GERM_NT = ~substring(germline, POSITION, POSITION)) %>%
    mutate_(MUTATED = ~(NT != GERM_NT & NT != "N" & NT != "-" & NT != "")) %>%
    mutate_(OBSERVED = ~(NT != "-" & NT != ""))
  return(pos_db)
}

# Find sequences carrying certain levels of mutation
#
# \code{mutationRangeSubset} determines the mutations in a \code{data.frame} of
# sequences and returns the subset of sequences that meet the given mutation
# count limits
# 
# @param  clip_db       A Change-O db data frame. See
#                       \link{findNovelAlleles} for a list of required
#                       columns.
# @param  germline      The germline to which all the sequences should be
#                       compared
# @param  pos_range     The range of positions within the sequences that should
#                       be analyzed for mutations
# @param  pos_range     The range of mutation counts that sequences can have
#                       and still be included
#
# @return  A data frame containing only the subset carrying the desired levels
# of mutation
#
mutationRangeSubset <- function(clip_db, germline, mut_range, pos_range){
  . = NULL
  pads = paste(rep("-", min(pos_range)-1), collapse="")
  clip_db$MUT_COUNT = clip_db$SEQUENCE_IMGT %>%
    substring(min(pos_range), max(pos_range)) %>%
    paste(pads, ., sep="") %>%
    getMutatedPositions(germline) %>%
    sapply(length)
  clip_db = clip_db %>%
    filter_(~MUT_COUNT %in% mut_range)
  return(clip_db)
}

# Find lower range of y-intercept confidence interval
#
# \code{findLowerY} finds the lower range of y-intercept confidence interval
# 
# @param    x         A vector of x values
# @param    y         A vector of y values
# @param    mut_min   The value where the the lowest mutation count should be
#                     found. See details.
# @param    alpha     The alpha cutoff the be used in constructing the
#                     confidence interval
#
# @details  If mut_min is 1, a y-intercept will be searched for at 0. If
# mut_min is above 1, then the "y-intercept" will be found at x = mut_min - 1.
#
# @return  A data frame containing only the subset carrying the desired levels
# of mutation
#
findLowerY = function(x, y, mut_min, alpha){
  y = y+1-mut_min
  lowerY = suppressWarnings(confint(lm(x ~ y),level=1-2*alpha)[[1]])
  return(lowerY)
}

# Enchanced substring extraction
#
# \code{superSubstring} is an enahnced version of \code{substring} in that
# it can find disjoint positions in one call.
#
# @param    string      a single string
# @param    positions   the positions to be extracted
#
# @return  a substring
# 
superSubstring = function(string, positions){
  if(length(string) != 1){ stop("Please submit only one string.") }
  chars = sapply(positions, function(x) substring(string, x, x))
  return(paste(chars, collapse=""))
}


# Layout multiple ggplots
#
# \code{multiplot} is a function provided by http://www.cookbook-r.com/ which
# allows for plotting multiple ggplot objects in one panel.
# 
# @param    ...       ggplot2 object(s)
# @param    plotlist  a list alternative to ...
# @param    file      an unused parameter, but present in the provided function
# @param    cols      Number of columns in layout
# @param    layout    A matrix specifying the layout. If present, 'cols' is
#                     ignored.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
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

