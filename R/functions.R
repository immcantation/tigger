
# The TIgGER Trifecta -----------------------------------------------------

#' Find novel alleles from repertoire sequencing data
#'
#' \code{findNovelAlleles} analyzes mutation patterns in sequences thought to
#' align to each germline allele in order to determine which positions
#' might be polymorphic.
#' 
#' The TIgGER allele-finding algorithm, briefly, works as follows:
#' Mutations are determined through comparison to the provided germline.
#' Mutation frequency at each *position* is determined as a function of
#' *sequence-wide* mutation counts. Polymorphic positions exhibit a high
#' mutation frequency despite sequence-wide mutation count. False positive of
#' potential novel alleles resulting from clonally-related sequences are guarded
#' against by ensuring that sequences perfectly matching the potential novel
#' allele utilize a wide range of combinations of J gene and junction length.
#' 
#' @param    data           a \code{data.frame} in Change-O format. See details.
#' @param    germline_db    a vector of named nucleotide germline sequences
#'                          matching the V calls in \code{data}.
#' @param    v_call         name of the column in \code{data} with V allele calls. 
#'                          Default is V_CALL.                                                    
#' @param    germline_min   the minimum number of sequences that must have a
#'                          particular germline allele call for the allele to
#'                          be analyzed
#' @param    auto_mutrange  if \code{TRUE}, the algorithm will attempt to
#'                          determine the appropriate mutation range
#'                          automatically using the mutation count of the most
#'                          common sequence assigned to each allele analyzed
#' @param    mut_range      the range of mutations that samples may carry and
#'                          be considered by the algorithm
#' @param    pos_range      the range of IMGT-numbered positions that should be
#'                          considered by the algorithm
#' @param    alpha          the alpha value used for determining whether the 
#'                          fit y-intercept is greater than the \code{y_intercept}
#'                          threshold
#' @param    y_intercept    the y-intercept threshold above which positions should be
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
#' @param    nproc          the number of processors to use
#'
#' @return
#' A \code{data.frame} with a row for each known allele analyzed.
#' Besides metadata on the the parameters used in the search, each row will have
#' either a note as to where the polymorphism-finding algorithm exited or a
#' nucleotide sequence for the predicted novel allele, along with columns providing
#' additional evidence.
#' 
#' The output contains the following columns:
#' \itemize{
#'   \item \code{GERMLINE_CALL}: The input (uncorrected) V call.
#'   \item \code{NOTE}: Comments regarding the inferrence.
#'   \item \code{POLYMORPHISM_CALL}: The novel allele call.
#'   \item \code{NT_SUBSTITUTIONS}: Mutations identified in the novel allele, relative
#'         to the reference germline (\code{GERMLINE_CALL})
#'   \item \code{NOVEL_IMGT}: The novel allele sequence.
#'   \item \code{NOVEL_IMGT_COUNT}:  The number of times the sequence \code{NOVEL_IMGT} 
#'         is found in the input data. Considers the subsequence of \code{NOVEL_IMGT} 
#'         in the \code{pos_range}.
#'   \item \code{NOVEL_IMGT_UNIQUE_J}: Number of distinct J calls associated to \code{NOVEL_IMGT} 
#'         in the input data. Considers the subsequence of \code{NOVEL_IMGT} in the \code{pos_range}.       
#'   \item \code{NOVEL_IMGT_UNIQUE_CDR3}: Number of distinct CDR3 sequences associated
#'         with \code{NOVEL_IMGT} in the input data. Considers the subsequence of \code{NOVEL_IMGT} 
#'         in the \code{pos_range}.                                              
#'   \item \code{PERFECT_MATCH_COUNT}: Final number of sequences retained to call the new 
#'         allele. These are unique sequences that have V segments that perfectly match 
#'         the predicted germline in the \code{pos_range}.
#'   \item \code{PERFECT_MATCH_FREQ}: \code{PERFECT_MATCH_COUNT / GERMLINE_CALL_COUNT}
#'   \item \code{GERMLINE_CALL_COUNT}: The number of sequences with the \code{GERMLINE_CALL} 
#'         in the input data that were initially considered for the analysis.
#'   \item \code{GERMLINE_CALL_FREQ}: The fraction of sequences with the \code{GERMLINE_CALL} 
#'         in the input data initially considered for the analysis.              
#'   \item \code{GERMLINE_IMGT}: Germline sequence for \code{GERMLINE_CALL}.
#'   \item \code{GERMLINE_IMGT_COUNT}:  The number of times the \code{GERMLINE_IMGT} 
#'         sequence is found in the input data.
#'   \item \code{MUT_MIN}: Minimum mutation considered by the algorithm.
#'   \item \code{MUT_MAX}: Maximum mutation considered by the algorithm.
#'   \item \code{MUT_PASS_COUNT}: Number of sequences in the mutation range.
#'   \item \code{POS_MIN}: First position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{POS_MAX}: Last position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{Y_INTERCEPT}: The y-intercept above which positions were considered 
#'         potentially polymorphic.
#'   \item \code{Y_INTERCEPT_PASS}: Number of positions that pass the \code{Y_INTERCEPT} threshold.
#'   \item \code{SNP_PASS}: Number of sequences that pass the \code{Y_INTERCEPT} threshold and are
#'         within the desired nucleotide range (\code{min_seqs}).
#'   \item \code{UNMUTATED_COUNT}: Number of unmutated sequences.
#'   \item \code{UNMUTATED_FREQ}: Number of unmutated sequences over \code{GERMLINE_IMGT_COUNT}.
#'   \item \code{UNMUTATED_SNP_J_GENE_LENGTH_COUNT}: Number of distinct combinations
#'         of SNP, J gene, and junction length.     
#'   \item \code{SNP_MIN_SEQS_J_MAX_PASS}: Number of SNPs that pass both the \code{min_seqs} 
#'         and \code{j_max} thresholds.
#'   \item \code{ALPHA}: Significance threshold to be used when constructing the 
#'         confidence interval for the y-intercept.
#'   \item \code{MIN_SEQS}: Input \code{min_seqs}. The minimum number of total sequences 
#'         (within the desired mutational range and nucleotide range) required 
#'         for the samples to be considered.
#'   \item \code{J_MAX}: Input \code{j_max}. The maximum fraction of sequences perfectly 
#'         aligning to a potential novel allele that are allowed to utilize to a particular 
#'         combination of junction length and J gene.
#'   \item \code{MIN_FRAC}: Input \code{min_frac}. The minimum fraction of sequences that must
#'         have usable nucleotides in a given position for that position to be considered.
#' }
#' 
#' The following comments can appear in the \code{NOTE} column:
#' 
#' \itemize{
#'   \item \emph{Novel allele found}: A novel allele was detected.
#'   \item \emph{Plurality sequence too rare}: No sequence is frequent enough to pass 
#'         the J test (\code{j_max}).
#'   \item \emph{A J-junction combination is too prevalent}: Not enough J diversity (\code{j_max}).
#'   \item \emph{No positions pass y-intercept test}: No positions above \code{y_intercept}.
#'   \item \emph{Insufficient sequences in desired mutational range}: 
#'         \code{mut_range} and \code{pos_range}.
#'   \item \emph{Not enough sequences}: Not enough sequences in the desired mutational 
#'         range and nucleotide range (\code{min_seqs}).
#'   \item \emph{No unmutated versions of novel allele found}: All observed variants of the 
#'         allele are mutated.
#' }
#' 
#' @seealso \link{plotNovel} to visualize the data supporting any
#' novel alleles hypothesized to be present in the data and
#' \link{inferGenotype} to determine if the novel alleles are frequent
#' enought to be included in the subject's genotype.
#' 
#' @examples
#' \donttest{
#' # Find novel alleles and return relevant data
#' novel <- findNovelAlleles(SampleDb, GermlineIGHV)
#' }
#' 
#' @export
findNovelAlleles <- function(data, germline_db,
                             v_call="V_CALL",
                             germline_min=200,
                             min_seqs=50,
                             auto_mutrange=TRUE,
                             mut_range=1:10,
                             pos_range=1:312,
                             y_intercept=0.125,
                             alpha=0.05,
                             j_max=0.15,
                             min_frac=0.75,
                             nproc=1) {
    . = idx = NULL
    
    # Keep only the db columns needed
    data <- data %>% 
        dplyr::select_('SEQUENCE_IMGT', v_call, 'J_CALL', 'JUNCTION_LENGTH', 'JUNCTION')
    
    # Keep only the columns we need and clean up the sequences
    missing = c("SEQUENCE_IMGT", v_call, "J_CALL", "JUNCTION_LENGTH") %>%
        setdiff(colnames(data))
    if (length(missing) != 0) {
        stop("Could not find required columns in the input data:\n  ",
             paste(missing, collapse="\n  "))
    }
    empty_junctions = sum(data$JUNCTION_LENGTH == 0, na.rm=TRUE)
    if (empty_junctions > 0) {
        stop(empty_junctions, " sequences have junction ", "length of zero. ",
             "Please remove these sequences.")
    }
    germlines = cleanSeqs(germline_db)
    names(germlines) = getAllele(names(germlines), first=FALSE, strip_d=FALSE)
    data$SEQUENCE_IMGT = cleanSeqs(data$SEQUENCE_IMGT)
    
    
    # Find which rows' calls contain which germline alleles
    cutoff =
        ifelse(germline_min < 1, round(nrow(data)*germline_min), germline_min)
    allele_groups = sapply(names(germlines), grep, data[[v_call]], fixed=TRUE,
                           simplify=FALSE)
    names(allele_groups) = names(germlines)
    allele_groups = allele_groups[sapply(allele_groups, length) >= cutoff]
    if(length(allele_groups) == 0){
        stop_message <- paste("Not enough sample sequences were assigned to any germline:\n",
                              " (1) germline_min is too large or\n",
                              " (2) sequences names don't match germlines.")
        stop(stop_message)
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
        foreach::registerDoSEQ()
    } else {
        #cluster_type = ifelse(Sys.info()['sysname'] == "Windows", "PSOCK", "FORK")
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        parallel::clusterExport(cluster, list("allele_groups",
                                              "germlines",
                                              "data",
                                              "min_seqs",
                                              "auto_mutrange",
                                              "mut_range",
                                              "pos_range",
                                              "y_intercept",
                                              "alpha",
                                              "j_max",
                                              "germline_min",
                                              "min_frac",
                                              "findLowerY",
                                              "mutationRangeSubset",
                                              "positionMutations",
                                              "superSubstring"),
                                envir=environment())
        doParallel::registerDoParallel(cluster)
    }
    
    out_list <- foreach(idx=iterators::icount(length(allele_groups))) %dopar% {
        # out_list <- lapply(1:length(allele_groups), function(idx) {  
        gc() 
        # message(paste0("idx=",idx))
        # Subset of data being analyzed
        allele_name = names(allele_groups)[idx]
        germline = germlines[allele_name]
        indicies = allele_groups[[allele_name]]
        db_subset = data[indicies, ]
        
        # If mutrange is auto, find most popular mutation count and start from there
        gpm = db_subset %>%
            dplyr::mutate_(V_CALL = ~allele_name) %>%
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
                                  NT_SUBSTITUTIONS=NA,
                                  NOVEL_IMGT = NA,
                                  NOVEL_IMGT_COUNT=NA,
                                  NOVEL_IMGT_UNIQUE_J=NA,
                                  NOVEL_IMGT_UNIQUE_CDR3=NA,
                                  PERFECT_MATCH_COUNT = NA,
                                  PERFECT_MATCH_FREQ = NA,                              
                                  GERMLINE_CALL_COUNT = length(indicies),
                                  GERMLINE_CALL_FREQ = round(length(indicies)/nrow(data), 3),
                                  MUT_MIN = NA,
                                  MUT_MAX = NA,
                                  MUT_PASS_COUNT=NA,
                                  GERMLINE_IMGT = as.character(germline),
                                  GERMLINE_IMGT_COUNT=NA,
                                  POS_MIN = min(pos_range),
                                  POS_MAX = max(pos_range),
                                  Y_INTERCEPT = y_intercept,
                                  Y_INTERCEPT_PASS = NA,
                                  SNP_PASS=NA,
                                  UNMUTATED_COUNT=NA,
                                  UNMUTATED_FREQ=NA,
                                  UNMUTATED_SNP_J_GENE_LENGTH_COUNT=NA,
                                  SNP_MIN_SEQS_J_MAX_PASS=NA,
                                  ALPHA = alpha,
                                  MIN_SEQS = min_seqs,
                                  J_MAX = j_max,
                                  MIN_FRAC = min_frac,
                                  stringsAsFactors = FALSE)
        for (mut_min in rev(mut_mins)) {
            gc()
            # message(paste0("|-- mut_min=",mut_min))
            if (mut_min == rev(mut_mins)[1]){
                df_run = df_run_empty
            } else {
                df_run = dplyr::bind_rows(df_run_empty, df_run)
            }
            mut_max = mut_min + diff(range(mut_range))
            df_run$MUT_MIN[1] = mut_min
            df_run$MUT_MAX[1] = mut_max
            
            # If no sequence is frequent enough to pass the J test, give up now
            if(nrow(gpm) < 1) {
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
            df_run$MUT_PASS_COUNT[1] <- nrow(db_subset_mm)
            
            if(nrow(db_subset_mm) < min_seqs){
                df_run$NOTE[1] = paste0("Insufficient sequences (",nrow(db_subset_mm),") in desired mutational range.")
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
                dplyr::group_by_(~POSITION) %>%
                dplyr::mutate_(PASS = ~mean(OBSERVED) >= min_frac) %>%
                dplyr::group_by_(~MUT_COUNT, ~POSITION) %>%
                dplyr::summarise_(POS_MUT_RATE = ~ mean(MUTATED)*unique(PASS) ) %>% 
                dplyr::ungroup()   
            
            rm(pos_db)
            gc()
            
            # Calculate y intercepts, find which pass the test
            pass_y = pos_muts %>%
                dplyr::group_by_(~POSITION) %>%
                dplyr::summarise_(Y_INT_MIN = ~findLowerY(POS_MUT_RATE, MUT_COUNT,
                                                          mut_min, alpha)) %>%
                dplyr::filter_(~Y_INT_MIN > y_intercept)
            
            df_run$Y_INTERCEPT_PASS[1] <- nrow(pass_y)
            
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
                dplyr::group_by(1:n()) %>%
                dplyr::mutate_(SNP_STRING = ~superSubstring(SEQUENCE_IMGT,
                                                            pass_y$POSITION)) %>%
                dplyr::filter_(~SNP_STRING != gl_substring) %>%
                dplyr::group_by_(~SNP_STRING) %>%
                dplyr::mutate_(STRING_COUNT = ~n()) %>%
                dplyr::filter_(~STRING_COUNT >= min_seqs)
            
            df_run$SNP_PASS[1] <- nrow(db_y_subset_mm)
            
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
                dplyr::filter_(~MUT_COUNT_MINUS_SUBSTRING == 0)
            
            df_run$UNMUTATED_COUNT[1] <- nrow(db_y_summary0)
            
            db_y_summary0 <- db_y_summary0 %>%
                dplyr::mutate_(J_GENE = ~getGene(J_CALL)) %>%
                dplyr::group_by_(~SNP_STRING, ~J_GENE, ~JUNCTION_LENGTH) %>%
                dplyr::summarise_(COUNT = ~n())
            
            df_run$UNMUTATED_SNP_J_GENE_LENGTH_COUNT[1] <- nrow(db_y_summary0)
            
            db_y_summary0 <- db_y_summary0 %>%
                dplyr::group_by_(~SNP_STRING) %>%
                dplyr::mutate_(FRACTION = ~COUNT/sum(COUNT)) %>%
                dplyr::summarise_(TOTAL_COUNT = ~sum(COUNT), MAX_FRAC = ~max(FRACTION))
            
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
            
            # db_y_summary = db_y_summary0 %>%
            #   filter_(~TOTAL_COUNT >= min_seqs & MAX_FRAC <= j_max)
            
            min_seqs_pass <- db_y_summary0$TOTAL_COUNT >= min_seqs
            j_max_pass <- db_y_summary0$MAX_FRAC <= j_max
            
            db_y_summary <- db_y_summary0[min_seqs_pass & j_max_pass, , drop=FALSE]
            
            df_run$SNP_MIN_SEQS_J_MAX_PASS[1] <- nrow(db_y_summary)
            
            if(nrow(db_y_summary) < 1){
                msg <- c(NA, NA)
                names(msg) <- c("j_max", "min_seqs")
                
                if (sum(min_seqs_pass) == 0) {
                    msg['min_seqs'] <- paste0("Not enough sequences (maximum total count is ",
                                              max(db_y_summary0$TOTAL_COUNT),
                                              ").")
                }
                
                if (sum(j_max_pass) == 0) {
                    msg['j_max'] <- paste0("A J-junction combination is too prevalent (",
                                           round(100*max(db_y_summary0$MAX_FRAC),1),"% of sequences).")
                }
                
                msg <- paste(na.omit(msg), collapse=" and ")
                df_run$NOTE[1] = paste("Position(s) passed y-intercept (",
                                       paste(pass_y$POSITION, collapse = ","),
                                       ") but ",
                                       msg,".", sep="")
                df_run$PERFECT_MATCH_COUNT[1] = max(db_y_summary0$TOTAL_COUNT)
                df_run$PERFECT_MATCH_FREQ[1] <- df_run$PERFECT_MATCH_COUNT[1]/df_run$GERMLINE_CALL_COUNT[1]
                if(mut_mins[1] == mut_min){
                    return(df_run)
                } else {
                    next
                }
            }
            
            germ_nts = unlist(strsplit(gl_substring,""))
            for (r in 1:nrow(db_y_summary)) {
                if (r > 1){
                    df_run = dplyr::bind_rows(df_run[1,], df_run)
                }
                # Create the new germline
                snp_nts = unlist(strsplit(db_y_summary$SNP_STRING[r],""))
                remain_mut = db_y_summary$SNP_STRING[r] %>%
                    getMutatedPositions(gl_substring) %>%
                    unlist() %>%
                    unique()
                germ = insertPolymorphisms(germline, pass_y$POSITION, snp_nts)
                is_known_allele <- germ == germlines
                if (sum(is_known_allele) == 0 ) {
                    names(germ) = mapply(paste, germ_nts[remain_mut],
                                         pass_y$POSITION[remain_mut],
                                         snp_nts[remain_mut], sep="") %>%
                        paste(collapse="_") %>%
                        paste(names(germline), ., sep="_")
                } else {
                    # If the match is with duplicated sequences in the reference germlines,
                    # use the first
                    known_allele_names <- sortAlleles(names(germlines)[is_known_allele],
                                                      method="position")
                    names(germ) = known_allele_names[1]
                }
                # Save the new germline to our data frame               
                df_run$POLYMORPHISM_CALL[1] = names(germ)
                df_run$NOVEL_IMGT[1] =  as.character(germ)
                df_run$PERFECT_MATCH_COUNT[1] = db_y_summary$TOTAL_COUNT[r]
                df_run$PERFECT_MATCH_FREQ[1] <- df_run$PERFECT_MATCH_COUNT[1]/df_run$GERMLINE_CALL_COUNT[1]
                df_run$NOTE[1] = "Novel allele found!"
            }
            
        } # end for each starting mutation counts
        return(df_run)
        
    } # end foreach allele
    
    if(nproc > 1) { stopCluster(cluster) }
    out_df <- dplyr::bind_rows(out_list)
    getMuSpec <- function(poly_call, germ_call) {
        sapply(1:length(poly_call), function(i){
            p <- gsub(germ_call[i], "", poly_call[i], fixed = T)
            p <- strsplit(p,"_")[[1]][-1]
            m <- gsub("([[:alpha:]])([[:digit:]]*)([[:alpha:]])", "\\2\\1>\\3", p)
            paste(m, collapse=",")
        })
    }
    
    # The number of records in the sequence dataset matching 
    # each exact NOVEL_IMGT sequence
    getDbMatch <- function(novel_imgt) {
        sapply(novel_imgt, function(n) {
            n <- substr(n, min(pos_range), max(pos_range))
            sum(grepl(gsub("[-\\.]","",n),
                      gsub("[-\\.]","",data$SEQUENCE_IMGT)))
        })
    }
    
    # The number of distinct J in the sequence dataset associated 
    # with the exact NOVEL_IMGT sequence
    getNumJ <- function(novel_imgt) {
        sapply(novel_imgt, function(n) {
            n <- substr(n, min(pos_range), max(pos_range))
            imgt_idx <- grepl(gsub("[-\\.]","",n),
                              gsub("[-\\.]","",data$SEQUENCE_IMGT))
            length(unique(getGene(data[['J_CALL']][imgt_idx])))
        })
    }
    
    
    # The number of distinct CDR3 in the sequence dataset associated 
    # with the exact NOVEL_IMGT sequence
    getNumCDR3 <- function(novel_imgt) {
        sapply(novel_imgt, function(n) {
            n <- substr(n, min(pos_range), max(pos_range))
            imgt_idx <- grepl(gsub("[-\\.]","",n),
                              gsub("[-\\.]","",data$SEQUENCE_IMGT))
            seq <- data[['JUNCTION']][imgt_idx]
            seq <- substr(seq, 4, stringi::stri_length(seq) - 3)
            length(unique(seq))
        })
    }
    
    idx <- which(!is.na(out_df$NOVEL_IMGT))
    if (length(idx)>0) {
        out_df$NT_SUBSTITUTIONS[idx] <- getMuSpec(out_df$POLYMORPHISM_CALL[idx],
                                                  out_df$GERMLINE_CALL[idx])
        out_df$NOVEL_IMGT_COUNT[idx] <- getDbMatch(out_df$NOVEL_IMGT[idx])
        out_df$NOVEL_IMGT_UNIQUE_J[idx] <- getNumJ(out_df$NOVEL_IMGT[idx])
        if ("JUNCTION" %in% colnames(data)) {
            out_df$NOVEL_IMGT_UNIQUE_CDR3[idx] <- getNumCDR3(out_df$NOVEL_IMGT[idx])
        }
    }
    out_df$GERMLINE_IMGT_COUNT <- getDbMatch(out_df$GERMLINE_IMGT)
    out_df$UNMUTATED_FREQ = out_df$UNMUTATED_COUNT/out_df$GERMLINE_CALL_COUNT
    rm(data)
    gc()
    
    return(out_df)
}

#' Select rows containing novel alleles
#' 
#' \code{selectNovel} takes the result from \link{findNovelAlleles} and
#' selects only the rows containing unique, novel alleles.
#' 
#' @details
#' If, for instance, subject has in his genome \code{IGHV1-2*02} and a novel 
#' allele equally close to \code{IGHV1-2*02} and \code{IGHV1-2*05}, the novel allele may be
#' detected by analyzing sequences that best align to either of these alleles.
#' If \code{keep_alleles} is \code{TRUE}, both polymorphic allele calls will
#' be retained. In the case that multiple mutation ranges are checked for the
#' same allele, only one mutation range will be kept in the output.
#' 
#' @param   novel           a \code{data.frame} of the type returned by
#'                          \link{findNovelAlleles}.
#' @param   keep_alleles    a \code{logical} indicating if different alleles
#'                          leading to the same novel sequence should be kept.
#'                          See Details.
#'                        
#' @return  A \code{data.frame} containing only unique, novel alleles (if any)
#' that were in the input.
#' 
#' @examples
#' novel <- selectNovel(SampleNovel)
#' 
#' @export
selectNovel <- function(novel, keep_alleles=FALSE) {
    # Remove non-novel rows
    novel = filter_(novel, ~!is.na(NOVEL_IMGT))
    
    if (keep_alleles) {
        novel = novel %>% 
            group_by_(~GERMLINE_CALL)
    }
    novel_set = novel %>%
        distinct_(~NOVEL_IMGT, .keep_all=TRUE) %>%
        ungroup()
    
    return(novel_set)
}

#' Visualize evidence of novel V alleles
#'
#' \code{plotNovel} is be used to visualize the evidence of any novel V
#' alleles found using \link{findNovelAlleles}. It can also be used to
#' visualize the results for alleles that did
#' 
#' @details
#' The first panel in the plot shows, for all sequences which align to a particular 
#' germline allele, the mutation frequency at each postion along the aligned 
#' sequece as a function of the sequence-wide mutation. Sequences that pass 
#' the novel allele test are colored red, while sequences that don't pass
#' the test are colored yellow. The second panel shows the nucleotide usage at the 
#' positions as a function of sequence-wide mutation count.
#' 
#' To avoid cases where a clonal expansion might lead to a false positive, tigger examines
#' the combinations of J gene and junction length among sequences which perfectly 
#' match the proposed germline allele.
#' 
#' @param    data           a \code{data.frame} in Change-O format. See
#'                          \link{findNovelAlleles} for details.
#' @param    novel_row      a single row from a data frame as output by
#'                          \link{findNovelAlleles} that contains a
#'                          polymorphism-containing germline allele
#' @param    v_call         name of the column in \code{data} with V allele
#'                          calls. Default is "V_CALL".
#' @param    ncol           number of columns to use when laying out the plots  
#' 
#' @examples
#' # Plot the evidence for the first (and only) novel allele in the example data
#' novel <- selectNovel(SampleNovel)
#' plotNovel(SampleDb, novel[1, ])
#' 
#' @export
plotNovel <- function(data, novel_row, v_call="V_CALL", ncol=1) {
    . = NULL
    
    # Use the data frame
    if(length(novel_row) > 0) {
        if(is.data.frame(novel_row) & nrow(novel_row) == 1) {
            pos_range = novel_row$POS_MIN:novel_row$POS_MAX
            germline = novel_row$GERMLINE_IMGT
            names(germline) = novel_row$GERMLINE_CALL
            mut_range = novel_row$MUT_MIN[1]:novel_row$MUT_MAX[1]
            novel_imgt = novel_row$NOVEL_IMGT
            names(novel_imgt) = novel_row$POLYMORPHISM_CALL
            min_frac = novel_row$MIN_FRAC
            note = novel_row$NOTE
        } else {
            stop("novel_row is not a data frame with only one row.")
        }
    }
    
    germline = cleanSeqs(germline)
    data$SEQUENCE_IMGT = cleanSeqs(data$SEQUENCE_IMGT)
    
    # Extract sequences assigned to the germline, determine which
    # have an appropriate range of mutations, and find the mutation
    # frequency of each position
    db_subset = data %>%
        select_(~SEQUENCE_IMGT, v_call, ~J_CALL, ~JUNCTION_LENGTH) %>%
        filter_(~grepl(names(germline),  data[[v_call]], fixed=TRUE))
    pos_db = db_subset %>%  
        mutationRangeSubset(germline, mut_range, pos_range)
    if (nrow(pos_db) == 0) {
        warning(paste0("Insufficient sequences (",nrow(pos_db),") in desired mutational range."))
        return (invisible(NULL))
    }
    pos_db <- pos_db %>%
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
    # Parse the note to find positions that passed y intercept if no novel found
    if(length(pass_y) == 0 & grepl("Position\\(s\\) passed y-intercept", note)){
        pass_y = note %>% gsub("Position\\(s\\) passed y-intercept \\(", "", .) %>%
            gsub("\\).*", "", .) %>% strsplit(",") %>% unlist %>% as.numeric
        p_y_f = sapply(pass_y, function (x) substring(germline, x, x))
        p_y_t = gsub(".", "?", p_y_f)
    }
    
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
    if (nrow(db_subset) == 0) {
        warning(paste0("Insufficient sequences (",nrow(db_subset),") with MUT_COUNT_NOVEL == 0."))
        return (invisible(NULL))
    }
    db_subset$JUNCTION_LENGTH = db_subset$JUNCTION_LENGTH %>%
        factor(levels=min(db_subset$JUNCTION_LENGTH):max(db_subset$JUNCTION_LENGTH))
    pos_muts$Polymorphic = pos_muts$Polymorphic %>%
        factor(levels = c("False", "True"))
    pos_db$NT = pos_db$NT %>%
        factor(levels = names(DNA_COLORS))
    pos_muts$GERMLINE = names(germline)
    
    # MAKE THE FIRST PLOT
    if(!is.na(novel_imgt)){
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
    } else{
        POLYCOLORS = setNames(DNA_COLORS[c(4,2)], c("False", "True"))
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
            guides(color = guide_legend("Passed y-intercept test",
                                        ncol = 2, reverse = TRUE))
    }
    # MAKE THE SECOND PLOT
    p2_data = mutate_(filter_(pos_db, ~POSITION %in% pass_y),
                      POSITION = ~to_from[as.character(POSITION)])
    if (nrow(p2_data)) {
        p2 = ggplot(p2_data,
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
    } else {
        p2_data = mutate_(filter_(pos_db,
                                  ~POSITION %in% names(which.max(table(pos_db$POSITION)))),
                          POSITION = ~"No positions pass y-intercept test.")
        p2 = ggplot(p2_data, aes_(~factor(MUT_COUNT))) +
            geom_bar(width=0.9) +
            facet_grid(POSITION ~ .) +
            xlab("Mutation Count (Sequence)") + ylab("Sequence Count") +
            theme_bw() +
            theme(legend.position=c(1,1), legend.justification=c(1,1),
                  legend.background=element_rect(fill = "transparent"))
    }
    # MAKE THE THIRD PLOT
    p3 = ggplot(db_subset, aes_(~JUNCTION_LENGTH, fill=~factor(J_GENE))) +
        geom_bar(width=0.9) +
        guides(fill = guide_legend("J Gene", ncol = 2)) +
        xlab("Junction Length") + ylab("Unmutated Sequence Count") +
        theme_bw() +
        theme(legend.position=c(1,1), legend.justification=c(1,1),
              legend.background=element_rect(fill = "transparent"))
    
    p2_height = length(unique(p2_data$POSITION))
    if (p2_height>1) { p2_height = 0.5 * p2_height}
    heights = c(1, p2_height, 1)
    multiplot(p1, p2, p3, cols = ncol, heights=heights)      
}

#' Infer a subject-specific genotype using a frequency method
#'
#' \code{inferGenotype} infers an subject's genotype using a frequency method.
#' The genotype is inferred by finding the minimum number set of alleles that 
#' can explain the majority of each gene's calls. The most common allele of 
#' each gene is included in the genotype first, and the next most common allele 
#' is added until the desired fraction of alleles can be explained. In this 
#' way, mistaken allele calls (resulting from sequences which
#' by chance have been mutated to look like another allele) can be removed.
#' 
#' @details
#' Allele calls representing cases where multiple alleles have been
#' assigned to a single sample sequence are rare among unmutated
#' sequences but may result if nucleotides for certain positions are
#' not available. Calls containing multiple alleles are treated as
#' belonging to all groups. If \code{novel} is provided, all
#' sequences that are assigned to the same starting allele as any
#' novel germline allele will have the novel germline allele appended
#' to their assignent prior to searching for unmutated sequences.
#'           
#' @param    data                 a \code{data.frame} containing V allele
#'                                calls from a single subject. If
#'                                \code{find_unmutated} is \code{TRUE}, then
#'                                the sample IMGT-gapped V(D)J sequence should 
#' @param    germline_db          named vector of sequences containing the
#'                                germline sequences named in
#'                                \code{allele_calls}. Only required if
#'                                \code{find_unmutated} is \code{TRUE}.
#' @param    novel             an optional \code{data.frame} of the type
#'                                novel returned by
#'                                \link{findNovelAlleles} containing
#'                                germline sequences that will be utilized if
#'                                \code{find_unmutated} is \code{TRUE}. See
#'                                Details.
#' @param    v_call               column in \code{data} with V allele calls.
#'                                Default is \code{"V_CALL"}.                            
#'                                be provided in a column \code{"SEQUENCE_IMGT"}
#' @param    fraction_to_explain  the portion of each gene that must be
#'                                explained by the alleles that will be included
#'                                in the genotype.
#' @param    gene_cutoff          either a number of sequences or a fraction of
#'                                the length of \code{allele_calls} denoting the
#'                                minimum number of times a gene must be
#'                                observed in \code{allele_calls} to be included
#'                                in the genotype.
#' @param    find_unmutated       if \code{TRUE}, use \code{germline_db} to
#'                                find which samples are unmutated. Not needed
#'                                if \code{allele_calls} only represent
#'                                unmutated samples.
#' 
#' @return
#' A \code{data.frame} of alleles denoting the genotype of the subject containing 
#' the following columns:
#'           
#' \itemize{
#'   \item \code{GENE}: The gene name without allele.
#'   \item \code{ALLELES}: Comma separated list of alleles for the given \code{GENE}.
#'   \item \code{COUNTS}: Comma separated list of observed sequences for each 
#'         corresponding allele in the \code{ALLELES} list.
#'   \item \code{TOTAL}: The total count of observed sequences for the given \code{GENE}.
#'   \item \code{NOTE}: Any comments on the inferrence.
#' }
#'           
#' @note
#' This method works best with data derived from blood, where a large
#' portion of sequences are expected to be unmutated. Ideally, there
#' should be hundreds of allele calls per gene in the input.
#' 
#' @seealso \link{plotGenotype} for a colorful visualization and
#'          \link{genotypeFasta} to convert the genotype to nucleotide sequences.
#'          See \link{inferGenotypeBayesian} to infer a subject-specific genotype 
#'          using a Bayesian approach.
#' 
#' @examples
#' # Infer IGHV genotype, using only unmutated sequences, including novel alleles
#' inferGenotype(SampleDb, germline_db=GermlineIGHV, novel=SampleNovel,
#'               find_unmutated=TRUE)
#' 
#' @export
inferGenotype <- function(data, germline_db=NA, novel=NA, v_call="V_CALL", 
                          fraction_to_explain=0.875, gene_cutoff=1e-4, 
                          find_unmutated=TRUE) {
    
    . = NULL
    allele_calls = getAllele(data[[v_call]], first=FALSE, strip_d=FALSE)
    # Find the unmutated subset, if requested
    if(find_unmutated){
        if(is.na(germline_db[1])){
            stop("germline_db needed if find_unmutated is TRUE")
        }
        if(!is.null(nrow(novel))){
            novel = filter_(novel, ~!is.na(POLYMORPHISM_CALL)) %>%
                select_(~GERMLINE_CALL, ~POLYMORPHISM_CALL, ~NOVEL_IMGT)
            if(nrow(novel) > 0){
                # Extract novel alleles if any and add them to germline_db
                novel_gl = novel$NOVEL_IMGT
                names(novel_gl) = novel$POLYMORPHISM_CALL
                germline_db = c(germline_db, novel_gl)
                # Add the novel allele calls to allele calls of the same starting allele
                for(r in 1:nrow(novel)){
                    ind = grep(novel$GERMLINE_CALL[r], allele_calls, fixed=TRUE)
                    allele_calls[ind] = allele_calls[ind] %>%
                        sapply(paste, novel$POLYMORPHISM_CALL[r], sep=",")
                }
            }
        }
        # Find unmutated sequences
        allele_calls = findUnmutatedCalls(allele_calls,
                                          as.character(data$SEQUENCE_IMGT),
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
    for (g in GENE) {
        # Keep only the part of the allele calls that uses the gene being analyzed
        ac = allele_calls[gene_groups[[g]]] %>%
            strsplit(",") %>%
            lapply(function(x) x[grep(paste(g, "\\*", sep=""), x)]) %>%
            sapply(paste, collapse=",")
        target = ceiling(fraction_to_explain*length(ac)) # how many we need to explain
        t_ac = table(ac) # table of allele calls
        potentials = unique(unlist(strsplit(names(t_ac),","))) # potential alleles
        # One allele? Easy!
        if (length(potentials) == 1 | length(t_ac) == 1) {
            genotype[genotype[,"GENE"]==g,"ALLELES"] = gsub("[^d\\*]*[d\\*]","",potentials )[1]
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
    if (find_unmutated == TRUE) {
        seqs = genotypeFasta(geno, germline_db)
        dist_mat = seqs %>%
            sapply(function(x) sapply((getMutatedPositions(seqs, x)), length)) %>%
            as.matrix
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
#' @param    genotype     a \code{data.frame} of alleles denoting a genotype, 
#'                        as returned by \link{inferGenotype}.
#' @param    facet_by     a column name in \code{genotype} to facet the plot by. 
#'                        If \code{NULL}, then do not facet the plot. 
#' @param    gene_sort    a string defining the method to use when sorting alleles.
#'                        If \code{"name"} then sort in lexicographic order. If
#'                        \code{"position"} then sort by position in the locus, as
#'                        determined by the final two numbers in the gene name.
#' @param    text_size    the point size of the plotted text.
#' @param    silent       if \code{TRUE} do not draw the plot and just return the ggplot
#'                        object; if \code{FALSE} draw the plot.
#' @param    ...          additional arguments to pass to ggplot2::theme.
#' 
#' @return  A ggplot object defining the plot.
#' 
#' @seealso \link{inferGenotype}
#' 
#' @examples
#' # Plot genotype
#' plotGenotype(SampleGenotype)
#' 
#' # Facet by subject
#' genotype_a <- genotype_b <- SampleGenotype
#' genotype_a$SUBJECT <- "A"
#' genotype_b$SUBJECT <- "B"
#' geno_sub <- rbind(genotype_a, genotype_b)
#' plotGenotype(geno_sub, facet_by="SUBJECT", gene_sort="pos")
#' 
#' @export
plotGenotype <- function(genotype, facet_by=NULL, gene_sort=c("name", "position"), 
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
#' @param    genotype     a \code{data.frame} of alleles denoting a genotype, 
#'                        as returned by \link{inferGenotype}.
#' @param    germline_db  a vector of named nucleotide germline sequences
#'                        matching the alleles detailed in \code{genotype}.
#' @param    novel        an optional \code{data.frame} containing putative
#'                        novel alleeles of the type returned by
#'                        \link{findNovelAlleles}.
#' 
#' @return   A named vector of strings containing the germline nucleotide
#'           sequences of the alleles in the provided genotype.
#' 
#' @seealso \link{inferGenotype}
#' 
#' @examples
#' # Find the sequences that correspond to the genotype
#' genotype_db <- genotypeFasta(SampleGenotype, GermlineIGHV, SampleNovel)
#' 
#' @export
genotypeFasta <- function(genotype, germline_db, novel=NA){
    if(!is.null(nrow(novel))){
        # Extract novel alleles if any and add them to germline_db
        novel = filter_(novel, ~!is.na(POLYMORPHISM_CALL)) %>%
            select_(~GERMLINE_CALL, ~POLYMORPHISM_CALL, ~NOVEL_IMGT)
        if(nrow(novel) > 0){
            novel_gl = novel$NOVEL_IMGT
            names(novel_gl) = novel$POLYMORPHISM_CALL
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
#' @details
#' In order to save time, initial gene assignments are preserved and
#' the allele calls are chosen from among those provided in \code{genotype_db},
#' based on a simple alignment to the sample sequence.
#' 
#' @param    data          a \code{data.frame} containing V allele calls from a
#'                         single subject and the sample IMGT-gapped V(D)J sequences under
#'                         \code{"SEQUENCE_IMGT"}.
#' @param    genotype_db   a vector of named nucleotide germline sequences
#'                         matching the calls detailed in \code{allele_calls}
#'                         and personalized to the subject
#' @param    v_call        name of the column in \code{data} with V allele
#'                         calls. Default is \code{"V_CALL"}.                 
#' @param    method        the method to be used when realigning sequences to
#'                         the genotype_db sequences. Currently, only \code{"hammming"}
#'                         (for Hamming distance) is implemented.
#' @param    path          directory containing the tool used in the
#'                         realignment method, if needed. Hamming distance does
#'                         not require a path to a tool.
#' @param    keep_gene     a string indicating if the gene (\code{"gene"}), 
#'                         family (\code{"family"}) or complete repertoire
#'                         (\code{"repertoire"}) assignments should be performed. 
#'                         Use of \code{"gene"} increases speed by minimizing required number of 
#'                         alignments, as gene level assignments will be maintained when possible.
#' 
#' @return   A modifed input \code{data.frame} containing the best allele call from 
#'           among the sequences listed in \code{genotype_db} in the 
#'           \code{V_CALL_GENOTYPED} column.
#' 
#' @examples
#' # Extract the database sequences that correspond to the genotype
#' genotype_db <- genotypeFasta(SampleGenotype, GermlineIGHV, novel=SampleNovel)
#' 
#' # Use the personlized genotype to determine corrected allele assignments
#' output_db <- reassignAlleles(SampleDb, genotype_db)
#' 
#' @export
reassignAlleles <- function(data, genotype_db, v_call="V_CALL",
                            method="hamming", path=NA,
                            keep_gene=c("gene", "family", "repertoire")){
    # Check arguments    
    keep_gene <- match.arg(keep_gene)
    
    # Extract data subset and prepare output vector
    v_sequences = as.character(data$SEQUENCE_IMGT)
    v_calls = getAllele(data[[v_call]], first=FALSE, strip_d=FALSE)
    v_call_genotyped = rep("", length(v_calls))
    
    if (keep_gene == "gene") { 
        v = getGene(v_calls, first = TRUE, strip_d=FALSE)
        geno = getGene(names(genotype_db),strip_d=TRUE)
        names(geno) = names(genotype_db)
    } else if (keep_gene == "family") {
        v <- getFamily(v_calls, first = TRUE, strip_d = FALSE)
        geno = getFamily(names(genotype_db),strip_d=TRUE)
        names(geno) = names(genotype_db)
    } else if (keep_gene == "repertoire") {
        v <- rep(v_call, length(v_calls))
        geno = rep(v_call, length(genotype_db))
        names(geno) = names(genotype_db)      
    } else {
        stop("Unknown keep_gene value: ", keep_gene)
    }
    
    # keep_gene == FALSE
    # Find which genotype genes/families are homozygous and assign those alleles first
    hetero = unique(geno[which(duplicated(geno))])
    homo = geno[!(geno %in% hetero)]
    homo_alleles = names(homo)
    names(homo_alleles) = homo
    homo_calls_i = which(v %in% homo)
    v_call_genotyped[homo_calls_i] = homo_alleles[v[homo_calls_i]]
    
    # Now realign the heterozygote sequences to each allele of that gene
    for (het in hetero){
        ind = which(v %in% het)
        if (length(ind) > 0){
            het_alleles = names(geno[which(geno == het)])
            het_seqs = genotype_db[het_alleles]
            if(method == "hamming"){
                dists = lapply(het_seqs, function(x)
                    sapply(getMutatedPositions(v_sequences[ind], x, match_instead=FALSE),
                           length))
                dist_mat = matrix(unlist(dists), ncol = length(het_seqs))
            } else {
                stop("Only Hamming distance is currently supported as a method.")
            }
            # The sapply-apply approach could become problematic when nrow(dist_mat)
            # is 1 and min(best_match) has multiple values, due to the fact that R 
            # does not always keep data structures unmutable
            # Explicitly specifying a list and subsequently keeping it as a list by
            # using lapply avoids that problem
            best_match = vector("list", length=nrow(dist_mat))
            for (i in 1:nrow(dist_mat)) {
                best_match[[i]] = which(dist_mat[i, ]==min(dist_mat[i, ]))
            }
            best_alleles = lapply(best_match, function(x) het_alleles[x])
            v_call_genotyped[ind] = unlist(lapply(best_alleles, paste, collapse=","))
        }
    }
    
    # Now realign the gene-not-in-genotype calls to every genotype allele
    hetero_calls_i = which(v %in% hetero)
    not_called = setdiff(1:length(v), c(homo_calls_i, hetero_calls_i))
    if(length(not_called)>1){
        if(method ==  "hamming"){
            dists = lapply(genotype_db, function(x)
                sapply(getMutatedPositions(v_sequences[not_called], x, match_instead=FALSE),
                       length))
            dist_mat = matrix(unlist(dists), ncol = length(genotype_db))
        } else {
            stop("Only Hamming distance is currently supported as a method.")
        }
        # The sapply-apply approach could become problematic when nrow(dist_mat)
        # is 1 and min(best_match) has multiple values, due to the fact that R 
        # does not always keep data structures unmutable
        # Explicitly specifying a list and subsequently keeping it as a list by
        # using lapply avoids that problem
        best_match = vector("list", length=nrow(dist_mat))
        for (i in 1:nrow(dist_mat)) {
            best_match[[i]] = which(dist_mat[i, ]==min(dist_mat[i, ]))
        }
        best_alleles = lapply(best_match, function(x) names(genotype_db[x]))
        v_call_genotyped[not_called] = unlist(lapply(best_alleles, paste, collapse=","))
    }
    
    if (all(v_call_genotyped == data[[v_call]])) {
        msg <- ("No allele assignment corrections made.") 
        if (all(v %in% homo) & length(hetero) > 0) {
            keep_opt <- eval(formals(reassignAlleles)$keep_gene)
            i <- match(keep_gene, keep_opt)
            rec_opt <- paste(keep_opt[(i+1):length(keep_opt)], collapse = ", ")
            msg <- paste(msg, "Consider setting keep_gene to one of:", rec_opt)
        }
        warning(msg)
    }
    
    data$V_CALL_GENOTYPED <- v_call_genotyped
    
    return(data)
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
#' seqs <- c("----GATA", "GAGAGAGA", "TANA")
#' ref <- "GATAGATA"
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
#' # Insert a mutation into a germline sequence
#' s2 <- s3 <- GermlineIGHV[1]
#' stringi::stri_sub(s2, 103, 103) <- "G"
#' stringi::stri_sub(s3, 107, 107) <- "C"
#' 
#' sample_seqs <- c(GermlineIGHV[2], s2, s3)
#' 
#' # Pretend that one sample sequence has received an ambiguous allele call
#' sample_alleles <- c(paste(names(GermlineIGHV[1:2]), collapse=","),
#'                     names(GermlineIGHV[2]),
#'                     names(GermlineIGHV[1]))
#' 
#' # Compare each sequence to its assigned germline(s) to determine the distance
#' getMutCount(sample_seqs, sample_alleles, GermlineIGHV)
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
#'                          where multiple calls are separated by a comma.
#' @param    germline_db    a vector of named nucleotide germline sequences
#' @param    sample_seqs    V(D)J-rearranged sample sequences matching the order
#'                          of the given \code{allele_calls}.
#' 
#' @return   A vector of strings containing the members of \code{allele_calls}
#'           that represent unmutated sequences.
#' 
#' @examples
#' # Find which of the sample alleles are unmutated
#' calls <- findUnmutatedCalls(SampleDb$V_CALL, SampleDb$SEQUENCE_IMGT, 
#'          germline_db=GermlineIGHV)
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

#' Find mutation counts for frequency sequences
#'
#' \code{getPopularMutationCount} determines which sequences occur frequently
#' for each V gene and returns the mutation count of those sequences.
#' 
#' @param  data          a \code{data.frame} in the Change-O format. See
#'                       \link{findNovelAlleles} for a list of required
#'                       columns.
#' @param  germline_db   A named list of IMGT-gapped germline sequences.
#' @param  gene_min      The portion of all unique sequences a gene must
#'                       constitute to avoid exclusion.
#' @param  seq_min       The number of copies of the V that must be present for
#'                       to avoid exclusion.
#' @param  seq_p_of_max  For each gene, fraction of the most common V sequence's
#'                       count that a sequence must meet to avoid exclusion.
#' @param  full_return   If \code{TRUE}, will return all \code{data} columns and
#'                       will include sequences with mutation count < 1.
#' 
#' @return  A data frame of genes that have a frequent sequence mutation count
#'          above 1.
#' 
#' @seealso \link{getMutatedPositions} can be used to find which positions
#'          of a set of sequences are mutated.
#' 
#' @examples
#' getPopularMutationCount(SampleDb, GermlineIGHV)
#' 
#' @export
getPopularMutationCount <- function(data, germline_db, gene_min = 1e-03,
                                    seq_min = 50, seq_p_of_max = 1/8,
                                    full_return = FALSE){
    modified_db = data %>%
        mutate_(V_GENE = ~getGene(V_CALL)) %>%
        group_by_(~V_GENE) %>%
        mutate_(V_GENE_N = ~n()) %>%
        group_by_(~1:n()) %>%
        mutate_(V_SEQUENCE_IMGT = ~substring(SEQUENCE_IMGT, 1, 312)) %>%
        # Count occurence of each unique IMGT-gapped V sequence
        group_by_(~V_GENE, ~V_SEQUENCE_IMGT) %>%
        mutate_(V_SEQUENCE_IMGT_N = ~n()) %>%
        # Determine count of most common sequence
        group_by_(~V_GENE) %>%
        mutate_(V_SEQUENCE_IMGT_N_MAX = ~max(V_SEQUENCE_IMGT_N)) %>%
        # Remove rare V genes, rare sequences, and sequences not making up a
        # sufficient proportion of sequences as compared to the most common
        ungroup %>%
        distinct_(~V_SEQUENCE_IMGT, .keep_all = TRUE) %>%
        filter_(~V_GENE_N >= (nrow(data)*gene_min)) %>%
        filter_(~V_SEQUENCE_IMGT_N >= seq_min) %>%
        mutate_(V_SEQUENCE_IMGT_P_MAX = ~V_SEQUENCE_IMGT_N/V_SEQUENCE_IMGT_N_MAX) %>%
        filter_(~V_SEQUENCE_IMGT_P_MAX >= seq_p_of_max)
    # Determine the mutation counts of the V sequences and append them to the db
    MUTATION_COUNT = getMutCount(modified_db$V_SEQUENCE_IMGT,
                                 modified_db$V_CALL,
                                 germline_db) %>% 
        sapply(function(x) min(unlist(x)))
    if (length(MUTATION_COUNT)==0){
        MUTATION_COUNT = integer(0)
    }
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
#' @param    sequence     starting nucletide sequence.
#' @param    positions    numeric vector of positions which to be changed.
#' @param    nucleotides  character vector of nucletides to which to change the
#'                        positions.
#'                        
#' @return   A sequence with the desired nucleotides in the provided locations.
#' 
#' @examples
#' insertPolymorphisms("HUGGED", c(1, 6, 2), c("T", "R", "I")) 
#' 
#' @export
insertPolymorphisms <- function(sequence, positions, nucleotides) {
    
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
#' @param    fasta_file       fasta-formatted file of immunoglobuling sequences.
#' @param    strip_down_name  if \code{TRUE}, will extract only the allele name
#'                            from the strings fasta file's sequence names.
#' @param    force_caps       if \code{TRUE}, will force nucleotides to
#'                            uppercase.
#'                            
#' @return   Named vector of strings respresenting Ig alleles.
#' 
#' @seealso  \link{writeFasta} to do the inverse.
#' 
#' @export
readIgFasta <- function(fasta_file, strip_down_name=TRUE, force_caps=TRUE) {
    all_char = readChar(fasta_file, file.info(fasta_file)$size)
    split_by_sequence = strsplit(all_char, "[ \t\r\n\v\f]?>")
    add_name_break = sapply(split_by_sequence, function(x) sub("[\r\n]",">",x))
    cleaned_up = sapply(add_name_break, function(x) gsub("[ \t\r\n\v\f]", "", x))
    broken_names = sapply(cleaned_up, strsplit, ">")
    
    seqs = sapply(broken_names, "[", 2)
    seq_names = sapply(broken_names, "[", 1)
    if(force_caps) { seqs = toupper(seqs) }
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
#' @param    file             the name of the output file.
#' @param    width            the number of characters to be printed per line.
#'                            if not between 1 and 255, width with be infinite.
#' @param    append           \code{logical} indicating if the output should be
#'                            appended to \code{file} instead of overwriting it
#' 
#' @return   A named vector of strings respresenting Ig alleles.
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
#' 
#' @param    allele_calls  a vector of strings respresenting IGHV allele names.
#' 
#' @return   Vector of strings respresenting updated IGHV allele names.
#' 
#' @note
#' IGMT has removed \code{IGHV2-5*10} and \code{IGHV2-5*07} as it has determined they
#' are actually alleles \code{02} and \code{04}, respectively. The updated allele 
#' names are based on IMGT release 201408-4.
#' 
#' @references
#' \enumerate{
#'   \item Xochelli et al. (2014) Immunoglobulin heavy variable (IGHV) genes
#'         and alleles: new entities, new names and implications for research and
#'         prognostication in chronic lymphocytic leukaemia. Immunogenetics. 67(1):61-6
#' }
#' 
#' @seealso Like \code{updateAlleleNames}, \link{sortAlleles} can help
#'          format a list of allele names.
#' 
#' @examples
#' # Create a vector that uses old gene/allele names.
#' alleles <- c("IGHV1-c*01", "IGHV1-f*02", "IGHV2-5*07")
#' 
#' # Update the alleles to the new names
#' updateAlleleNames(alleles)
#' 
#' @export
updateAlleleNames <- function(allele_calls) {
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
#' non-duplicated counterparts (e.g. \code{IGHV1-69D*01} comes after \code{IGHV1-69*01} 
#' but before \code{IGHV1-69*02}), and non-localized genes (e.g. \code{IGHV1-NL1*01}) 
#' come last within their gene family.
#' 
#' @param    allele_calls  a vector of strings respresenting Ig allele names.
#' @param    method        a string defining the method to use when sorting alleles.
#'                         If \code{"name"} then sort in lexicographic order. If
#'                         \code{"position"} then sort by position in the locus, as
#'                         determined by the final two numbers in the gene name.
#' @return   A sorted vector of strings respresenting Ig allele names.
#' 
#' @seealso Like \code{sortAlleles}, \link{updateAlleleNames} can help
#'          format a list of allele names.
#' 
#' @examples
#' # Create a list of allele names
#' alleles <- c("IGHV1-69D*01","IGHV1-69*01","IGHV1-2*01","IGHV1-69-2*01",
#'              "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
#'              "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")
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
        mutate_(GENE1 = ~gsub("[^-]+[-S]([^-\\*D]+).*","\\1",SUBMITTED_CALLS)) %>%
        mutate_(GENE1 = ~as.numeric(gsub("[^0-9]+", "99", GENE1))) %>%
        # If there is a second gene number, determine that, too
        mutate_(GENE2 = ~gsub("[^-]+[-S][^-]+-?","",GENE)) %>%
        mutate_(GENE2 = ~as.numeric(gsub("[^0-9]+", "99", GENE2))) %>%
        mutate_(ALLELE = ~getAllele(SUBMITTED_CALLS)) %>%      
        mutate_(ALLELE = ~(sub("[^\\*]+\\*|[^\\*]+$","",
                               ALLELE))) %>%
        mutate_(ALLELE = ~as.numeric(sub("_.+$","",
                                         ALLELE)))
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
#' \code{cleanSeqs} capitalizes nucleotides and replaces all characters 
#' besides \code{c("A", "C", "G", "T", "-", ".")} with \code{"N"}. 
#' 
#' @param    seqs  a vector of nucleotide sequences.
#' 
#' @return   A modified vector of nucleotide sequences.
#' 
#' @seealso \link{sortAlleles} and \link{updateAlleleNames} can
#'          help format a list of allele names.
#' 
#' @examples
#' # Clean messy nucleotide sequences
#' seqs <- c("AGAT.taa-GAG...ATA", "GATACAGTXXZZAGNNPPACA")
#' cleanSeqs(seqs)
#' 
#' @export
cleanSeqs <- function(seqs) {
    # . = NULL
    # seqs %>%
    #     toupper %>%
    #     gsub(".", "-", . , fixed = TRUE) %>%
    #     gsub("[^ACGT-]", "N", .) %>%
    #     return

    return (gsub("[^ACGT\\.\\-]", "N", toupper(seqs)))
}


# Private Functions -------------------------------------------------------

# Find muations-by-position compared to a germline
#
# \code{positionMutations} duplicates the rows of a data frame for each
# position to be analyzed and determines if each sample is mutated at that
# position
# 
# @param  data          a Change-O db data.frame. See
#                       \link{findNovelAlleles} for a list of required
#                       columns.
# @param  germline      the germline to which all the sequences should be
#                       compared
# @param  pos_range     the range of positions within the sequence for which
#                       the rows should be duplicated and checked for mutation
# 
# @return  A data frame with rows duplicated for all the positions to be
# analyzed and a column indicating whether the position is mutated in
# comparison to the germline
#
positionMutations <- function(data, germline, pos_range){
    . = NULL
    pos_db = pos_range %>%
        length() %>%
        rep("data", .) %>%
        paste(collapse=",") %>%
        paste("bind_rows(",., ")") %>%
        parse(text=.) %>%
        eval()
    pos_db$POSITION = c(sapply(pos_range, rep, nrow(data)))
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
# @param  data          a Change-O db data frame. See
#                       \link{findNovelAlleles} for a list of required
#                       columns.
# @param  germline      the germline to which all the sequences should be
#                       compared
# @param  pos_range     the range of positions within the sequences that should
#                       be analyzed for mutations
# @param  pos_range     the range of mutation counts that sequences can have
#                       and still be included
#
# @return
# A data.frame containing only the subset carrying the desired levels
# of mutation
#
mutationRangeSubset <- function(data, germline, mut_range, pos_range){
    . = NULL
    pads = paste(rep("-", min(pos_range)-1), collapse="")
    data$MUT_COUNT = data$SEQUENCE_IMGT %>%
        substring(min(pos_range), max(pos_range)) %>%
        paste(pads, ., sep="") %>%
        getMutatedPositions(germline) %>%
        sapply(length)
    data = data %>%
        filter_(~MUT_COUNT %in% mut_range)
    return(data)
}

# Find lower range of y-intercept confidence interval
#
# \code{findLowerY} finds the lower range of y-intercept confidence interval
# 
# @details  If mut_min is 1, a y-intercept will be searched for at 0. If
# mut_min is above 1, then the "y-intercept" will be found at x = mut_min - 1.
#
# @param    x         A vector of x values
# @param    y         A vector of y values
# @param    mut_min   The value where the the lowest mutation count should be
#                     found. See details.
# @param    alpha     The alpha cutoff the be used in constructing the
#                     confidence interval
#
# @return  A data frame containing only the subset carrying the desired levels
# of mutation
#
findLowerY = function(x, y, mut_min, alpha){
    y = y + 1 - mut_min
    lowerY = suppressWarnings(confint(lm(x ~ y), level=1 - 2*alpha)[[1]])
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
# @param    heights   A numeric vector A numeric vector or unit object 
#                     describing the heights of the rows in the layout. Will
#                     be passed to grid.layout. Default is all plots have 
#                     the same height.
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, heights=NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    ncol <- cols
    nrow <- ceiling(numPlots/cols)
    if (is.null(heights)) { heights = rep(1,nrow) }
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * nrow),
                         ncol = cols, nrow = nrow)
    }
    grob <- gridExtra::arrangeGrob(grobs=plots, 
                                   nrow=nrow, ncol=ncol, layout_matrix = layout,
                                   heights=heights)
    p <- ggplot() + 
        layer(data = data.frame(x = NA),
              stat = StatIdentity, 
              position = PositionIdentity, 
              # geom = GeomDrawGrob, 
              geom = GeomCustomAnn,
              inherit.aes = FALSE, 
              params = list(grob = grob, 
                            xmin = 0,
                            xmax = 1, 
                            ymin = 0, 
                            ymax = 1)) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))
    p
}
