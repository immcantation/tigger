# Find non triplet gaps in a nucleotide sequence
hasNonImgtGaps <- function (seq) {
    len <- ceiling(nchar(seq)/3)*3
    codons <- substring(seq, seq(1, len-2, 3), seq(3, len, 3))
    gaps_lengths <- nchar(gsub("[^\\.\\-]", "", codons))
    if (any(gaps_lengths %% 3 != 0)) {
        TRUE
    } else {
        FALSE
    }
}

# Compare two IMGT gapped sequences and find AA mutations
getMutatedAA <- function(ref_imgt, novel_imgt) {
    if (grepl("N", ref_imgt)) {
        stop("Unexpected N in ref_imgt")
    }     
    if (grepl("N", novel_imgt)) {
        stop("Unexpected N in novel_imgt")
    }          
    
    if (hasNonImgtGaps(ref_imgt)) {
        warning("Non IMGT gaps found in ref_imgt")
    }
    
    if (hasNonImgtGaps(novel_imgt)) {
        warning("Non IMGT gaps found in novel_imgt")
    }
    
    ref_imgt <- strsplit(alakazam::translateDNA(ref_imgt),"")[[1]]
    novel_imgt <- strsplit(alakazam::translateDNA(novel_imgt),"")[[1]]
    mutations <- c()
    diff_idx <- which(ref_imgt != novel_imgt)
    if (length(diff_idx)>0) {
        mutations <- paste0(diff_idx, ref_imgt[diff_idx],">",novel_imgt[diff_idx])
    }
    mutations
}


#' Generate evidence
#'
#' Generate evidence of the final novel V alleles.
#' 
#' @return a data.frame
#' 
#' @param gt  A table of alleles generated with \link{inferGenotype} denoting 
#'            the genotype of the subject.
#' @param nv  A data.frame returned by \link{findNovelAlleles}.
#' @param germline_nv A vector of named nucleotide germline sequences 
#'            matching the V calls in clip_db. This is an extended version of 
#'            \code{germline_db} in that includes the original reference datatabase
#'            and any novel alleles identified after running 
#'            \link{findNovelAlleles} and \link{inferGenotype}.
#' @param germline_db The original input germline database used to by
#'            \link{findNovelAlleles} to identify novel alleles in 
#'            \code{db}.
#' @param iteration_id Column name that identifies iterations from 
#'            if \link{itigger} was used.
#' @param fields  Column names of fields used to split the data in \link{itigger}
#' @return   Returns \code{gt} with additional columns providing supporting evidence
#'           for each inferred allele.
#'    \itemize{
#'       \item \code{FIELD_ID} Data subset identifier, defined with the input paramter \code{fields}.
#'       \item A variable number of columns, specified with the input parameter \code{fields}.
#'       \item \code{ITERATION} The iteration number.
#'       \item \code{POLYMORPHISM_CALL} The new allele call.
#'       \item \code{CLOSEST_REFERENCE} The closest reference gene and allele in the input germline 
#'                                database (\code{germline}).
#'       \item \code{NT_DIFF} Number of nucleotides that differ between the new allele and
#'                      the closest reference (\code{CLOSEST_REFERENCE}) in the input (\code{germline}).
#'       \item \code{NT_SUBSTITUTIONS} A \code{character} with specific nucleotide differences (e.g. 112G>A),
#'                               comma separated.
#'       \item \code{AA_DIFF} Number of aminoacids that differ between the new allele and the closest 
#'                      reference (\code{CLOSEST_REFERENCE}) in the input (\code{germline}).
#'       \item \code{AA_SUBSTITUTIONS} A \code{character} with specific aminoacid differences (e.g. 96A>N),
#'                               comma separated.
#'       \item \code{SEQUENCES} Number of sequences unambiguosly assigned to this allele.
#'       \item \code{UNMUTATED_SEQUENCES} Number of records with the unmutated new allele sequence.
#'       \item \code{UNMUTATED_FREQUENCY} Proportion of records with the unmutated new allele 
#'                                  sequence. \code{UNMUTATED_SEQUENCES}/\code{SEQUENCE}
#'       \item \code{ALLELIC_PERCENTAGE} Percentage at which this (unmutated) allele is observed in the 
#'                                 sequence dataset, compared  to other (unmutated) alleles.
#'       \item \code{UNIQUE_JS} Number of unique J sequences found associated with the new allele. The sequences
#'                        are those who have been unambiguously assigned the new alelle (\code{POLYMORPHISM_CALL})
#'       \item \code{UNIQUE_CDR3S} Number of unique CDR3s associated with the inferred allele
#'       \item \code{GERMLINE_CALL} See \link{findNovelAlleles}
#'       \item \code{MUT_MIN} See \link{findNovelAlleles}
#'       \item \code{MUT_MAX} See \link{findNovelAlleles}
#'       \item \code{POS_MIN} See \link{findNovelAlleles}
#'       \item \code{POS_MAX} See \link{findNovelAlleles}
#'       \item \code{Y_INTERCEPT} See \link{findNovelAlleles}
#'       \item \code{ALPHA} See \link{findNovelAlleles}
#'       \item \code{MIN_SEQS} See \link{findNovelAlleles}
#'       \item \code{J_MAX} See \link{findNovelAlleles}
#'       \item \code{MIN_FRAC} See \link{findNovelAlleles} 
#'       \item \code{NOVEL_IMGT} See \link{findNovelAlleles} 
#'       \item \code{CLOSEST_REFERENCE_IMGT} Sequence of the closest reference gene and 
#'                          allele in the input germline database (\code{germline}).
#'       \item \code{GERMLINE_IMGT} See \link{findNovelAlleles}
#'       \item \code{NOTE} See \link{findNovelAlleles}
#'    }     
#' @export
generateEvidence <- function(gt, nv, germline_nv, germline_input,
                             db, iteration_id=NULL, fields=NULL) {
    
    # Find closest reference
    .findClosestReference <- function(seq, allele_calls, ref_germ, 
                                      exclude_self=F) {
        closest <- getMutCount(seq,
                               paste(allele_calls, collapse=","),
                               ref_germ)
        min_dist <- min(unlist(closest))
        closest_idx <- which(unlist(closest) == min_dist)
        closest_names <- unique(allele_calls[closest_idx])
        if (exclude_self & names(seq) %in% closest_names) {
            warning("Excluding self")
            closest_names <- closest_names[closest_names!=names(seq)] # not self
        }        
        if (length(closest_names) > 1) {
            warning(paste0("More than one closest reference found for ", 
                           names(seq),": ", 
                           paste(closest_names, collapse=",")))
            # Keep the one with less mutated positions
            mut_pos_count <- sapply(gsub("[^_]","",closest_names), nchar)
            closest_names <- closest_names[mut_pos_count==min(mut_pos_count)]
            # Pick same length
            if (length(closest_names) > 1 ) {
                idx <- which(
                    sapply(ref_germ[closest_names],nchar) == nchar(ref_germ[names(seq)])
                )
                if (length(idx) > 0 ) {
                    closest_names <- closest_names[idx]   
                }
            }         
            # Pick same allele
            if (length(closest_names) > 1 ) {
                idx <- which(
                    getAllele(closest_names) == gsub("_.+", "", getAllele(names(seq)))
                )
                if (length(idx) > 0 ) {
                    closest_names <- closest_names[idx]
                }
            } 
            # Pick not duplicated
            if (length(closest_names) > 1 ) {
                idx <- !grepl("D\\*", closest_names)
                if (any(idx)) {
                    closest_names <- closest_names[idx]   
                }
            }            
            # If still more than one, err and TODO
            if (length(closest_names) > 1 ) {
                stop(paste0("Multiple closest reference found for ", 
                            names(seq),":\n", 
                            paste(closest_names, collapse=",")))
            } else {
                warning(paste0("Use: ", 
                               closest_names, 
                               " (less mutated positions, not D, same length, same allele)"))
            }
        }
        paste(closest_names, collapse=",")        
    }
    
    
    # Subset to novel alleles
    final_gt <- gt %>%
        dplyr::group_by(GENE) %>%
        dplyr::filter(duplicated(ALLELES) == FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            ALLELES=strsplit(as.character(ALLELES),","),
            COUNTS=strsplit(as.character(COUNTS),",")) %>%
        tidyr::unnest(ALLELES, COUNTS) %>%
        dplyr::rename(ALLELE=ALLELES) %>%
        dplyr::mutate(POLYMORPHISM_CALL=paste0(GENE,"*" ,ALLELE)) %>%
        dplyr::filter(POLYMORPHISM_CALL %in% nv$POLYMORPHISM_CALL)
    
    # Add info from nv
    final_gt <- merge(final_gt %>% 
                          dplyr::rename(NOTE_GT=NOTE), 
                      nv, 
                      by=c(iteration_id, fields, "POLYMORPHISM_CALL"))
    
    # Add message if the same novel img sequence found from
    # different starting alleles, these will be novel imgt sequences
    # with more than one polymorphism call
    final_gt <- final_gt %>%
        dplyr::group_by(NOVEL_IMGT) %>%
        dplyr::mutate(NUM_CALLS=length(unique(POLYMORPHISM_CALL))) %>%
        dplyr::ungroup()
    idx_mult <- which(final_gt$NUM_CALLS > 1)
    final_gt$NUM_CALLS <- NULL
    if (length(idx_mult)>0) {
        final_gt$NOTE_GT[idx_mult] <- paste(
            final_gt$NOTE_GT[idx_mult],
            " Found multiple polymorphism calls for the same NOVEL_IMGT.", 
            sep="")
    }
    
    
    if (nrow(final_gt)>0) {
        
        .addEvidence <- function(df, germline_nv, germline_input) { 
            polymorphism <- df[['POLYMORPHISM_CALL']]
            novel_imgt <- df[["NOVEL_IMGT"]]
            names(novel_imgt) <- polymorphism
            gene <- df[['GENE']]
            allele <- df[['ALLELE']]
            germline_call <- df[['GERMLINE_CALL']]
            this_germline <- germline_input
            V_CALL_GENOTYPED <- db[["V_CALL_GENOTYPED"]]
            
            SEQUENCES <- sum(V_CALL_GENOTYPED==polymorphism)
            df[["SEQUENCES"]] <- SEQUENCES
            closest_ref_input <- .findClosestReference(novel_imgt,
                                                       names(germline_input), 
                                                       germline_input,
                                                       exclude_self=F)
            closest_ref <- .findClosestReference(novel_imgt,
                                                 names(germline_nv), 
                                                 germline_nv,
                                                 exclude_self=F)
            
            if (getGene(closest_ref_input) != getGene(closest_ref)) {
                warning("closest reference gene difference")
            }
            
            if (closest_ref != polymorphism) {
                warning(paste0("closest reference allele (",
                               closest_ref
                               ,") different from POLYMORPHISM_CALL allele (",
                               polymorphism,")"))
            }
            
            ## TODO: this still not clear.
            ## Any diff using sequence_imgt instead of germline[[polymorphism]]?
            df[["CLOSEST_REFERENCE"]] <- closest_ref_input
            
            nt_diff <- unlist(getMutatedPositions(novel_imgt, germline_nv[[closest_ref_input]]))
            nt_diff_string <- ""
            if (length(nt_diff) > 0 ) {
                nt_diff_string <- paste(paste(
                    nt_diff, 
                    strsplit(germline_nv[[closest_ref_input]],"")[[1]][nt_diff], 
                    ">",
                    strsplit(germline_nv[[polymorphism]],"")[[1]][nt_diff],
                    sep=""), collapse=",")    
            } 
            
            df[["NT_DIFF"]] <- length(nt_diff)
            df[["NT_SUBSTITUTIONS"]] <- nt_diff_string
            
            diff_aa <- getMutatedAA(germline_nv[[closest_ref_input]], germline_nv[[polymorphism]])
            
            if (length(diff_aa)>0) {
                df[["AA_DIFF"]] <- length(diff_aa)
                df[["AA_SUBSTITUTIONS"]] <- paste(diff_aa,collapse=",")
            } else {
                df[["AA_DIFF"]] <- 0
                df[["AA_SUBSTITUTIONS"]] <- ""
            }
            
            df[["COUNTS"]] <- as.numeric(df[["COUNTS"]])
            df[["TOTAL"]] <- as.numeric(df[["TOTAL"]])
            df[["UNMUTATED_SEQUENCES"]] <- as.numeric(df[["COUNTS"]])
            df[["UNMUTATED_FREQUENCY"]] <- as.numeric(df[["COUNTS"]])/SEQUENCES
            
            df[["ALLELIC_PERCENTAGE"]] <- 100*df[["UNMUTATED_SEQUENCES"]]/as.numeric(df[["TOTAL"]])
            
            df[["UNIQUE_JS"]] <- db %>%
                dplyr::filter(V_CALL_GENOTYPED==polymorphism)  %>%
                dplyr::distinct(J_CALL) %>% nrow()
            df[["UNIQUE_CDR3S"]] <- db %>%
                dplyr::filter(V_CALL_GENOTYPED==polymorphism)  %>%
                dplyr::distinct(translateDNA(JUNCTION, trim=TRUE)) %>% 
                nrow()
            # Add closest germline
            df[["CLOSEST_REFERENCE_IMGT"]] <- cleanSeqs(germline_nv[[closest_ref_input]])
            
            data.frame(df, stringsAsFactors = F)
        }
        
        final_gt <- final_gt %>%
            dplyr::rowwise() %>%
            do(.addEvidence(., germline_nv=germline_nv, germline_input=germline_input))
        
    }
    final_gt
}
