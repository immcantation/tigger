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
#' \code{generateEvidence} builds a table of evidence metrics for the final novel V 
#' allele detection and genotyping inferrences.
#' 
#' @param data          a \code{data.frame} containing sequence data that has been
#'                      passed through \link{reassignAlleles} to correct the allele 
#'                      assignments.
#' @param novel         the \code{data.frame} returned by \link{findNovelAlleles}.
#' @param genotype      the \code{data.frame} of alleles generated with \link{inferGenotype} 
#'                      denoting the genotype of the subject.
#' @param genotype_db   a vector of named nucleotide germline sequences in the genotype.
#'                      Returned by \link{genotypeFasta}. 
#' @param germline_db   the original uncorrected germline database used to by
#'                      \link{findNovelAlleles} to identify novel alleles.
#' @param fields        character vector of column names used to split the data to 
#'                      identify novel alleles, if any. If \code{NULL} then the data is 
#'                      not divided by grouping variables.
#'  
#' @return
#' Returns the \code{genotype} input \code{data.frame} with the following additional columns 
#' providing supporting evidence for each inferred allele:
#' 
#' \itemize{
#'   \item \code{FIELD_ID}: Data subset identifier, defined with the input paramter \code{fields}.
#'   \item A variable number of columns, specified with the input parameter \code{fields}.
#'   \item \code{POLYMORPHISM_CALL}: The novel allele call.
#'   \item \code{NOVEL_IMGT}: The novel allele sequence.
#'   \item \code{CLOSEST_REFERENCE}: The closest reference gene and allele in 
#'         the \code{germline_db} database.
#'   \item \code{CLOSEST_REFERENCE_IMGT}: Sequence of the closest reference gene and 
#'         allele in the \code{germline_db} database.
#'   \item \code{GERMLINE_CALL}: The input (uncorrected) V call.
#'   \item \code{GERMLINE_IMGT}: Germline sequence for \code{GERMLINE_CALL}.
#'   \item \code{NT_DIFF}: Number of nucleotides that differ between the new allele and
#'         the closest reference (\code{CLOSEST_REFERENCE}) in the \code{germline_db} database.
#'   \item \code{NT_SUBSTITUTIONS}: A comma separated list of specific nucleotide 
#'         differences (e.g. \code{112G>A}) in the novel allele.
#'   \item \code{AA_DIFF}: Number of amino acids that differ between the new allele and the closest 
#'         reference (\code{CLOSEST_REFERENCE}) in the \code{germline_db} database.
#'   \item \code{AA_SUBSTITUTIONS}: A comma separated list with specific amino acid 
#'         differences (e.g. \code{96A>N}) in the novel allele.
#'   \item \code{SEQUENCES}: Number of sequences unambiguosly assigned to this allele.
#'   \item \code{UNMUTATED_SEQUENCES}: Number of records with the unmutated novel allele sequence.
#'   \item \code{UNMUTATED_FREQUENCY}: Proportion of records with the unmutated novel allele 
#'         sequence (\code{UNMUTATED_SEQUENCES / SEQUENCE}).
#'   \item \code{ALLELIC_PERCENTAGE}: Percentage at which the (unmutated) allele is observed 
#'         in the sequence dataset compared  to other (unmutated) alleles.
#'   \item \code{UNIQUE_JS}: Number of unique J sequences found associated with the 
#'         novel allele. The sequences are those who have been unambiguously assigned 
#'         to the novel allelle (\code{POLYMORPHISM_CALL}).
#'   \item \code{UNIQUE_CDR3S}: Number of unique CDR3s associated with the inferred allele.
#'   \item \code{MUT_MIN}: Minimum mutation considered by the algorithm.
#'   \item \code{MUT_MAX}: Maximum mutation considered by the algorithm.
#'   \item \code{POS_MIN}: First position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{POS_MAX}: Last position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{Y_INTERCEPT}: The y-intercept above which positions were considered 
#'         potentially polymorphic.
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
#'   \item \code{NOTE}: Comments regarding the novel allele inferrence.
#' }
#' 
#' @seealso 
#' See \link{findNovelAlleles}, \link{inferGenotype} and \link{genotypeFasta} 
#' for generating the required input.
#' 
#' @examples
#' \donttest{
#' # Generate input data
#' novel <- findNovelAlleles(SampleDb, GermlineIGHV)
#' genotype <- inferGenotype(SampleDb, find_unmutated=TRUE, germline_db=GermlineIGHV,
#'                           novel=novel)
#' genotype_db <- genotypeFasta(genotype, GermlineIGHV, novel)
#' data_db <- reassignAlleles(SampleDb, genotype_db)
#' 
#' # Assemble evidence table
#' evidence <- generateEvidence(data_db, novel, genotype, genotype_db, GermlineIGHV)
#' }
#' 
#' @export
generateEvidence <- function(data, novel, genotype, genotype_db, 
                             germline_db, fields=NULL) {
    # Visibility hack
    . <- NULL
    
    # Define set of sequences containing genotype and uncorrected calls
    germline_set <- c(germline_db[!names(germline_db) %in% names(genotype_db)], 
                      genotype_db)
    
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
    final_gt <- genotype %>%
        dplyr::group_by(.data$GENE) %>%
        dplyr::filter(!duplicated(.data$ALLELES)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(ALLELES=strsplit(as.character(.data$ALLELES), ","),
                      COUNTS=strsplit(as.character(.data$COUNTS), ",")) %>%
        tidyr::unnest(.data$ALLELES, .data$COUNTS) %>%
        dplyr::mutate(POLYMORPHISM_CALL=paste0(.data$GENE, "*" , .data$ALLELES)) %>%
        dplyr::filter(.data$POLYMORPHISM_CALL %in% novel$POLYMORPHISM_CALL)  %>%
        dplyr::rename(ALLELE="ALLELES")

        
    # Add info from novel
    final_gt <- dplyr::inner_join(dplyr::rename(final_gt, NOTE_GT="NOTE"), 
                                  novel, 
                                  by=c(fields, "POLYMORPHISM_CALL"))
    
    # Add message if the same novel img sequence found from
    # different starting alleles, these will be novel imgt sequences
    # with more than one polymorphism call
    final_gt <- final_gt %>%
        dplyr::group_by(.data$NOVEL_IMGT) %>%
        dplyr::mutate(NUM_CALLS=length(unique(.data$POLYMORPHISM_CALL))) %>%
        dplyr::ungroup()
    idx_mult <- which(final_gt$NUM_CALLS > 1)
    final_gt$NUM_CALLS <- NULL
    if (length(idx_mult) > 0) {
        final_gt$NOTE_GT[idx_mult] <- paste(
            final_gt$NOTE_GT[idx_mult],
            " Found multiple polymorphism calls for the same NOVEL_IMGT.", 
            sep="")
    }
    
    
    if (nrow(final_gt)>0) {
        
        .addEvidence <- function(df, germline_set, germline_db) { 
            polymorphism <- df[['POLYMORPHISM_CALL']]
            novel_imgt <- df[["NOVEL_IMGT"]]
            names(novel_imgt) <- polymorphism
            #gene <- df[['GENE']]
            #allele <- df[['ALLELE']]
            #germline_call <- df[['GERMLINE_CALL']]
            #this_germline <- germline_db
            v_call_genotyped <- data[["V_CALL_GENOTYPED"]]
            
            SEQUENCES <- sum(v_call_genotyped == polymorphism)
            df[["SEQUENCES"]] <- SEQUENCES
            closest_ref_input <- .findClosestReference(novel_imgt,
                                                       names(germline_db), 
                                                       germline_db,
                                                       exclude_self=F)
            closest_ref <- .findClosestReference(novel_imgt,
                                                 names(germline_set), 
                                                 germline_set,
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
            
            nt_diff <- unlist(getMutatedPositions(novel_imgt, germline_set[[closest_ref_input]]))
            nt_diff_string <- ""
            if (length(nt_diff) > 0 ) {
                nt_diff_string <- paste(paste(
                    nt_diff, 
                    strsplit(germline_set[[closest_ref_input]],"")[[1]][nt_diff], 
                    ">",
                    strsplit(germline_set[[polymorphism]],"")[[1]][nt_diff],
                    sep=""), collapse=",")    
            } 
            
            df[["NT_DIFF"]] <- length(nt_diff)
            df[["NT_SUBSTITUTIONS"]] <- nt_diff_string
            
            diff_aa <- getMutatedAA(germline_set[[closest_ref_input]], germline_set[[polymorphism]])
            
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
            
            df[["UNIQUE_JS"]] <- data %>%
                dplyr::filter(.data$V_CALL_GENOTYPED == polymorphism)  %>%
                dplyr::distinct(.data$J_CALL) %>% 
                nrow()
            df[["UNIQUE_CDR3S"]] <- data %>%
                dplyr::filter(.data$V_CALL_GENOTYPED == polymorphism)  %>%
                dplyr::distinct(translateDNA(.data$JUNCTION, trim=TRUE)) %>% 
                nrow()
            # Add closest germline
            df[["CLOSEST_REFERENCE_IMGT"]] <- cleanSeqs(germline_set[[closest_ref_input]])
            
            data.frame(df, stringsAsFactors=FALSE)
        }
        
        final_gt <- final_gt %>%
            dplyr::rowwise() %>%
            do(.addEvidence(., germline_set=germline_set, germline_db=germline_db)) %>%
            dplyr::mutate(NOTE=trimws(paste(.data$NOTE_GT, .data$NOTE, sep=" "))) %>%
            dplyr::select(-c("NOTE_GT"))
    }
    
    return(final_gt)
}
