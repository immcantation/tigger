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
        mutations <- paste0(diff_idx, ref_imgt[diff_idx],">",
                            replace(novel_imgt[diff_idx], is.na(novel_imgt[diff_idx]),"-"))
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
#' @param j_call        name of the column in \code{data} with J allele calls. 
#'                      Default is \code{j_call}.
#' @param junction      Junction region nucleotide sequence, which includes
#'                      the CDR3 and the two flanking conserved codons. Default
#'                      is \code{junction}.
#' @param fields        character vector of column names used to split the data to 
#'                      identify novel alleles, if any. If \code{NULL} then the data is 
#'                      not divided by grouping variables.
#'  
#' @return
#' Returns the \code{genotype} input \code{data.frame} with the following additional columns 
#' providing supporting evidence for each inferred allele:
#' 
#' \itemize{
#'   \item \code{field_id}: Data subset identifier, defined with the input paramter \code{fields}.
#'   \item A variable number of columns, specified with the input parameter \code{fields}.
#'   \item \code{polymorphism_call}: The novel allele call.
#'   \item \code{novel_imgt}: The novel allele sequence.
#'   \item \code{closest_reference}: The closest reference gene and allele in 
#'         the \code{germline_db} database.
#'   \item \code{closest_reference_imgt}: Sequence of the closest reference gene and 
#'         allele in the \code{germline_db} database.
#'   \item \code{germline_call}: The input (uncorrected) V call.
#'   \item \code{germline_imgt}: Germline sequence for \code{germline_call}.
#'   \item \code{nt_diff}: Number of nucleotides that differ between the new allele and
#'         the closest reference (\code{closest_reference}) in the \code{germline_db} database.
#'   \item \code{nt_substitutions}: A comma separated list of specific nucleotide 
#'         differences (e.g. \code{112G>A}) in the novel allele.
#'   \item \code{aa_diff}: Number of amino acids that differ between the new allele and the closest 
#'         reference (\code{closest_reference}) in the \code{germline_db} database.
#'   \item \code{aa_substitutions}: A comma separated list with specific amino acid 
#'         differences (e.g. \code{96A>N}) in the novel allele.
#'   \item \code{sequences}: Number of sequences unambiguosly assigned to this allele.
#'   \item \code{unmutated_sequences}: Number of records with the unmutated novel allele sequence.
#'   \item \code{unmutated_frequency}: Proportion of records with the unmutated novel allele 
#'         sequence (\code{unmutated_sequences / sequences}).
#'   \item \code{allelic_percentage}: Percentage at which the (unmutated) allele is observed 
#'         in the sequence dataset compared  to other (unmutated) alleles.
#'   \item \code{unique_js}: Number of unique J sequences found associated with the 
#'         novel allele. The sequences are those who have been unambiguously assigned 
#'         to the novel allelle (\code{polymorphism_call}).
#'   \item \code{unique_cdr3s}: Number of unique CDR3s associated with the inferred allele.
#'         The sequences are those who have been unambiguously assigned to the 
#'         novel allelle (polymorphism_call).
#'   \item \code{mut_min}: Minimum mutation considered by the algorithm.
#'   \item \code{mut_max}: Maximum mutation considered by the algorithm.
#'   \item \code{pos_min}: First position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{pos_max}: Last position of the sequence considered by the algorithm (IMGT numbering).
#'   \item \code{y_intercept}: The y-intercept above which positions were considered 
#'         potentially polymorphic.
#'   \item \code{alpha}: Significance threshold to be used when constructing the 
#'         confidence interval for the y-intercept.
#'   \item \code{min_seqs}: Input \code{min_seqs}. The minimum number of total sequences 
#'         (within the desired mutational range and nucleotide range) required 
#'         for the samples to be considered.
#'   \item \code{j_max}: Input \code{j_max}. The maximum fraction of sequences perfectly 
#'         aligning to a potential novel allele that are allowed to utilize to a particular 
#'         combination of junction length and J gene.
#'   \item \code{min_frac}: Input \code{min_frac}. The minimum fraction of sequences that must
#'         have usable nucleotides in a given position for that position to be considered.
#'   \item \code{note}: Comments regarding the novel allele inferrence.
#' }
#' 
#' @seealso 
#' See \link{findNovelAlleles}, \link{inferGenotype} and \link{genotypeFasta} 
#' for generating the required input.
#' 
#' @examples
#' \donttest{
#' # Generate input data
#' novel <- findNovelAlleles(AIRRDb, SampleGermlineIGHV,
#'     v_call="v_call", j_call="j_call", junction="junction", 
#'     junction_length="junction_length", seq="sequence_alignment")
#' genotype <- inferGenotype(AIRRDb, find_unmutated=TRUE, 
#'                           germline_db=SampleGermlineIGHV,
#'                           novel=novel,
#'                           v_call="v_call", seq="sequence_alignment")
#' genotype_db <- genotypeFasta(genotype, SampleGermlineIGHV, novel)
#' data_db <- reassignAlleles(AIRRDb, genotype_db, 
#' v_call="v_call", seq="sequence_alignment")
#' 
#' # Assemble evidence table
#' evidence <- generateEvidence(data_db, novel, genotype, 
#'                              genotype_db, SampleGermlineIGHV,
#'                              j_call = "j_call", 
#'                              junction = "junction")
#' }
#' 
#' @export
generateEvidence <- function(data, novel, genotype, genotype_db, 
                             germline_db, j_call="j_call", junction="junction",
                             fields=NULL) {
    # Visibility hack
    . <- NULL
    
    # Define set of sequences containing genotype and uncorrected calls
    germline_set <- c(germline_db[!names(germline_db) %in% names(genotype_db)], 
                      genotype_db)
    
    # Find closest reference
    .findClosestReference <- function(seq, allele_calls, ref_germ, 
                                      exclude_self=F, multiple=F) {
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
            if (length(closest_names) > 1 & multiple==FALSE) {
                msg <- paste0("Multiple closest reference found for ", 
                              names(seq),":\n", 
                              paste(closest_names, collapse=","))
                stop(msg)   
            } 
            warning(paste0("Use: ", 
                           paste(closest_names, collapse=","), 
                           " (less mutated positions, not D, same length, same allele)"))
            
        }
        closest_names
    }
    
    # Subset to novel alleles
    unnest_cols <- c("alleles", "counts")
    final_gt <- genotype %>%
        dplyr::group_by(.data$gene) %>%
        dplyr::filter(!duplicated(.data$alleles)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(alleles=strsplit(as.character(.data$alleles), ","),
                      counts=strsplit(as.character(.data$counts), ",")) %>%
        tidyr::unnest(cols=unnest_cols) %>%
        dplyr::mutate(polymorphism_call=paste0(.data$gene, "*" , .data$alleles)) %>%
        dplyr::filter(.data$polymorphism_call %in% novel$polymorphism_call)  %>%
        dplyr::rename(allele="alleles")

        
    # Add info from novel
    final_gt <- dplyr::inner_join(dplyr::rename(final_gt, note_gt="note"), 
                                  novel, 
                                  by=c(fields, "polymorphism_call"))
    
    # Add message if the same novel img sequence found from
    # different starting alleles, these will be novel imgt sequences
    # with more than one polymorphism call
    final_gt <- final_gt %>%
        dplyr::group_by(.data$novel_imgt) %>%
        dplyr::mutate(NUM_CALLS=length(unique(.data$polymorphism_call))) %>%
        dplyr::ungroup()
    idx_mult <- which(final_gt$NUM_CALLS > 1)
    final_gt$NUM_CALLS <- NULL
    if (length(idx_mult) > 0) {
        final_gt$note_gt[idx_mult] <- paste(
            final_gt$note_gt[idx_mult],
            " Found multiple polymorphism calls for the same novel_imgt.", 
            sep="")
    }
    
    
    if (nrow(final_gt)>0) {
        
        .addEvidence <- function(df, germline_set, germline_db) { 
            polymorphism <- df[['polymorphism_call']]
            novel_imgt <- df[["novel_imgt"]]
            names(novel_imgt) <- polymorphism

            v_call_genotyped <- data[["v_call_genotyped"]]
            
            sequences <- sum(v_call_genotyped == polymorphism)
            df[["sequences"]] <- sequences
            closest_ref_input <- .findClosestReference(novel_imgt,
                                                       names(germline_db), 
                                                       germline_db,
                                                       exclude_self=F)
            closest_ref <- .findClosestReference(novel_imgt,
                                                 names(germline_set), 
                                                 germline_set,
                                                 exclude_self=F, multiple=T)
            
            if (all(getGene(closest_ref_input) != getGene(closest_ref))) {
                warning("closest reference gene difference")
            }
            
            if (all(closest_ref != polymorphism)) {
                warning(paste0("closest reference allele (",
                               closest_ref
                               ,") different from polymorphism_call allele (",
                               polymorphism,")"))
            }
            
            df[["closest_reference"]] <- closest_ref_input
            
            nt_diff <- unlist(getMutatedPositions(novel_imgt, germline_set[[closest_ref_input]]))
            nt_diff_string <- ""
            if (nchar(novel_imgt) < nchar(germline_set)[[closest_ref_input]]) {
                nt_diff <- c(nt_diff, (nchar(novel_imgt)+1):nchar(germline_set[[closest_ref_input]]))
            }
            if (length(nt_diff) > 0 ) {
                ref_nt <- strsplit(germline_set[[closest_ref_input]],"")[[1]][nt_diff]
                novel_nt <- strsplit(germline_set[[polymorphism]],"")[[1]][nt_diff]
                nt_diff_string <- paste(paste(
                    nt_diff, 
                    ref_nt, 
                    ">",
                    replace(novel_nt, is.na(novel_nt), "-"),
                    sep=""), collapse=",")    
            } 
            
            df[["nt_diff"]] <- length(nt_diff)
            df[["nt_substitutions"]] <- nt_diff_string
            
            diff_aa <- getMutatedAA(germline_set[[closest_ref_input]], germline_set[[polymorphism]])
            
            if (length(diff_aa)>0) {
                df[["aa_diff"]] <- length(diff_aa)
                df[["aa_substitutions"]] <- paste(diff_aa,collapse=",")
            } else {
                df[["aa_diff"]] <- 0
                df[["aa_substitutions"]] <- ""
            }
            
            df[["counts"]] <- as.numeric(df[["counts"]])
            df[["total"]] <- as.numeric(df[["total"]])
            df[["unmutated_sequences"]] <- as.numeric(df[["counts"]])
            df[["unmutated_frequency"]] <- as.numeric(df[["counts"]])/sequences
            
            df[["allelic_percentage"]] <- 100*df[["unmutated_sequences"]]/as.numeric(df[["total"]])
            
            if (sequences > 0) {
                df[["unique_js"]] <- data %>%
                    dplyr::filter(.data$v_call_genotyped == polymorphism)  %>%
                    dplyr::distinct(.data[[j_call]]) %>% 
                    nrow()
                df[["unique_cdr3s"]] <- data %>%
                    dplyr::filter(.data$v_call_genotyped == polymorphism)  %>%
                    dplyr::distinct(translateDNA(.data[[junction]], trim=TRUE)) %>% 
                    nrow()
            } else {
                df[["unique_js"]]  <- NA
                df[["unique_cdr3s"]] <- NA
            }

            # Add closest germline
            df[["closest_reference_imgt"]] <- cleanSeqs(germline_set[[closest_ref_input]])
            
            data.frame(df, stringsAsFactors=FALSE)
        }
        
        final_gt <- final_gt %>%
            dplyr::rowwise() %>%
            do(.addEvidence(., germline_set=germline_set, germline_db=germline_db)) %>%
            dplyr::mutate(note=trimws(paste(.data$note_gt, .data$note, sep=" "))) %>%
            dplyr::select(-c("note_gt"))
    }
    
    return(final_gt)
}
