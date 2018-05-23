#' Iterative TIgGER 
#' 
#'  
#' Iteratively run the TIgGER pipeline to find novel alleles, infer subject-specific genotypes
#' and use it to correct allele calls. 
#' This wrapper function will run \link{findNovelAlleles}, \link{inferGenotype}
#' and \link{reassignAlleles} while new germlines are found and until the maximun 
#' number of iterations specified by the user is reached. With each iteration, 
#' \code{v_call} and \code{germline} are updated so that the next iteration will
#' use as input the corrected \code{v_call} (\code{V_CALL_GENOTYPED}) and the
#' \code{germline} database that includes the new alleles germlines.
#'
#' \code{itigger} 
#' 
#' @param    db             a \code{data.frame} in Change-O format. See details.
#' @param    germline       a vector of named nucleotide germline sequences
#'                          matching the V calls in \code{db}.
#' @param    v_call         name of the column in db with V allele calls. 
#'                          Default is V_CALL.   
#' @param    fields         name of the column(s) in \code{db} that will be 
#'                          used to split \code{db} in subsets to be 
#'                          analyzed independently.         
#' @param    verbose        whether to keep the results for all iterations (TRUE)
#'                          or just the last one (FALSE)       
#' @param    max.iter       maximum number of iterations
#' @param    ...            additional arguments for \link{findNovelAlleles}
#' 
#' @return Returns a named list with the fields:
#' \itemize{
#'    \item \code{db} The input \code{db} with the additional column
#'                    \code{V_CALL_GENOTYPED}, which has the the correctd \code{V_CALL}.
#'                    If \code{verbose} was set to \code{TRUE}, there will be as many
#'                    \code{V_CALL_GENOTYPED_X} (where X is a number) columns as iterations.
#'    \item \code{fields} The input \code{fields}
#'    \item \code{v_call} The input \code{v_call}
#'    \item \code{nv} A \code{data.frame} with the results of \code{findNovelAlleles} obtained 
#'                    after all iterations.
#'    \item \code{gt} A \code{data.frame} with the results of \code{inferGenotype} obtained 
#'                    after all iterations.
#'    \item \code{new_germlines} A \code{character} vector with the new germlines.
#'    \item \code{germline} A \code{character} with all germlines, old (input) and new.
#'    \item \code{summary} A \code{data.frame} with information on each inferred allele in 
#'                         the last iteration.
#'    \itemize{
#'       \item \code{FIELD_ID} Data subset identifier, defined with the input paramter \code{fields}.
#'       \item A variable number of columns, specified with the input parameter \code{fields}.
#'       \item{ITERATION} The iteration number.
#'       \item{POLYMORPHISM_CALL} The new allele call.
#'       \item{CLOSEST_REFERENCE} The closest reference gene and allele in the input germline 
#'                                database (\code{germline}).
#'       \item{NT_DIFF} Number of nucleotides that differ between the new allele and
#'                      the closest reference (\code{CLOSEST_REFERENCE}) in the input (\code{germline}).
#'       \item{NT_SUBSTITUTIONS} A \code{character) with specific nucleotide differences (e.g. G112A),
#'                               comma separated.
#'       \item{AA_DIFF} Number of aminoacids that differ between the new allele and the closest 
#'                      reference (\code{CLOSEST_REFERENCE}) in the input (\code{germline}).
#'       \item{AA_SUBSTITUTIONS} A \code{character) with specific aminoacid differences (e.g. A96N),
#'                               comma separated.
#'       \item{SEQUENCES} Number of sequences unambiguosly assigned to this allele.
#'       \item{UNMUTATED_SEQUENCES} Number of records with the unmutated new allele sequence.
#'       \item{UNMUTATED_FREQUENCY} Proportion of records with the unmutated new allele 
#'                                  sequence. \code{UNMUTATED_SEQUENCES}/\code{SEQUENCE}
#'       \item{ALLELIC_PERCENTAGE} Percentage at which this (unmutated) allele is observed in the 
#'                                 sequence dataset, compared  to other (unmutated) alleles.
#'       \item{UNIQUE_JS} Number of unique J sequences found associated with the new allele. The sequences
#'                        are those who have been unambiguously assigned the new alelle (\item{POLYMORPHISM_CALL})
#'       \item{UNIQUE_CDR3S}
#'       \item{GERMLINE_CALL}
#'       \item{MUT_MIN}
#'       \item{MUT_MAX}
#'       \item{POS_MIN}
#'       \item{POS_MAX}
#'       \item{Y_INTERCEPT}
#'       \item{ALPHA}
#'       \item{MIN_SEQS}
#'       \item{J_MAX}
#'       \item{MIN_FRAC}
#'       \item{NOVEL_IMGT}
#'       \item{CLOSEST_GERMLINE_IMGT}
#'       \item{GERMLINE_IMGT}
#'       \item{NOTE}
#'    }
#' }
#'         
#' @seealso \link{findNovelAlleles},  \link{inferGenotype}, \link{reassignAlleles}, and \link{plotTigger}.
#' @examples
#' \dontrun{
#' data(sample_db)
#' data(germline_ighv)

#' # Find novel alleles and return relevant data
#' novel_alleles <- itigger(sample_db, germline_ighv, max.iter=Inf)
#' novel_alleles$summary
#' }
#' @export
itigger <- function(db, germline, 
                   v_call="V_CALL",
                   fields=NULL, 
                   max.iter=1, 
                   verbose=TRUE,
                   ...) {
    gt_cols <- grepl("V_CALL_GENOTYPED(_*)", colnames(db))
    if ( any(gt_cols) ) {
        warning(paste0("Column ", paste(colnames(db)[gt_cols], sep=", "), " removed\n"))
        db[,gt_cols] <- NULL
    }
    
    db[['FIELD_ID']] <- db %>%
        dplyr::group_by_(.dots=fields) %>%
        dplyr::group_indices()
    db <- db %>% dplyr::ungroup()
    
    FIELD_ID_label <- db %>%
        dplyr::select_(.dots=c(fields, "FIELD_ID")) %>%
        dplyr::distinct() %>%
        data.frame()
    
    foundAlleles <- list()
    
    for (idx in FIELD_ID_label[["FIELD_ID"]]) {
        
        label <- merge(list(FIELD_ID=idx), FIELD_ID_label)
        message(paste(c("\n>",gsub("^,","",
                                   paste(paste(names(label), 
                                               label, sep="="), 
                                         collapse="-")), collapse=", ")))

        # (re)set germline and v_call
        germline_idx <- germline
        v_call_idx <- v_call
        
        db_idx <- db %>%
            dplyr::filter(FIELD_ID==idx) %>%
            data.frame()
        
        all_nv <- all_gt <- list()
        
        i <- 0
        while (i < max.iter) {
            i <- i+1
            i_char <- as.character(i)
            message(paste0("   Iteration ",i, " of ", max.iter))
            message(paste0("   |- v_call: ",v_call_idx))
            
            nv <- tryCatch( findNovelAlleles(db_idx, germline_db=germline_idx,
                                             v_call=v_call_idx,...),
                            error=function(e){
                                m <- gsub("\\\n","", geterrmessage())
                                message(e)
                                message("")
                                return(data.frame( NOTE=m,stringsAsFactors = FALSE))
                            })
            
            # save nv even if no new alleles
            # this implies nv can be one iteration ahead of genotyping 
            if (verbose) {
                all_nv[[i_char]] <- nv
            } else {
                all_nv[[i_char]] <- nv
            }
            selected <- selectNovel(nv)
            
            message(paste0("   |- selectNovel: ", nrow(selected)," rows"))
            if (nrow(selected) > 0) {
            
                message("     ... infer genotype")
                gt <- inferGenotype(db_idx, germline_db=germline_idx, 
                                    novel_df=nv, v_call=v_call_idx)
                if (verbose) {
                    all_gt[[i_char]] <- gt
                } else {
                    all_gt[["1"]] <- gt
                }
                genotype_seqs <- genotypeFasta(gt, germline_idx, nv)
                novel_alleles_found <- any(genotype_seqs %in% germline_idx == F)
                
                if (novel_alleles_found) {
                    
                    message("     ... reassign alleles")
                    new_v_call <- reassignAlleles(db_idx, genotype_seqs, 
                                                  v_call=v_call_idx)[['V_CALL_GENOTYPED']]
                    ## append new calls to db
                    if (verbose) {
                        v_call_idx <- paste0("V_CALL_GENOTYPED_",i)
                    } else {
                        v_call_idx <- "V_CALL_GENOTYPED"
                    }
                    db_idx[[v_call_idx]] <- new_v_call
                    germline_idx <- c(germline_idx,
                                     genotype_seqs[names(genotype_seqs) %in% names(germline_idx) == F])

                } else {
                    message(paste0("Done. Reason: no new genotypes in iteration ", i))
                    break
                }
            } else {
                message(paste0("Done. Reason: empty selectNovel in iteration ", i))
                break
            }
        }  # end iteration loop    
        
        # If verbose=TRUE, add column V_CALL_GEBOTYPED with last iteration
        genotype_cols <- grep("V_CALL_GENOTYPED_", colnames(db_idx))
        if (length(genotype_cols)>0) {
            iterations <- as.numeric(gsub(".*_","",colnames(db_idx)[genotype_cols]))
            last_genotype_col <- genotype_cols[which.max(iterations)]
            db_idx[["V_CALL_GENOTYPED"]] <- db_idx[[colnames(db_idx)[last_genotype_col]]]
        } 
        
        new_germlines <- names(germline_idx) %in% names(germline) == FALSE
        num_new_germlines <- sum(new_germlines)
        new_germlines <- germline_idx[new_germlines]
        message(paste0(" |- ",num_new_germlines, " new germlines in field ", label))
        all_nv <- bind_rows(all_nv, .id="ITERATION")
        if (nrow(all_nv)>0) {
            all_nv <- cbind(label, all_nv)
        }
        all_gt <- bind_rows(all_gt, .id="ITERATION")
        if (nrow(all_gt)>0) {
            all_gt <- cbind(label, all_gt)
        }
        foundAlleles[[as.character(idx)]] <- list(
            db=db_idx,
            nv=all_nv,
            gt=all_gt,
            germline=germline_idx,
            new_germlines=new_germlines 
        )
        
    } # end fields loop
    
    # Final unique germlines, considering all fields
    all_germ_names <- unlist(lapply(lapply(foundAlleles, '[[', "germline"), names))
    all_germ <- unlist(lapply(foundAlleles, '[[', "germline"), use.names = F)
    names(all_germ) <- all_germ_names
    all_germ <- all_germ[match(unique(names(all_germ)), names(all_germ))]
    
    new_germlines <- all_germ[names(all_germ) %in% names(germline)== FALSE]
    message(paste0("\nTOTAL distinct new germlines after genotyping: ", length(new_germlines)))
    
    all_nv <- bind_rows(lapply(foundAlleles, '[[', "nv"))
    all_gt <- bind_rows(lapply(foundAlleles, '[[', "gt"))
    
    # TODO: check this when 0 novel alleles found
    # Subset to last iteration with novel alleles
    final_gt <- all_gt %>%
        dplyr::group_by_(.dots=c("FIELD_ID")) %>%
        dplyr::filter(duplicated(ALLELES) == FALSE) %>%
        dplyr::filter(ITERATION==max(ITERATION)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            ALLELES=strsplit(as.character(ALLELES),","),
            COUNTS=strsplit(as.character(COUNTS),",")) %>%
        tidyr::unnest(ALLELES, COUNTS) %>%
        dplyr::filter(grepl("_",ALLELES)) %>%
        dplyr::rename(ALLELE=ALLELES) %>%
        dplyr::mutate(POLYMORPHISM_CALL=paste0(GENE,"*" ,ALLELE))
    
    final_gt <- merge(final_gt %>% 
                          dplyr::rename(NOTE_GT=NOTE), 
                      all_nv, 
                      by=c("FIELD_ID", fields, "ITERATION", "POLYMORPHISM_CALL")) 
    
    # Find closest reference in original germline
    findClosestReference <- function(seq, allele_calls) {
        closest <- getMutCount(seq,
                               paste(allele_calls, collapse=","),
                               all_germ)
        min_dist <- min(unlist(closest))
        closest_idx <- which(unlist(closest) == min_dist)
        if (length(closest_idx) > 1) {
            stop("More than one closest reference found")
        }
        paste(unique(allele_calls[closest_idx]), collapse=",")

    }
    
    
    # Add additional evidence to final table
    addEvidence <- function(dfr) {
        for (i in 1:nrow(dfr)) {
            this_field <- as.character(dfr[["FIELD_ID"]][i])
            polymorphism <- dfr[['POLYMORPHISM_CALL']][i]
            novel_imgt <- dfr[["NOVEL_IMGT"]][i]
            gene <- dfr[['GENE']][i]
            allele <- dfr[['ALLELE']][i]
            germline_call <- dfr[['GERMLINE_CALL']][i]
            this_germline <- foundAlleles[[this_field]][['germline']]
            V_CALL_GENOTYPED <- foundAlleles[[this_field]][['db']][["V_CALL_GENOTYPED"]]
            
            SEQUENCES <- sum(V_CALL_GENOTYPED==polymorphism)
            dfr[["SEQUENCES"]][i] <- SEQUENCES

            closest_ref_input <- findClosestReference(novel_imgt,
                                                names(germline))
            closest_ref <- findClosestReference(novel_imgt,
                                                names(all_germ))
            
            if (getGene(closest_ref_input) != getGene(closest_ref)) {
                warning("closest reference gene difference")
            }
            
            if (closest_ref != polymorphism) {
                warning(paste0("closest reference (",
                               getAllele(closest_ref)
                               ,") different from POLYMORPHISM_CALL (",
                               getAllele(polymorphism),")"))
            }
            
            ## TODO: this still not clear.
            ## Any diff using sequence_imgt instead of germline[[polymorphism]]?
            dfr[["CLOSEST_REFERENCE"]][i] <- closest_ref_input
            
            nt_diff <- unlist(getMutatedPositions(novel_imgt, all_germ[[closest_ref_input]]))
            nt_diff_string <- paste(paste(
                nt_diff, 
                strsplit(all_germ[[closest_ref_input]],"")[[1]][nt_diff], 
                ">",
                strsplit(all_germ[[polymorphism]],"")[[1]][nt_diff],
                sep=""), collapse=",")
                
            dfr[["NT_DIFF"]][i] <- length(nt_diff)
            dfr[["NT_SUBSTITUTIONS"]][i] <- nt_diff_string
            
            aa_substitutions <- calcObservedMutations(
                                        all_germ[[polymorphism]], 
                                        all_germ[[closest_ref_input]],
                                        regionDefinition = NULL,
                                        mutationDefinition = NULL,
                                        returnRaw = T)
            
            pos_R <- aa_substitutions$pos %>%
                dplyr::filter(R > 0) %>%
                dplyr::select(position) %>%
                c()
            
            dfr[["AA_DIFF"]][i] <- length(pos_R$position)

            
            poly_aa <- strsplit(translateDNA(all_germ[[polymorphism]]),"")[[1]]
            germ_aa <- strsplit(translateDNA(all_germ[[closest_ref_input]]),"")[[1]]
            
            dfr[["AA_SUBSTITUTIONS"]][i] <- paste(paste(
                                                    pos_R$position/3, 
                                                    germ_aa[pos_R$position/3], 
                                                    ">",
                                                    poly_aa[pos_R$position/3],
                                                    sep=""), collapse=",")
            dfr[["UNMUTATED_SEQUENCES"]][i] <- as.numeric(dfr[["COUNTS"]][i])
            dfr[["UNMUTATED_FREQUENCY"]][i] <- as.numeric(dfr[["COUNTS"]][i])/SEQUENCES
            
            dfr[["ALLELIC_PERCENTAGE"]][i] <- 100*dfr[["UNMUTATED_SEQUENCES"]][i]/as.numeric(dfr[["TOTAL"]][i])
            
            dfr[["UNIQUE_JS"]][i] <- foundAlleles[[this_field]]$db %>%
                dplyr::filter(V_CALL_GENOTYPED==polymorphism)  %>%
                dplyr::distinct(J_CALL) %>% nrow()
            dfr[["UNIQUE_CDR3S"]][i] <- foundAlleles[[this_field]]$db %>%
                dplyr::filter(V_CALL_GENOTYPED==polymorphism)  %>%
                dplyr::distinct(translateDNA(JUNCTION, trim=TRUE)) %>% 
                nrow()
            # Add closest germline
            dfr[["CLOSEST_GERMLINE_IMGT"]][i] <- cleanSeqs(all_germ[[closest_ref_input]])
        }
        dfr
    }
    
    final_gt <- addEvidence(final_gt)
    
    ## Filter, reorder,...
    final_gt <- final_gt %>%
        dplyr::select_(.dots=c(
            "FIELD_ID", fields, "ITERATION", "POLYMORPHISM_CALL",
            "CLOSEST_REFERENCE",
            "NT_DIFF", "NT_SUBSTITUTIONS",
            "AA_DIFF", "AA_SUBSTITUTIONS",
            "SEQUENCES", "UNMUTATED_SEQUENCES", "UNMUTATED_FREQUENCY",
            "ALLELIC_PERCENTAGE",
            "UNIQUE_JS", "UNIQUE_CDR3S",
            "GERMLINE_CALL", "MUT_MIN", "MUT_MAX", "POS_MIN", "POS_MAX",
            "Y_INTERCEPT", "ALPHA", "MIN_SEQS", "J_MAX", "MIN_FRAC",
            "NOVEL_IMGT", "CLOSEST_GERMLINE_IMGT", "GERMLINE_IMGT",
            "NOTE")
        ) 

    list(
         db=bind_rows(lapply(foundAlleles, '[[', "db")),
         fields=fields,
         v_call=v_call,
         nv=all_nv,
         gt=all_gt,
         new_germlines=new_germlines,
         germline=lapply(foundAlleles, '[[', "germline"),
         summary=final_gt
    )        
}


#' Visualize genotypes and evidence of novel V alleles
#'
#' \code{plotTigger} takes the output of \link{itigger} and uses
#' \link{findNovelAlleles} and \link{plotGenotype} to visualize genotypes
#' and evidence of the final novel V alleles (\code{tigger_list$summary}).
#' 
#' @return a list of length 2, with elements named \code{polymorphisms} and
#'         \code{genotypes}. \code{polymorphisms} contains a list of figures 
#'         generated by \link{plotNovel}. The names in the list start with a 
#'         number that matches the row from \code{tigger_list$summary} used to make
#'         the figure. The following number and characters are delimit groupings
#'         generated with \code{tigger_list$fields}. 
#'         \code{genotypes} contains a list of figures generated 
#'         with \link{plotGenotype}. As before, the names are made 
#'         from the groupings created with \code{tigger_list$fields}. 
#' 
#' @param    tigger_list    a \code{list} generated with \link{itigger}
#' 
#' @examples
#' \dontrun{
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Find novel alleles and return relevant data
#' novel_alleles <- itigger(sample_db, germline_ighv, max.iter=2)
#' # Plot the evidence for the first (and only) novel allele in the example data
#' tigger_plots <- plotTigger(novel_alleles)
#' # plots is a list of 2 elements, with one plot each in this example.
#' tigger:::multiplot(tigger_plots$polymorphisms[[1]][[1]], tigger_plots$genotypes[[1]], cols=2)
#' }
#' @export
plotTigger <- function(tigger_list) {
    # For each field
    # subset to last iteration 
    # then visualize evidence
    plot_list <- list(polymorphisms=list(), 
                      genotypes=list())

    gt <- tigger_list$gt %>%
        dplyr::group_by_(.dots=c("FIELD_ID")) %>%
        dplyr::filter(ITERATION==max(ITERATION))
    
    if (nrow(gt) > 0 ){
        FIELD_IDs <- unique(gt[,c("FIELD_ID", tigger_list$fields), drop=F])
        plot_list[["genotypes"]] <- lapply(1:nrow(FIELD_IDs), function(i) {
            plotGenotype(merge(FIELD_IDs[i,], gt, by=c("FIELD_ID", tigger_list$fields)), silent=T)
        })
        names(plot_list[["genotypes"]]) <- apply(FIELD_IDs, 1, paste, collapse="_")
    }
    
    
    nv <- tigger_list$summary %>%
        dplyr::mutate(ROW_ID=1:n(),
                      LABEL=paste(c(FIELD_ID, tigger_list$fields), collapse="_"))
    if (nrow(nv) > 0 ){
        for (i in 1:nrow(nv)) {
            this_field <- nv[['LABEL']][i]
            this_db <- merge(nv[i,], tigger_list[['db']], by="FIELD_ID")
            row_id <- as.character(nv[["ROW_ID"]][i])
            plot_list[["polymorphisms"]][[this_field]][[row_id]] <- plotNovel(this_db,
                      nv[i,], v_call=tigger_list$v_call)
        }
     
    }
    
    plot_list
}
