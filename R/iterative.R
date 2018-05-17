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
#' @seealso \link{findNovelAlleles},  \link{inferGenotype} and  \link{reassignAlleles}
#' @examples
#' \dontrun{
#' data(sample_db)
#' data(germline_ighv)

#' # Find novel alleles and return relevant data
#' novel_alleles <- itiger(sample_db, germline_ighv, max.iter=Inf)
#' novel_alleles$final
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
                all_nv[[i]] <- nv
            } else {
                all_nv[[1]] <- nv
            }
            selected <- selectNovel(nv)
            
            message(paste0("   |- selectNovel: ", nrow(selected)," rows"))
            if (nrow(selected) > 0) {
            
                message("     ... infer genotype")
                gt <- inferGenotype(db_idx, germline_db=germline_idx, 
                                    novel_df=nv, v_call=v_call_idx)
                if (verbose) {
                    all_gt[[i]] <- gt
                } else {
                    all_gt[[1]] <- gt
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
    
    final_gt <- all_gt %>%
        dplyr::group_by_(.dots=c("FIELD_ID")) %>%
        dplyr::filter(ITERATION==max(ITERATION)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            ALLELES=strsplit(as.character(ALLELES),","),
            COUNTS=strsplit(as.character(COUNTS),",")) %>%
        tidyr::unnest(ALLELES, COUNTS) %>%
        dplyr::filter(grepl("_",ALLELES)) %>%
        dplyr::rename(ALLELE=ALLELES) %>%
        dplyr::mutate(POLYMORPHISM_CALL=paste0(GENE,"*" ,ALLELE))
    
    # add additional evidence to final table
    addEvidence <- function(dfr) {
        for (i in 1:nrow(dfr)) {
            this_field <- as.character(dfr[["FIELD_ID"]][i])
            polymorphism <- dfr[['POLYMORPHISM_CALL']][i]
            gene <- dfr[['GENE']][i]
            allele <- dfr[['ALLELE']]
            germline <- foundAlleles[[this_field]][['germline']]
            V_CALL_GENOTYPED <- foundAlleles[[this_field]][['db']][["V_CALL_GENOTYPED"]]
            SEQUENCES_AMB <- sum(grepl(polymorphism, V_CALL_GENOTYPED, fixed = T))
            dfr[["SEQUENCES_AMB"]][i] <- SEQUENCES_AMB
            SEQUENCES <- sum(V_CALL_GENOTYPED==polymorphism)
            GENE_COUNT <- sum(grepl(gene, V_CALL_GENOTYPED, fixed = T))
            ALLELIC_PERCENTAGE <- 100*SEQUENCES/GENE_COUNT
            ALLELIC_AMB_PERCENTAGE <- 100*SEQUENCES_AMB/GENE_COUNT  
            dfr[["SEQUENCES"]][i] <- SEQUENCES
            dfr[["SEQUENCES_AMB"]][i] <- SEQUENCES_AMB
            TOTAL_AMB <- sum(grepl(gene, V_CALL_GENOTYPED, fixed = T))    
            dfr[["TOTAL_AMB"]][i] <- TOTAL_AMB
            # TODO: check this works if closest is a new allele
            closest_ref <- strsplit(polymorphism, "_")[[1]][1]
            NT_DIFF <- length(strsplit(polymorphism, "_")[[1]])-1
            dfr[["NT_DIFF"]][i] <- NT_DIFF
            ## TODO: calcObservedMutations
            dfr[["AA_DIF"]][i] <- sum(calcObservedMutations(germline[[polymorphism]], 
                                   germline[[closest_ref]],
                                   regionDefinition = NULL,
                                   mutationDefinition = NULL))
            dfr[["AA_SUBSTITUTIONS"]][i] <- calcObservedMutations(germline[[polymorphism]], 
                                                         germline[[closest_ref]],
                                                         regionDefinition = NULL,
                                                         mutationDefinition = NULL)[["SEQ_S"]]
        }
        dfr
    }
    final_gt <- addEvidence(final_gt)
    final_gt <- merge(final_gt %>% 
              dplyr::rename(NOTE_GT=NOTE), 
          all_nv, 
          by=c("FIELD_ID", fields, "ITERATION", "POLYMORPHISM_CALL")) 
    
    list(
         db=bind_rows(lapply(foundAlleles, '[[', "db")),
         fields=fields,
         v_call=v_call,
         # db_id=digest(db),
         nv=all_nv,
         gt=all_gt,
         new_germlines=new_germlines,
         germline=lapply(foundAlleles, '[[', "germline"),
         final=final_gt
    )        
}


#' Visualize genotypes and evidence of novel V alleles
#'
#' \code{plotTigger} takes the output of \link{itiger} and uses
#' \link{findNovelAlleles} and \link{plotGenotype} to visualize genotypes
#' and evidence of the final novel V alleles (\code{tigger_list$final}).
#' 
#' @return a list of length 2, with elements named \code{polymorphisms} and
#'         \code{genotypes}. \code{polymorphisms} contains a list of figures 
#'         generated by \link{plotNovel}. The names in the list start with a 
#'         number that matches the row from \code{tigger_list$final} used to make
#'         the figure. The following number and characters are delimit groupings
#'         generated with \code{tigger_list$fields}. 
#'         \code{genotypes} contains a list of figures generated 
#'         with \link{plotGenotype}. As before, the names are made 
#'         from the groupings created with \code{tigger_list$fields}. 
#' 
#' @param    tigger_list    a \code{list} generated with \link{itiger}
#' 
#' @examples
#' \dontrun{
#' # Load example data and germlines
#' data(sample_db)
#' data(germline_ighv)
#' 
#' # Find novel alleles and return relevant data
#' novel_alleles <- itigger(sample_db, germline_ighv)
#' # Plot the evidence for the first (and only) novel allele in the example data
#' plots <- plotTigger(novel_alleles)
#' # plots is a list of 2 elements, with one plot each (in this example)
#' tigger:::multiplot(plots$polymorphisms[[1]], plots$genotypes[[1]], cols=2)
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
    
    
    nv <- tigger_list$final %>%
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
