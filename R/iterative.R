#' Iterative TIgGER 
#' 
#'  
#' Iteratively run the TIgGER pipeline to find novel alleles, infer subject-specific genotypes
#' and use it to correct allele calls. 
#' 
#' The wrapper function \code{itigger} will run \link{findNovelAlleles}, \link{inferGenotype}
#' and \link{reassignAlleles} while new germlines are found and until the maximun 
#' number of iterations specified by the user is reached. With each iteration, 
#' \code{v_call} and \code{germline} are updated so that the next iteration will
#' use as input the corrected \code{v_call} (\code{V_CALL_GENOTYPED}) and the
#' \code{germline} database that includes the new alleles germlines.
#'
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
#'    \item \code{summary} A \code{data.frame} generated with \code{generateEvidence} 
#'                    providing supporting evidence on each inferred allele.
#' }
#'         
#' @seealso \link{findNovelAlleles},  \link{inferGenotype}, \link{reassignAlleles}, 
#'          \link{generateEvidence} and \link{plotTigger}.
#' @examples
#' \dontrun{
#' data(sample_db)
#' data(germline_ighv)

#' # Find novel alleles and return relevant data
#' novel_alleles <- itigger(sample_db, germline_ighv, max.iter=5)
#' novel_alleles$summary
#' }
#' @export
itigger <- function(db, germline, 
                   v_call="V_CALL",
                   fields=NULL, 
                   max.iter=1, 
                   verbose=TRUE,
                   keep_gene="repertoire",
                   ...) {
    
    warning("\n\nThis function is under active development and not fully tested.\n\n", 
            immediate.=TRUE)
    
    gt_cols <- grepl("V_CALL_GENOTYPED(_*)", colnames(db))
    if ( any(gt_cols) ) {
        warning(paste0("Column ", paste(colnames(db)[gt_cols], sep=", "), " removed\n"),
                immediate.=TRUE)
        db[,gt_cols] <- NULL
    }
    
    db[['FIELD_ID']] <- db %>%
        dplyr::group_by_(.dots=fields) %>%
        dplyr::group_indices() %>%
        as.character()
    
    db <- db %>%
        dplyr::ungroup()
    
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
        
        all_nv <- all_gt <- all_germline_input <- list()
        genotyped_alleles <- c()
        
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
                                return(data.frame(NOVEL_IMGT=NA,NOTE=m,stringsAsFactors = FALSE))
                            })
            
            # save nv and gt even if no new alleles
            if (verbose) {
                all_nv[[i_char]] <- nv
            } else {
                all_nv[["1"]] <- nv
            }
            
            # save input germlines, used by generateEvidence
            all_germline_input[[i_char]] <- germline_idx
            
            if (ncol(nv)>2) {
                message("     ... infer genotype")
                gt <- inferGenotype(db_idx, germline_db=germline_idx, 
                                    novel_df=nv, v_call=v_call_idx)
            } else {
                gt <- data.frame(stringsAsFactors = F)
            }
            if (verbose) {
                all_gt[[i_char]] <- gt
            } else {
                all_gt[["1"]] <- gt
            }
            
            selected <- selectNovel(nv)
            
            message(paste0("   |- selectNovel: ", nrow(selected)," rows"))
            if (nrow(selected) > 0) {
            
                genotype_seqs <- genotypeFasta(gt, germline_idx, nv)
                # novel_alleles_found <- any(genotype_seqs %in% germline_idx == F)
                novel_alleles_found <- any(names(genotype_seqs) %in% genotyped_alleles == F)
                # update seen genotyped alleles
                genotyped_alleles <- unique(c(genotyped_alleles, names(genotype_seqs)))
                
                if (novel_alleles_found) {
                    
                    message("     ... reassign alleles")
                    new_v_call <- reassignAlleles(db_idx, genotype_seqs, 
                                                  v_call=v_call_idx, 
                                                  keep_gene = keep_gene)[['V_CALL_GENOTYPED']]
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
        
        # If verbose=TRUE, add column V_CALL_GENOTYPED with last iteration
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
            new_germlines=new_germlines,
            germline_input=all_germline_input
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

    # Arrange by ITERATION for generateEvidence to keep the first time
    # the novel allele was found
    final_gt <- all_gt %>%
        dplyr::arrange(as.numeric(ITERATION)) %>%
        dplyr::group_by_(.dots=c("FIELD_ID",fields, "GENE")) %>%
        dplyr::filter(duplicated(ALLELES) == FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            ALLELES=strsplit(as.character(ALLELES),","),
            COUNTS=strsplit(as.character(COUNTS),",")) %>%
        tidyr::unnest(ALLELES, COUNTS) %>%
        dplyr::mutate(POLYMORPHISM_CALL=paste0(GENE,"*" ,ALLELES)) %>%
        dplyr::filter(POLYMORPHISM_CALL %in% all_nv$POLYMORPHISM_CALL)
    
    genotyped_db <- bind_rows(lapply(foundAlleles, '[[', "db"))

    # Group by field and iteration
    # (the input germline changes with each iteration, as
    # it is expanded with new alleles found)
    # then generate evidence
    final_gt$FIELD_ITERATION_ID <- final_gt %>%
        dplyr::group_by(FIELD_ID, ITERATION) %>%
        dplyr::group_indices()
    final_gt <- bind_rows(lapply(
        unique(final_gt$FIELD_ITERATION_ID), 
        function(f_i_id) {
            this_gt <- final_gt %>%
                dplyr::filter(FIELD_ITERATION_ID==f_i_id)
            this_field <- this_gt[['FIELD_ID']][1]
            this_iteration <- this_gt[['ITERATION']][1]
            this_germline <- foundAlleles[[this_field]][['germline_input']][[this_iteration]]
            generateEvidence(this_gt, 
                             all_nv %>%
                                 dplyr::filter(FIELD_ID==this_field & ITERATION==this_iteration),
                             germline_nv=all_germ,
                             germline_input = this_germline,
                             db = genotyped_db %>%
                                 dplyr::filter(FIELD_ID==this_field),
                             iteration_id = "ITERATION", fields=c("FIELD_ID", fields))
        })) %>%
        dplyr::select(-FIELD_ITERATION_ID)
    
    list(
         db=genotyped_db,
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
    # For each field and allele
    # subset to first iteration 
    # then visualize evidence
    plot_list <- list(polymorphisms=list(), 
                      genotypes=list(),
                      all_genotypes=NA)

    gt <- tigger_list$gt %>%
        dplyr::group_by_(.dots=c("FIELD_ID", "GENE", "ALLELES")) %>%
        dplyr::filter(ITERATION==min(ITERATION)) %>%
        dplyr::ungroup()
    
    if (nrow(gt) > 0 ){
        FIELD_IDs <- unique(gt[,c("FIELD_ID", tigger_list$fields), drop=F])
        plot_list[["genotypes"]] <- lapply(1:nrow(FIELD_IDs), function(i) {
            plotGenotype(merge(FIELD_IDs[i,], gt, by=c("FIELD_ID", tigger_list$fields)), silent=T)
        })
        names(plot_list[["genotypes"]]) <- apply(FIELD_IDs, 1, paste, collapse="_")
        plot_list[["all_genotypes"]] <- plotGenotype(gt, facet_by = "FIELD_ID", silent=T)
    }
    
    nv <- tigger_list$summary
    
    if (nrow(nv) > 0 ){

        nv <- nv %>%
            dplyr::mutate(ROW_ID=1:n(),
                          LABEL=paste(c(FIELD_ID, tigger_list$fields), collapse="_"),
                          FIELD_ID=as.character(FIELD_ID))
        
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
