#' Infer a subject-specific genotype using a Bayesian approach
#'
#' \code{inferGenotypeBayesian} infers an subject's genotype by applying a Bayesian framework
#' with a Dirichlet prior for the multinomial distribution. Up to four distinct alleles are
#' allowed in an individual’s genotype. Four likelihood distributions were generated by
#' empirically fitting three high coverage genotypes from three individuals
#' (Laserson and Vigneault et al, 2014). A posterior probability is calculated for the
#' four most common alleles. The certainty of the highest probability model was
#' calculated using a Bayes factor (the most likely model divided by second-most likely model).
#' The larger the Bayes factor (K), the greater the certainty in the model.
#'
#' @details
#' Allele calls representing cases where multiple alleles have been
#' assigned to a single sample sequence are rare among unmutated
#' sequences but may result if nucleotides for certain positions are
#' not available. Calls containing multiple alleles are treated as
#' belonging to all groups. If \code{novel} is provided, all
#' sequences that are assigned to the same starting allele as any
#' novel germline allele will have the novel germline allele appended
#' to their assignment prior to searching for unmutated sequences.
#'
#' @param    data            a \code{data.frame} containing V allele
#'                           calls from a single subject. If \code{find_unmutated}
#'                           is \code{TRUE}, then the sample IMGT-gapped V(D)J sequence
#'                           should be provided in column \code{sequence_alignment}
#' @param    v_call          column in \code{data} with V allele calls.
#'                           Default is \code{"v_call"}.
#' @param    seq             name of the column in \code{data} with the
#'                           aligned, IMGT-numbered, V(D)J nucleotide sequence.
#'                           Default is \code{"sequence_alignment"}.
#' @param    find_unmutated  if \code{TRUE}, use \code{germline_db} to
#'                           find which samples are unmutated. Not needed
#'                           if \code{allele_calls} only represent
#'                           unmutated samples.
#' @param    germline_db     named vector of sequences containing the
#'                           germline sequences named in \code{allele_calls}.
#'                           Only required if \code{find_unmutated} is \code{TRUE}.
#' @param    novel           an optional \code{data.frame} of the type
#'                           novel returned by \link{findNovelAlleles} containing
#'                           germline sequences that will be utilized if
#'                           \code{find_unmutated} is \code{TRUE}. See Details.
#' @param    priors          a numeric vector of priors for the multinomial distribution.
#'                           The \code{priors} vector must be nine values that defined
#'                           the priors for the heterozygous (two allele),
#'                           trizygous (three allele), and quadrozygous (four allele)
#'                           distributions. The first two values of \code{priors} define
#'                           the prior for the heterozygous case, the next three values are for
#'                           the trizygous case, and the final four values are for the
#'                           quadrozygous case. Each set of priors should sum to one.
#'                           Note, each distribution prior is actually defined internally
#'                           by set of four numbers, with the unspecified final values
#'                           assigned to \code{0}; e.g., the heterozygous case is
#'                           \code{c(priors[1], priors[2], 0, 0)}. The prior for the
#'                           homozygous distribution is fixed at \code{c(1, 0, 0, 0)}.
#'
#' @return
#' A \code{data.frame} of alleles denoting the genotype of the subject with the log10
#' of the likelihood of each model and the log10 of the Bayes factor. The output
#' contains the following columns:
#'
#' \itemize{
#'   \item \code{gene}: The gene name without allele.
#'   \item \code{alleles}: Comma separated list of alleles for the given \code{gene}.
#'   \item \code{counts}: Comma separated list of observed sequences for each
#'         corresponding allele in the \code{alleles} list.
#'   \item \code{total}: The total count of observed sequences for the given \code{gene}.
#'   \item \code{note}: Any comments on the inference.
#'   \item \code{kh}: log10 likelihood that the \code{gene} is homozygous.
#'   \item \code{kd}: log10 likelihood that the \code{gene} is heterozygous.
#'   \item \code{kt}: log10 likelihood that the \code{gene} is trizygous
#'   \item \code{kq}: log10 likelihood that the \code{gene} is quadrozygous.
#'   \item \code{k_diff}: log10 ratio of the highest to second-highest zygosity likelihoods.
#' }
#'
#' @note
#' This method works best with data derived from blood, where a large
#' portion of sequences are expected to be unmutated. Ideally, there
#' should be hundreds of allele calls per gene in the input.
#'
#' @seealso \link{plotGenotype} for a colorful visualization and
#'          \link{genotypeFasta} to convert the genotype to nucleotide sequences.
#'          See \link{inferGenotype} to infer a subject-specific genotype using
#'          a frequency method
#'
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of
#'          vaccine-induced immune responses. PNAS. 2014 111(13):4928-33.
#' }
#'
#' @examples
#' # Infer IGHV genotype, using only unmutated sequences, including novel alleles
#' inferGenotypeBayesian(AIRRDb, germline_db=SampleGermlineIGHV, novel=SampleNovel,
#'                       find_unmutated=TRUE, v_call="v_call", seq="sequence_alignment")
#'
#' @export
inferGenotypeBayesian <- function(data, germline_db=NA, novel=NA,
                                  v_call="v_call", seq="sequence_alignment",
                                  find_unmutated=TRUE,
                                  priors=c(0.6, 0.4, 0.4, 0.35, 0.25, 0.25, 0.25, 0.25, 0.25)){
    # Visibility hack
    . <- NULL

    allele_calls <- getAllele(data[[v_call]], first=FALSE, strip_d=FALSE)
    # Find the unmutated subset, if requested
    if(find_unmutated){
        if(is.na(germline_db[1])){
            stop("germline_db needed if find_unmutated is TRUE")
        }
        if (!is.null(nrow(novel))) {
            novel <- novel %>%
                filter(!is.na(!!rlang::sym("polymorphism_call"))) %>%
                select(!!!rlang::syms(c("germline_call", "polymorphism_call", "novel_imgt")))
            if(nrow(novel) > 0){
                # Extract novel alleles if any and add them to germline_db
                novel_gl <- novel$novel_imgt
                names(novel_gl) <- novel$polymorphism_call
                germline_db <- c(germline_db, novel_gl)
                # Add the novel allele calls to allele calls of the same starting allele
                for(r in 1:nrow(novel)){
                    ind <- grep(novel$germline_call[r], allele_calls, fixed=TRUE)
                    allele_calls[ind] <- allele_calls[ind] %>%
                        sapply(paste, novel$polymorphism_call[r], sep=",")
                }
            }
        }
        # Find unmutated sequences
        allele_calls <- findUnmutatedCalls(allele_calls,
                                          as.character(data[[seq]]),
                                          germline_db)
        if(length(allele_calls) == 0){
            stop("No unmutated sequences found! Set 'find_unmutated' to 'FALSE'.")
        }
    }

    # Find which rows' calls contain which genes
    gene_regex <- allele_calls %>% strsplit(",") %>% unlist() %>%
        getGene(strip_d=FALSE) %>%  unique() %>% paste("\\*", sep="")
    gene_groups <- sapply(gene_regex, grep, allele_calls, simplify=FALSE)
    names(gene_groups) <- gsub("\\*", "", gene_regex, fixed=TRUE)
    gene_groups <- gene_groups[sortAlleles(names(gene_groups))]

    # Make a table to store the resulting genotype
    gene <- names(gene_groups)
    #   alleles = counts = note = rep("", length(gene))
    #   total = sapply(gene_groups, length)
    #   genotype = cbind(gene, alleles, counts, total, note)
    alleles <- counts <- kh <- kd <- kt <- kq <- k_diff <- note <- rep("", length(gene))
    total <- sapply(gene_groups, length)
    genotype <- cbind(gene, alleles, counts, total, note, kh, kd, kt, kq, k_diff)

    # For each gene, find which alleles to include
    for (g in gene){
        # Keep only the part of the allele calls that uses the gene being analyzed
        ac <- allele_calls[gene_groups[[g]]] %>%
            strsplit(",") %>%
            lapply(function(x) x[grep(paste(g, "\\*", sep=""), x)]) %>%
            sapply(paste, collapse=",")
        t_ac <- table(ac) # table of allele calls
        potentials <- unique(unlist(strsplit(names(t_ac),","))) # potential alleles

        regexpotentials <- paste(gsub("\\*","\\\\*", potentials),"$",sep="")
        regexpotentials <-
            paste(regexpotentials,gsub("\\$",",",regexpotentials),sep="|")
        tmat <-
            sapply(regexpotentials, function(x) grepl(x, names(t_ac),fixed=FALSE))

        if (length(potentials) == 1 | length(t_ac) == 1){
            seqs_expl <- t(as.data.frame(apply(t(as.matrix(tmat)), 2, function(x) x *
                                                  t_ac)))
            rownames(seqs_expl) <- names(t_ac)[1]
        }else{
            seqs_expl <- as.data.frame(apply(tmat, 2, function(x) x *
                                                t_ac))
        }
        #       seqs_expl = as.data.frame(apply(tmat, 2, function(x) x*t_ac))
        colnames(seqs_expl) <- potentials
        # Add low (fake) counts
        sapply(colnames(seqs_expl), function(x){if(sum(rownames(seqs_expl) %in% paste(x)) == 0){
            seqs_expl <<- rbind(seqs_expl,rep(0,ncol(seqs_expl)));
            rownames(seqs_expl)[nrow(seqs_expl)] <<- paste(x)
            seqs_expl[rownames(seqs_expl) %in% paste(x),paste(x)] <<- 0.01

        }})

        # Build ratio dependent allele count distribution of multi assigned reads
        seqs_expl_single <- seqs_expl[grep(',',rownames(seqs_expl),invert = T),]

        seqs_expl_multi <- seqs_expl[grep(',',rownames(seqs_expl),invert = F),]
        if(is.null(nrow(seqs_expl_multi))){
            seqs_expl_multi <- t(as.data.frame(seqs_expl_multi))
            rownames(seqs_expl_multi) <- grep(',',rownames(seqs_expl),invert = F,value = T)
        }

        if(!is.null(nrow(seqs_expl_single))  && nrow(seqs_expl_single) !=0 && nrow(seqs_expl_single) != nrow(seqs_expl)){
            if(nrow(seqs_expl_multi)>1){
                seqs_expl_multi <- seqs_expl_multi[order(nchar(row.names(seqs_expl_multi))),]
            }
            sapply(1:nrow(seqs_expl_multi),function(x){
                genes <- unlist(strsplit(row.names(seqs_expl_multi)[x],','));
                counts <- seqs_expl_single[rownames(seqs_expl_single) %in% genes,genes]
                counts <- colSums(counts)
                counts_to_distribute <- seqs_expl_multi[x,genes]

                new_counts <- counts+((counts_to_distribute*counts)/sum(counts))
                for(i in 1:length(new_counts)){
                    gene_tmp <- names(new_counts)[i]
                    seqs_expl_single[rownames(seqs_expl_single) %in% gene_tmp,gene_tmp] <<- new_counts[i]
                }
            })
        }

        # Cycle through the table, including alleles to explain more sequences,
        # until we explain enough sequences
        #included = counts = character(0)
        #tot_expl = 0

        seqs_expl <- if(is.null(nrow(seqs_expl_single)) || nrow(seqs_expl_single) ==0 ){seqs_expl}else{seqs_expl_single}
        seqs_expl <- round(seqs_expl)
        if(sum(rowSums(seqs_expl) == 0 ) != 0){
            seqs_expl <- seqs_expl[rowSums(seqs_expl)!= 0, ]
        }

        allele_tot <- sort(apply(seqs_expl, 2, sum),decreasing=TRUE)
        len=min(length(allele_tot),4);
        #print(priors)
        probs <-get_probabilites_with_priors(sort(c(allele_tot,rep(0,4-len)),decreasing = T)[1:4],priors = priors)
        probs[probs==-Inf] <- -1000
        names(probs) <- c('H','D','T','Q')

        k <- sort(as.numeric(probs),decreasing = T);

        probs<-c(probs,k[1]-k[2])
        names(probs)[5] <- "k_diff"

        genotype[genotype[, "gene"] == g, "alleles"] <- paste(gsub("[^d\\*]*[d\\*]",
                                                                  "", names(allele_tot)[1:len]), collapse = ",")
        genotype[genotype[, "gene"] == g, "counts"] <- paste(as.numeric(allele_tot)[1:len],
                                                            collapse = ",")
        genotype[genotype[, "gene"] == g, "kh"] <- probs[1];
        genotype[genotype[, "gene"] == g, "kd"] <- probs[2];
        genotype[genotype[, "gene"] == g, "kt"] <- probs[3];
        genotype[genotype[, "gene"] == g, "kq"] <- probs[4];
        genotype[genotype[, "gene"] == g, "k_diff"] <- probs[5];
        #     }

    }


    geno <- as.data.frame(genotype, stringsAsFactors = FALSE)

    # Check for indistinguishable calls
    if(find_unmutated == TRUE){
        seqs <- genotypeFasta(geno, germline_db)
        dist_mat <- seqs %>%
            sapply(function(x) sapply((getMutatedPositions(seqs, x)), length))
        rownames(dist_mat) <- colnames(dist_mat)
        for (i in 1:nrow(dist_mat)){ dist_mat[i,i] = NA }
        same <- which(dist_mat == 0, arr.ind=TRUE)
        if (nrow(same) > 0 ) {
            for (r in 1:nrow(same)) {
                inds <- as.vector(same[r,])
                geno[getGene(rownames(dist_mat)[inds][1]),]$note <-
                    paste(rownames(dist_mat)[inds], collapse=" and ") %>%
                    paste("Cannot distinguish", .)
            }
        }
    }
    rownames(geno) <- NULL
    return(geno)
}


# Calculate models likelihood
#
#
# @param    X      a vector of counts
# @param    alpha_dirichlet      alpha parameter for dirichlet distribution
# @param    epsilon    epsilon
# @param    priors      a vector of priors
#
# @return  log10 of the likelihoods
get_probabilites_with_priors <- function(X, alpha_dirichlet=c(0.5,0.5,0.5,0.5)*2,
                                         epsilon=0.01,
                                         priors=c(0.5,0.5,0.33,0.33,0.33,0.25,0.25,0.25,0.25)){
    ## Hypotheses
    X<-sort(X,decreasing=TRUE)

    H1<-c(1,0,0,0)
    H2<-c(priors[1],priors[2],0,0)
    H3<-c(priors[3],priors[4],priors[5],0)
    H4<-c(priors[6],priors[7],priors[8],priors[9])

    E1<-ddirichlet((H1+epsilon)/sum(H1+epsilon),alpha_dirichlet+X)
    E2<-ddirichlet((H2+epsilon)/sum(H2+epsilon),alpha_dirichlet+X)
    E3<-ddirichlet((H3+epsilon)/sum(H3+epsilon),alpha_dirichlet+X)
    E4<-ddirichlet((H4+epsilon)/sum(H4+epsilon),alpha_dirichlet+X)



    while(sort(c(E1,E2,E3,E4),decreasing=TRUE)[2] == 0 ){

        X <- X/10
        E1<-ddirichlet((H1+epsilon)/sum(H1+epsilon),alpha_dirichlet+X)
        E2<-ddirichlet((H2+epsilon)/sum(H2+epsilon),alpha_dirichlet+X)
        E3<-ddirichlet((H3+epsilon)/sum(H3+epsilon),alpha_dirichlet+X)
        E4<-ddirichlet((H4+epsilon)/sum(H4+epsilon),alpha_dirichlet+X)

    }

    return(log10(c(E1,E2,E3,E4)))
}
