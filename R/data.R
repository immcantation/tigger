#' Human IGHV germlines
#'
#' A \code{character} vector of all 344 human IGHV germline gene segment alleles
#' in IMGT Gene-db release 201408-4.
#'
#' @name germline_ighv
#' @docType data
#' @format Values correspond to IMGT-gaped nuceltoide sequences (with
#' nucleotides capitalized and gaps represented by ".") while names correspond
#' to stripped-down IMGT allele names (e.g. "IGHV1-18*01").
#' 
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
NULL


#' Example human Rep-Seq data
#'
#' Example VDJ-rearranged immunoglobulin Rep-Seq sequences derived from a single
#' individual (PGP1), sequenced on the Roche 454 platform, and thought by
#' IMGT/V-QUEST to utilize IGHV1 family alleles.
#'
#' @name sample_db
#' @docType data
#' @format A \code{data.frame} where rows correspond to unique VDJ sequences and
#' columns include:
#' \itemize{
#'   \item IMGT-gapped nucleotide sequence (\code{"SEQUENCE_IMGT"})
#'   \item IMGT/V-QUEST allele calls (\code{"V_CALL"}, \code{"D_CALL"}, and
#'     \code{"J_CALL"})
#'   \item Junction length (\code{"JUNCTION_LENGTH"})
#' }
#' 
#' @references Gadala-Maria \emph{et al}. (2015) Automated analysis of
#' high-throughput B cell sequencing data reveals a high frequency of novel
#' immunoglobulin V gene segment alleles. \emph{PNAS}. 112(8):E862-70.
#' @keywords data
NULL

#' Example of Analyzed Rep-Seq data
#'
#' Example VDJ-rearranged immunoglobulin Rep-Seq sequences derived from a single
#' individual (PGP1), sequenced on the Roche 454 platform, and thought by
#' IMGT/V-QUEST to utilize IGHV1 family alleles, as processed by
#' \link{findNovelAlleles}.
#'
#' @name novel_df
#' @docType data
#' @format A \code{data.frame} where rows correspond to alleles checked for
#' polymorphisms and columns give results as well as paramaters used to run
#' the test.
#' 
#' @references Gadala-Maria \emph{et al}. (2015) Automated analysis of
#' high-throughput B cell sequencing data reveals a high frequency of novel
#' immunoglobulin V gene segment alleles. \emph{PNAS}. 112(8):E862-70.
#' @keywords data
NULL

#' Example of an Inferred Genotype
#'
#' Example VDJ-rearranged immunoglobulin Rep-Seq sequences derived from a single
#' individual (PGP1), sequenced on the Roche 454 platform, and thought by
#' IMGT/V-QUEST to utilize IGHV1 family alleles, as processed by
#' \link{findNovelAlleles} and \link{inferGenotype}
#'
#' @name genotype
#' @docType data
#' @format A \code{data.frame} where rows correspond to genes carried by an
#' individual and columns lists the alleles of those genes and their counts.
#' 
#' @references Gadala-Maria \emph{et al}. (2015) Automated analysis of
#' high-throughput B cell sequencing data reveals a high frequency of novel
#' immunoglobulin V gene segment alleles. \emph{PNAS}. 112(8):E862-70.
#' @keywords data
NULL