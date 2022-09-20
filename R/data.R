#' Human IGHV germlines
#'
#' A \code{character} vector of all human IGHV germline gene segment alleles
#' in IMGT/GENE-DB (2019-06-01, 372 alleles). 
#' See IMGT data updates: https://www.imgt.org/IMGTgenedbdoc/dataupdates.html.
#'
#' @name GermlineIGHV
#' @docType data
#' @format Values correspond to IMGT-gaped nuceltoide sequences (with
#' nucleotides capitalized and gaps represented by ".") while names correspond
#' to stripped-down IMGT allele names (e.g. "IGHV1-18*01").
#' 
#' @references 
#' \enumerate{
#'   \item Xochelli, et al. (2014) Immunoglobulin heavy variable (IGHV) genes and 
#'         alleles: new entities, new names and implications for research and 
#'         prognostication in chronic lymphocytic leukaemia. Immunogenetics. 67(1):61-6.
#' }
#' 
#' @keywords data
NULL

#' Example Human IGHV germlines
#'
#' A \code{character} vector of all 344 human IGHV germline gene segment alleles
#' in IMGT/GENE-DB release 2014-08-4.
#'
#' @name SampleGermlineIGHV
#' @docType data
#' @format Values correspond to IMGT-gaped nuceltoide sequences (with
#' nucleotides capitalized and gaps represented by ".") while names correspond
#' to stripped-down IMGT allele names (e.g. "IGHV1-18*01").
#' 
#' @references 
#' \enumerate{
#'   \item Xochelli, et al. (2014) Immunoglobulin heavy variable (IGHV) genes and 
#'         alleles: new entities, new names and implications for research and 
#'         prognostication in chronic lymphocytic leukaemia. Immunogenetics. 67(1):61-6.
#' }
#' 
#' @keywords data
NULL

#' Example human immune repertoire data
#'
#' A \code{data.frame} of example V(D)J immunoglobulin sequences derived from a 
#' single individual (PGP1), sequenced on the Roche 454 platform, and assigned by
#' IMGT/HighV-QUEST to IGHV1 family alleles.
#'
#' @name SampleDb
#' @docType data
#' @format A \code{data.frame} where rows correspond to unique V(D)J sequences and
#' columns include:
#' \itemize{
#'   \item \code{"SEQUENCE_IMGT"}: IMGT-gapped V(D)J nucleotide sequence.
#'   \item \code{"V_CALL"}: IMGT/HighV-QUEST V segment allele calls.
#'   \item \code{"D_CALL"}: IMGT/HighV-QUEST D segment allele calls.
#'   \item \code{"J_CALL"}: IMGT/HighV-QUEST J segment allele calls.
#'   \item \code{"JUNCTION_LENGTH"}: Junction region length.
#' }
#' @seealso See \link{AIRRDb} for an AIRR formatted version of \code{SampleDb}.
#' @references
#' \enumerate{
#'   \item Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
#'         sequencing data reveals a high frequency of novel immunoglobulin V gene 
#'         segment alleles. PNAS. 112(8):E862-70.
#' }
#' 
#' @keywords data
NULL

#' Example human immune repertoire data
#'
#' A \code{data.frame} of example V(D)J immunoglobulin sequences derived from a 
#' single individual (PGP1), sequenced on the Roche 454 platform, and assigned by
#' IMGT/HighV-QUEST to IGHV1 family alleles.
#'
#' @name AIRRDb
#' @docType data
#' @format A \code{data.frame} where rows correspond to unique V(D)J sequences and
#' columns include:
#' \itemize{
#'   \item \code{"sequence_alignment"}: IMGT-gapped V(D)J nucleotide sequence.
#'   \item \code{"v_call"}: IMGT/HighV-QUEST V segment allele calls.
#'   \item \code{"d_call"}: IMGT/HighV-QUEST D segment allele calls.
#'   \item \code{"j_call"}: IMGT/HighV-QUEST J segment allele calls.
#'   \item \code{"junction_length"}: Junction region length.
#' }
#' @seealso See \link{SampleDb} for Change-O formatted version of \code{AIRRDb}.
#' @references
#' \enumerate{
#'   \item Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
#'         sequencing data reveals a high frequency of novel immunoglobulin V gene 
#'         segment alleles. PNAS. 112(8):E862-70.
#' }
#' 
#' @keywords data
NULL

#' Example novel allele detection results
#'
#' A \code{data.frame} of novel allele detection results from \link{findNovelAlleles}. 
#' Source data was a collection of V(D)J immunoglobulin sequences derived from a single
#' individual (PGP1), sequenced on the Roche 454 platform, and assigned by
#' IMGT/HighV-QUEST to IGHV1 family alleles.
#'
#' @name SampleNovel
#' @docType data
#' @format A \code{data.frame} where rows correspond to alleles checked for
#' polymorphisms and columns give results as well as paramaters used to run
#' the test.
#' 
#' @seealso See \link{findNovelAlleles} for detailed column descriptions.
#' 
#' @references
#' \enumerate{
#'   \item Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
#'         sequencing data reveals a high frequency of novel immunoglobulin V gene 
#'         segment alleles. PNAS. 112(8):E862-70.
#' }
#' 
#' @keywords data
NULL

#' Example genotype inferrence results
#'
#' A \code{data.frame} of genotype inference results from \link{inferGenotype}
#' after novel allele detection via \link{findNovelAlleles}.
#' Source data was a collection of V(D)J immunoglobulin sequences derived from a single
#' individual (PGP1), sequenced on the Roche 454 platform, and assigned by
#' IMGT/HighV-QUEST to IGHV1 family alleles.
#' 
#' @name SampleGenotype
#' @docType data
#' @format A \code{data.frame} where rows correspond to genes carried by an
#' individual and columns lists the alleles of those genes and their counts.
#' 
#' @seealso See \link{inferGenotype} for detailed column descriptions.
#' 
#' @references
#' \enumerate{
#'   \item Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
#'         sequencing data reveals a high frequency of novel immunoglobulin V gene 
#'         segment alleles. PNAS. 112(8):E862-70.
#' }
#' 
#' @keywords data
NULL