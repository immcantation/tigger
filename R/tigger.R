# Project documentation for tigger

#' tigger
#' 
#' Here we provide a \strong{T}ool for \strong{I}mmuno\strong{g}lobulin
#' \strong{G}enotype \strong{E}lucidation via \strong{R}ep-Seq (TIgGER). 
#' TIgGER inferrs the set of Ig alleles carried by an
#' individual (including any novel alleles) and then uses this set of alleles to
#' correct the initial assignments given to sample sequences by existing tools.
#' 
#' @details
#' Immunoglobulin repertoire sequencing (AIRR-Seq, Rep-Seq) data is currently the
#' subject of much study. A key step in analyzing these data involves assigning
#' the closest known V(D)J germline alleles to the (often somatically mutated)
#' sample sequences using a tool such as IMGT/HighV-QUEST. However, if the
#' sample utilizes alleles not in the germline database used for alignment, this
#' step will fail. Additionally, this alignment has an associated error rate of
#' ~5%, notably among sequences carrying a large number of somatic
#' mutations. The purpose of TIgGER is to address these issues.
#' 
#' @section  Allele detection and genotyping:
#' \itemize{
#'   \item  \link{findNovelAlleles}:       Detect novel alleles.
#'   \item  \link{plotNovel}:              Plot evidence of novel alleles.
#'   \item  \link{inferGenotype}:          Infer an Ig genotype using a frequency approach.
#'   \item  \link{inferGenotypeBayesian}:  Infer an Ig genotype using a Bayesian approach.
#'   \item  \link{plotGenotype}:           A colorful genotype visualization.
#'   \item  \link{genotypeFasta}:          Convert a genotype to sequences.
#'   \item  \link{reassignAlleles}:        Correct allele calls.
#'   \item  \link{generateEvidence}:       Generate evidence for the genotype and 
#'                                         allele detection inferrence.
#' }
#' 
#' @section  Mutation handling:
#' \itemize{
#'   \item  \link{getMutatedPositions}:      Find mutation locations.
#'   \item  \link{getMutCount}:              Find distance from germline.
#'   \item  \link{findUnmutatedCalls}:       Subset unmutated sequences.
#'   \item  \link{getPopularMutationCount}:  Find most common sequence's
#'                                           mutation count.
#'   \item  \link{insertPolymorphisms}:      Insert SNPs into a sequence.
#' }
#' 
#' @section  Input, output and formatting:
#' \itemize{
#'   \item  \link{readIgFasta}:        Read a fasta file of Ig sequences.
#'   \item  \link{updateAlleleNames}:  Correct outdated allele names.
#'   \item  \link{sortAlleles}:        Sort allele names intelligently.
#'   \item  \link{cleanSeqs}:          Standardize sequence format.
#' }
#' 
#' @name        tigger
#' @docType     package
#' @references
#' \enumerate{
#'   \item Gadala-Maria, et al. (2015) Automated analysis of high-throughput B cell 
#'         sequencing data reveals a high frequency of novel immunoglobulin V gene 
#'         segment alleles. PNAS. 112(8):E862-70.
#' }
#' 
#' @import      ggplot2
#' @importFrom  alakazam    getAllele getGene getFamily translateDNA DNA_COLORS checkColumns
#' @importFrom  doParallel  registerDoParallel
#' @importFrom  dplyr       do n desc %>%
#'                          glimpse distinct group_indices
#'                          as_data_frame data_frame
#'                          bind_cols bind_rows combine inner_join
#'                          filter select arrange
#'                          group_by ungroup
#'                          mutate pull rename slice
#'                          summarise transmute
#' @importFrom  foreach     foreach %dopar% registerDoSEQ
#' @importFrom  graphics    plot
#' @importFrom  gridExtra   arrangeGrob
#' @importFrom  gtools      ddirichlet
#' @importFrom  iterators   icount
#' @importFrom  lazyeval    interp
#' @importFrom  parallel    clusterEvalQ clusterExport makeCluster stopCluster
#' @importFrom  rlang       .data := sym syms
#' @importFrom  stats       na.omit setNames ecdf sd cor cov median mad
#'                          confint lm
#' @importFrom  stringi     stri_length stri_detect_fixed stri_replace_all_regex
#'                          stri_sub stri_sub<- stri_trans_toupper
#' @importFrom  tidyr       gather spread unnest
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- paste("As of v1.0.0 the AIRR Rearrangement schema is now the default file format.",
                 "A description of the standard is available at https://docs.airr-community.org.",
                 "The legacy Change-O format is supported through arguments to each function",
                 "that allow the input column names to be explicitly defined.",
                 sep="\n")
    packageStartupMessage(msg)
}
