# Project documentation for tigger
# 
# @author     Daniel Gadala-Maria
# @copyright  Copyright 2016 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Unported
# @version    0.2.7
# @date       2016.07.24


#' tigger
#' 
#' Here we provide a \strong{T}ool for \strong{I}mmuno\strong{g}lobulin
#' \strong{G}enotype \strong{E}lucidation via
#' \strong{R}ep-Seq (TIgGER). TIgGER inferrs the set of Ig alleles carried by an
#' individual (including any novel alleles) and then uses this set of alleles to
#' correct the initial assignments given to sample sequences by existing tools.
#' 
#' @details Immunoglobulin Repertoire-Sequencing (Rep-Seq) data is currently the
#' subject of much study. A key step in analyzing these data involves assigning
#' the closest known V(D)J germline alleles to the (often somatically mutated)
#' sample sequences using a tool such as IMGT/HighV-QUEST. However, if the
#' sample utilizes alleles not in the germline database used for alignment, this
#' step will fail. Additionally, this alignment has an associated error rate of
#' ~5 percent, notably among sequences carrying a large number of somatic
#' mutations. The purpose of TIgGER is to address these issues.
#' 
#' @section  Core tigger functions:
#' \itemize{
#'   \item  \link{findNovelAlleles}:   Detect novel alleles
#'   \item  \link{plotNovel}:          Plot evidence of novel alleles
#'   \item  \link{inferGenotype}:      Infer an Ig genotype
#'   \item  \link{plotGenotype}:       A colorful genotype visualization
#'   \item  \link{genotypeFasta}:      Convert a genotype to sequences
#'   \item  \link{reassignAlleles}:    Correct allele calls
#' }
#' 
#' @section  Mutation-related functions:
#' \itemize{
#'   \item  \link{getMutatedPositions}:      Find mutation locations
#'   \item  \link{getMutCount}:              Find distance from germline
#'   \item  \link{findUnmutatedCalls}:       Subset unmutated sequences
#'   \item  \link{getPopularMutationCount}:  Find most common sequence's
#'                                           mutation count
#'   \item  \link{insertPolymorphisms}:      Insert SNPs into a sequence
#' }
#' 
#' @section  Input and formatting:
#' \itemize{
#'   \item  \link{readIgFasta}:        Read a fasta file of Ig sequences
#'   \item  \link{updateAlleleNames}:  Correct outdated allele names
#'   \item  \link{sortAlleles}:        Sort allele names intelligently
#'   \item  \link{cleanSeqs}:          Standardize sequence format
#' }
#' 
#' @name        tigger
#' @docType     package
#' @references  Gadala-Maria \emph{et al}. (2015) Automated analysis of
#' high-throughput B cell sequencing data reveals a high frequency of novel
#' immunoglobulin V gene segment alleles. \emph{PNAS}. 112(8):E862-70.
#' 
#' @import      alakazam 
#' @import      shazam
#' @import      doParallel
#' @importFrom  graphics    plot
#' @importFrom  stats       na.omit setNames ecdf sd cor cov median mad
#'                          confint lm
#' @importFrom  tidyr       gather gather_ spread spread_
#' @importFrom  dplyr       do n desc %>%
#'                          glimpse distinct distinct_
#'                          as_data_frame data_frame data_frame_
#'                          bind_cols bind_rows combine
#'                          filter filter_ select select_ arrange arrange_
#'                          group_by group_by_ ungroup
#'                          mutate mutate_ transmute transmute_
#'                          rename rename_ summarise summarise_
#'                          slice slice_
#' @import      foreach
#' @import      ggplot2
#' @importFrom  grid        grid.layout grid.newpage pushViewport viewport
#' @import      iterators
#' @import      parallel
#' 
NULL
