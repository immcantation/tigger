# Project documentation for tigger
# 
# @author     Daniel Gadala-Maria
# @copyright  Copyright 2016 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Unported
# @version    0.2.3.999
# @date       2016.04.13


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
#'   \item  \code{\link{findNovelAlleles}}:   Detect novel alleles
#'   \item  \code{\link{plotNovel}}:          Plot evidence of novel alleles
#'   \item  \code{\link{inferGenotype}}:      Infer an Ig genotype
#'   \item  \code{\link{plotGenotype}}:       A colorful genotype visualization
#'   \item  \code{\link{genotypeFasta}}:      Convert a genotype to sequences
#'   \item  \code{\link{reassignAlleles}}:    Correct allele calls
#' }
#' 
#' @section  Mutation-related functions:
#' \itemize{
#'   \item  \code{\link{getMutatedPositions}}:      Find mutation locations
#'   \item  \code{\link{getMutCount}}:              Find distance from germline
#'   \item  \code{\link{findUnmutatedCalls}}:       Subset unmutated sequences
#'   \item  \code{\link{getPopularMutationCount}}:  Find most common sequence's
#'                                                  mutation count
#'   \item  \code{\link{insertPolymorphisms}}:      Insert SNPs into a sequence
#' }
#' 
#' @section  Input and formatting:
#' \itemize{
#'   \item  \code{\link{readIgFasta}}:        Read a fasta file of Ig sequences
#'   \item  \code{\link{updateAlleleNames}}:  Correct outdated allele names
#'   \item  \code{\link{sortAlleles}}:        Sort allele names intelligently
#'   \item  \code{\link{cleanSeqs}}:          Standardize sequence format
#' }
#' 
#' @name        tigger
#' @docType     package
#' @references  Gadala-Maria \emph{et al}. (2015) Automated analysis of
#' high-throughput B cell sequencing data reveals a high frequency of novel
#' immunoglobulin V gene segment alleles. \emph{PNAS}. 112(8):E862-70.
#' 
#' @import alakazam 
#' @import shazam
#' @import doParallel
#' @import tidyr
#' @import dplyr
#' @import foreach
#' @import ggplot2
#' @import grid
#' @import iterators
#' @import parallel
#' 
NULL
