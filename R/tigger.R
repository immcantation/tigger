# Project documentation for tigger
# 
# @author     Daniel Gadala-Maria
# @copyright  Copyright 2015 Kleinstein Lab, Yale University. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    2.0
# @date       2015.03.04


#' tigger
#' 
#' Here we provide a *T*ool for *I*mmuno*g*lobulin *G*enotype *E*lucidation via
#' *R*ep-Seq (TIgGER). TIgGER inferrs the set of Ig alleles carried by an
#' individual (including any novel alleles) and then uses this set of alleles to
#' correct the initial assignments given to sample sequences by existing tools.
#' 
#' @details Immunoglobulin Repertoire-Sequencing (Rep-Seq) data is currently the
#' subject of much study. A key step in analyzing these data involves assigning
#' the closest known V(D)J germline alleles to the (often somatically mutated)
#' sample sequences using a tool such as IMGT/HighV-QUEST. However, if the
#' sample utilizes alleles not in the germline database used for alignment, this
#' step will fail. Additionally, this alignment has an associated error rate of
#' ~5%, notably among sequences carrying a large number of somatic mutations.
#' The purpose of TIgGER is to address these issues.
#' 
#' @name        tigger
#' @docType     package
#' @references  Gadala-Maria D, Yaari G, Uduman M, Kleinstein SH (2015) Automated analysis of high-throughput B cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles. PNAS. 112(8):E862-70 
#' 
#' @import dplyr
#' 
NULL
