# TIgGER #

High-throughput sequencing of B cell immunoglobulin receptors is providing unprecedented insight into adaptive immunity. A key step in analyzing these data involves assignment of the germline V, D and J gene segment alleles that comprise each immunoglobulin sequence by matching them against a database of known V(D)J alleles. However, this process will fail for sequences that utilize previously undetected alleles, whose frequency in the population is unclear.

**TIgGER is a computational method that significantly improves V(D)J allele assignments by first determining the complete set of gene segments carried by an individual (including novel alleles) from V(D)J-rearrange sequences. TIgGER can then infer a subject's genotype from these sequences, and use this genotype to correct the initial V(D)J allele assignments.**

The application of TIgGER identifies a surprisingly high frequency of novel alleles, highlighting the critical need for this approach. To cite TIgGER, please use:

[Gadala-Maria D, Yaari G, Uduman M, Kleinstein SH (2015) "Automated analysis of high-throughput B cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles." *PNAS* 112(8), E862-E870](http://www.pnas.org/content/112/8/E862.abstract)


###Core Abilities###
* Detecting novel alleles
* Inferring a subject's genotype
* Correcting preliminary allele calls

###Required Input###
* A table of V(D)J-rearranged sequences from a single individual, with columns containing the following:
    * V(D)J sequences (in IMGT-gapped format)
    * Names of preliminary V allele calls
    * Name of preliminary J allele calls
    * Length of the junction region 
* Germline Ig sequences in IMGT-gapped fasta format

The former can be created through the use of [IMGT/HighV-QUEST](http://www.imgt.org/) and [Change-O CLT](http://clip.med.yale.edu/changeo/download.php).

### Requirements ###

Software             | Link
---------------------|-------------------------------------------
R Studio (IDE)       | https://www.rstudio.com/
alakazam (R package) | https://bitbucket.org/kleinstein/alakazam/
shm (R package)      | https://bitbucket.org/kleinstein/shm
dplyr (R package)    | https://cran.rstudio.com/web/packages/dplyr/
ggplot2 (R package)  | https://cran.rstudio.com/web/packages/ggplot2/
grid (R package)     | https://cran.rstudio.com/src/contrib/Archive/grid/


### Build Instructions ###

Install build dependencies:
```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

Building with Rstudio:

-  Build -> Configure Build Tools
-  Check use devtools option
-  Check use roxygen option
-  Select configure roxygen options and check everything.
-  Build -> Build and Reload

Building from the R console:

```R
library(roxygen2)
library(devtools)
document()
install_deps()
build(vignettes=FALSE)
install()
```

####Step-by-step Usage Example####
Please see the [TIgGER vignette](http://clip.med.yale.edu/tigger/Tigger-Vignette.pdf).

### Contact ###

For help, questions, or suggestions, please contact [daniel.gadala-maria@yale.edu](daniel.gadala-maria@yale.edu)