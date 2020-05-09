# Imports
library(tigger)
library(alakazam)
library(profvis)

#### Load example data ####

data("AIRRDb")
data("GermlineIGHV")

#### Find novel alleles ####

profvis({
    nv <- findNovelAlleles(AIRRDb, GermlineIGHV,nproc=6)
})
