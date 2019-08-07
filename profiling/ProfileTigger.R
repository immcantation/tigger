# Imports
library(tigger)
library(alakazam)
library(profvis)

#### Load example data ####

data("airrDb")
data("GermlineIGHV")

#### Find novel alleles ####

profvis({
    nv <- findNovelAlleles(airrDb, GermlineIGHV,nproc=6)
})
