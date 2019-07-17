# Imports
library(tigger)
library(alakazam)
library(profvis)

#### Load example data ####

data("SampleDb")
data("GermlineIGHV")

#### Find novel alleles ####

profvis({
    nv <- findNovelAlleles(SampleDb, GermlineIGHV,nproc=6)
})
