dir.create(path = 'data',showWarnings = FALSE)
dir.create(path = 'figures', showWarnings = FALSE)
# install.packages("tidyverse")
# install.packages("BiocManger")
# BiocManager::install('ggbio')

# = 
num_vec <-  c(1, 4, 9, 818)
char_vec <-  c('a', 'rintu','kutum')

c(1,4,'a')
# = 


# list function
# 
ls()

#----------------
# S3 class
# S4 class

rm(list=ls())
alu_whole <- readr::read_tsv(
  'data/alu_whole',
  col_names = FALSE
)

View(alu_whole)

jaspar_mef2c <- readr::read_tsv(
  'data/Jaspar_MEF2C',
  col_names = FALSE
)

View(jaspar_mef2c)
#--------- chr 1
library(ggbio)
p.ideo <- ggbio::Ideogram(genome = "hg19")

library(GenomicRanges)

p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))


# IRanges(start = 50324340, end = 50324355)
# https://bioconductor.org/packages/release/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf

#' Home work
#' Convert the jaspar_mef2c S3 data.frame {X2, X3 column} to IRanges object?
#' 



