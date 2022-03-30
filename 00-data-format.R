dir.create(path = 'data',showWarnings = FALSE)
dir.create(path = 'figures', showWarnings = FALSE)
# install.packages("tidyverse")
# install.packages("BiocManger")
# BiocManager::install('ggbio')

# = 
num_vec <-  c(1,4,9,818)
char_vec <-  c('a', 'rintu','kutum')

c(1,4,'a')
# = 


# list function
# 
ls()

#----------------
rm(list=ls())
alu_whole <- readr::read_tsv(
  'data/alu_whole',
  col_names = FALSE
)

jaspar_mef2c <- readr::read_tsv(
  'data/Jaspar_MEF2C',
  col_names = FALSE
)


