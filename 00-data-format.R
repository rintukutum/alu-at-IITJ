dir.create(path = 'data',showWarnings = FALSE)
dir.create(path = 'figures', showWarnings = FALSE)
# install.packages("tidyverse")
# install.packages("BiocManger")
# BiocManager::install('ggbio')
# https://www.openintro.org/

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
#--------------- STEP 01
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
colnames(jaspar_mef2c) <- c(
  'chr','start','end','info','dot','strand'
)
View(jaspar_mef2c)
#--------- chr 1
library(ggbio)
p.ideo <- ggbio::Ideogram(genome = "hg19")

library(GenomicRanges)

p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))

#' Reference for S4 class
#' http://adv-r.had.co.nz/S4.html
#' 
#  IRanges(start = 50324340, end = 50324355)
#  https://bioconductor.org/packages/release/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf

#' Home work
#' Convert the jaspar_mef2c S3 data.frame {X2, X3 column} to IRanges object?
#' 

IRanges(start = 1, end = 10)

# vector
jm_start <- jaspar_mef2c$start
jm_end <- jaspar_mef2c$end

jm_irange <- IRanges(start = jm_start, end = jm_end)

#' https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html

# extarct chr information from data.frame
jm_chrs <- jaspar_mef2c$chr

jm_strand <- jaspar_mef2c$strand

jm_hg <- 'hg19'

jm_grange <- GRanges(
  seqnames = Rle(jm_chrs),
  ranges = jm_irange,
  strand = strand(jm_strand)
)
save(jm_grange,jaspar_mef2c,
     file = 'data/jm.RData')

#--------------- STEP 02
rm(list=ls())
alu_whole <- readr::read_tsv(
  'data/alu_whole',
  col_names = FALSE
)
#' review of the ucsc.score for alu
colnames(alu_whole) <- c(
  'chr','alu.start','alu.end','alu.type',
  'ucsc.score','strand'
)
alu_grange <- GRanges(
  seqnames = Rle(alu_whole$chr),
  ranges = IRanges(start = alu_whole$alu.start,
                   end = alu_whole$alu.end),
  strand = strand(Rle(alu_whole$strand))
)

load('data/jm.RData')


# barplot(sort(table(alu_whole$chr),
#              decreasing = TRUE)[1:10],
#              horiz = TRUE,
#         las=2)

idx_chr1 <- alu_whole$chr == 'chr1'
alu_chr1 <- alu_whole[idx_chr1,]

# barplot(sort(table(alu_chr$alu.type)),
#         horiz = TRUE,
#         las=2)

chr1_alu_types <- unique(alu_chr1$alu.type)

n_chr1_alu_types <- length(chr1_alu_types)

loop <- 1:n_chr1_alu_types

# for(index in 1:5){
#   #print(index)
# }
par(mfrow=c(2,2))
for(index in 1:4){
  alu_type <- alu_types[index]
  idx_alu <- alu_chr$alu.type == alu_type
  alu_score <- alu_chr[idx_alu,]
  boxplot(alu_score$ucsc.score~alu_score$strand,
          ylab='Score',main= alu_type,
          xlab='Strand')
}

library(ggplot2)
ggplot(data = alu_chr1,
       aes(x = strand,y=log10(ucsc.score))) +
  geom_boxplot() +
  facet_wrap(facets = 'alu.type')

