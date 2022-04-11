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

idx.rm <- grep('\\_ ',alu_whole$chr)

alu_whole_chrs <- alu_whole[-idx.rm,]
alu_grange <- GRanges(
  seqnames = Rle(alu_whole_chrs$chr),
  ranges = IRanges(start = alu_whole_chrs$alu.start,
                   end = alu_whole_chrs$alu.end),
  strand = strand(Rle(alu_whole_chrs$strand))
)

load('data/jm.RData')

table(alu_whole_chrs[,c('chr','strand')])



# barplot(sort(table(alu_whole$chr),
#              decreasing = TRUE)[1:10],
#              horiz = TRUE,
#         las=2)
alu_whole <- alu_whole_chrs
idx_chr1 <- alu_whole$chr == 'chr1'
alu_chr1 <- alu_whole[idx_chr1,]

barplot(sort(table(alu_chr1$alu.type)),
        horiz = TRUE,
        las=2)

chr1_alu_types <- unique(alu_chr1$alu.type)

n_chr1_alu_types <- length(chr1_alu_types)

loop <- 1:n_chr1_alu_types

# for(index in 1:5){
#   #print(index)
# }
par(mfrow=c(2,2))
for(index in 1:4){
  alu_type <- chr1_alu_types[index]
  idx_alu <- alu_chr1$alu.type == alu_type
  alu_score <- alu_chr1[idx_alu,]
  #-------
  boxplot(
    alu_score$ucsc.score ~ alu_score$strand,
    ylab ='Score',
    main = alu_type,
    xlab ='Strand',
    col  = c('red','skyblue')
  )
}

library(ggplot2)
ggplot(data = alu_chr1,
       aes(x = strand,
           y = log10(ucsc.score))) +
  geom_boxplot() +
  facet_wrap(facets = 'alu.type')
#-----------
# distribution of alu
alu_type_chr <- data.frame(table(alu_whole[,c('chr','alu.type','strand')]))
ggplot(alu_type_chr,aes(x=chr,y=alu.type)) +
  geom_tile(aes(fill=log10(Freq+1))) +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(
    angle=45,vjust = 1,hjust = 1))

library(Homo.sapiens)
library(ggbio)
library(biovizBase)
# p <- ggbio() + circle(
#   jm_grange,geom='rect', aes(color = strand)) +
#   circle(hg19sub, geom = "ideo", fill = "gray70") +
#   circle(hg19sub, geom = "scale", size = 2) + 
#   circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
# 
# #-------
# pos_strand <- jm_grange[jm_grange@strand@values == '+',]
# neg_strand <- jm_grange[jm_grange@strand@values == '-',]
# p <- ggbio() + 
#   circle(pos_strand,geom='rect', color = 'steelblue') +
#   circle(neg_strand,geom='rect', color = 'red') +
#   circle(hg19sub, geom = "ideo", fill = "gray70") +
#   circle(hg19sub, geom = "scale", size = 2) + 
#   circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)

#--------
jaspar_mef2c <- readr::read_tsv(
  'data/Jaspar_MEF2C',
  col_names = FALSE
)
colnames(jaspar_mef2c) <- c(
  'chr','start','end','info','dot','strand'
)
jaspar_mef2c <- jaspar_mef2c[-grep('X',jaspar_mef2c$chr),]
# vector
jm_start <- jaspar_mef2c$start
jm_end <- jaspar_mef2c$end
jm_irange <- IRanges(start = jm_start, end = jm_end)
jm_chrs <- jaspar_mef2c$chr
jm_chrs <- gsub('chr','',jm_chrs)
jm_strand <- jaspar_mef2c$strand

jm_g <- GRanges(
  seqnames = Rle(jm_chrs),
  ranges = jm_irange,
  strand = strand(jm_strand)
)
genome(jm_g) <- 'hg19'
pos_s <- jm_g[jm_g@strand@values == '+',]
pos_s <- jm_g[jm_g@strand@values == '-',]
p <- ggbio() + 
  circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) + 
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3) +
  circle(pos_s, geom='rect', color = 'steelblue') +
  circle(pos_s, geom='rect', color = 'red')
