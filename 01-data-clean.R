#' ref https://ggplot2-book.org/index.html
#' 
#' ref https://jokergoo.github.io/circlize_book/book/
#' 
rm(list=ls())
datasets <- list.files(path = 'data/', pattern = 'DATASET', full.names = TRUE)

alu_mef2c <- readr::read_csv(datasets[5])
alu_mef2c_exons <- readr::read_csv(datasets[1])
alu_mef2c_ig <- readr::read_csv(datasets[2])
alu_mef2c_intron <- readr::read_csv(datasets[3])
alu_mef2c_promo <- readr::read_csv(datasets[4])
total <- nrow(alu_mef2c_exons) + nrow(alu_mef2c_ig) +
  nrow(alu_mef2c_intron) + nrow(alu_mef2c_promo)

alu_freq <- c(
  'exons' = nrow(alu_mef2c_exons),
  'intergenic' = nrow(alu_mef2c_ig),
  'intron' = nrow(alu_mef2c_intron),
  'promoter' = nrow(alu_mef2c_promo))

alu.prop <- round(c('exons' = nrow(alu_mef2c_exons)/total,
  'intergenic' = nrow(alu_mef2c_ig)/total,
  'intron' = nrow(alu_mef2c_intron)/total,
  'promoter' = nrow(alu_mef2c_promo)/total
  )*100,2)

pie(alu.prop)

#-------
non_alu_mef2c <- readr::read_csv(datasets[10])
non_alu_mef2c_exons <- readr::read_csv(datasets[6])
non_alu_mef2c_ig <- readr::read_csv(datasets[7])
non_alu_mef2c_intron <- readr::read_csv(datasets[8])
non_alu_mef2c_promo <- readr::read_csv(datasets[9])
non_total <- nrow(non_alu_mef2c_exons) + nrow(non_alu_mef2c_ig) +
  nrow(non_alu_mef2c_intron) + nrow(non_alu_mef2c_promo)


non_alu_freq <- c('exons' = nrow(non_alu_mef2c_exons),
              'intergenic' = nrow(non_alu_mef2c_ig),
              'intron' = nrow(non_alu_mef2c_intron),
              'promoter' = nrow(non_alu_mef2c_promo))
non_alu.prop <- round(c('exons' = nrow(non_alu_mef2c_exons)/non_total,
                    'intergenic' = nrow(non_alu_mef2c_ig)/non_total,
                    'intron' = nrow(non_alu_mef2c_intron)/non_total,
                    'promoter' = nrow(non_alu_mef2c_promo)/non_total
)*100,2)

pdf('figures/Alu-MEF2C.pdf', width = 5, height = 5)
pie(alu.prop, main = 'Alu with MEF2C')
dev.off()

png('figures/Alu-MEF2C.png', 
    width = 800, height = 700, 
    res = 150)
pie(alu.prop, main = 'Alu with MEF2C')
dev.off()

pie(alu.prop, main = 'Alu with MEF2C')
dev.off()

pdf('figures/Alu-MEF2C.pdf',width = 5, height = 5)
pie(non_alu.prop, main = 'Non-Alu with MEF2C')
dev.off()

names(non_alu.prop) <- paste0(names(non_alu.prop),'(',non_alu.prop,'%)')
png('figures/non-Alu-MEF2C.png', 
    width = 800, height = 700, 
    res = 150)
pie(non_alu.prop, main = 'Non-Alu with MEF2C')
dev.off()

prop.test(alu_freq, non_alu_freq)$p.value

xx <- data.frame(alu_freq,non_alu_freq)
colnames(xx) <- c('Alu MEF2C','Non Alu MEF2C')
xx$genome.feature <- rownames(xx)

xx_m <- reshape2::melt(xx)
colnames(xx_m)[1:2] <- c('to', 'from')

chordDiagram(xx_m)