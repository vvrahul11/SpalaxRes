plotit <- function (vector, title, genes){
  barplot(vector
          , horiz = TRUE
          , col= rev(rep(topo.colors(2), each = 10))
          , main= title
          , xlab = "log2Foldchange"
          #, ylab = "Genes"
          , names.arg = genes
          , cex.names = 1.5
          , cex.axis = 2
          , font = 2
          #, font.axis = "Bold"
          , las = 1
          , font.lab = 2
          ,  cex.lab = 1.8)
legend("topright"
         , inset=.09
         , title="Direction of gene regulation"
         , c("Down", "Up")
         , fill=topo.colors(2)
         , cex=0.9)
}

topgenes <- read.csv("C:/projects/Final Rat Trans Project Result2023/blindmolerat_rnaseqanalysis-2024/top10genesWork/Top20genesUpandDownAll.csv")
s.wound <- topgenes[, c(1, 2, 3)]
row.names(s.wound) = s.wound$S.Wound_Genes

m.wound <- topgenes[, c(4, 5, 6)]
row.names(m.wound) = m.wound$M.Wound_Genes

s.scab <- topgenes[, c(7, 8, 9)]
row.names(s.scab) = s.scab$S.Scab_Genes

m.scab <- topgenes[, c(10, 11, 12)]
row.names(m.scab) = m.scab$M.Scab_Genes

plotit(s.wound[,2],
       "Splalax wound up- and down-regulated genes"
       , levels(droplevels(s.wound$S.Wound_Genes)))

plotit(m.wound[,2],
       "Mouse wound up- and down-regulated genes"
       , levels(droplevels(m.wound$M.Wound_Genes)))

plotit(s.scab[,2], "Splalax scab up- and down-regulated genes"
          , levels(droplevels(s.wound$S.Wound_Genes)))


plotit(m.scab[,2], "Mouse scab up- and down-regulated genes"
          , levels(droplevels(m.scab$M.Scab_Genes)))

