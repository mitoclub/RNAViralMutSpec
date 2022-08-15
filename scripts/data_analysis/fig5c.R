library(seqinr)
library(wordcloud)
library(tm)
setwd("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/")
pdf("new_data/figures/37.Corr_FrffU_Gainers.pdf", height = 21.32, width = 20.66)
virusff = read.csv('new_data/data_obtained/23.VirusesFF.csv')
virusff = virusff[virusff$VirType == 'ssRnaPlus',]
viruslg = read.csv('new_data/data_obtained/23.VirusesLG.csv')

viruses = merge(virusff, viruslg, by = "VirSpeciesCommonName")  

plot(viruses$FourFoldFrT, viruses$GainFr, pch = 19, col = "lightblue", xlab='Four Fold Frequency U', ylab='Gainers Frequency', xlim = c(0.15,0.8), ylim = c(0.1,0.25), cex.lab = 1.7, cex.axis = 1.7)

# Regression line
abline(lm(viruses$GainFr ~ viruses$FourFoldFrT), col = "red", lwd = 3)

# Pearson correlation
text(paste("Correlation:", round(cor(viruses$FourFoldFrT, viruses$GainFr, method='spearman'), 2)), x = 0.3, y = 0.2, cex=2.5)
par(new=TRUE)
textplot(viruses$FourFoldFrT, viruses$GainFr,viruses$VirSpeciesCommonName, xlab='', ylab='', cex=1.7, col = viruses$Col, xlim = c(0.15,0.8), ylim = c(0.1,0.25), xaxt = "n", yaxt = "n")

dev.off()