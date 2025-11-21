library(seqinr)
library(ggplot2)
library(tune)
require(RColorBrewer)
library('scales')
library("ggExtra")
#library(wordcloud)
#library(tm)
setwd("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/")
#pdf("new_data/figures/39.emmadata.GLFR.pdf", height = 20, width = 20)
png("new_data/figures/39.emmadata.GLFR.png", height = 1500, width = 1500)

emma_data = read.csv("new_data/data/FullData.txt", sep='\t', header= TRUE)
emma_data = emma_data[emma_data$Strand == 'positive ',]
emma_data$col = 'black'
emma_data[emma_data$Order == 'Nidovirales',]$col = 'deeppink4'
emma_data[emma_data$Family == 'Coronaviridae',]$col = 'red'

#emma_data$all = rowSums(emma_data[,6:25])
#emma_data = subset(emma_data, all > quantile(emma_data$all, 0.1) & all < quantile(emma_data$all, 0.9))

myplotfr <- function(data) {
  loosers = c('P', 'A', 'G', 'R', 'H', 'Q', 'D', 'E', 'T', 'S')
  gainers = c('F', 'I', 'M', 'Y')
  data$losersfr = rowSums(data[,loosers]) / rowSums(data[,6:25])
  data$gainerfr = rowSums(data[,gainers]) / rowSums(data[,6:25])
  data$all = rowSums(data[,6:25])
  
  #data = subset(data, all > quantile(data$all, 0.05) & all < quantile(data$all, 0.95))

  COL = c()
  for (rownum in 1:length(row.names(data))){
    if (data$Family[rownum] == 'Coronaviridae'){
      COL = append(COL, 'red', after = length(COL))
    }else if (data$Order[rownum] == 'Nidovirales'){
      COL = append(COL, 'green', after = length(COL))
    }else{
      COL = append(COL, 'black', after = length(COL))
    }
  }

  # plot = ggplot(data = data[data$col == 'black',], aes(x=losersfr, y=gainerfr)) + geom_point(data = data[data$col == 'black',], alpha = 1/25, size=2, aes(x=losersfr, y=gainerfr, colour = 'other'), colour='black') +
  # geom_point(data = data[data$col == 'deeppink4',], alpha = 1, size=3, aes(x=losersfr, y=gainerfr, colour = 'Nidovirales'), colour='deeppink4') +
  # geom_point(data = data[data$col == 'red',], size=3, alpha=1, aes(x=losersfr, y=gainerfr, colour = 'Coronaviridae'), colour='red')  + 
  #   labs(x = "Losers Frequency", y="Gainers Frequency") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                                                                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=20),
  #                                                                                                                                          axis.title=element_text(size=40,face="bold"))
  data$sep = ''
  data[data$col == 'black',]$sep = 'other +RNA viruses'
  data[data$col == 'red',]$sep = 'Coronaviridae'
  data[data$col == 'deeppink4',]$sep = 'other Nidovirales'
  plot = ggplot(data = data, aes(x=losersfr, y=gainerfr, colour=factor(sep))) + scale_color_manual(values=c("#FF0000", "#000000", "#008000")) + geom_point(size = ifelse(data$col == 'black', 2, 4), aes(x=losersfr, y=gainerfr, alpha = factor(sep))) +
    scale_alpha_manual(values = c(1, 0.02, 1)) + 
    labs(x = "Losers Frequency", y="Gainers Frequency") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=20),
                                                                axis.title=element_text(size=40,face="bold")) +
    theme(legend.title=element_blank())+
    theme(legend.text = element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    theme(legend.position = c(0.85, 0.9))+
    scale_x_continuous(limits = c(0.35, 0.65)) +
    scale_y_continuous(limits = c(0.05, 0.35))
  margolt = ggMarginal(plot, groupColour = TRUE, groupFill = TRUE)
  
  # hist(data[data$Order != 'Nidovirales',]$all, breaks = 30, col=rgb(1,0,0,0.5) )
  # 
  # # Second with add=T to plot on top
  # hist(data[data$Order == 'Nidovirales',]$all, breaks = 30, col=rgb(0,0,1,0.5))
  margolt
  
}

myplotfr(emma_data)


dev.off()

##############
