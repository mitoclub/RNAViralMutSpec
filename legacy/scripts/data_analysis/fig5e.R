library(seqinr)
library(ggplot2)
library(tune)
require(RColorBrewer)
library('scales')
#library(wordcloud)
#library(tm)
setwd("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/")
pdf("new_data/figures/37.emmadata.GL.pdf", height = 20, width = 20)

emma_data = read.csv("new_data/data/FullData.txt", sep='\t', header= TRUE)
emma_data = emma_data[emma_data$Strand == 'positive ',]
emma_data$col = 'black'
emma_data[emma_data$Order == 'Nidovirales',]$col = 'deeppink4'
emma_data[emma_data$Family == 'Coronaviridae',]$col = 'red'


Matrix = as.matrix(emma_data[,6:25])

PCA <- prcomp(Matrix, scale = TRUE)

FinalMatrix = data.frame(emma_data$Order, emma_data$Family, emma_data$Species, predict(PCA)[, 1],predict(PCA)[, 2],predict(PCA)[, 3])
FinalMatrix$col = emma_data$col
Min = min(c(FinalMatrix$predict.PCA....1.,FinalMatrix$predict.PCA....2.))
Max = max(c(FinalMatrix$predict.PCA....1.,FinalMatrix$predict.PCA....2.))


PCbiplot <- function(PC, FM, x="PC1", y="PC2") {
  min_a = min(FM$predict.PCA....1.)
  if (min(FM$predict.PCA....2.) < min_a){
    min_a = min(FM$predict.PCA....2.)
  }
  max_a = max(FM$predict.PCA....1.)
  if (max(FM$predict.PCA....2.) > max_a){
    max_a = max(FM$predict.PCA....2.)
  }
  pov <- PCA$sdev^2/sum(PCA$sdev^2)
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)

  loosers = c('P', 'A', 'G', 'R', 'H', 'Q', 'D', 'E', 'T', 'S')
  gainers = c('F', 'I', 'M', 'Y')
  colvec = c()
  for (aa in row.names(PC$rotation)){
    if (aa %in% loosers){
      colvec = append(colvec, 'blue', after = length(colvec))
    }else if (aa %in% gainers){
      colvec = append(colvec, 'red', after = length(colvec))
    }else{
      colvec = append(colvec, 'black', after = length(colvec))
    }
  }
  
  for (aa_num in 1:length(row.names(PC$rotation))){
    row.names(PC$rotation)[aa_num] = aaa(row.names(PC$rotation)[aa_num])
  }
  
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  newpc_l = subset(datapc[datapc$varnames %in% aaa(loosers),], select = -c(varnames))
  newpc_g = subset(datapc[datapc$varnames %in% aaa(gainers),], select = -c(varnames))
  los_row = colSums(newpc_l)/nrow(newpc_l)
  los_row = data.frame(as.list(los_row))
  los_row$varnames = 'LOSERS'
  gain_row = colSums(newpc_g)/nrow(newpc_g)
  gain_row = data.frame(as.list(gain_row))
  gain_row$varnames = 'GAINERS'
  datapc = rbind(datapc, gain_row)
  colvec = append(colvec, 'red4', after = length(colvec))
  datapc = rbind(datapc, los_row)
  colvec = append(colvec, 'blue4', after = length(colvec))
  
  data$col = emma_data$col
  plot <- ggplot(data = FM[FM$col == 'black',], aes(x=predict.PCA....1., y=predict.PCA....2.)) + 
    #stat_density_2d(aes(alpha = ..level..), geom = "polygon")+
    #scale_fill_distiller(palette='Greys', direction=1)+
    geom_point(data = FM[FM$col == 'black',], alpha = 1/25, size=2, aes(x=predict.PCA....1., y=predict.PCA....2., colour = 'other +RNA viruses'), colour='black') + 
    geom_point(data = FM[FM$col == 'deeppink4',], alpha = 1, size=3, aes(x=predict.PCA....1., y=predict.PCA....2., colour = 'other Nidovirales'), colour='deeppink4') + 
   geom_point(data = FM[FM$col == 'red',], size=3, alpha=1, aes(x=predict.PCA....1., y=predict.PCA....2., colour = 'Coronaviridae'), colour='red') + 
    labs(x = sprintf("PC1 (%s)", percent(pov[1], accuracy = 0.01)), y=sprintf("PC2 (%s)",percent(pov[2], accuracy = 0.01)))
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + geom_text(data=datapc, aes(x=v1*1.1, y=v2*1.1, label=varnames), size = 8, vjust=1, color=colvec)
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.5,"cm")), alpha=0.8, color=colvec, size=1.5)
  plot <- plot + coord_equal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold")) +
    scale_x_continuous(sec.axis = dup_axis(), limits = c(-5, 5)) +
    scale_y_continuous(sec.axis = dup_axis(), limits = c(-5, 5))+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=15))
  plot <- plot + tune::coord_obs_pred()
  plot
}

PCbiplot(PCA, FinalMatrix)

dev.off()

##############
