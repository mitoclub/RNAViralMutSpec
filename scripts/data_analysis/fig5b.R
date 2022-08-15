rm(list=ls(all=TRUE))
library(seqinr)
setwd("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/")
pdf("new_data/figures/24.CodonUsageParser.PCA.RealPca_LG_plus.pdf", height = 21.36, width = 20.7)
#### read summary file with all viruses
Sum = read.csv("new_data/data/CodonUsage/NcbiVirusHumanRna2.csv", header= TRUE)
for (species_segment in unique(Sum[Sum$Segment != '',]$Species)){
  accesions_for_segm = (Sum[Sum$Species == species_segment,]$Accession)
  full_file = read.csv(sprintf('new_data/data/CodonUsage/AllFromSergey2/%s.csv', accesions_for_segm[1]),sep = ',')
  names(full_file)[1] = c('Codons')
  if (ncol(full_file) > 2){
    full_file$segment_count = rowSums(full_file[, c(2:ncol(full_file))])
    full_file = full_file[,c('Codons', 'segment_count')]
  }else{
    names(full_file)[2] = c('segment_count')
  }
  for (accesions_num in (2:length(accesions_for_segm))){
    new_segment = read.csv(sprintf('new_data/data/CodonUsage/AllFromSergey2/%s.csv', accesions_for_segm[accesions_num]),sep = ',')
    names(new_segment)[1] = c('Codons')
    if (ncol(new_segment) > 2){
      new_segment[,sprintf('segment_count_%s', accesions_num)] = rowSums(new_segment[, c(2:ncol(new_segment))])
      new_segment = new_segment[,c('Codons', sprintf('segment_count_%s', accesions_num))]
    }else{
      names(new_segment)[2] = sprintf('segment_count_%s', accesions_num)
    }
    full_file = merge(full_file,new_segment, by = 'Codons', all = TRUE)
    full_file$segment_count = full_file$segment_count + new_segment[, sprintf('segment_count_%s', accesions_num)]
  }
  full_file = full_file[,c('Codons', 'segment_count')]
  write.csv(full_file, sprintf('new_data/data/CodonUsage/without_duplicates2/%s.csv', accesions_for_segm[1]))
}
NoDubDir = "new_data/data/CodonUsage/without_duplicates2/"
VecOfFilesCompared = list.files(path = NoDubDir);

Sum = Sum[!duplicated(Sum$Species),]
for (accession in Sum$Accession){
  if (paste(accession, '.csv', '') %in% VecOfFilesCompared){
    print('Уже есть такой файл')
  }else{
    file = read.csv(sprintf('new_data/data/CodonUsage/AllFromSergey2/%s.csv', accession),sep = ',')
    write.csv(file, sprintf('new_data/data/CodonUsage/without_duplicates2/%s.csv', accession))
  }
}
# merge this file with Sergeys files

VecRnaPlus = c('ssRNA(+)')
VecRnaMinus = c('ssRNA(-)')

InputDir = "new_data/data/CodonUsage/without_duplicates2/"
VecOfFiles = list.files(path = InputDir); length(VecOfFiles) 
Final=data.frame()

for (i in 1:length(VecOfFiles))
{ # i = 16   
  if (file.info(paste(InputDir,VecOfFiles[i],sep = ''))$size > 0)
  {
    #Vir = read.csv(paste(InputDir,VecOfFiles[i],sep = ''))  # }
    
    CU = read.csv(paste(InputDir,VecOfFiles[i],sep = ''), row.names=1)
    names(CU)[1] = c('Codon')
    CU$NumbTotal = apply(as.matrix(CU[,seq(2,ncol(CU))]),1,FUN = sum)
    VirRefSeqId = gsub(".csv",'',VecOfFiles[i])
    VirSpeciesCommonName = Sum[Sum$Accession == VirRefSeqId,]$Species
    VirSpeciesGenBankTitle = Sum[Sum$Accession == VirRefSeqId,]$GenBank_Title
    Molecule_type = Sum[Sum$Accession == VirRefSeqId,]$Molecule_type
    VirFamily = Sum[Sum$Accession == VirRefSeqId,]$Family
    VirType = 'NA'
    if (Molecule_type %in% VecRnaPlus) {VirType = 'ssRnaPlus'
    }else if (Molecule_type %in% VecRnaMinus) {VirType = 'ssRnaMinus'
    }else {VirType = 'DS/Ambi'}
    
    Title = paste(VirRefSeqId,VirSpeciesCommonName,VirSpeciesGenBankTitle,VirFamily,VirType,sep = '\n')
    
    CU$Amino = ''
    for (j in 1:nrow(CU))
    { # j = 1
      CU$Amino[j] = aaa(translate(s2c(as.character(CU$Codon[j]))))
    }
    
    CU$names = paste(CU$Amino,CU$Codon, sep = ' : ')
    
    CU$Codon  <- factor(CU$Codon , levels = c("AAG","AAA","AAC","AAT","ACG","ACA","ACC","ACT","AGG","AGA","AGC","AGT","ATG","ATA","ATC","ATT","CAG","CAA","CAC","CAT","CCG","CCA","CCC","CCT","CGG","CGA","CGC","CGT","CTG","CTA","CTC","CTT","GAG","GAA","GAC","GAT","GCG","GCA","GCC","GCT","GGG","GGA","GGC","GGT","GTG","GTA","GTC","GTT","TAG","TAA","TAC","TAT","TCG","TCA","TCC","TCT","TGG","TGA","TGC","TGT","TTG","TTA","TTC","TTT"))
    CU = CU[order(CU$Codon),]
    ColVec=c(rgb(0,0,1,0.25),rgb(0,0,1,0.75),rgb(1,0,0,0.25),rgb(1,0,0,0.75))
    # par(mfrow=c(2,1))
    # barplot(CU$NumbTotal, names = (CU$names), las = 2, col = ColVec, main = Title) # dev.off()
    
    #### relative bias:
    AaVec = unique(CU$Amino); length(AaVec)
    ExtractThird <- function(x) {unlist(strsplit(x,''))[3]}
    CU$ThirdNuc = apply(as.matrix(as.character(CU$Codon)),1,FUN=ExtractThird)
    
    GainersGegen = c('TTC','TTT','TTG','TTA','ATC','ATT','ATA','TAC','TAT','TAG','TAA');
    LoosersGegen = c('CCC','CCA','CCT','CCG','GCC','GCA','GCT','GCG','CGC','CGA','CGT','CGG','GGC','GGA','GGG','GGT','ACC','ACG','ACA','ACT','AGC','AGA','AGT','AGG','CAC','CAT','CAA','CAG','GAC','GAA','GAG','GAT');
    GainLooseCU = CU
    GainLooseCU$type = ''
    GainLooseCU[GainLooseCU$Codon %in% GainersGegen,]$type = 'Gainer'
    GainLooseCU[GainLooseCU$Codon %in% LoosersGegen,]$type = 'Looser'
    
    GainLooseCUAgg = aggregate(GainLooseCU$NumbTotal, by = list(GainLooseCU$type), FUN = sum)
    
    names(GainLooseCUAgg)=c('codon_type','number')
    
    
    GainLooseCUAgg$Fr = GainLooseCUAgg$number/sum(GainLooseCUAgg$number)
    
    GainFr = GainLooseCUAgg[GainLooseCUAgg$codon_type == 'Gainer',]$Fr
    LoosFr = GainLooseCUAgg[GainLooseCUAgg$codon_type == 'Looser',]$Fr
  }  
  OneLine = c(GainFr,LoosFr,VirRefSeqId,VirSpeciesCommonName,VirFamily,VirType)
  Final = rbind(Final,OneLine)
  names(Final) = c('GainFr','LoosFr','VirRefSeqId','VirSpeciesCommonName','VirFamily','VirType')
}

Final = Final[Final$VirType == 'ssRnaMinus' | Final$VirType == 'ssRnaPlus',]
#### ANALYSIS:
Final$GainFr = as.numeric(Final$GainFr)
Final$LoosFr = as.numeric(Final$LoosFr)

FinalAgg = aggregate(list(Final$GainFr,Final$LoosFr), by = list(Final$VirSpeciesCommonName,Final$VirFamily,Final$VirType), FUN = mean)
names(FinalAgg)=c('VirSpeciesCommonName','VirFamily','VirType','GainFr','LoosFr')
FinalAgg = FinalAgg[FinalAgg$VirType == 'ssRnaPlus',]

FinalAgg$Col = 'black'
FinalAgg[FinalAgg$VirFamily == 'Coronaviridae',]$Col = 'red'

write.csv(x = FinalAgg, 'new_data/data_obtained/23.VirusesLG.csv')

plot(FinalAgg[FinalAgg$VirType == 'ssRnaPlus',]$LoosFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus',]$GainFr, col = 'red', xlab = 'Loosers Frequency', ylab = 'Gainers Frequency', pch = 16, xlim = c(0.35,0.6), ylim = c(0.1,0.25), cex.lab = 1.7, cex.axis = 1.7)
par(new=TRUE)
textplot(FinalAgg[FinalAgg$VirType == 'ssRnaPlus',]$LoosFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus',]$GainFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus',]$VirSpeciesCommonName, xlab = '', ylab = '', cex=2, xlim = c(0.35,0.6), ylim = c(0.1,0.25), col = FinalAgg$Col, xaxt = "n", yaxt = "n")
#textplot(FinalAgg[FinalAgg$VirType == 'ssRnaPlus'& FinalAgg$Col == '',]$LoosFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus'& FinalAgg$Col == '',]$GainFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus'& FinalAgg$Col == '',]$VirSpeciesCommonName, xlab = 'Loosers Frequency', ylab = 'Gainers Frequency', cex=2, col= 'black', xlim = c(0.3,0.6), ylim = c(0.1,0.3))
#par(new=TRUE)
#textplot(FinalAgg[FinalAgg$VirType == 'ssRnaPlus'& FinalAgg$Col == 'SARS',]$LoosFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus'& FinalAgg$Col == 'SARS',]$GainFr,FinalAgg[FinalAgg$VirType == 'ssRnaPlus'& FinalAgg$Col == 'SARS',]$VirSpeciesCommonName, xlab = 'Loosers Frequency', ylab = 'Gainers Frequency', cex=2, col= 'red', xlim = c(0.3,0.6), ylim = c(0.1,0.3))

#aes(colour = factor(VirType)), size = 3
#LGplot <- ggplot(FinalAgg, aes(x= LoosFr, y = GainFr)) + 
#  geom_point() + ggtitle("The shares of gainers and losers among virus codons") +
#  xlab("Loosers Frequency") + 
#  ylab("Gainers Frequency") + 
#  theme_classic() + theme(text = element_text(size = 40)) + theme(axis.text = element_text(size = 20)) 

### geom_label_repel
#LGplot + 
#  geom_label_repel(size = 8, aes(label = VirSpeciesCommonName))

dev.off()
