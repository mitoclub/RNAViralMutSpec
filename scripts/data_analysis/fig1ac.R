#Painting 12 and 192 component Mutation Spectrum
library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)
library(stringr)
library(seqinr)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
'%not_in%' <- Negate('%in%')
rm(list=ls(all=TRUE))

#Mutations from GISAID import
#mut_df = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/data/mutation_dists.filtered.splits.csv", row.names=1)
# U filtered data with codons and Aa
setwd("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/")



mut_df = read.csv('new_data/data/U_mut_df_filtred_split.csv', row.names=1)

GainersGegen = c('UUC','UUU','UUG','UUA','AUC','AUU','AUA','UAC','UAU','UAG','UAA');
LoosersGegen = c('CCC','CCA','CCU','CCG','GCC','GCA','GCU','GCG','CGC','CGA','CGU','CGG','GGC','GGA','GGG','GGU','ACC','ACG','ACA','ACU','AGC','AGA','AGU','AGG','CAC','CAU','CAA','CAG','GAC','GAA','GAG','GAU');

Inter_U = c('CUU', 'CUC', 'CUA', 'CUG', 'UCU', 'UCC', 'UCA', 'UCG', 'GUU', 'GUC', 'GUA', 'GUG', 'UGU', 'UGC', 'UGA', 'UGG') 
GainersGegen_u = GainersGegen %>% str_replace_all(c("T" = "U"))
LoosersGegen_u = LoosersGegen %>% str_replace_all(c("T" = "U"))

mut_df = mut_df[mut_df$AaSub == 'NS',]
 
nrow(mut_df[(mut_df$RefCodon %in% LoosersGegen_u) & (mut_df$AltCodon %in% GainersGegen_u),]) + nrow(mut_df[(mut_df$RefCodon %in% Inter_U) & (mut_df$AltCodon %in% GainersGegen_u),]) + nrow(mut_df[(mut_df$RefCodon %in% LoosersGegen_u) & (mut_df$AltCodon %in% Inter_U),])

nrow(mut_df[(mut_df$RefCodon %in% GainersGegen_u) & (mut_df$AltCodon %in% LoosersGegen_u),]) + nrow(mut_df[(mut_df$RefCodon %in% GainersGegen_u) & (mut_df$AltCodon %in% Inter_U),]) + nrow(mut_df[(mut_df$RefCodon %in% Inter_U) & (mut_df$AltCodon %in% LoosersGegen_u),])
#u_mut_df = read.csv('new_data/data/mutation_dists.filtered.splits.csv', row.names=1)

#reading reference data (Wuhan genome NC_045512.2) for normalization
ideal_table = read.csv("new_data/data/U_ideal_aaa.csv", row.names=1)

### Here I add to df some information about codons and AA and save it
#mut_df = subset(u_mut_df, select = -c(mut_id) )
# mut_df = u_mut_df
# mut_df$pos = mut_df$pos+1
# 
# names(mut_df)[1] <- 'Pos'
# mut_df$NeighL = ''
# mut_df$NeighR = ''
# 
# #
# mut_df[mut_df$Pos != 1 & mut_df$Pos != 2 & mut_df$Pos != 29902 & mut_df$Pos != 29903,]$NeighL = toupper(str_split_fixed(mut_df[mut_df$Pos != 1 & mut_df$Pos != 2 & mut_df$Pos != 29902 & mut_df$Pos != 29903,]$parent_nucl_context, '', 5)[,2])
# #
# mut_df[mut_df$Pos != 1 & mut_df$Pos != 2 & mut_df$Pos != 29902 & mut_df$Pos != 29903,]$NeighR = toupper(str_split_fixed(mut_df[mut_df$Pos != 1 & mut_df$Pos != 2 & mut_df$Pos != 29902 & mut_df$Pos != 29903,]$parent_nucl_context, '', 5)[,4])
# #
# ann_gens = ideal_table[, c(1,3,4,5,8)]
# ann_gens = distinct(ann_gens)
# nrow(mut_df)
# sapply(mut_df, class)
# sapply(ann_gens, class)
# mut_df$Pos = as.integer(mut_df$Pos)
# mut_df = merge(x = mut_df, y = ann_gens, by = "Pos", all.x=TRUE)
# nrow(mut_df)
# #
# mut_df$RefCodon = ''
# mut_df$AltCodon = ''
# mut_df[mut_df$NucInCodon == 1,]$RefCodon = toupper(paste(str_split_fixed(mut_df[mut_df$NucInCodon == 1,]$parent_nucl_context, '', 5)[,3],str_split_fixed(mut_df[mut_df$NucInCodon == 1,]$parent_nucl_context, '', 5)[,4],str_split_fixed(mut_df[mut_df$NucInCodon == 1,]$parent_nucl_context, '', 5)[,5],sep=''))
# mut_df[mut_df$NucInCodon == 1,]$AltCodon = toupper(paste(str_split_fixed(mut_df[mut_df$NucInCodon == 1,]$child_nucl_context, '', 5)[,3],str_split_fixed(mut_df[mut_df$NucInCodon == 1,]$child_nucl_context, '', 5)[,4],str_split_fixed(mut_df[mut_df$NucInCodon == 1,]$child_nucl_context, '', 5)[,5],sep=''))
# #
# mut_df[mut_df$NucInCodon == 2,]$RefCodon = toupper(paste(str_split_fixed(mut_df[mut_df$NucInCodon == 2,]$parent_nucl_context, '', 5)[,2],str_split_fixed(mut_df[mut_df$NucInCodon == 2,]$parent_nucl_context, '', 5)[,3],str_split_fixed(mut_df[mut_df$NucInCodon == 2,]$parent_nucl_context, '', 5)[,4],sep=''))
# mut_df[mut_df$NucInCodon == 2,]$AltCodon = toupper(paste(str_split_fixed(mut_df[mut_df$NucInCodon == 2,]$child_nucl_context, '', 5)[,2],str_split_fixed(mut_df[mut_df$NucInCodon == 2,]$child_nucl_context, '', 5)[,3],str_split_fixed(mut_df[mut_df$NucInCodon == 2,]$child_nucl_context, '', 5)[,4],sep=''))
# 
# mut_df[mut_df$NucInCodon == 3,]$RefCodon = toupper(paste(str_split_fixed(mut_df[mut_df$NucInCodon == 3,]$parent_nucl_context, '', 5)[,1],str_split_fixed(mut_df[mut_df$NucInCodon == 3,]$parent_nucl_context, '', 5)[,2],str_split_fixed(mut_df[mut_df$NucInCodon == 3,]$parent_nucl_context, '', 5)[,3],sep=''))
# mut_df[mut_df$NucInCodon == 3,]$AltCodon = toupper(paste(str_split_fixed(mut_df[mut_df$NucInCodon == 3,]$child_nucl_context, '', 5)[,1],str_split_fixed(mut_df[mut_df$NucInCodon == 3,]$child_nucl_context, '', 5)[,2],str_split_fixed(mut_df[mut_df$NucInCodon == 3,]$child_nucl_context, '', 5)[,3],sep=''))
# #
# mut_df$RefAa = ''
# mut_df$AltAa = ''
# ref_string = toString(mut_df[mut_df$GenType == 'translated',]$RefCodon)
# ref_string = str_replace_all(ref_string, ', ', '')
# ref_string = s2c(ref_string)
# mut_df[mut_df$GenType == 'translated',]$RefAa = aaa(translate(seq = ref_string))
# 
# alt_string = toString(mut_df[mut_df$GenType == 'translated',]$AltCodon)
# alt_string = str_replace_all(alt_string, ', ', '')
# alt_string = s2c(alt_string)
# mut_df[mut_df$GenType == 'translated',]$AltAa = aaa(translate(seq = alt_string))
# 
# mut_df[mut_df$RefCodon == 'CTT' | mut_df$RefCodon == 'CTA' | mut_df$RefCodon == 'CTG' | mut_df$RefCodon == 'CTC',]$RefAa = 'Leu_CT'
# mut_df[mut_df$RefCodon == 'TTA' | mut_df$RefCodon == 'TTG',]$RefAa = 'Leu_TT'
# 
# mut_df[mut_df$RefCodon == 'TCT' | mut_df$RefCodon == 'TCA' | mut_df$RefCodon == 'TCG' | mut_df$RefCodon == 'TCC',]$RefAa = 'Ser_TC'
# mut_df[mut_df$RefCodon == 'AGT' | mut_df$RefCodon == 'AGC',]$RefAa = 'Ser_AG'
# 
# mut_df[mut_df$RefCodon == 'CGC' | mut_df$RefCodon == 'CGA' | mut_df$RefCodon == 'CGT' | mut_df$RefCodon == 'CGG',]$RefAa = 'Arg_CG'
# mut_df[mut_df$RefCodon == 'AGG' | mut_df$RefCodon == 'AGA',]$RefAa = 'Arg_AG'
# 
# mut_df[mut_df$RefCodon == 'TAG' | mut_df$RefCodon == 'TAA',]$RefAa = '*_TA'
# mut_df[mut_df$RefCodon == 'TGA',]$RefAa = '*_TG'
# 
# mut_df[mut_df$AltCodon == 'CTT' | mut_df$AltCodon == 'CTA' | mut_df$AltCodon == 'CTG' | mut_df$AltCodon == 'CTC',]$AltAa = 'Leu_CT'
# mut_df[mut_df$AltCodon == 'TTA' | mut_df$AltCodon == 'TTG',]$AltAa = 'Leu_TT'
# 
# mut_df[mut_df$AltCodon == 'TCT' | mut_df$AltCodon == 'TCA' | mut_df$AltCodon == 'TCG' | mut_df$AltCodon == 'TCC',]$AltAa = 'Ser_TC'
# mut_df[mut_df$AltCodon == 'AGT' | mut_df$AltCodon == 'AGC',]$AltAa = 'Ser_AG'
# 
# mut_df[mut_df$AltCodon == 'CGC' | mut_df$AltCodon == 'CGA' | mut_df$AltCodon == 'CGT' | mut_df$AltCodon == 'CGG',]$AltAa = 'Arg_CG'
# mut_df[mut_df$AltCodon == 'AGG' | mut_df$AltCodon == 'AGA',]$AltAa = 'Arg_AG'
# 
# mut_df[mut_df$AltCodon == 'TAG' | mut_df$AltCodon == 'TAA',]$AltAa = '*_TA'
# mut_df[mut_df$AltCodon == 'TGA',]$AltAa = '*_TG'
# 
# mut_df$AaSub = ''
# mut_df[mut_df$GenType == 'translated',]$AaSub = ifelse(mut_df[mut_df$GenType == 'translated',]$RefAa == mut_df[mut_df$GenType == 'translated',]$AltAa, 'S', 'NS')
# #
# names(mut_df)[2] <- 'RefNuc'
# names(mut_df)[3] <- 'AltNuc'
# mut_df$RefNuc <- gsub('T', 'U', mut_df$RefNuc)
# mut_df$AltNuc <- gsub('T', 'U', mut_df$AltNuc)
# mut_df$NeighL <- gsub('T', 'U', mut_df$NeighL)
# mut_df$NeighR <- gsub('T', 'U', mut_df$NeighR)
# mut_df$RefCodon <- gsub('T', 'U', mut_df$RefCodon)
# mut_df$RefCodon <- gsub('T', 'U', mut_df$RefCodon)
# mut_df$RefCodon <- gsub('T', 'U', mut_df$RefCodon)
# mut_df$AltCodon <- gsub('T', 'U', mut_df$AltCodon)
# mut_df$AltCodon <- gsub('T', 'U', mut_df$AltCodon)
# mut_df$AltCodon <- gsub('T', 'U', mut_df$AltCodon)
# mut_df$RefAa <- gsub('T', 'U', mut_df$RefAa)
# mut_df$RefAa <- gsub('T', 'U', mut_df$RefAa)
# mut_df$AltAa <- gsub('T', 'U', mut_df$AltAa)
# mut_df$AltAa <- gsub('T', 'U', mut_df$AltAa)
# #
# mut_df[mut_df$RefCodon == 'ACC' | mut_df$RefCodon == 'ACG' | mut_df$RefCodon == 'ACA' | mut_df$RefCodon == 'ACU',]$RefAa = 'Thr'
# mut_df[mut_df$AltCodon == 'ACC' | mut_df$AltCodon == 'ACG' | mut_df$AltCodon == 'ACA' | mut_df$AltCodon == 'ACU',]$AltAa = 'Thr'
# 
# mut_df[mut_df$RefCodon == 'UGG',]$RefAa = 'Trp'
# mut_df[mut_df$AltCodon == 'UGG',]$AltAa = 'Trp'
# 
# mut_df[mut_df$RefCodon == 'UAU' | mut_df$RefCodon == 'UAC',]$RefAa = 'Tyr'
# mut_df[mut_df$AltCodon == 'UAU' | mut_df$AltCodon == 'UAC',]$AltAa = 'Tyr'
# 
# 
# #
# 
# write.csv(x = mut_df, 'new_data/data/U_mut_df_filtred_split.csv')

###


# all_codons = mut_df[mut_df$GenType=='translated',]
# #mut_df = mut_df[mut_df$Pos != 1 & mut_df$Pos != 2 & mut_df$Pos != 29902 & mut_df$Pos != 29903,]
# all_codons$FromWithNeigh = paste(all_codons$NeighL,all_codons$RefNuc,all_codons$NeighR, sep='')
# all_codons$ToWithNeigh = paste(all_codons$NeighL,all_codons$AltNuc,all_codons$NeighR, sep='')
# all_codons$MutSpecWithNeigh = paste(all_codons$FromWithNeigh, all_codons$ToWithNeigh, sep = '>')
# all_mut_with_neigh_list = unique(all_codons$MutSpecWithNeigh)
# 
# 
# # Now I am selecting necessary data (similar for each dataframe, for exapmle if I need four fold mutation, I choose them from Gisaid data and from RefSeq for normalization)
# 
# ff_aa = c('Gly', 'Pro', 'Ser_UC', 'Ala', 'Thr', 'Leu_CU', 'Val', 'Arg_CG')
# 
# # FF
# #mut_df = mut_df[(mut_df$RefAa %in% ff_aa) & (mut_df$NucInCodon == 3),]
# #ideal_table = ideal_table[(ideal_table$RefAa %in% ff_aa) & (ideal_table$NucInCodon == 3),]
# 
# #nuc_all_usage = as.data.frame(table(ideal_table['RefNuc']) / 3)
# #names(nuc_all_usage)=c('Nuc', 'Count')
# #write.csv(x = nuc_all_usage, 'new_data/data_obtained/Nuc_FF_usage.csv')
# 
# #Syn
# mut_df = mut_df[mut_df$AaSub=='S',]
# ideal_table = ideal_table[ideal_table$AaSub=='S',]


draw_mutspec = function(df, ann, n_mutspec = 12, all_mut_list=all_mut_with_neigh_list, prefix = 'All'){
  
  all_codons = df[df$GenType=='translated',]
  #mut_df = mut_df[mut_df$Pos != 1 & mut_df$Pos != 2 & mut_df$Pos != 29902 & mut_df$Pos != 29903,]
  all_codons$FromWithNeigh = paste(all_codons$NeighL,all_codons$RefNuc,all_codons$NeighR, sep='')
  all_codons$ToWithNeigh = paste(all_codons$NeighL,all_codons$AltNuc,all_codons$NeighR, sep='')
  all_codons$MutSpecWithNeigh = paste(all_codons$FromWithNeigh, all_codons$ToWithNeigh, sep = '>')
  all_mut_with_neigh_list = unique(all_codons$MutSpecWithNeigh)
  
  ff_aa = c('Gly', 'Pro', 'Ser_UC', 'Ala', 'Thr', 'Leu_CU', 'Val', 'Arg_CG')
  if (prefix == 'Syn'){
    df = df[df$AaSub=='S',]
    ann = ann[ann$AaSub=='S',]
  }else if (prefix == 'FF'){
    df = df[(df$RefAa %in% ff_aa) & (df$NucInCodon == 3),]
    ann = ann[(ann$RefAa %in% ff_aa) & (ann$NucInCodon == 3),]
  }else{
    df = df
    ann = ann
  }
  coul <- brewer.pal(12, "Paired")
  mutations = c('A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G')
  if(n_mutspec == 12){
    mut_list = data.frame(mutations)
    names(mut_list)=c('NucSubst')
    ann$MutExp <- 1
    ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
    Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
    names(Exp12Comp) = c('NucSubst','ExpFr')
    # Here are expected mutation from RefSeq, we will normalize Gisaid data on it
    
    df$MutObs = 1
    df$NucSubst = paste(df$RefNuc,df$AltNuc,sep='>')
    Obs12Comp = aggregate(df$MutObs, by = list(df$NucSubst), FUN = sum);
    names(Obs12Comp) = c('NucSubst','ObsFr')
    # Here are all mutations from Gisaid data
    
    MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
    MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
    MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
    MutSpec12Comp$MutSpec = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr  #Normalization on all possible same mutations in RefSeq
    MutSpec12Comp$MutSpec[is.na(MutSpec12Comp$MutSpec)] <- 0
    write.csv(x = MutSpec12Comp, sprintf('new_data/data_obtained/07.%s_MutSpec_12comp.csv', prefix))
    
    #now plotting
    #png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/new_data/figures/07.FF_MutSpec_12comp.png')
    #barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for FF nucleotides', cex.names = 0.7, col=coul)
    #dev.off()
    fig <- ggplot(MutSpec12Comp,aes(x = NucSubst, y = MutSpec, ymin = 0, fill = NucSubst))+
      geom_bar(stat = "identity", width = 0.6)+ 
      #coord_flip() + 
      xlab("Mutations") +
      ylab("Substitution rate (observed / expected)")+
      theme(legend.position="none")+
      #theme(axis.text.x = element_text(size=3.2, angle = 90))+
      ggtitle(sprintf('Covid Mut spec 12 components for %s nucleotides', prefix)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      theme(axis.title=element_text(size=20)) +
      theme(axis.text=element_text(size=20))+ 
      scale_y_continuous(expand = c(0,0))
    #theme(axis.text.y = element_text(hjust = 3.4))+
    #geom_text(aes(label=NucSubst), size = 1.4, angle = 90)
    fig + scale_fill_brewer(palette = "Paired")
    fig %>% ggplotly
    ggsave(sprintf('new_data/figures/07.%s_MutSpec_12comp.pdf', prefix), plot = fig, device = 'pdf',  height = 10 , width = 12)
  }else if(n_mutspec == 192){
    mut_list = all_mut_list
    mut_list = data.frame(mut_list)
    names(mut_list)=c('NucSubst')
    
    ann$MutExp = 1
    ann$FromWithNeigh = paste(ann$NeighL,ann$RefNuc,ann$NeighR, sep='')
    ann$ToWithNeigh = paste(ann$NeighL,ann$AltNuc,ann$NeighR, sep='')
    ann$MutSpecWithNeigh = paste(ann$FromWithNeigh, ann$ToWithNeigh, sep = '>')
    Exp192Comp = aggregate(ann$MutExp, by = list(ann$MutSpecWithNeigh), FUN = sum);
    names(Exp192Comp) = c('NucSubst','ExpFr')
    ann$Colour = paste(ann$RefNuc, ann$AltNuc, sep = '>')
    Exp_with_colour = aggregate(ann$Colour, by = list(ann$MutSpecWithNeigh), FUN=unique)
    names(Exp_with_colour) = c('NucSubst','Colour')
    Exp192Comp = merge(Exp_with_colour, Exp192Comp, by = 'NucSubst', all.x = TRUE)
    Exp192Comp = Exp192Comp[order(Exp192Comp$Colour),]
    rownames(Exp192Comp) <- 1:nrow(Exp192Comp)
    Exp192Comp$Colour = factor(Exp192Comp$Colour)
    
    df$MutObs = 1
    df$FromWithNeigh = paste(df$NeighL,df$RefNuc,df$NeighR, sep='')
    df$ToWithNeigh = paste(df$NeighL,df$AltNuc,df$NeighR, sep='')
    df$MutSpecWithNeigh = paste(df$FromWithNeigh, df$ToWithNeigh, sep = '>')
    Obs192Comp = aggregate(df$MutObs, by = list(df$MutSpecWithNeigh), FUN = sum);
    names(Obs192Comp) = c('NucSubst','ObsFr')
    ThreeNuc = rep(c('A','C','G','U'),times=48)
    FirstNuc = rep(c('A','C','G','U'),times=12)
    
    MutSpec192Comp = merge(mut_list,Obs192Comp, by = 'NucSubst', all.x = TRUE)
    MutSpec192Comp = merge(MutSpec192Comp,Exp192Comp, by = 'NucSubst', all.x = TRUE)
    
    MutSpec192Comp$Colour = paste(str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,2],str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,6],sep='>')
    
    MutSpec192Comp$ObsToExp = MutSpec192Comp$ObsFr/MutSpec192Comp$ExpFr
    MutSpec192Comp$ObsToExp[is.na(MutSpec192Comp$ObsToExp)] <- 0
    #MutSpec192Comp$MutSpec = MutSpec192Comp$ObsToExp/sum(MutSpec192Comp$ObsToExp)
    MutSpec192Comp$Context = paste(str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,1], 'x', str_split_fixed(MutSpec192Comp$NucSubst, '', 7)[,3],sep='')
    write.csv(x = MutSpec192Comp, sprintf('new_data/data_obtained/07.%s_MutSpec_192comp.csv', prefix))
    
    fig <- ggplot(MutSpec192Comp,aes(x = reorder(NucSubst, Colour, FUN = unique, order = is.ordered(NucSubst)), y = ObsToExp, ymin = 0, fill = Colour))+
      geom_bar(stat = "identity", width = 1, color = "black")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      #coord_flip() + 
      xlab("Mutations") +
      ylab("Substitution rate (observed / expected)")+
      theme(legend.position="none")+
      #theme(axis.text.x = element_text(size=3.2, angle = 90))+
      ggtitle(sprintf('Covid Mut spec 192 components for %s nucleotides', prefix)) +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      theme(axis.title=element_text(size=20)) +
      theme(axis.text=element_text(size=20)) +
      #annotate(geom = "text", x = seq_len(nrow(MutSpec192Comp)), y = -5, label = ThreeNuc, size = 2)+
      #annotate(geom = "text", x = 2.5 + 4 * (0:47), y = -10, label = FirstNuc, size = 3)+
      geom_vline(xintercept = seq(16.5, 192, by=16))+
      #annotate(geom = "text", x = 4.5 + 4 * (0:47), y = -5, label = rep(c('|'),times=48), size = 3)+
      annotate(geom = "text", x = 8.5 + 16 * (0:11), y = -5, label = mutations, size = 5) +
      theme_bw() + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme(#axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank()
            #panel.grid.minor = element_blank(),
            #panel.border = element_blank(),
            #panel.background = element_blank()
            )
      #+
      #annotate(geom = "text", x = -5, y = -4, label = 'Third Nucleotide', size = 2)+
      #annotate(geom = "text", x = -5, y = -10, label = 'First Nucleotide', size = 3)
      #theme(axis.text.y = element_text(hjust = 3.4))+
      #geom_text(aes(label=NucSubst), size = 1.4, angle = 90)
    fig %>% ggplotly
    options(repr.plot.width=20, repr.plot.height=8)
    ggsave(sprintf('new_data/figures/07.%s_MutSpec_192comp.pdf', prefix), plot = fig, device = 'pdf', width=17, height=10,)
  }
}

draw_mutspec(mut_df, ideal_table, n_mutspec = 192, prefix='FF')


ideal_table = ideal_table[ideal_table$GenType=='translated',]
ideal_table$codon_count=1
codon_ref = aggregate(ideal_table$codon_count, by = list(ideal_table$RefCodon), FUN = sum);
names(codon_ref) = c('Codons','Number_of_codons')
codon_ref$Number_of_codons = codon_ref$Number_of_codons/9
codon_ref$Number_of_codons = round(codon_ref$Number_of_codons)
write.csv(codon_ref, 'new_data/data_obtained/Ref_Codon_count.csv')




# Mut spec by date
# mutations = c('A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G')
# mut_spec_date = data.frame(mutations)
# names(mut_spec_date)=c('NucSubst')
# ann = ideal_table
# for (i in unique(mut_df$date_split)){
#   df = mut_df[mut_df$date_split == i,]
#   mut_list = data.frame(mutations)
#   names(mut_list)=c('NucSubst')
#   ann$MutExp <- 1
#   ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
#   Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
#   names(Exp12Comp) = c('NucSubst','ExpFr')
#   # Here are expected mutation from RefSeq, we will normalize Gisaid data on it
# 
#   df$MutObs = 1
#   df$NucSubst = paste(df$RefNuc,df$AltNuc,sep='>')
#   Obs12Comp = aggregate(df$MutObs, by = list(df$NucSubst), FUN = sum);
#   names(Obs12Comp) = c('NucSubst','ObsFr')
#   # Here are all mutations from Gisaid data
# 
#   MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
#   MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
#   MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
#   MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr  #Normalization on all possible same mutations in RefSeq
#   MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
#   MutSpec12Comp[toString(i)] = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
#   MutSpec12Comp = MutSpec12Comp[,c('NucSubst', toString(i))]
#   mut_spec_date = merge(mut_spec_date,MutSpec12Comp, by = 'NucSubst', all.x = TRUE)
# }
# 
# 
# mut_spec_date = mut_spec_date[,c('NucSubst', '0', '1', '2', '3', '4', '5', '6', '7', '8')]
# names(mut_spec_date) = c('NucSubst','2019-12_-_2020-04', '2020-05_-_2020-07','2020-08_-_2020-10','2020-11_-_2020-12','2021-01_-_2021-02','2021-03_-_2021-04','2021-05_-_2021-06','2021-07_-_2021-08','2021-09_-_2021-10')
# 
# molten <- melt(setDT(mut_spec_date), id.vars = c("NucSubst"))
# write.csv(x = molten, 'new_data/data_obtained/MutChange_12spec_FF_normed_byone.csv')
# fig = ggplot(molten, aes(x = variable, y = as.numeric(value), group = NucSubst, colour = NucSubst)) +
#   theme(axis.text.x = element_text(angle = 90))+
#   xlab("Dates") +
#   ylab("Share Of Mutations")+
#   ggtitle('Mutation change over time')+
#   geom_line() + 
#   geom_point()
# ggsave('new_data/figures/MutChange_12spec_normed_one.svg', plot = fig, device = 'svg',  limitsize = F)
