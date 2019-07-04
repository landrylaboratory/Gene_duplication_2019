#######################################################################################
#                            Script to analyse our PCA                                #
#######################################################################################

#Colonies flagged as irregular by gitter (as S or S, C flags) or that did not grow on the
#last diploid selection step or on DMSO medium were considered as missing data. 
#We considered only bait-prey pairs with at least four replicates and used the median of 
#colony sizes as PCA signal. The data was finally filtered based on the completeness of 
#paralogous pairs so we could test HMs and HETs systematically. 

rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(vioplot)
library(circlize)
library(VennDiagram)
library(mixtools)
library(magrittr)

#to set the working directory
setwd("dir")

df.med=read.table("output/PCA_med_data_2019_02.tab", sep="\t", header=T) #1328

#Filter SSD with MSA sequence similarity:
#If sequence similarity < 20 %, paralogs need to be in the same phylome (from phylomDB)
#to be selected

MSA_seq_ident_under_20 <- read.table("data/MSA_seq_ident_under_20.txt", sep="\t", header=T)
MSA_seq_ident_under_20 <- filter(MSA_seq_ident_under_20, Duplication=="SSD")
MSA_seq_ident_under_20$pair<- apply(cbind(as.character(MSA_seq_ident_under_20$P1), 
                             as.character(MSA_seq_ident_under_20$P2)), 
                       1, function(x) paste(sort(x), collapse="."))

MSA_seq_ident_under_20 <- select(MSA_seq_ident_under_20, pair, Same_phylome)
df.med <- left_join(df.med, MSA_seq_ident_under_20, by='pair')
df.med <- filter(df.med, is.na(Same_phylome) | Same_phylome==1)


#To check zize distribution
colplot=c("#00A44A", "#28F040", "#0083B3", "#00F5FF", "#003B63",
          "#7D0A90", "#CB10EA",  "#C98EF6", 
          "#8C0011", "#FF011B", "#FF7605", "#FDEF33")


df.med %$% plot(density(med_MTX[type=="HM" & Duplication=="SSD"]), col=colplot[1], 
                              main = "Size of colonies 30C : density depending dimers", 
                xlim = c(5,11), ylim = c(0, 0.6))
df.med %$% lines(density(med_MTX[type=="HET" & Duplication=="SSD"]), col=colplot[3])
df.med %$% lines(density(med_MTX[type=="HM" & Duplication=="WGD"]), col=colplot[11])
df.med %$% lines(density(med_MTX[type=="HET" & Duplication=="WGD"]), col=colplot[12])
legend("topright", c("HM_SSD", "HET_SSD", "HM_WGD", "HET_WGD"), 
       lty=c(1,1), lwd=c(2.0,2.0), col=c(colplot[1], colplot[3], colplot[11], colplot[12]))



##############################################################################################################
# Determination of a threshold of size colony which represents the presence of protein - protein interaction #                         #
##############################################################################################################

prob.type=c()

plot_normalmixEM <- function(data)
{
  plot(density(data$med_MTX), main=paste("Density - ", dupli, sep=""))
  abline(v=16, col="gray")
}

pdf("without_low_qual/normalmixEM_for_threshold_2019_02.pdf")
for (dupli in unique(df.med$Duplication)) 
{
  dupli.subset=filter(df.med, Duplication==dupli)
  plot_normalmixEM(dupli.subset)
  # Exclude values that are less than 5, to have 2 distributions
  df.med.sup5=subset(df.med, df.med$med_MTX > 5)
  
  for (type in unique(df.med.sup5$type)) 
  {
    type.subset=df.med.sup5[df.med.sup5$type == type,]
    test=normalmixEM(type.subset$med_MTX)
    
    plot(test, whichplots=2, main2 = paste("Density Curves", dupli, type, sep=" - "))
    lines(density(type.subset$med_MTX), lty=2, lwd=2)
    
    if (test$mu[1] < test$mu[2])
      # mu: a k-vector of initial med_MTXs for the mean parameters
    {
      type.subset2=cbind(dupli, type, test$mu[1], test$sigma[1])
      # sigma: a k-vector of initial med_MTXs for the standard deviation parameters
      type.subset2=as.data.frame(type.subset2)
      names(type.subset2)=c("Duplication", "type", "mean.type", "sd.type")
    }
    
    if (test$mu[2] < test$mu[1]) 
    {
      type.subset2=cbind(dupli, type, test$mu[2], test$sigma[2])
      type.subset2=as.data.frame(type.subset2)
      names(type.subset2)=c("Duplication", "type", "mean.type", "sd.type")
    }
    prob.type=rbind(prob.type, type.subset2)
  }
}
dev.off()


head(prob.type)
nrow(prob.type)

prob.type$Duplication=as.character(prob.type$Duplication)
prob.type$type=as.character(prob.type$type)
prob.type$mean.type=as.numeric(as.character(prob.type$mean.type))
prob.type$sd.type=as.numeric(as.character(prob.type$sd.type))

prob.type$seuil=prob.type$mean.type + 2.5*(prob.type$sd.type)
prob.type

df.med.seuil <- left_join(df.med, prob.type, by=c("Duplication", "type"))

#z-score
df.med.seuil %<>% mutate (Zscore = (med_MTX - mean.type) / sd.type)

#to determine if there are interaction or not
df.med.seuil %<>% mutate (interaction = ifelse(Zscore>2.5, 1, 0))
(nrow(filter(df.med.seuil,interaction==1))/nrow(df.med.seuil))*100 #37.59% of ppi tested are positive
write.table(df.med.seuil, file="without_low_qual/PCA_med_seuil_data_2019_06.tab", sep="\t",  quote=F, row.names=F)

#selection of complete pairs: with both HM tested and at least 1 HET
df.med.paralogs <- df.med.seuil

df.nbHET.tot = df.med.paralogs %>% group_by(pair) %>% filter(type=="HET") %>% 
  summarise(nbHET.tot = n()) %>% as.data.frame()

df.nbHM.tot = df.med.paralogs  %>% group_by(pair) %>% filter(type=="HM") %>% 
  summarise(nbHM.tot = n()) %>% as.data.frame()
df.ntype = full_join(df.nbHM.tot, df.nbHET.tot, by="pair")
nrow(df.ntype) #294 pairs

df.med.paralogs  = left_join(df.med.paralogs , df.ntype, by="pair")
df.compl.pairs = filter(df.med.paralogs, nbHM.tot==2 & !is.na(nbHET.tot))
nrow(df.compl.pairs) #1072
nrow(df.compl.pairs) / 4 #268 pairs fully tested
nrow(filter(df.compl.pairs, nbHET.tot==1)) #0


#####################################################################################################
#              Looking for heterodimers found in a direction and not in the other                   #  
#               (P1.DHFR[1-2] - P2.DHFR[3] but not P1.DHFR[3] - P2.DHFR[1-2] or the contrary)       #
#####################################################################################################
#If we observe HET in a direction and not in an other (ex: P1.prey-P2.bait shows PPI but not P2.prey-P1.bait),
#it could potentially significate that one DHFR insertion was faild (here P2.prey or P1.bait)
#in this case, we can't observe HM of the bad strain

#Here, we check if the case with only one direction showing HET ppi are significatively more associated
#with only one homodimer observation

#check in pairs showing interactions
df.nbHET.inter = df.med.paralogs %>% group_by(pair) %>% filter(interaction==1) %>% 
  filter(type=="HET") %>% summarise(nbHET.inter = n()) %>% as.data.frame()

only1HET <- filter(df.nbHET.inter, nbHET.inter==1)

(nrow(only1HET)/nrow(df.nbHET.inter))*100 
#16.98% with only one of 2 directions showing ppi detection

#####################################################################################################
#                     Determine motif interaction pattern for each pairs                           #
#####################################################################################################
#motif NI : no interaction
#motif 1HM : one of the 2 paralogs forms HM / no HET
#motif 2HM : both paralogs form HM / no HET
#motif HET : HET / no HM
#motif 1HM.HET : one of the 2 paralogs forms HM / HET
#motif 2HM.HET : both paralogs form HM / HET

P <- colsplit(df.compl.pairs$pair, "[.]", names=c("P1", "P2"))
df.compl.pairs <- cbind(df.compl.pairs, P)

df.compl.pairs %<>% group_by(Duplication, pair, Mat_A, Mat_alpha) %>% 
  mutate(ppi_type=ifelse(ppi==paste(c(P1,P1), collapse="."), "HM1.PCA",
                  ifelse(ppi==paste(c(P2,P2), collapse="."), "HM2.PCA",
                  ifelse(ppi==paste(c(P1,P2), collapse="."), "HET.P1P2",
                  ifelse(ppi==paste(c(P2,P1), collapse="."), "HET.P2P1", "pb"))))) %>% as.data.frame()

interac <- df.compl.pairs %>% select(.,Duplication, pair, P1, P2, ppi_type, interaction)
interac %<>% spread(., ppi_type, interaction)

medMTX <- df.compl.pairs %>% select(.,Duplication, pair, P1, P2, ppi_type, med_MTX)
medMTX %<>% spread(., ppi_type, med_MTX)
colnames(medMTX) <- c("Duplication", "pair", "P1", "P2", 
                      "med.MTX.HET.P1P2", "med.MTX.HET.P2P1", "med.MTX.HM.P1", "med.MTX.HM.P2")

rep <- df.compl.pairs %>% select(.,Duplication, pair, P1, P2, ppi_type, replicats)
rep %<>% spread(., ppi_type, replicats)
colnames(rep) <- c("Duplication", "pair", "P1", "P2", 
                      "n.rep.HET.P1P2", "n.rep.HET.P2P1", "n.rep.HM.P1", "n.rep.HM.P2")

df.med.nbrI <- left_join(rep, medMTX, by=c("Duplication", "pair", "P1", "P2"))
df.med.nbrI <-  left_join(df.med.nbrI, interac, by=c("Duplication", "pair", "P1", "P2")) %>%
  mutate(nbHET.inter=HET.P1P2 + HET.P2P1,
         nbHM.inter=HM1.PCA + HM2.PCA) %>% as.data.frame()


##simple motifs
#noI : no interaction
#HM : one of the 2 paralogs forms HM / no HET or both paralogs form HM / no HET
#HET : HET / no HM
#HM&HET : one of the 2 paralogs forms HM / HET or both paralogs form HM / HET


df.motif <- df.med.nbrI %>% group_by(Duplication, pair) %>% 
  mutate(motif.categories.PCA = if (nbHM.inter==0 & nbHET.inter==0) "NI"
            else if ((nbHM.inter==1 & nbHET.inter==0) | (nbHM.inter==2 & nbHET.inter==0)) "HM"
            else if (nbHM.inter==0 & (nbHET.inter==2 | nbHET.inter==1)) "HET"
            else if ((nbHM.inter==1 & (nbHET.inter==2 | nbHET.inter==1)) | (nbHM.inter==2 & (nbHET.inter==2 | nbHET.inter==1))) "HM&HET"
            else "pb") %>% 
  mutate(motif.number.PCA = if (nbHM.inter==0 & nbHET.inter==0) 'NI'
         else if (nbHM.inter==1 & nbHET.inter==0) '1HM'
           else if (nbHM.inter==2 & nbHET.inter==0) '2HM'
             else if (nbHM.inter==0 & nbHET.inter>0) 'HET'
               else if (nbHM.inter==1 & nbHET.inter>0) '1HM.HET'
                 else if (nbHM.inter==2 & nbHET.inter>0) '2HM.HET') %>% as.data.frame()

write.table(df.motif, file="without_low_qual/PCA_motif_interaction_2019_06_phylom.tab", sep="\t",  quote=F, row.names=F)
