rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(reshape2)
library(cowplot)
library(stringr)

setwd("/Users/axellemarchant/Documents/postdoc_Landry/AMarchant_2016-2019/papier_AMarchant_2019")

WGD <- read.table("data/WGD.csv", header = F, sep=";")
WGD$Duplication = "WGD"

SSD <- read.table("data/duplication_SDS_1paire.txt", header = F, sep="\t")
SSD$Duplication = "SSD"

df <- rbind(WGD, SSD)
colnames(df) <- c("P1", "P2", "Duplication")

df$pair <- apply(cbind(as.character(df$P1), 
                          as.character(df$P2)), 
                    1, function(x) paste(sort(x), collapse="."))

#Organization of P1 P2 by alphabet order
df <- select(df, Duplication, pair)
P <- colsplit(df$pair, "[.]", names=c("P1", "P2"))
df <- cbind(df, P)

bg <- read.table("data/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.166.tab2.txt", 
                 header = T,fill=T, quote="", sep="\t",stringsAsFactors=FALSE, comment.char = "")


#consider only physical interactions 
bg.p = bg %>% filter(Experimental.System.Type =="physical")

#look at the techniques to make sure we do not use things that do not reflect ppis
unique(bg.p$Experimental.System)

#remove techniques that are not ppis
bg.p <- filter(bg.p, Experimental.System!="Co-fractionation" & 
                 Experimental.System!="Affinity Capture-RNA" & 
                 Experimental.System!="Protein-RNA" & 
                 Experimental.System!="Co-purification" &
                 Experimental.System!="Proximity Label-MS" & 
                 Experimental.System!="Co-localization")

#reduce the sice of the dataframe so it is easier to handle
bg.p %<>% select(Systematic.Name.Interactor.A,Systematic.Name.Interactor.B,
                     Throughput, Experimental.System,Pubmed.ID)

#create a unique id per ppi
bg.p$ppi <- apply(cbind(as.character(bg.p$Systematic.Name.Interactor.A), 
                            as.character(bg.p$Systematic.Name.Interactor.B)), 
                      1, function(x) paste(sort(x), collapse="."))

bg.p %<>% mutate(direct = ifelse((Experimental.System=="PCA" |
                                     Experimental.System=="Two-hybrid" |
                                     Experimental.System=="Biochemical Activity" |
                                     Experimental.System=="Protein-peptide"), 1, 0),
                   in.vivo = ifelse((Experimental.System=="Affinity Capture-Luminescence" |
                                      Experimental.System=="Affinity Capture-MS" |
                                      Experimental.System=="Affinity Capture-Western" | #not sure
                                      Experimental.System=="PCA" | 
                                      Experimental.System=="Two-hybrid" |
                                      Experimental.System=="FRET"), 1, 0),
                   regul.nat = ifelse((Experimental.System=="PCA" | Experimental.System=="FRET"), 1, 0),
                   HM = ifelse(Systematic.Name.Interactor.A==Systematic.Name.Interactor.B, 1, 0),
                   H.throughput =ifelse((Throughput=="High Throughput" | Throughput=="High Throughput|Low Throughput"), 1, 0))

#get list of studies in which there could be HM
#considering cases where there is at least one HM reported
#code if ppi a HM
#group by study
#count HMs over all PPIs

bg.p$Pubmed.ID <- as.character(bg.p$Pubmed.ID)

bg.p_hm_stats  <- bg.p %>% group_by(Pubmed.ID) %>%
                  summarise(tot_hm = sum(HM), tot_ppi=n()) %>%
                  mutate(percent_hm = tot_hm/tot_ppi)

hist(bg.p_hm_stats$percent_hm)

#get a list of study where there is at least on hm. will be use below.

see_hm <- unique(bg.p_hm_stats$Pubmed.ID[bg.p_hm_stats$percent_hm>0])

#for each study, did HM were be tested
bg.p %<>% group_by(Pubmed.ID) %>% mutate(see_hm_in_study = ifelse((sum(HM) > 0), 1, 0)) %>% as.data.frame()

#for each ppi, test if could be seen as hm or not. 

#score each PPI as NA,0,1
#NA, could not be seen, consider missing data
#0, was not seen but could have (see conditions below)
#1, was seen 

#record data in data_Frame uploaded above. Add three columns to code biogrid data
#put NA as default value.
bg.p %<>% mutate(ExpRef = paste(Experimental.System,"(",Pubmed.ID, ")", sep="")) 

bg.p.HET <- filter(bg.p, Systematic.Name.Interactor.A!=Systematic.Name.Interactor.B)

df.HET <- left_join(df, bg.p.HET, by=c("pair"="ppi"))
df.HET <- filter(df.HET, !is.na(ExpRef))

bg.p.HM <- filter(bg.p, Systematic.Name.Interactor.A==Systematic.Name.Interactor.B)
write.table(bg.p.HM, file="output/bg.p.HM_2019_02.csv", sep="\t",  quote=F, row.names=F)

#Selection of WGD and SSD pairs
df.HM1 <- left_join(df, bg.p.HM, by=c("P1"="Systematic.Name.Interactor.A"))
df.HM1 <- filter(df.HM1, !is.na(ExpRef))
df.HM2 <- left_join(df, bg.p.HM, by=c("P2"="Systematic.Name.Interactor.A"))
df.HM2 <- filter(df.HM2, !is.na(ExpRef))

studies.bait <- unique(select(bg.p, ExpRef, Systematic.Name.Interactor.A))
studies.bait$P1.bait.in.study <- 1
studies.prey <- unique(select(bg.p, ExpRef, Systematic.Name.Interactor.B))
studies.prey$P1.prey.in.study <- 1

df.HET <- left_join(df.HET, studies.bait, by=c("ExpRef" , "P1"="Systematic.Name.Interactor.A"))
df.HET <- left_join(df.HET, studies.prey, by=c("ExpRef" , "P1"="Systematic.Name.Interactor.B"))
df.HM1 <- left_join(df.HM1, studies.bait, by=c("ExpRef" , "P1"="Systematic.Name.Interactor.A"))
df.HM1 <- left_join(df.HM1, studies.prey, by=c("ExpRef" , "P1"="Systematic.Name.Interactor.B"))
df.HM2 <- left_join(df.HM2, studies.bait, by=c("ExpRef" , "P1"="Systematic.Name.Interactor.A"))
df.HM2 <- left_join(df.HM2, studies.prey, by=c("ExpRef" , "P1"="Systematic.Name.Interactor.B"))

studies.bait <- rename(studies.bait, "P2.bait.in.study"="P1.bait.in.study")
studies.prey <- rename(studies.prey, "P2.prey.in.study"="P1.prey.in.study")

df.HET <- left_join(df.HET, studies.bait, by=c("ExpRef" , "P2"="Systematic.Name.Interactor.A"))
df.HET <- left_join(df.HET, studies.prey, by=c("ExpRef" , "P2"="Systematic.Name.Interactor.B"))
df.HM1 <- left_join(df.HM1, studies.bait, by=c("ExpRef" , "P2"="Systematic.Name.Interactor.A"))
df.HM1 <- left_join(df.HM1, studies.prey, by=c("ExpRef" , "P2"="Systematic.Name.Interactor.B"))
df.HM2 <- left_join(df.HM2, studies.bait, by=c("ExpRef" , "P2"="Systematic.Name.Interactor.A"))
df.HM2 <- left_join(df.HM2, studies.prey, by=c("ExpRef" , "P2"="Systematic.Name.Interactor.B"))

df.HET$HET <- 1
df.HM1$HM1 <- 1
df.HM2$HM2 <- 1

df.HET <- select(df.HET, Duplication, pair, P1, P2, ExpRef, direct, in.vivo, regul.nat, H.throughput, 
                 P1.bait.in.study, P1.prey.in.study, P2.bait.in.study, P2.prey.in.study, HET, see_hm_in_study)
df.HM1 <- select(df.HM1, Duplication, pair, P1, P2, ExpRef, direct, in.vivo, regul.nat, H.throughput, 
                 P1.bait.in.study, P1.prey.in.study, P2.bait.in.study, P2.prey.in.study, HM1)
df.HM2 <- select(df.HM2, Duplication, pair, P1, P2, ExpRef, direct, in.vivo, regul.nat, H.throughput, 
                 P1.bait.in.study, P1.prey.in.study, P2.bait.in.study, P2.prey.in.study, HM2)


#Concatenate HM1, HM2 and HET data per pair
df.tot <- full_join(df.HM1, df.HM2, by=c("Duplication", "pair", "P1", "P2", "ExpRef", "direct", "in.vivo", "regul.nat", "H.throughput", 
                                         "P1.bait.in.study", "P1.prey.in.study", "P2.bait.in.study", "P2.prey.in.study"))
df.tot <- full_join(df.tot, df.HET, by=c("Duplication", "pair", "P1", "P2", "ExpRef", "direct", "in.vivo", "regul.nat", "H.throughput", 
                                         "P1.bait.in.study", "P1.prey.in.study", "P2.bait.in.study", "P2.prey.in.study"))

df.tot %<>% mutate(see_hm_in_study= ifelse(is.na(see_hm_in_study) & (HM1==1 | HM2==1),1,
                                    ifelse(is.na(see_hm_in_study) & (HM1==0 & HM2==0), 0, see_hm_in_study)))

df.tot %<>% mutate(HM1.2 = ifelse((is.na(P1.bait.in.study) |  is.na(P1.prey.in.study)), NA,
                                 ifelse(see_hm_in_study==0, NA, 
                                 ifelse((is.na(HM1) & P1.bait.in.study==1  & P1.prey.in.study==1 & see_hm_in_study==1), 0, 
                                 ifelse(HM1==1, 1, 'pb')))),
                   HM2.2 = ifelse((is.na(P2.bait.in.study) |  is.na(P2.prey.in.study)), NA,
                                  ifelse(see_hm_in_study==0, NA, 
                                  ifelse((is.na(HM2) & P2.bait.in.study==1  & P2.prey.in.study==1 & see_hm_in_study==1), 0, 
                                  ifelse(HM2==1, 1, 'pb')))),
                   HET.2 = ifelse((is.na(P1.bait.in.study) & is.na(P1.prey.in.study)) | (is.na(P2.bait.in.study) & is.na(P2.prey.in.study)), NA,
                           ifelse(((is.na(P1.bait.in.study) & is.na(P2.bait.in.study)) | (is.na(P1.prey.in.study) & is.na(P2.prey.in.study))), NA,
                           ifelse(((!is.na(P1.bait.in.study) & !is.na(P2.prey.in.study)) | (!is.na(P2.bait.in.study) & !is.na(P1.prey.in.study))) & is.na(HET), 0,
                           ifelse(HET==1, 1, 'pb')))))


#Summarise all data per pair
df.tot.sum <- df.tot %>% group_by(Duplication, pair, P1, P2) %>% summarise(timesrep.bg.HM1 = sum(HM1.2[!is.na(HM1.2)]),
                                        list.bg.HM1.Pubmed.id = paste(sort(unique(c(as.character(ExpRef[HM1.2==1])))), collapse=";"),
                                        HM1bg = ifelse(sum(HM1.2[!is.na(HM1.2)]) > 0, 1, 0),
                                        timesrep.bg.HM2 = sum(HM2.2[!is.na(HM2.2)]),
                                        HM2bg = ifelse(sum(HM2.2[!is.na(HM2.2)]) > 0, 1, 0),
                                        list.bg.HM2.Pubmed.id = paste(sort(unique(c(as.character(ExpRef[HM2.2==1])))), collapse=";"),
                                        timesrep.bg.HET = sum(HET.2[!is.na(HET.2)]),
                                        HETbg = ifelse(sum(HET.2[!is.na(HET.2)]) > 0, 1, 0),
                                        list.bg.HET.Pubmed.id = paste(sort(unique(c(as.character(ExpRef[HET.2==1])))), collapse=";")) %>% as.data.frame()

df.tot.sum$list.bg.HM1.Pubmed.id[df.tot.sum$list.bg.HM1.Pubmed.id==''] <- 'NA'
df.tot.sum$list.bg.HM2.Pubmed.id[df.tot.sum$list.bg.HM2.Pubmed.id==''] <- 'NA'
df.tot.sum$list.bg.HET.Pubmed.id[df.tot.sum$list.bg.HET.Pubmed.idT==''] <- 'NA'

df.tot.sum %<>% as.data.frame()



#################################################
#     Complete with Stynen et al., 2018 data      #
################################################

Michnick <- read.table("data/MedianValuesHomoMichnick_dupli_statue_2019_02_22_AM.tab", header = T, sep="\t")
Michnick <- select(Michnick, orf, interaction)

prot_remove <- read.table("data/Prot_removed.csv", header=T)
colnames(prot_remove) <- "orf"
prot_remove$remove_MatA <- "remove"
Michnick <- left_join(Michnick, prot_remove, by="orf")
Michnick <- filter(Michnick, is.na(remove_MatA)) %>% select(., orf, interaction)


colnames(Michnick) <- c("P1", "HM1.Stynen")
df.Mich.P1 <- left_join(df, Michnick, by="P1")
colnames(Michnick) <- c("P2", "HM2.Stynen")
df.Mich.P2 <- left_join(df, Michnick, by="P2")
df.Mich <- full_join(df.Mich.P1, df.Mich.P2, by=c("Duplication", "pair", "P1", "P2"))
df.tot.sum.Mich <- full_join(df.tot.sum, df.Mich, by=c("Duplication", "pair", "P1", "P2"))
df.tot.sum.Mich %<>% group_by(pair) %>% mutate(HM1bg.S = ifelse(((is.na(HM1bg) & !is.na(HM1.Stynen)) | (HM1bg==0 & !is.na(HM1.Stynen) & HM1.Stynen==1)), HM1.Stynen, HM1bg))
df.tot.sum.Mich %<>% group_by(pair) %>% mutate(HM2bg.S = ifelse(((is.na(HM2bg) & !is.na(HM2.Stynen)) | (HM2bg==0 & !is.na(HM2.Stynen) & HM2.Stynen==1)), HM2.Stynen, HM2bg))

######################################
#       Complete with PDB data       #
######################################

PDB <- read.table("data/paralogs_PDB_structures.txt", header=T)
PDB %<>% select(.,-P1_unspecific, -P2_unspecific)
colnames(PDB) <- c("P1", "P2", "Duplication", "comp.list.PDB.monomer.P1",                                       
                   "comp.list.PDB.monomer.P2", "comp.list.PDB.HM1", 
                   "comp.list.PDB.HM2", "comp.list.PDB.HET")
PDB$pair <- apply(cbind(as.character(PDB$P1), 
                        as.character(PDB$P2)), 
                  1, function(x) paste(sort(x), collapse="."))

#when we have no observation. 0 to say 'unlikely' NA to say 'no data'
#unlikely to see a HM if we have a monomer
#unlikely to see HET if we have the 2 HM observed or the 2 monomers observed
PDB %<>% mutate(HM1.PDB = ifelse(is.na(comp.list.PDB.HM1) & is.na(comp.list.PDB.monomer.P1), NA,
                          ifelse(is.na(comp.list.PDB.HM1) & !is.na(comp.list.PDB.monomer.P1), 0, 
                          ifelse(!is.na(comp.list.PDB.HM1), 1, 'pb'))),
                HM2.PDB = ifelse(is.na(comp.list.PDB.HM2) & is.na(comp.list.PDB.monomer.P2), NA,
                          ifelse(is.na(comp.list.PDB.HM2) & !is.na(comp.list.PDB.monomer.P2), 0, 
                          ifelse(!is.na(comp.list.PDB.HM2), 1, 'pb'))),
                HET.PDB = ifelse(is.na(comp.list.PDB.HET)  & is.na(comp.list.PDB.HM1) & is.na(comp.list.PDB.monomer.P1), NA,
                          ifelse(is.na(comp.list.PDB.HET)  & is.na(comp.list.PDB.HM2) & is.na(comp.list.PDB.monomer.P2), NA,
                          ifelse(is.na(comp.list.PDB.HET) & (!is.na(comp.list.PDB.HM1) | !is.na(comp.list.PDB.monomer.P1)) & 
                                (!is.na(comp.list.PDB.HM2) | !is.na(comp.list.PDB.monomer.P2)), 0,
                          ifelse(!is.na(comp.list.PDB.HET), 1, 'pb')))))


df.tot.sum.Mich.PDB <- full_join(df.tot.sum.Mich, PDB, by=c("Duplication", "pair", "P1", "P2"))
df.tot.sum.Mich.PDB$HM1.PDB <- as.numeric(df.tot.sum.Mich.PDB$HM1.PDB)
df.tot.sum.Mich.PDB$HM2.PDB <- as.numeric(df.tot.sum.Mich.PDB$HM2.PDB)
df.tot.sum.Mich.PDB$HET.PDB <- as.numeric(df.tot.sum.Mich.PDB$HET.PDB)
df.tot.sum.Mich.PDB$HM1bg.S <- as.numeric(df.tot.sum.Mich.PDB$HM1bg.S)
df.tot.sum.Mich.PDB$HM2bg.S <- as.numeric(df.tot.sum.Mich.PDB$HM2bg.S)
df.tot.sum.Mich.PDB$HETbg <- as.numeric(df.tot.sum.Mich.PDB$HETbg)
df.tot.sum.Mich.PDB$HET.PDB <- as.numeric(df.tot.sum.Mich.PDB$HET.PDB)

df.tot.sum.Mich.PDB %<>% group_by(pair) %>% mutate(HM1bg.S.PDB = ifelse(((is.na(HM1bg.S) & !is.na(HM1.PDB)) | (HM1bg.S==0 & !is.na(HM1.PDB) & HM1.PDB==1)), HM1.PDB, HM1bg.S))
df.tot.sum.Mich.PDB %<>% group_by(pair) %>% mutate(HM2bg.S.PDB = ifelse(((is.na(HM2bg.S) & !is.na(HM2.PDB)) | (HM2bg.S==0 & !is.na(HM2.PDB) & HM2.PDB==1)), HM2.PDB, HM2bg.S))
df.tot.sum.Mich.PDB %<>% group_by(pair) %>% mutate(HETbg.PDB = ifelse(((is.na(HETbg) & !is.na(HET.PDB)) | (HETbg==0 & !is.na(HET.PDB) & HET.PDB==1)), HET.PDB, HETbg))
df.tot.sum.Mich.PDB %<>% as.data.frame()

###############################################
#      Complete with Kim et al., 2019  data     #
###############################################
Kim <- read.table("data/orf_tested_by_Kim_et_al_2019.txt", header = T, sep="\t")
Kim.falsePos <- read.table("data/fasle_positiv_Kim_et_al_2019.txt", header = T, sep="\t")
Kim.falsePos <- select(Kim.falsePos, ORF)
Kim.falsePos$false <- 'yes'
#remove potential false +
Kim <- left_join(Kim, Kim.falsePos, by="ORF")
Kim <- filter(Kim, is.na(false))

HM.Kim <- read.table("data/HM_Kim_et_al_2019.txt", header = T, sep="\t")
HM.Kim <- select(HM.Kim, ORF)
HM.Kim$HM.Kim <- 1
Kim <- left_join(Kim, HM.Kim, by="ORF")
Kim$HM.Kim[is.na(Kim$HM.Kim)] <- 0
Kim <- select(Kim, ORF, HM.Kim)
colnames(Kim) <- c("P1", "HM1.Kim")

df.Kim <- left_join(df, Kim, by="P1")
colnames(Kim) <- c("P2", "HM2.Kim")
df.Kim <- left_join(df.Kim, Kim, by="P2")

df.tot.sum.Mich.PDB.Kim <- full_join(df.tot.sum.Mich.PDB, df.Kim, by=c("Duplication", "pair", "P1", "P2"))

##################################################################################
#     Final  result of previously reported HM and HET for duplicated proteins   #
#################################################################################
df.tot.sum.Mich.PDB.Kim %<>% group_by(pair) %>% mutate(HM1bg.S.PDB.Kim = ifelse(((is.na(HM1bg.S.PDB) & !is.na(HM1.Kim)) | (HM1bg.S.PDB==0 & !is.na(HM1.Kim) & HM1.Kim==1)), HM1.Kim, HM1bg.S.PDB))
df.tot.sum.Mich.PDB.Kim %<>% group_by(pair) %>% mutate(HM2bg.S.PDB.Kim = ifelse(((is.na(HM2bg.S.PDB) & !is.na(HM2.Kim)) | (HM2bg.S.PDB==0 & !is.na(HM2.Kim) & HM2.Kim==1)), HM2.Kim, HM2bg.S.PDB))
df.tot.sum.Mich.PDB.Kim %<>% as.data.frame()


#####################################################
#     Select complete pairs for both HM and HET     #
#####################################################
complete_pairs <- filter(df.tot.sum.Mich.PDB.Kim, !is.na(HM1bg.S.PDB.Kim) & !is.na(HM2bg.S.PDB.Kim) & !is.na(HETbg.PDB)) %>% as.data.frame() #498

#Motifs of interactions observed based on previously reported data
complete_pairs %<>% group_by(pair) %>% 
  mutate(motif.categories.bg.S.PDB.Kim=ifelse((HM1bg.S.PDB.Kim==1 | HM2bg.S.PDB.Kim==1) & HETbg.PDB==1, "HM&HET",
                          ifelse((HM1bg.S.PDB.Kim==1 | HM2bg.S.PDB.Kim==1) & HETbg.PDB==0, "HM",
                          ifelse(HM1bg.S.PDB.Kim==0 & HM2bg.S.PDB.Kim==0 & HETbg.PDB==1, "HET",
                          ifelse(HM1bg.S.PDB.Kim==0 & HM2bg.S.PDB.Kim==0 & HETbg.PDB==0, "NI", "PB"))))) %>% as.data.frame()


complete_pairs %<>% group_by(pair) %>% 
  mutate(motif.number.bg.S.PDB.Kim=ifelse(HM1bg.S.PDB.Kim==0 & HM2bg.S.PDB.Kim==0 & HETbg.PDB==0, 1,
                           ifelse(((HM1bg.S.PDB.Kim==1 & HM2bg.S.PDB.Kim==0) | (HM1bg.S.PDB.Kim==0 & HM2bg.S.PDB.Kim==1)) & HETbg.PDB==0, 2,
                           ifelse(HM1bg.S.PDB.Kim==1 & HM2bg.S.PDB.Kim==1 & HETbg.PDB==0, 3,
                           ifelse(HM1bg.S.PDB.Kim==0 & HM2bg.S.PDB.Kim==0 & HETbg.PDB==1, 4,
                           ifelse(((HM1bg.S.PDB.Kim==1 & HM2bg.S.PDB.Kim==0) | (HM1bg.S.PDB.Kim==0 & HM2bg.S.PDB.Kim==1)) & HETbg.PDB==1, 5,
                           ifelse(HM1bg.S.PDB.Kim==1 & HM2bg.S.PDB.Kim==1 & HETbg.PDB==1, 6, "pb")))))))


##########################################
#     Complete with data of our PCA      #
##########################################

our_PCA <- read.table("output/PCA_motif_interaction_2019_03.tab", sep="\t", header = T)

our_PCA <- full_join(our_PCA, complete_pairs, by=c("Duplication", "pair", "P1", "P2")) #622

#remove proteins interacting with everythings (defined as false positive in Tarassov et al., 2008) 
#and not repertoried in biogrid
prot_remove <- read.table("data/Prot_removed.csv", header=T)
colnames(prot_remove) <- "P1"
prot_remove$P1.Tarassov.removed <- "remove"
our_PCA <- left_join(our_PCA, prot_remove, by="P1")
colnames(prot_remove) <- c("P2", "P2.Tarassov.removed")
our_PCA <- left_join(our_PCA, prot_remove, by="P2") #622 ok
nrow(filter(our_PCA, P1.Tarassov.removed == 'remove' | P2.Tarassov.removed == 'remove')) #123

our_PCA %<>% mutate(P1.to.remove = ifelse((P1.Tarassov.removed=='remove' & (is.na(HM1bg.S.PDB.Kim) | is.na(HETbg.PDB))), 1, 0),
                    P2.to.remove = ifelse((P2.Tarassov.removed=='remove' & (is.na(HM2bg.S.PDB.Kim) | is.na(HETbg.PDB))), 1, 0))
nrow(filter(our_PCA, P1.to.remove==1 |  P2.to.remove==1)) #27

our_PCA %<>% filter(., (is.na(P1.to.remove) | P1.to.remove==0) &  (is.na(P2.to.remove) | P2.to.remove==0)) #595

#infere complete motif with both our PCA and previously reported interactions
nrow(filter(our_PCA, is.na(motif.categories.PCA))) #309
nrow(filter(our_PCA, is.na(motif.categories.bg.S.PDB.Kim))) #97
nrow(filter(our_PCA, !is.na(motif.categories.PCA) & !is.na(motif.categories.bg.S.PDB.Kim))) #189
filter(our_PCA, !is.na(motif.categories.PCA) & !is.na(motif.categories.bg.S.PDB.Kim)) %>% 
  filter(., motif.categories.PCA==motif.categories.bg.S.PDB.Kim) %>% nrow(.) #85


our_PCA$HM1bg.S.PDB.Kim <- as.numeric(our_PCA$HM1bg.S.PDB.Kim)
our_PCA %<>% group_by(pair) %>% mutate(HM1.final = ifelse((is.na(P1.Tarassov.removed) & !is.na(HM1.PCA) & !is.na(HM1bg.S.PDB.Kim) & (HM1.PCA==1 | HM1bg.S.PDB.Kim==1)), 1, 
                                                   ifelse((is.na(P1.Tarassov.removed) & !is.na(HM1.PCA) & !is.na(HM1bg.S.PDB.Kim) & HM1.PCA==0 & HM1bg.S.PDB.Kim==0), 0,
                                                   ifelse(is.na(P1.Tarassov.removed) & is.na(HM1.PCA), HM1bg.S.PDB.Kim,
                                                   ifelse(is.na(P1.Tarassov.removed) & is.na(HM1bg.S.PDB), HM1.PCA, 
                                                   ifelse(P1.Tarassov.removed== "remove" & is.na(HM1bg) & !is.na(HM1.Kim), HM1.Kim,
                                                   ifelse(P1.Tarassov.removed== "remove" & !is.na(HM1bg) & is.na(HM1.Kim), HM1bg,
                                                   ifelse(P1.Tarassov.removed== "remove" & !is.na(HM1bg) & !is.na(HM1.Kim) & (HM1bg+HM1.Kim) > 0, 1,
                                                   ifelse(P1.Tarassov.removed== "remove" & !is.na(HM1bg) & !is.na(HM1.Kim) & (HM1bg+HM1.Kim) == 0, 0, "pb"))))))))) %>% as.data.frame()

our_PCA$HM2bg.S.PDB.Kim <- as.numeric(our_PCA$HM2bg.S.PDB.Kim)
our_PCA %<>% group_by(pair) %>% mutate(HM2.final = ifelse((is.na(P2.Tarassov.removed) & !is.na(HM2.PCA) & !is.na(HM2bg.S.PDB.Kim) & (HM2.PCA==1 | HM2bg.S.PDB.Kim==1)), 1, 
                                                   ifelse((is.na(P2.Tarassov.removed) &!is.na(HM2.PCA) & !is.na(HM2bg.S.PDB.Kim) & HM2.PCA==0 & HM2bg.S.PDB.Kim==0), 0,
                                                   ifelse(is.na(P2.Tarassov.removed) & is.na(HM2.PCA), HM2bg.S.PDB.Kim,
                                                   ifelse(is.na(P2.Tarassov.removed) & is.na(HM2bg.S.PDB.Kim), HM2.PCA, 
                                                   ifelse(P2.Tarassov.removed== "remove" & is.na(HM2bg) & !is.na(HM2.Kim), HM2.Kim,
                                                   ifelse(P2.Tarassov.removed== "remove" & !is.na(HM2bg) & is.na(HM2.Kim), HM2bg,
                                                   ifelse(P2.Tarassov.removed== "remove" & !is.na(HM2bg) & !is.na(HM2.Kim) & (HM2bg+HM2.Kim) > 0, 1,
                                                   ifelse(P2.Tarassov.removed== "remove" & !is.na(HM2bg) & !is.na(HM2.Kim) & (HM2bg+HM2.Kim) == 0, 0, "pb"))))))))) %>% as.data.frame()


our_PCA %<>% group_by(pair) %>% mutate(HET.ourPCA = ifelse((!is.na(HET.P1P2) & !is.na(HET.P2P1) & (HET.P1P2==1 | HET.P2P1==1)), 1, 
                                                   ifelse((!is.na(HET.P1P2) & !is.na(HET.P2P1) & HET.P1P2==0 & HET.P2P1==0), 0,
                                                   ifelse(is.na(HET.P1P2) & is.na(HET.P2P1), NA, "pb")))) %>% as.data.frame()

our_PCA$HETbg.PDB <- as.numeric(our_PCA$HETbg.PDB)
our_PCA %<>% group_by(pair) %>% mutate(HET.final = ifelse((is.na(P1.Tarassov.removed) & is.na(P2.Tarassov.removed) & !is.na(HET.ourPCA) & !is.na(HETbg.PDB) & (HET.ourPCA==1 | HETbg.PDB==1)), 1, 
                                                   ifelse((is.na(P1.Tarassov.removed) & is.na(P2.Tarassov.removed) & !is.na(HET.ourPCA) & !is.na(HETbg.PDB) & HET.ourPCA==0 & HETbg.PDB==0), 0,
                                                   ifelse(is.na(P1.Tarassov.removed) & is.na(P2.Tarassov.removed) & is.na(HET.ourPCA), HETbg.PDB,
                                                   ifelse(is.na(P1.Tarassov.removed) & is.na(P2.Tarassov.removed) & is.na(HETbg.PDB), HET.ourPCA, 
                                                   ifelse(P1.Tarassov.removed== "remove" | P2.Tarassov.removed== "remove",  HETbg.PDB, "pb")))))) %>% as.data.frame()

head(our_PCA)


#Motif interaction with previously reported and our PCA data
our_PCA %<>% group_by(pair) %>% 
  mutate(motif.categories=ifelse((HM1.final==1 | HM2.final==1) & HET.final==1, "HM&HET",
                               ifelse((HM1.final==1 | HM2.final==1) & HET.final==0, "HM",
                               ifelse(HM1.final==0 & HM2.final==0 & HET.final==1, "HET",
                               ifelse(HM1.final==0 & HM2.final==0 & HET.final==0, "NI", "PB"))))) %>% as.data.frame()


our_PCA %<>% group_by(pair) %>% 
  mutate(motif.number=ifelse(HM1.final==0 & HM2.final==0 & HET.final==0, 1,
                             ifelse(((HM1.final==1 & HM2.final==0) | (HM1.final==0 & HM2.final==1)) & HET.final==0, 2,
                             ifelse(HM1.final==1 & HM2.final==1 & HET.final==0, 3,
                             ifelse(HM1.final==0 & HM2.final==0 & HET.final==1, 4,
                             ifelse(((HM1.final==1 & HM2.final==0) | (HM1.final==0 & HM2.final==1)) & HET.final==1, 5,
                             ifelse(HM1.final==1 & HM2.final==1 & HET.final==1, 6, "pb"))))))) %>% as.data.frame()
nrow(our_PCA) #595

write.table(our_PCA, file="output/PCA_completed_by_Biog_Sty_2019_03.csv", sep="\t",  quote=F, row.names=F)


#comparisons between previously reported interactions and our PCA


#HM and HET never found before:
#HM
nrow(filter(our_PCA, HM1.PCA==1 & (is.na(HM1bg.S.PDB.Kim) | HM1bg.S.PDB.Kim==0))) #33
nrow(filter(our_PCA, HM2.PCA==1 & (is.na(HM2bg.S.PDB.Kim) | HM2bg.S.PDB.Kim==0))) #21
#tot HM new = 55
nrow(filter(our_PCA, (HET.P1P2==1 | HET.P2P1==1) & (is.na(HETbg.PDB) | HETbg.PDB==0))) #22
#tot new interaction = 77

#HM and HET found before and tested but not detected here:
nrow(filter(our_PCA, HM1.PCA==0 & HM1bg.S.PDB.Kim==1)) #42
nrow(filter(our_PCA, HM2.PCA==0 & HM2bg.S.PDB.Kim==1)) #48
#tot HM not detected = 90
nrow(filter(our_PCA, HET.P1P2==0 & HET.P2P1== 0 & HETbg.PDB==1)) #49
#tot 136

