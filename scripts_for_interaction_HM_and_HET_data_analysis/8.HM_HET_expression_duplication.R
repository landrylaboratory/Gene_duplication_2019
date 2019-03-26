rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

###########################################################
#        our PCA completed by Biogrid and Stynley         #
###########################################################

summary_table <- read.table("output/PCA_completed_by_Biog_Sty_2019_03.csv", sep="\t", header = T)

#nrow = 595 OK

#####################################
#    Compare with other databases  #
####################################
#reported in other studies and not detected in our PCA
#HMs
nrow(filter(summary_table, HM1.PCA==0 & HM1bg.S.PDB.Kim==1)) +
nrow(filter(summary_table, HM2.PCA==0 & HM2bg.S.PDB.Kim==1))
#HETs
nrow(filter(summary_table, HET.ourPCA==0 & (HETbg==1 | HET.PDB==1)))
#reported in other studies and not tested in our PCA
#HETs
nrow(filter(summary_table, is.na(HET.ourPCA) & (HETbg==1 | HET.PDB==1)))
#HMs
nrow(filter(summary_table, is.na(HM2.PCA) & HM2bg.S.PDB.Kim==1)) +
nrow(filter(summary_table, is.na(HM1.PCA) & HM1bg.S.PDB.Kim==1))
#reported in other studies and detected in our PCA
#HMs
nrow(filter(summary_table, HM1.PCA==1 & (HM1bg.S.PDB.Kim==1 | is.na(HM1bg.S.PDB.Kim)))) +
nrow(filter(summary_table, HM2.PCA==1 & (HM2bg.S.PDB.Kim==1 | is.na(HM2bg.S.PDB.Kim))))
#HETs
nrow(filter(summary_table, HET.ourPCA==1 & (HETbg==1 | HET.PDB==1)))
#never reported in other studies and detected in our PCA
nrow(filter(summary_table, HM1.PCA==1 & (HM1bg.S.PDB.Kim==0 | is.na(HM1bg.S.PDB.Kim)))) +
nrow(filter(summary_table, HM2.PCA==1 & (HM2bg.S.PDB.Kim==0 | is.na(HM2bg.S.PDB.Kim))))
nrow(filter(summary_table, HET.ourPCA==1 & (HETbg==0 | is.na(HETbg)) & (HET.PDB==0 | is.na(HET.PDB))))



#SFig1 --> comparison with biogrid
HM1.PCA <- select(summary_table, med.MTX.HM.P1, timesrep.bg.HM1)
colnames(HM1.PCA) <- c("median.size.PCA", "N.report.in.biogrid")
HM2.PCA <- select(summary_table, med.MTX.HM.P2, timesrep.bg.HM2)
colnames(HM2.PCA) <- c("median.size.PCA", "N.report.in.biogrid")
HET_P1P2 <- select(summary_table, med.MTX.HET.P1P2, timesrep.bg.HET)
colnames(HET_P1P2) <- c("median.size.PCA", "N.report.in.biogrid")
Interact <- rbind(HM1.PCA, HM2.PCA, HET_P1P2)
Interact$N.report.in.biogrid <- as.numeric(as.character(Interact$N.report.in.biogrid))
Interact$N.report.in.biogrid[is.na(Interact$N.report.in.biogrid)] <- 0
Interact$N.report.in.biogrid[Interact$N.report.in.biogrid>=3] <- "3+"
Interact$N.report.in.biogrid <- as.character(Interact$N.report.in.biogrid)

write.table(Interact, file="output/comp_Biogrid.csv", sep="\t",  quote=F, row.names=F)

#comp Michnick
Michnick <- read.table("data/MedianValuesHomoMichnick_dupli_statue_2019_02_22_AM.tab", sep="\t", header=T)
Michnick <- select(Michnick,orf, med)
colnames(Michnick) <- c("P1", "med.MTX.HM.P1.Stynen")
summary_table <- left_join(summary_table, Michnick, by="P1")
colnames(Michnick) <- c("P2", "med.MTX.HM.P2.Stynen")
summary_table <- left_join(summary_table, Michnick, by="P2")

HMP1 <- filter(summary_table, !is.na(med.MTX.HM.P1.Stynen)) %>% 
  select(., med.MTX.HM.P1, med.MTX.HM.P1.Stynen) %>% as.data.frame()
colnames(HMP1) <- c("med_MTX_our_PCA", "med_MTX_Stynen")
HMP2 <- filter(summary_table, !is.na(med.MTX.HM.P2.Stynen)) %>% 
  select(., med.MTX.HM.P2, med.MTX.HM.P2.Stynen) %>% as.data.frame()
colnames(HMP2) <- c("med_MTX_our_PCA", "med_MTX_Stynen")
HM_Stynen_our_PCA <- rbind(HMP1, HMP2)

cor.test(HM_Stynen_our_PCA$med_MTX_our_PCA, HM_Stynen_our_PCA$med_MTX_Stynen, method = "spearman") 
#R=0.5784954   p-value < 2.2e-16

write.table(HM_Stynen_our_PCA, file="output/comp_Stynen.csv", sep="\t",  quote=F, row.names=F)


#comp Tarassov
Tarassov <- read.table("data/Tarassov_1153878s4.csv", sep="\t", header=T)
Tarassov %<>% group_by(MATa.orf_name, MATalpha.orf_name) %>% 
  mutate(interaction = if (intensity >= 23000 & z.score >= 2.4 & PPV. >= 97.7) 1
         else 0) %>% as.data.frame()
Tarassov <-  Tarassov %>% select(MATa.orf_name, MATalpha.orf_name, intensity)
colnames(Tarassov) <- c("P1", "P2", "intensity_Tarassov")

Tarassov$P1 <- as.character(Tarassov$P1)
Tarassov$P2 <- as.character(Tarassov$P2)
Tarassov_HM <- filter(Tarassov, P1==P2) %>% 
select(., P1, intensity_Tarassov)
colnames(Tarassov_HM) <- c("P1", "intensity.HMP1.Tarassov")
summary_table <- left_join(summary_table, Tarassov_HM, by="P1")
colnames(Tarassov_HM) <- c("P2", "intensity.HMP2.Tarassov")
summary_table <- left_join(summary_table, Tarassov_HM, by="P2")

Tarassov$pair <- apply(cbind(as.character(Tarassov$P1), 
                        as.character(Tarassov$P2)), 
                  1, function(x) paste(sort(x), collapse="."))
Tarassov <- Tarassov %>% group_by(P1, P2, intensity_Tarassov, pair) %>%
              mutate(sens = ifelse(pair==(paste(c(P1,P2), collapse=".")), "P1P2", 
                            ifelse(pair==(paste(c(P2,P1), collapse=".")), "P2P1", "pb"))) %>% ungroup()
TarassovP1P2 <- filter(Tarassov, sens=="P1P2") 
TarassovP1P2 <- select(TarassovP1P2,-P1, -P2, -sens)
colnames(TarassovP1P2) <- c("intensity.HETP1P2.Tarassov", "pair")
TarassovP2P1 <- filter(Tarassov, sens=="P2P1") %>% select(.,-P1, -P2, -sens)
colnames(TarassovP2P1) <- c("intensity.HETP2P1.Tarassov", "pair")

summary_table <- left_join(summary_table, TarassovP1P2, by="pair")
summary_table <- left_join(summary_table, TarassovP2P1, by="pair")

HMP1 <- filter(summary_table, !is.na(intensity.HMP1.Tarassov)) %>% 
  select(., med.MTX.HM.P1, intensity.HMP1.Tarassov) %>% as.data.frame()
colnames(HMP1) <- c("med_MTX_our_PCA", "intensity_Tarassov")
HMP2 <- filter(summary_table, !is.na(intensity.HMP2.Tarassov)) %>% 
  select(., med.MTX.HM.P2, intensity.HMP2.Tarassov) %>% as.data.frame()
colnames(HMP2) <- c("med_MTX_our_PCA", "intensity_Tarassov")
HETP1P2 <- filter(summary_table, !is.na(intensity.HETP1P2.Tarassov)) %>% 
  select(.,med.MTX.HET.P1P2, intensity.HETP1P2.Tarassov) %>% as.data.frame()
colnames(HETP1P2) <- c("med_MTX_our_PCA", "intensity_Tarassov")
HETP2P1 <- filter(summary_table, !is.na(intensity.HETP2P1.Tarassov)) %>% 
  select(.,med.MTX.HET.P2P1, intensity.HETP2P1.Tarassov) %>% as.data.frame()
colnames(HETP2P1) <- c("med_MTX_our_PCA", "intensity_Tarassov")
HM_Tarassov_our_PCA <- rbind(HMP1, HMP2, HETP1P2, HETP2P1)

cor.test(HM_Tarassov_our_PCA$med_MTX_our_PCA, HM_Tarassov_our_PCA$intensity_Tarassov, method = "spearman") 
#R=0.5551074    p-value = 2.911e-08

write.table(HM_Tarassov_our_PCA, file="output/comp_Tarassov.csv", sep="\t",  quote=F, row.names=F)

### --> FigS1

###############################################################
### $ Global transcriptomic / proteomic expression analysis ###
###############################################################
RNAseq <- read.table("data/count_table.csv", header=T, sep="\t")
RNAseq %<>% group_by(gene) %>% mutate(log10_med_DMSO = log10(median(DMSO1, DMSO2, DMSO2)+1),
                                      log10_med_MTX = log10(median(MTX1, MTX2, MTX3)+1)) %>% as.data.frame()
RNAseq <- select(RNAseq, gene, log10_med_DMSO, log10_med_MTX)

colnames(RNAseq) <- c("P1", "log10.RNAseq.DMSO.P1", "log10.RNAseq.MTX.P1")
summary_table <- left_join(summary_table, RNAseq, by = "P1")
colnames(RNAseq) <- c("P2", "log10.RNAseq.DMSO.P2", "log10.RNAseq.MTX.P2")
summary_table <- left_join(summary_table, RNAseq, by = "P2")

paxdb <- read.table("data/PaxDB_4932-WHOLE_ORGANISM-integrated_AM.txt") %>% select(V2, V3)
colnames(paxdb) <- c("P1", "log10.abundance.paxDB.P1")
paxdb$log10.abundance.paxDB.P1 = log10(paxdb$log10.abundance.paxDB.P1)
paxdb <- select(paxdb, P1, log10.abundance.paxDB.P1)
summary_table <- left_join(summary_table, paxdb, by = "P1")
colnames(paxdb) <- c("P2", "log10.abundance.paxDB.P2")
summary_table <- left_join(summary_table, paxdb, by = "P2")

#comparison RNAseq with paxDB

P1 <- select(summary_table, log10.abundance.paxDB.P1, log10.RNAseq.DMSO.P1, log10.RNAseq.MTX.P1)
colnames(P1) <- c("log10_abundance_paxDB", "log10_RNAseq_DMSO", "log10_RNAseq_MTX")
P2 <- select(summary_table, log10.abundance.paxDB.P2, log10.RNAseq.DMSO.P2, log10.RNAseq.MTX.P2)
colnames(P2) <- c("log10_abundance_paxDB", "log10_RNAseq_DMSO", "log10_RNAseq_MTX")
P1P2 <- rbind(P1, P2)
cor.test(P1P2$log10_abundance_paxDB, P1P2$log10_RNAseq_DMSO, method = "spearman") 
#p-value < 2.2e-16 R =0.6658133 
cor.test(P1P2$log10_abundance_paxDB, P1P2$log10_RNAseq_MTX, method = "spearman") 
#p-value < 2.2e-16 R =0.5771592 
colnames(P1P2) <- c("log10_abundance_paxDB", "DMSO", "MTX")
P1P2 <- gather(P1P2, media, log10_RNAseq, DMSO:MTX)


#expression versus HM PCA signal
comp_E_PCA_P1 <- select(summary_table, P1, med.MTX.HM.P1, log10.RNAseq.MTX.P1)
colnames(comp_E_PCA_P1) <- c("paralog", "med_PCA_size_MTX", "log10_RNAseq_MTX")
comp_E_PCA_P2 <- select(summary_table, P2, med.MTX.HM.P2, log10.RNAseq.MTX.P2)
colnames(comp_E_PCA_P2) <- c("paralog", "med_PCA_size_MTX", "log10_RNAseq_MTX")
comp_E_PCA <- rbind(comp_E_PCA_P1, comp_E_PCA_P2)
cor.test(comp_E_PCA$med_PCA_size_MTX, comp_E_PCA$log10_RNAseq_MTX, method = "spearman")
#R=0.2905069 , p-value = 1.379e-12


#HET PCA signal versus the abundance of the less expressed paralog
HM1 <- select(summary_table, pair, HM1.PCA, HM1bg.S.PDB.Kim, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
HM2 <- select(summary_table, pair, HM2.PCA, HM2bg.S.PDB.Kim, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM1) <- c("pair", "PCA", "known", "log10.RNAseq.MTX", "log10.abundance.paxDB")
colnames(HM2) <- c("pair", "PCA", "known", "log10.RNAseq.MTX", "log10.abundance.paxDB")
HM <- rbind(HM1, HM2)

HET <- select(summary_table, pair, HET.ourPCA, HETbg.PDB, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
minHET <- HET %>% group_by(pair, HET.ourPCA, HETbg.PDB) %>% summarise(log10.RNAseq.MTX = min(log10.RNAseq.MTX.P1, log10.RNAseq.MTX.P2),
                                                                   log10.abundance.paxDB =  min(log10.abundance.paxDB.P1, log10.abundance.paxDB.P2)) %>% as.data.frame()

colnames(minHET) <- c("pair", "PCA", "known", "log10.RNAseq.MTX", "log10.abundance.paxDB")
HM_minHET <- rbind(HM, minHET)

#--> As previously observed, we found a correlation between PCA signal and expression level, 
#both at the level of mRNA and protein abundance:
cor.test(HM_minHET$PCA, HM_minHET$log10.RNAseq.MTX, method = "spearman")
#r = 0.3159929  p-value < 2.2e-16
cor.test(HM_minHET$PCA, HM_minHET$log10.abundance.paxDB, method = "spearman")
#r = 0.4531822 p-value < 2.2e-16


#Focusing only on HMs:
HM_minHET_known <- filter(HM_minHET, known==1)
cor.test(HM_minHET_known$PCA, HM_minHET_known$log10.RNAseq.MTX, method = "spearman")
#r = 0.2709791 p-value = 2.248e-07
cor.test(HM_minHET_known$PCA, HM_minHET_known$log10.abundance.paxDB, method = "spearman")
#r=0.2942679 p-value < 1.676e-08

### $ Transcriptomic expression in MTX and homomers detection
#expectation of HM detection depending on expression
HM1.PCA <- select(summary_table, P1, HM1.PCA, log10.RNAseq.MTX.P1, HM1bg.S.PDB.Kim)
HM2.PCA <- select(summary_table, P2, HM2.PCA, log10.RNAseq.MTX.P2, HM2bg.S.PDB.Kim)
colnames(HM1.PCA) <- c("paralog", "HM", "log10_RNAseq_MTX", "HMbg.S.PDB.Kim")
colnames(HM2.PCA) <- c("paralog", "HM", "log10_RNAseq_MTX", "HMbg.S.PDB.Kim")
HM <- rbind(HM1.PCA, HM2.PCA)

HM <- HM %>% .[complete.cases(.), ]

HM <- filter(HM, HMbg.S.PDB.Kim==1)

lines <- (ksmooth(HM$log10_RNAseq_MTX, HM$HM, "normal", 
                  bandwidth = 0.8, x.points = seq(0, 5, by=0.001))) %>% as.data.frame()
#FigS3
HM1_rep_undec <- filter(summary_table, HM1.PCA==0 & HM1bg.S.PDB.Kim==1)
HM1_rep_undec <- select(HM1_rep_undec, P1, HM1.PCA, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
colnames(HM1_rep_undec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")
HM2_rep_undec <- filter(summary_table, HM2.PCA==0 & HM2bg.S.PDB.Kim==1)
HM2_rep_undec <- select(HM2_rep_undec, P2, HM2.PCA, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM2_rep_undec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")

HM_rep_undec <- rbind(HM1_rep_undec, HM2_rep_undec)

HM1_rep_dec <- filter(summary_table, HM1.PCA==1 & HM1bg.S.PDB.Kim==1)
HM1_rep_dec <- select(HM1_rep_dec, P1, HM1.PCA, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
colnames(HM1_rep_dec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")
HM2_rep_dec <- filter(summary_table, HM2.PCA==1 & HM2bg.S.PDB.Kim==1)
HM2_rep_dec <- select(HM2_rep_dec, P2, HM2.PCA, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM2_rep_dec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")

HM_rep_dec <- rbind(HM1_rep_dec, HM2_rep_dec)
wilcox.test(HM_rep_undec$log10_RNAseq_MTX, HM_rep_dec$log10_RNAseq_MTX)
#p-value = 1.42e-05
wilcox.test(HM_rep_undec$log10_abundance_paxDB, HM_rep_dec$log10_abundance_paxDB)
#p-value = 1.651e-05


P_HM_sachant_Exp <- ksmooth(HM$log10_RNAseq_MTX, HM$HM, "normal", bandwidth = 0.8, x.points = seq(0, 5, by=0.001))
P_HM_sachant_Exp <- P_HM_sachant_Exp %>% as.data.frame()
colnames(P_HM_sachant_Exp) <- c("log10_RNAseq_MTX", "Proba_HM")

quantile(HM$log10_RNAseq_MTX, c(0.5, 0.9), na.rm=T)
#2.541579 3.674452

filter(P_HM_sachant_Exp, log10_RNAseq_MTX==2.541)
#  log10_RNAseq_MTX  Proba_HM
#1            2.541 0.5844266
filter(P_HM_sachant_Exp, log10_RNAseq_MTX==3.674)
#log10_RNAseq_MTX Proba_HM
#1            3.674 0.7974532
wilcox.test(HM$log10_RNAseq_MTX~HM$HM)
#p-value = 1.42e-05

expP1P2 <- dplyr::select(summary_table, pair, log10.RNAseq.MTX.P1, log10.RNAseq.MTX.P2)

expP1P2 <- expP1P2 %>% group_by(pair) %>% mutate(log10.RNAseq.MTX.P1_round = round(log10.RNAseq.MTX.P1, digits = 3),
                                                 log10.RNAseq.MTX.P2_round = round(log10.RNAseq.MTX.P2, digits = 3))
P_HM_sachant_Exp <- P_HM_sachant_Exp %>% as.data.frame()
colnames(P_HM_sachant_Exp) <- c("log10.RNAseq.MTX.P1_round", "P_HM_Exp_P1")
expP1P2 <- left_join(expP1P2, P_HM_sachant_Exp, by="log10.RNAseq.MTX.P1_round")
colnames(P_HM_sachant_Exp) <- c("log10.RNAseq.MTX.P2_round", "P_HM_Exp_P2")
expP1P2 <- left_join(expP1P2, P_HM_sachant_Exp, by="log10.RNAseq.MTX.P2_round")
expP1P2 <- dplyr::select(expP1P2, pair, P_HM_Exp_P1, P_HM_Exp_P2)

summary_table <- left_join(summary_table, expP1P2, by="pair")


#motif3_or_5 <- filter(summary_table, motif.number.PCA==3 | motif.number.PCA ==5)
#nrow(motif3_or_5) #35
#nrow(filter(motif3_or_5, (HM1.PCA==0 & log10.RNAseq.MTX.P1 < log10.RNAseq.MTX.P2) | 
#              (HM2.PCA==0 & log10.RNAseq.MTX.P2 < log10.RNAseq.MTX.P1))) #18
#nrow(filter(motif3_or_5, (HM1.PCA==0 & log10.abundance.paxDB.P1 < log10.abundance.paxDB.P2) | 
#              (HM2.PCA==0 & log10.abundance.paxDB.P2 < log10.abundance.paxDB.P1))) #20


#"However, no significant difference was observed between P(HM|Exp) of the paralogs 
#with HM and P(HM|Exp) of the paralogs without HM"
#motif3_or_5HM <- motif3_or_5 %>% group_by(pair) %>% 
#  summarise(PHM=ifelse(HM1.PCA==1, "P1", "P2"),
#            P_HM_Exp=ifelse(HM1.PCA==1, P_HM_Exp_P1 , P_HM_Exp_P2))
#motif3_or_5noHM <- motif3_or_5 %>% group_by(pair) %>% 
#  summarise(PHM=ifelse(HM1.PCA==0, "P1", "P2"),
#            P_HM_Exp=ifelse(HM1.PCA==0, P_HM_Exp_P1 , P_HM_Exp_P2))
#chisq.test(motif3_or_5HM$P_HM_Exp, motif3_or_5noHM$P_HM_Exp)
#p-value = 0.2371

#HM reported in BG, tested by PCA but not detected are less expressed?
HM1.BG.noPCA <- filter(summary_table, HM1.PCA==0 & HM1bg.S.PDB.Kim==1) %>% select(., P1, HM1.PCA, HM1bg.S.PDB.Kim, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
HM2.BG.noPCA <- filter(summary_table, HM2.PCA==0 & HM2bg.S.PDB.Kim==1) %>% select(., P2, HM2.PCA, HM2bg.S.PDB.Kim, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM1.BG.noPCA) <- c("orf", "PCA", "HMbg.S.PDB.Kim", "log10_RNAseq_MTX", "log10_abundance_paxDB")
colnames(HM2.BG.noPCA) <- c("orf", "PCA", "HMbg.S.PDB.Kim", "log10_RNAseq_MTX", "log10_abundance_paxDB")
HM.BG.noPCA <- rbind(HM1.BG.noPCA, HM2.BG.noPCA)

HM1.BG.PCA <- filter(summary_table, HM1.PCA==1 & HM1bg.S.PDB.Kim==1) %>% select(., P1, HM1.PCA, HM1bg.S.PDB.Kim, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
HM2.BG.PCA <- filter(summary_table, HM2.PCA==1 & HM2bg.S.PDB.Kim==1) %>% select(., P2, HM2.PCA, HM2bg.S.PDB.Kim, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM1.BG.PCA) <- c("orf", "PCA", "HMbg.S.PDB.Kim", "log10_RNAseq_MTX", "log10_abundance_paxDB")
colnames(HM2.BG.PCA) <- c("orf", "PCA", "HMbg.S.PDB.Kim", "log10_RNAseq_MTX", "log10_abundance_paxDB")
HM.BG.PCA <- rbind(HM1.BG.PCA, HM2.BG.PCA)

summary(HM.BG.noPCA$log10_RNAseq_MTX)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.978   2.301   2.290   2.687   3.686 
summary(HM.BG.PCA$log10_RNAseq_MTX)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   2.307   2.713   2.734   3.238   4.316 
wilcox.test(HM.BG.noPCA$log10_RNAseq_MTX, HM.BG.PCA$log10_RNAseq_MTX)
#p-value = 1.42e-05

summary(HM.BG.noPCA$log10_abundance_paxDB)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1139  0.9360  1.5587  1.4918  1.9015  3.0806 
summary(HM.BG.PCA$log10_abundance_paxDB)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.1203  1.4930  1.8982  1.9859  2.4783  3.7280 
wilcox.test(HM.BG.noPCA$log10_abundance_paxDB, HM.BG.PCA$log10_abundance_paxDB)
#p-value = 1.651e-05

quantile(c(summary_table$log10.RNAseq.MTX_P2, summary_table$log10.RNAseq.MTX.P1), na.rm=T)
#0%      25%      50%      75%     100% 
#0.000000 1.982271 2.430558 2.847573 5.364649 

##
#number of cases with only one HM and showing a higher expression for the paralog forming a HM
HM1 <- summary_table %>% mutate(diff_expP2P1 = log10.RNAseq.MTX.P2-log10.RNAseq.MTX.P1) %>% 
  filter(., nbHM.inter==1)
nrow(HM1) #77
nrow(filter(HM1, (HM1.PCA==1 & diff_expP2P1<0) | (HM2.PCA==1 & diff_expP2P1>0))) #52 --> #52/77 = 0.68
HM1 <- HM1 %>% group_by(pair) %>% summarise(exp_HM_positif = ifelse(HM1.PCA==1, log10.RNAseq.MTX.P1, log10.RNAseq.MTX.P2),
                                            exp_HM_negatif = ifelse(HM1.PCA==0, log10.RNAseq.MTX.P1, log10.RNAseq.MTX.P2))

wilcox.test(HM1$exp_HM_positif, HM1$exp_HM_negatif)
#p-value = 0.0001104

###############################################################
##        Study factor influencing interaction motifs        ##
###############################################################

############################################################
#                         Duplication                      #
############################################################
# Figure 2. C : frequency of interaction motifs in function of Duplication
chisq.test(summary_table$Duplication, summary_table$motif.categories)
#p-value = 4.044e-06

tot <- summary_table %>% group_by(Duplication) %>% filter(., !is.na(motif.categories)) %>%  
  summarise(tot = length(pair))
n <- summary_table %>% group_by(Duplication, motif.categories) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.categories)) %>% as.data.frame()

freq <- full_join(n, tot, by="Duplication") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n

#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

diff_SSD_WGD_per_motif.simple <- freq %>% group_by(motif.categories) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()


#result
#motif.categories p.value_fisherTest
#1              HET       3.836343e-01
#2               HM       2.529560e-05
#3           HM&HET       1.630493e-06
#4               NI       2.860273e-01

#--> Fig2E

# Separating pairs with only one HM to those with 2 HMS coded with motif number
chisq.test(summary_table$Duplication, summary_table$motif.number)
#--> FigS8

####################################################################################################
#   Compare onhologs from Homeologization (Homeologs) and from tetraploidization (True Onhologs)    #
####################################################################################################

Gabaldon <- read.table("data/Gabaldon.csv", header = T, sep=";")
Gabaldon$pair <- apply(cbind(as.character(Gabaldon$Seed.protein), 
                             as.character(Gabaldon$Ohnolog)), 
                       1, function(x) paste(sort(x), collapse="."))
Gabaldon <- unique(select(Gabaldon, pair, Inferred.topology))

#4 pairs were found belonging to pre-KLE and real WGD depending of the sens seed.protein / onholog,
#we should remove them
Gabaldon %>% group_by(pair) %>% summarise(nbr=n()) %>% filter(., nbr>1)
# A tibble: 4 x 2
#pair              nbr
#1 YCL040W.YDR516C     2
#2 YER039C.YGL225W     2
#3 YER064C.YIL056W     2
#4 YLR293C.YOR185C     2

Gabaldon <- Gabaldon %>% filter(., pair!= "YCL040W.YDR516C" &
                            pair!= "YER039C.YGL225W" &
                            pair!= "YER064C.YIL056W" &
                            pair!= "YLR293C.YOR185C")


#Inferred.topology A - A = same parent = True ohnologs
#Inferred.topology A - AB = two different parents = Homeologs

#Test compare HM and HM&HET in function of Inffered.topology

Gabaldon %<>% group_by(pair) %>% 
  mutate(Origin.of.WGDs = ifelse(Inferred.topology=="A – A","True_ohnologs","Homeologs")) %>%
  select(., pair, Origin.of.WGDs)

summary_table <- left_join(summary_table, Gabaldon, by='pair')

tot <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs) %>% filter(., !is.na(motif.categories)) %>%  
  summarise(tot = length(pair))
n <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs, motif.categories) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.categories)) %>% as.data.frame()

freq <- full_join(n, tot, by="Origin.of.WGDs") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n


#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))
diff_Inffered.topo_per_motif <- freq %>% group_by(motif.categories) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
#result:
#motif.categories p.value_fisherTest
#1              HET         0.40946337
#2               HM         0.04951537
#3           HM&HET         0.01095495
#4               NI         1.00000000


freq$Origin.of.WGDs <- factor(freq$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))


#--> Fig2E2 


# Separating pairs with only one HM to those with 2 HMS coded with motif number
tot <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs) %>% filter(., !is.na(motif.number)) %>%  
  summarise(tot = length(pair))
n <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs, motif.number) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.number)) %>% as.data.frame()

freq <- full_join(n, tot, by="Origin.of.WGDs") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n


#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

diff_Inffered.topo_per_motif <- freq %>% group_by(motif.number) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
#motif.number p.value_fisherTest
#1            1         1.00000000
#2            2         0.04222921
#3            3         1.00000000
#4            4         0.40946337
#5            5         0.22081352
#6            6         0.10468644


############################################################
#              HET and age of duplication HM               #
############################################################
SSD_age <- read.table("output/SSD_corrected_age_2019_02_AM.csv", header = T)
SSD_age <- select(SSD_age, pair, age_dupli2)
colnames(SSD_age) <- c("pair", "age.group")
summary_table <- left_join(summary_table, SSD_age, by="pair")
summary_table$age.group[summary_table$Origin.of.WGDs=="True_ohnologs"] <- 1
summary_table$age.group[summary_table$Origin.of.WGDs=="Homeologs"] <- 2

#FigS6B Age groups --> ratio HM&HET/HM
HM <- filter(summary_table, motif.categories=="HM" & !is.na(age.group)) %>% 
  group_by(age.group, Duplication) %>%  summarise(HM=length(motif.categories)) %>% as.data.frame()
HM_HET <-  filter(summary_table, motif.categories=="HM&HET" & !is.na(age.group)) %>% 
  group_by(age.group, Duplication) %>%  summarise("HM_HET"=length(motif.categories)) %>% as.data.frame()
HM$age.group <- as.character(HM$age.group)
HM_HET$age.group <- as.character(HM_HET$age.group)
HM$HM <- as.numeric(as.character(HM$HM))
HM_HET$HM_HET <- as.numeric(as.character(HM_HET$HM_HET))
ratio_HM_HET <- left_join(HM, HM_HET, by=c("age.group", "Duplication")) 
ratio_HM_HET$HM_HET[is.na(ratio_HM_HET$HM_HET)] <- 0
ratio_HM_HET <- ratio_HM_HET %>% mutate(ratioHETsurHM=(HM_HET/HM))
ratio_HM_HET <- ratio_HM_HET %>% group_by(age.group, Duplication) %>%
  mutate(effectif=paste(c(HM_HET,HM), collapse = "/"))

#not enougth SSD with HM&HET or HM observed per age group


#############################################
#                 Complexes                 #
#############################################

#Import information of involvment of proteins into different complexes
all <- read.table("data/dscercomplexes_combo (1).tsv", header=T,comment.char="@", sep="\t", fill=T,quote = "\"")

all <- all[,c(2,3)]
names(all) <- c("complex", "subunits")

#split to multiple lines
test <- all %>% 
  mutate(orfs = strsplit(as.character(subunits), ",")) %>% 
  unnest(orfs) %>%
  filter(!is.na(orfs), orfs !="") %>% 
  group_by(complex) %>%
  mutate(size_comp = n()) %>%
  ungroup()

test <- test %>% mutate(orfs = str_replace_all(orfs , " ", "")) 

comp_per_gene <- test %>%  filter(size_comp>2) %>% 
  group_by(orfs) %>%
  summarise(list.of.cplx = paste(sort(unique(c(as.character(complex)))), collapse=";")) %>%
  ungroup()

summary_table <- left_join(summary_table, comp_per_gene, by = c("P1"="orfs")) %>%
  left_join(., comp_per_gene, by=c("P2"="orfs"))

summary_table <- dplyr::rename(summary_table, list.of.cplx.P1 = list.of.cplx.x)
summary_table <- dplyr::rename(summary_table, list.of.cplx.P2 = list.of.cplx.y)

#check if P1 and P2 are in same complex

same_comp <- function(col1,col2){
  
  if (col1 !="" & col2 !="" & col1 !="NA" & col2 !="NA" & !is.na(col1) & !is.na(col2)) {
    a <- unlist(strsplit(col1, split=";"))
    b <- unlist(strsplit(col2, split=";"))
    res = intersect(a,b)
    resl =length(res)
    
    if (resl>0) {
      return(resl)
    } 
    else {
      return(0)
    }
  }
  
  else{
    return(NA)
  }
}
same_comp <- Vectorize(same_comp)

summary_table$list.of.cplx.P1 = as.character(summary_table$list.of.cplx.P1)
summary_table$list.of.cplx.P2 = as.character(summary_table$list.of.cplx.P2)

summary_table <- summary_table %>% 
  mutate(sim.complex= same_comp(list.of.cplx.P1,list.of.cplx.P2)) %>%
  mutate(shared.complex = ifelse(sim.complex>0,1,0))

data.frame(table(summary_table$motif.categories, summary_table$shared.complex))

data.frame(table(summary_table$motif.categories))

############################################################
#                  Function similarity                     #
############################################################


go <- read.table("data/go_slim_mapping.tab", sep="\t", header=F)
names(go) <- c("orfid", "geneid", "code", "cat", "desc", "gocode", "status")

#Some missing data are replaced by terms
#cellular_component, biological_process, molecular_function
#We removed these terms because this does not make sense
test <- data.frame(table(go$desc))
tofilter <- c("cellular_component","not_yet_annotated", "biological_process","molecular_function", "other")

go <- go %>% mutate(desc=ifelse(desc %in% tofilter, "NA", as.character(desc)))
go$orfid = as.character(go$orfid)

#some genes have annotations on multiple lines. Need to create a list of those annotations

GOFsum <- go %>% filter(cat=="F") %>% 
  dplyr::group_by(orfid) %>% 
  dplyr::summarise(molfunc = paste(sort(unique(c(as.character(desc)))), collapse=";")) %>%
  ungroup()

GOPsum <- go %>% filter(cat=="P") %>% 
  dplyr::group_by(orfid) %>% 
  dplyr::summarise(cellproc = paste(sort(unique(c(as.character(desc)))), collapse=";")) %>%
  ungroup()

GOCsum <- go %>% filter(cat=="C") %>% 
  dplyr::group_by(orfid) %>% 
  dplyr::summarise(cellcomp = paste(sort(unique(c(as.character(desc)))), collapse=";")) %>%
  ungroup()

#merge
GO_stuff <- full_join(GOFsum, GOPsum, by=c("orfid"="orfid")) %>%
  full_join(., GOCsum, by=c("orfid"="orfid")) 

#add this data to dataframe
summary_table <- left_join(summary_table, GO_stuff, by = c("P1"="orfid")) %>%
  left_join(., GO_stuff, by=c("P2"="orfid"))

summary_table <- dplyr::rename(summary_table, mol.fct.P1 = molfunc.x)
summary_table <- dplyr::rename(summary_table, mol.fct.P2 = molfunc.y)
summary_table <- dplyr::rename(summary_table, bio.proc.P1 = cellproc.x)
summary_table <- dplyr::rename(summary_table, bio.proc.P2 = cellproc.y)
summary_table <- dplyr::rename(summary_table, cell.comp.P1 = cellcomp.x)
summary_table <- dplyr::rename(summary_table, cell.comp.P2 = cellcomp.y)

#write a function that can calculate a similarity for go terms
#take the intersection over the union to measure a score
similarity <- function(col1,col2){
  
  if (!is.na(col1) & !is.na(col2) & col1 !="" & col2 !="") {
    a <- unlist(strsplit(col1, split=";"))
    b <- unlist(strsplit(col2, split=";"))
    
    res = intersect(a,b)
    resl = length(res)
    tot = length(union(a,b))
    
    if (resl>0) {
      return(resl/tot)
    } 
    else {
      return(0)
    }
  }
  
  else{
    return(NA)
  }
}
similarity <- Vectorize(similarity)


#the analsysis below show that heteromers seem to be more similar functionally than non heteromers
summary_table <- summary_table %>%
  mutate(sim.mol.fct.P1P2= similarity(mol.fct.P1, mol.fct.P2),
         sim.bio.proc.P1P2= similarity(bio.proc.P1, bio.proc.P2),
         sim.cell.comp.P1P2= similarity(cell.comp.P1, cell.comp.P2)) 

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}

summary_table$int <- interaction(summary_table$motif.categories, summary_table$Duplication)



#look at phenotype file
phe <- read.table("data/phenotype_data.tab", sep="\t", header=F, fill=T,as.is = T,quote = "\"")
names(phe) <- c("id", "type", "name", "sgdID","ref", "exp_type", "mutant_type","allele", "back", "pheno", "chem", "cond",
                "details", "reporter")


#filter to keep relevant data only
phe <- phe %>% filter(back=="S288C")

#check what type of mutations we have
data.frame(table(phe$mutant_type))
#we could focus on null because this is where most data is
phe <- phe %>% filter(mutant_type=="null")
#filter to keep only paralogs and see how much information there is

pset <- c(as.character(summary_table$P1), as.character(summary_table$P2))
phe <- phe %>% filter(id %in% pset)
# look at the information on 
head(data.frame(table(phe$pheno)))

#code information as to whether viable, non-viable, haploinsufficient or haplosufficient
phe <- phe %>% mutate(haploins = ifelse(pheno=="haploinsufficient", 1,
                                  ifelse(pheno=="haploproficient", 0, "NA")),
                viability = ifelse(pheno=="viable", 1,
                                   ifelse(pheno=="inviable", 0, "NA")))

#keep information with specific phenotypes, those that contain :

phe <- phe %>% mutate(phenotype_cat = ifelse(str_detect(pheno, ":"), paste(pheno,chem, cond, sep="-"), NA)) %>%
  mutate(phenotype_cat  = str_replace_all(phenotype_cat, ";", ","))  

head(data.frame(table(phe$phenotype_cat)))

temp_phe <- phe %>% filter(!is.na(phenotype_cat)) %>% dplyr::group_by(id) %>% 
  dplyr:: summarise(list_pheno = paste(sort(unique(c(as.character(phenotype_cat)))), collapse=";")) %>%
  ungroup()

#add to summary_table data.frame
summary_table <- left_join(summary_table, temp_phe, by = c("P1"="id")) %>%
  left_join(., temp_phe, by=c("P2"="id"))

summary_table$list_pheno.x[is.na(summary_table$list_pheno.x)]<- "NA"
summary_table$list_pheno.y[is.na(summary_table$list_pheno.y)]<- "NA"

summary_table <- summary_table %>% dplyr::mutate(sim.pheno.P1P2= similarity(list_pheno.x,list_pheno.y)) 

resPhe <- summary_table  %>% filter(motif.categories !="NI" & motif.categories !="HET") %>%
  dplyr::group_by(Duplication) %>% 
  do(w = wilcox.test(sim.pheno.P1P2~motif.categories, data=., paired=FALSE)) %>% 
  dplyr::summarise(Duplication, Wilcox = w$p.value)


summary_table <- dplyr::rename(summary_table, list.pheno.P1 = list_pheno.x)
summary_table <- dplyr::rename(summary_table, list.pheno.P2 = list_pheno.y)


####test for similarity of genetic interactions####
#hypothesis is that heteromers would maintain more similarity

#some cleaning to do here
sp <- fread("data/cc_NN.txt", sep="\t", header=T) 
##col names and row names are not unique so need to make them unique and then switch
spt <-sp
#remove first colum and second row
sp <- sp[-1,]
sp <- sp[,-1]

#save rown names and col names
rown <- sp$V2
sp <- sp[,-1]
coln <- names(sp)


rown <- data.frame(cbind(seq(1,length(rown),1), rown))
names(rown) <- c("row_num", "orf1")
rown$orf1 <- as.character(rown$orf1)

coln<- data.frame(cbind(seq(1,length(coln),1), coln))
names(coln) <- c("col_num", "orf2")
coln$orf2 <- as.character(coln$orf2)

names(sp) <- as.character(coln$col_num)

#transform the matrix into list of pairs of genes
sp[lower.tri(sp, diag = TRUE ) ]  <- NA
sp$row_num <- as.character(rownames(sp))

sp.molten <- melt( sp, na.rm= TRUE, id.vars="row_num",
                   value.name="corrgi", variable.name="col_num" )

sp.molten <- left_join(sp.molten,rown,by = c("row_num" = "row_num"))
sp.molten <-left_join(sp.molten,coln,by = c("col_num" = "col_num"))

#consider only the genes that are in the paralog set so 
#we consider only the genes we studied
pset <- c(as.character(summary_table$P1), as.character(summary_table$P2))
sp.molten$orf1 = as.character(sp.molten$orf1)
sp.molten$orf2 = as.character(sp.molten$orf2)

sp.molten %<>% dplyr::mutate(orf1_in = ifelse(orf1 %in% pset,1,0), 
                             orf2_in = ifelse(orf2 %in% pset,1,0)) %>%
  filter(orf1_in==1, orf2_in==1)

#create unique id pair to use as unique ID to merge with summary_table
sp.molten$pairID <- apply(cbind(as.character(sp.molten$orf1), 
                                as.character(sp.molten$orf2)), 
                          1, function(x) paste(sort(x), collapse="."))

sp.molten$corrgi <- as.numeric(sp.molten$corrgi)
#some pairs have been measured several times. Take the median
mean_co_gi <- sp.molten %>% dplyr::group_by(pairID) %>% 
  dplyr::summarise(med.gi.cor=median(corrgi))


#combine with summary_table

summary_table <- left_join(summary_table, mean_co_gi, by=c("pair"="pairID"))

Dupli <- select(summary_table, pair, Duplication)
test <- left_join(mean_co_gi, Dupli,  by=c("pairID"='pair'))
test <- test %>% dplyr::mutate(type_pair = ifelse(Duplication=="SSD", 'SSD',
                                         ifelse(Duplication=="WGD", "WGD","Not_in_duplicated_set")))

test$type_pair[is.na(test$type_pair)] <- "NoDup"

pairwise.wilcox.test(test$med.gi.cor,as.factor(test$type_pair))

#paralogs have more similar genetic interaction profiles than non paralog pairs (among the same set of genes)
#this is expected
#now compare pairs of het versus non het

resGI <- summary_table  %>% filter(motif.categories !="NI" & motif.categories !="HET") %>%
  dplyr:: group_by(Duplication) %>% 
  do(w = t.test(med.gi.cor~motif.categories, data=., paired=FALSE)) %>% 
  dplyr::summarise(Duplication, TTest = w$p.value)

######Localization
loc <- read.table("data/localization_data.csv", sep = "\t", header=T, row.names = NULL)
loc <- select(loc, yORF, localization.summary)
colnames(loc) <- c("P1", "localization.P1")
summary_table <- left_join(summary_table, loc, by="P1")
colnames(loc) <- c("P2", "localization.P2")
summary_table <- left_join(summary_table, loc, by="P2")


summary_table$localization.P1[summary_table$localization.P1==""] <- NA
summary_table$localization.P2[summary_table$localization.P2==""] <- NA

#localization separate with , --> change the function similarity
similarity <- function(col1,col2){
  
  if (!is.na(col1) & !is.na(col2) & col1 !="" & col2 !="") {
    a <- unlist(strsplit(col1, split=","))
    b <- unlist(strsplit(col2, split=","))
    
    res = intersect(a,b)
    resl = length(res)
    tot = length(union(a,b))
    
    if (resl>0) {
      return(resl/tot)
    } 
    else {
      return(0)
    }
  }
  
  else{
    return(NA)
  }
}
similarity <- Vectorize(similarity)

#summary_table <- select(summary_table, Duplication, pair, motif.categories, localization.P1, localization.P2)
summary_table$localization.P1 <- as.character(summary_table$localization.P1)
summary_table$localization.P2 <- as.character(summary_table$localization.P2)

summary_table <- summary_table %>% dplyr::mutate(sim.loc= similarity(localization.P1, localization.P2)) 


########## transcription factor 
factor.transc <- read.table("data/RegulationTwoColumnTable_Documented_2013927.csv", header = F, sep=";")
colnames(factor.transc) <- c("TF", "orf")
factor.transc$TF <- as.character(factor.transc$TF)
factor.transc$orf  <- as.character(factor.transc$orf)

stand.to.syst <- read.table("data/stand_to_syst_name.csv", header = F, sep="\t")
stand.to.syst <- unique(stand.to.syst)
colnames(stand.to.syst) <- c("stand", "syst")
stand.to.syst$stand <- as.character(stand.to.syst$stand)
stand.to.syst$syst <- as.character(stand.to.syst$syst)

factor.transc <- left_join(factor.transc, stand.to.syst, by=c("orf"="stand"))
factor.transc <- select(factor.transc, TF, syst)
colnames(factor.transc) <- c("TF.P1", "P1")

summary_tablef <- select(summary_table, pair, P1, P2)
summary_tablef <- left_join(summary_tablef, factor.transc, by="P1")
colnames(factor.transc) <- c("TF.P2", "P2")
summary_tablef <- left_join(summary_tablef, factor.transc, by="P2")
summary_tablef2 <- summary_tablef %>%  dplyr::group_by(pair, P1, P2) %>% dplyr::summarise(tot.TF.P1 = paste(sort(unique(TF.P1)), collapse = ","),
                                                                    tot.TF.P2 = paste(sort(unique(TF.P2)), collapse = ","))

similarity <- function(col1,col2){
  
  if (!is.na(col1) & !is.na(col2) & col1 !="" & col2 !="") {
    a <- unlist(strsplit(col1, split=","))
    b <- unlist(strsplit(col2, split=","))
    
    res = intersect(a,b)
    resl = length(res)
    tot = length(union(a,b))
    
    if (resl>0) {
      return(resl/tot)
    } 
    else {
      return(0)
    }
  }
  
  else{
    return(NA)
  }
}
similarity <- Vectorize(similarity)


summary_tablef2$tot.TF.P1 <- as.character(summary_tablef2$tot.TF.P1)
summary_tablef2$tot.TF.P2 <- as.character(summary_tablef2$tot.TF.P2)

summary_tablef2 <- summary_tablef2 %>% dplyr::mutate(sim.TF= similarity(tot.TF.P1, tot.TF.P2)) 

summary_table <- left_join(summary_table, summary_tablef2, by=c("pair", "P1", "P2"))

resTF <- summary_table %>% filter(motif.categories =="HM" | motif.categories =="HM&HET") %>%
  dplyr::group_by(Duplication) %>% 
  do(w = wilcox.test(sim.F~motif.categories, data=., paired=FALSE)) %>% 
  dplyr::summarise(Duplication, Wilcox = w$p.value)

summary_table$int <- interaction(summary_table$motif.categories, summary_table$Duplication)

############################################################
#                  Percentage identity                    #
############################################################

sequences <- read.table("data/SGD_seq_S288C_ref_orf_trans_2015_01_13_AM.csv", sep="\t", header=F)
colnames(sequences) <- c("P1", "P1_seq")

pairs <- select(summary_table, pair, P1, P2)

pairs_seq <- left_join(pairs, sequences, by="P1")
colnames(sequences) <- c("P2", "P2_seq")
pairs_seq <- left_join(pairs_seq, sequences, by="P2")

write.table(pairs_seq, file="output/pairs_seq_2019_02_AM.csv", sep="\t",  quote=F, row.names=F)


library(Biostrings)
pairs_seq_align <- pairs_seq %>% group_by(pair, P1, P2) %>% 
  summarise (pid = pid((pairwiseAlignment(pattern = P1_seq, subject = P2_seq)), type="PID1")) %>% ungroup()
write.table(pairs_seq_align, file="output/pairs_perc_ident_2019_02_AM.csv", sep="\t",  quote=F, row.names=F)
pairs_seq_align <- select(pairs_seq_align, pair, pid)

summary_table <- left_join(summary_table, pairs_seq_align, by="pair")

#--> Fig2F


##########################################################################
#                 Co-expression of paralogs and HET                      #
##########################################################################
### Separation of paralogs in space or time can solve the constraints imposed by molecular pleiotropy

#Figure 6
expression_matrix <- read.table("data/Ihmels_expression_matrix.csv", sep="\t", header=F)
orfs <- read.table("data/Ihmels_orfs.csv", sep="\t", header=F)
expression <- cbind(orfs, expression_matrix)

pairs <- select(summary_table, pair, P1, P2)
names(expression)[1] <- "P1"
abond_P1 <- left_join(pairs, expression, by="P1")
abond_P1$paralog = "P1"
names(expression)[1] <- "P2"
abond_P2 <- left_join(pairs, expression, by="P2")
abond_P2$paralog = "P2"

abond_correl <- rbind(abond_P1, abond_P2)
abond_correl <- gather(abond_correl, condition, expression, V1:V1706)

abond_correl <- abond_correl %>% dplyr::group_by(pair, P1, P2) %>% 
  dplyr::summarise(expression.correl.coeff=cor(expression[paralog=="P1"], expression[paralog=="P2"], method="pearson", use="na.or.complete")) %>% as.data.frame()
abond_correl <- select(abond_correl, pair, expression.correl.coeff)

summary_table <- left_join(summary_table, abond_correl, by="pair")


#correlation coexpression vs PI
SSD <- filter(summary_table, Duplication=="SSD")
WGD <- filter(summary_table, Duplication=="WGD")
cor.test(SSD$expression.correl.coeff, SSD$pid, method = "spearman") #R=0.31131   p-value = 2.663e-08
cor.test(WGD$expression.correl.coeff, WGD$pid, method = "spearman") #R=0.2163993  p-value = 0.000271

#correlation with sequence identity and correlation of expression, 
#we need to test if each of them has an effect on its own
#only HM and HM&HET
fromHM <- filter(summary_table, motif.categories=="HM" | motif.categories=="HM&HET")

glm1<-glm(fromHM$motif.categories ~ 
            fromHM$expression.correl.coeff+
            fromHM$pid, family=binomial)

summary(glm1)
anova(glm1,test="LRT") 


SSD <- filter(fromHM, Duplication=="SSD" & (motif.categories=="HM" | motif.categories=="HM&HET"))
glm1<-glm(SSD$motif.categories ~ 
            SSD$expression.correl.coeff+
            SSD$pid, family=binomial)

summary(glm1)
anova(glm1,test="LRT") 

#                            summary_table Deviance Resid. summary_table Resid. Dev  Pr(>Chi)    
#NULL                                          236     289.38              
#SSD$expression.correl.coeff  1   17.935       235     271.44 2.286e-05 ***
#SSD$pid                      1   13.014       234     258.43 0.0003091 ***


WGD <- filter(fromHM, Duplication=="WGD" & (motif.categories=="HM" | motif.categories=="HM&HET"))

glm1<-glm(WGD$motif.categories ~ 
            WGD$expression.correl.coeff+
            WGD$pid, family=binomial)

summary(glm1)
anova(glm1,test="LRT") 

#summary_table Deviance Resid. summary_table Resid. Dev  Pr(>Chi)    
#NULL                                          217     301.55              
#WGD$expression.correl.coeff  1   1.9756       216     299.58    0.1599    
#WGD$pid                      1  18.6102       215     280.97 1.604e-05 ***



Homeologs <- filter(WGD, Origin.of.WGDs=="Homeologs" & (motif.categories=="HM" | motif.categories=="HM&HET"))

glm1<-glm(Homeologs$motif.categories ~ 
            Homeologs$expression.correl.coeff+
            Homeologs$pid, family=binomial)

summary(glm1)
anova(glm1,test="LRT") 


GD <- filter(WGD, Origin.of.WGDs=="True_ohnologs" & (motif.categories=="HM" | motif.categories=="HM&HET"))

glm1<-glm(GD$motif.categories ~ 
            GD$expression.correl.coeff+
            GD$pid, family=binomial)

summary(glm1)
anova(glm1,test="LRT") 
#summary_table Deviance Resid. summary_table Resid. Dev Pr(>Chi)  
#NULL                                          62     78.742           
#GD$expression.correl.coeff  1   0.5333        61     78.209  0.46524  
#GD$pid                      1   3.1689        60     75.040  0.07505 .



##########################################
###             Dependency            ####
##########################################
dep <- read.table("data/paralog_dep.csv", header=T, sep= "\t")
dep$pair <- apply(cbind(as.character(dep$P1), 
                                                as.character(dep$P2)), 
                                   1, function(x) paste(sort(x), collapse="."))
dep <- select(dep, pair, dependency)
summary_table <- left_join(summary_table, dep, by='pair')
dep <- filter(summary_table, !is.na(dependency))
dep %<>% mutate(nHM = HM1.final+HM2.final)
ct_hom <- table(dep$dependency, dep$motif.categories)

#HET HM HM&HET NI
#Dep     0  2     16  0
#Indep   0 11     12  2

nrow(filter(dep, dependency=="Dep", motif.categories=="HM&HET", nHM==2))
nrow(filter(dep, dependency=="Dep", motif.categories=="HM&HET", nHM==1))
nrow(filter(dep, dependency=="Indep", motif.categories=="HM&HET", nHM==2))
nrow(filter(dep, dependency=="Indep", motif.categories=="HM&HET", nHM==1))

#Similarity of fct and dep?
wilcox.test(dep$sim.bio.proc.P1P2[dep$dependency=="Indep"], dep$sim.bio.proc.P1P2[dep$dependency=="Dep"])
#p-value = 0.5458
wilcox.test(dep$sim.mol.fct.P1P2[dep$dependency=="Indep"], dep$sim.mol.fct.P1P2[dep$dependency=="Dep"])
#p-value = 0.3021
wilcox.test(dep$mean_gi_cor[dep$dependency=="Indep"], dep$mean_gi_cor[dep$dependency=="Dep"])
#p-value = 0.1957
wilcox.test(dep$sim.pheno.P1P2[dep$dependency=="Indep"], dep$sim.pheno.P1P2[dep$dependency=="Dep"])
#p-value = 0.4135

Deluna <- read.table("data/response_paralog_deletion_Deluna2010.csv", header=T, sep= "\t")
syst.name <- read.table("data/resp_del_paralog.tsv", header=F, sep= "\t")

colnames(syst.name) <- c("Syst.name.GFP", "X1.GFP")
Deluna <- left_join(Deluna, syst.name, by="X1.GFP")

colnames(syst.name) <- c("Syst.name.Δx2", "Δx2")
Deluna <- left_join(Deluna, syst.name, by="Δx2")

filter(Deluna, is.na(Syst.name.GFP))
Deluna$Syst.name.GFP <- as.character(Deluna$Syst.name.GFP)
Deluna %<>% group_by(ID) %>% mutate(Syst.name.GFP = ifelse(is.na(Syst.name.GFP) & X1.GFP!="RHR2" & X1.GFP!="HOR2", X1.GFP,
                                                    ifelse(is.na(Syst.name.GFP) & X1.GFP=="RHR2", "YIL053W", 
                                                    ifelse(is.na(Syst.name.GFP) & X1.GFP=="HOR2", "YER062C", Syst.name.GFP))))
filter(Deluna, is.na(Syst.name.Δx2))
Deluna$Syst.name.Δx2 <- as.character(Deluna$Syst.name.Δx2)
Deluna %<>% group_by(ID) %>% mutate(Syst.name.Δx2 = ifelse(is.na(Syst.name.Δx2) & Δx2!="RHR2" & Δx2!="HOR2", Δx2,
                                                    ifelse(is.na(Syst.name.Δx2) & Δx2=="RHR2", "YIL053W", 
                                                    ifelse(is.na(Syst.name.Δx2) & Δx2=="HOR2", "YER062C", Syst.name.Δx2))))
Deluna <- Deluna %>% ungroup()

Deluna$pair <- apply(cbind(as.character(Deluna$Syst.name.GFP), 
                        as.character(Deluna$Syst.name.Δx2)), 
                  1, function(x) paste(sort(x), collapse="."))

Signif <- read.table("data/response_paralog_deletion_Deluna2010_signif.csv", header=F, sep= "\t")
colnames(Signif) <- c("X1.GFP", "response")

Deluna <- left_join(Deluna, Signif, by="X1.GFP")
Deluna$response <- as.character(Deluna$response)
Deluna$response[is.na(Deluna$response)] <- "no_rep"

Deluna <- select(Deluna, pair, X1.GFP, Δx2, response)
Deluna %<>% group_by(pair) %>% summarise(Response = ifelse (paste(sort(response), collapse=".") == "no_rep.no_rep", "no_rep",
                                                    ifelse(paste(sort(response), collapse=".") =="neg.no_rep", "neg",
                                                    ifelse(paste(sort(response), collapse=".") =="no_rep.pos", "pos",
                                                    ifelse(paste(sort(response), collapse=".") =="pos.pos", "pos",
                                                    ifelse(paste(sort(response), collapse=".") =="neg.neg", "neg", 
                                                           paste(sort(response), collapse=".")))))))



Deluna <- left_join(Deluna, summary_table, by="pair")
Deluna <- select(Deluna, pair, P1, P2, HM1.final, HM2.final, HET.final, motif.categories, log10.RNAseq.MTX.P1, log10.RNAseq.MTX.P2, Response)
ct_hom <- table(Deluna$Response, Deluna$motif.categories)

#HET HM HM&HET NI
#neg      1  1      2  0
#no_rep   4 33     39  4
#pos      0  6      9  1

write.table(summary_table, file="output/TableS1.csv", sep="\t",  quote=F, row.names=F)
