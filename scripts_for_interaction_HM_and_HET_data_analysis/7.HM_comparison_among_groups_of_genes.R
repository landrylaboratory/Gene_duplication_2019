rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)

##################################################################################
#########Enrichment of homomers among proteins coding by duplicated genes#########
##################################################################################
dff.S.bg.Kim.PCA <- read.table("output/HM.data.csv", sep="\t", header=T)

#we find that about 31% of yeast proteins form HMs:
nrow(filter(dff.S.bg.Kim.PCA, HM.bg.kim.S.PCA==1)) / nrow(dff.S.bg.Kim.PCA)

#keeping successive SSDs

dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para== "ssd-successiv"] <- "ssd"
ct_hom <- table(droplevels(dff.S.bg.Kim.PCA)$type_para, dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA)

dbg <- data.frame(class=c("No HM S", "HM S", 
                          "No HM SSD", "HM SSD",
                          "No HM 2D", "HM 2D",
                          "No HM WGD", "HM WGD"),
                  count=c(ct_hom[1,1], ct_hom[1,2],
                          ct_hom[2,1],ct_hom[2,2],
                          ct_hom[3,1],ct_hom[3,2],
                          ct_hom[4,1],ct_hom[4,2]))
ct_hom <-t(ct_hom)
chisq.test(ct_hom) #p-value < 2.2e-16

fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA, as.factor(dff.S.bg.Kim.PCA$type_para), workspace=2e8)
#p-value < 2.2e-16
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"]), workspace=2e8)
#p-value < 2.2e-16
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"]), workspace=2e8)
#p-value =  9.875e-06
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"]), workspace=2e8)
#p-value = 2.438e-13

tot <- sum(dbg$count)
dbg <- dbg %>% separate(class, c("HMstatue", "dupli"), "HM ")
dbg$HMstatue[(dbg$HMstatue)==""] <- "Yes"
dbg$HMstatue[(dbg$HMstatue)=="No "] <- "No"
dbg <- dbg %>% group_by(dupli) %>% mutate(sum_per_dupli=sum(count)) %>% as.data.frame()
dbg_percent <- filter(dbg, HMstatue=="Yes") %>% group_by(dupli) %>% 
  summarise(percentHM = (count/sum_per_dupli)*100) %>% as.data.frame()

genome <- filter(dbg, HMstatue=="Yes") %>% group_by(HMstatue) %>% 
  summarise(sum_per_genome = (sum(count)/tot)*100) %>% as.data.frame()
genome$HMstatue[(genome$HMstatue)=="Yes"] <- "Total"
colnames(genome) <- c("dupli", "percentHM")
dbg_percent <- rbind(genome, dbg_percent)

dbg_percent$n <- c(1,5,2,3,4)
dbg_percent$dupli <- factor(dbg_percent$dupli, levels = dbg_percent$dupli[order(dbg_percent$n)])  

###comparison expression WGD / SSD
RNAseq <- read.table("data/count_table.csv", header=T, sep="\t")
RNAseq %<>% group_by(gene) %>% mutate(log10_med_DMSO = log10(median(DMSO1, DMSO2, DMSO2)+1),
                                      log10_med_MTX = log10(median(MTX1, MTX2, MTX3)+1)) %>% as.data.frame()
RNAseq <- select(RNAseq, gene, log10_med_MTX)
colnames(RNAseq) <- c("orf", "log10_RNA_med_MTX")
dff.RNA <- left_join(dff.S.bg.Kim.PCA, RNAseq, by="orf")
dff.RNA <- filter(dff.RNA, !is.na(log10_RNA_med_MTX))
dff.RNA.ssd <- filter(dff.RNA, type_para=="ssd")
dff.RNA.wgd <- filter(dff.RNA, type_para=="wgd")
wilcox.test(dff.RNA.ssd$log10_RNA_med_MTX, dff.RNA.wgd$log10_RNA_med_MTX)
#p-value = 1.851e-11
dff.RNA.S <- filter(dff.RNA, type_para=="S")
wilcox.test(dff.RNA.wgd$log10_RNA_med_MTX,dff.RNA.S$log10_RNA_med_MTX)
#p-value < 2.2e-16
wilcox.test(dff.RNA.ssd$log10_RNA_med_MTX,dff.RNA.S$log10_RNA_med_MTX)
#p-value < 2.2e-16

cor.test(dff.RNA$med, dff.RNA$log10_RNA_med_MTX, method = "spearman")
#p-value = 0.01066 R=0.04228025

###With HM yet reported in BG
dff.S.bg.Kim.PCA <- left_join(dff.S.bg.Kim.PCA, RNAseq, by="orf")
dff.S.bg.Kim.HMreport <- filter(dff.S.bg.Kim.PCA, n.study>0)
cor.test(dff.S.bg.Kim.HMreport$med, dff.S.bg.Kim.HMreport$log10_RNA_med_MTX, method = "spearman")
#p-value = 0.05399 R=0.04590053 

###With paxdb
paxdb <- read.table("data/PaxDB_4932-WHOLE_ORGANISM-integrated_AM.txt") %>% select(V2, V3)
colnames(paxdb) <- c("orf", "paxdb")

dff.S.bg.Kim.PCA <- left_join(dff.S.bg.Kim.PCA, paxdb, by="orf")
cor.test(dff.S.bg.Kim.PCA$med, dff.S.bg.Kim.PCA$paxdb, method = "spearman")
#p-value = 1.06e-09 R=0.1016517

dff.S.bg.Kim.HMreport <- filter(dff.S.bg.Kim.PCA, n.study>0)
cor.test(dff.S.bg.Kim.HMreport$med, dff.S.bg.Kim.HMreport$paxdb, method = "spearman")
#p-value = 2.297e-07 R=0.1228353 

dff.paxdb.wgd <- filter(dff.S.bg.Kim.PCA, type_para=="wgd")
dff.paxdb.ssd <- filter(dff.S.bg.Kim.PCA, type_para=="ssd")
dff.paxdb.S <- filter(dff.S.bg.Kim.PCA, type_para=="S")
wilcox.test(dff.paxdb.ssd$paxdb, dff.paxdb.S$paxdb) #p-value < 2.2e-16
wilcox.test(dff.paxdb.S$paxdb, dff.paxdb.wgd$paxdb) #p-value < 2.2e-16
wilcox.test(dff.paxdb.ssd$paxdb, dff.paxdb.wgd$paxdb) #p-value =  0.04574

#glm
dff.S.bg.Kim.PCA <- dff.S.bg.Kim.PCA %>% mutate(Duplication.binaire = ifelse(type_para=="S", 0, 1))


glm1<-glm(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA ~ 
            dff.S.bg.Kim.PCA$Duplication.binaire+
            dff.S.bg.Kim.PCA$log10_RNA_med_MTX, family=binomial)
summary(glm1)
anova(glm1,test="LRT") 


#                                     Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#NULL                                                  6070     7673.8              
#dff.S.bg.Kim.PCA$Duplication.binaire  1   104.44      6069     7569.4 < 2.2e-16 ***
#dff.S.bg.Kim.PCA$log10_RNA_med_MTX    1   413.23      6068     7156.2 < 2.2e-16 ***


glm2<-glm(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA ~ 
            dff.S.bg.Kim.PCA$Duplication.binaire+
            dff.S.bg.Kim.PCA$paxdb, family=binomial)
summary(glm2)
anova(glm2,test="LRT") 

#Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#NULL                                                  5982     7606.3              
#dff.S.bg.Kim.PCA$Duplication.binaire  1   89.592      5981     7516.7 < 2.2e-16 ***
#dff.S.bg.Kim.PCA$paxdb                1   64.516      5980     7452.2 9.574e-16 ***

write.table(dff.S.bg.Kim.PCA, file="output/HM_expression.csv", sep="\t",  quote=F, row.names=F)
