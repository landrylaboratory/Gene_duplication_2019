###################################################################################
#       Comparison of HM proportion among singletons and duplicates proteins      #
###################################################################################
rm(list=ls())

#to set the working directory
setwd("dir")

library(dplyr)
library(tidyr)
library(stringr)
library(jtools)

dff.S.bg.Kim.PCA <- read.table("without_low_qual/HM.data.csv", sep="\t", header=T)

#we find that about 32% of yeast proteins tested form HMs:
nrow(filter(dff.S.bg.Kim.PCA, HM.bg.kim.S.PCA==1)) / nrow(dff.S.bg.Kim.PCA)*100
#32.18

#keeping successive SSDs

dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para== "ssd-successiv"] <- "ssd"
ct_hom <- table(droplevels(dff.S.bg.Kim.PCA)$type_para, dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA)

#Are they differences of HM formation depending of the duplication
dbg <- data.frame(class=c("No HM S", "HM S", 
                          "No HM SSD", "HM SSD",
                          "No HM 2D", "HM 2D",
                          "No HM WGD", "HM WGD"),
                  count=c(ct_hom[1,1], ct_hom[1,2],
                          ct_hom[2,1],ct_hom[2,2],
                          ct_hom[3,1],ct_hom[3,2],
                          ct_hom[4,1],ct_hom[4,2]))
ct_hom <-t(ct_hom)
chisq.test(ct_hom) #p-value < 2.2e-16 --> yes

fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA, as.factor(dff.S.bg.Kim.PCA$type_para), workspace=2e8)
#p-value < 2.2e-16

# Differences of HM formation between singletons (S) and small scale duplicats (ssd)?
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"]), workspace=2e8)
#p-value < 2.2e-16 --> yes

# Differences of HM formation between singletons (S) and whole genome duplicats (wgd)?
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"]), workspace=2e8)
#p-value =  1.353e-05 --> yes

# Differences of HM formation between singletons (S) and paralogs which were duplicated by both duplication mechanims ssd and wgd?
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"]), workspace=2e8)
#p-value =  4.998e-07 --> yes


#Check for proportion of HM amond SSD, WGD, double duplication (2D) and singletons (S)
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
dbg_percent
#dupli percentHM n
#1 Total  32.22951 1
#2    2D  43.90244 5
#3     S  24.99008 2
#4   SSD  38.47664 3
#5   WGD  32.71889 4

###comparison of the expression between S, WGD and SSD
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
#p-value = 1.326e-13
dff.RNA.S <- filter(dff.RNA, type_para=="S")
wilcox.test(dff.RNA.wgd$log10_RNA_med_MTX,dff.RNA.S$log10_RNA_med_MTX)
#p-value < 2.2e-16
wilcox.test(dff.RNA.ssd$log10_RNA_med_MTX,dff.RNA.S$log10_RNA_med_MTX)
#p-value < 2.2e-16


###Does the expression of a gene impact the report of HM in BG
dff.S.bg.Kim.PCA <- left_join(dff.S.bg.Kim.PCA, RNAseq, by="orf")
dff.S.bg.Kim.HMreport <- filter(dff.S.bg.Kim.PCA, n.study>0)
cor.test(dff.S.bg.Kim.HMreport$med, dff.S.bg.Kim.HMreport$log10_RNA_med_MTX, method = "spearman")
#p-value = 0.0521 R=0.04650405 


###Same analysis with protein expression (using paxdb database)
paxdb <- read.table("data/PaxDB_4932-WHOLE_ORGANISM-integrated_AM.txt") %>% select(V2, V3)
colnames(paxdb) <- c("orf", "paxdb")

dff.S.bg.Kim.PCA <- left_join(dff.S.bg.Kim.PCA, paxdb, by="orf")
cor.test(dff.S.bg.Kim.PCA$med, dff.S.bg.Kim.PCA$paxdb, method = "spearman")
#p-value = 1.404e-09 R=0.1011999 

dff.S.bg.Kim.HMreport <- filter(dff.S.bg.Kim.PCA, n.study>0)
cor.test(dff.S.bg.Kim.HMreport$med, dff.S.bg.Kim.HMreport$paxdb, method = "spearman")
#p-value = 2.185e-07 R=0.1236839 

dff.paxdb.wgd <- filter(dff.S.bg.Kim.PCA, type_para=="wgd")
dff.paxdb.ssd <- filter(dff.S.bg.Kim.PCA, type_para=="ssd")
dff.paxdb.S <- filter(dff.S.bg.Kim.PCA, type_para=="S")
wilcox.test(dff.paxdb.ssd$paxdb, dff.paxdb.S$paxdb) #p-value < 2.2e-16
wilcox.test(dff.paxdb.S$paxdb, dff.paxdb.wgd$paxdb) #p-value < 2.2e-16
wilcox.test(dff.paxdb.ssd$paxdb, dff.paxdb.wgd$paxdb) #p-value =  0.00789

#If paralogs more expressed than singletons, which parameter impact the HM formation?
#The duplication or the expression --> glm (see Table S5)
dff.S.bg.Kim.PCA <- dff.S.bg.Kim.PCA %>% mutate(Duplication.binaire = ifelse(type_para=="S", 0, 1))


glm1<-glm(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA ~ 
            dff.S.bg.Kim.PCA$Duplication.binaire+
            dff.S.bg.Kim.PCA$log10_RNA_med_MTX, family=binomial)
summary(glm1)
export_summs(glm1)
anova(glm1,test="LRT") 


glm1<-glm(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA ~ 
            dff.S.bg.Kim.PCA$Duplication.binaire+
            dff.S.bg.Kim.PCA$paxdb, family=binomial)
summary(glm1)
export_summs(glm1)
anova(glm1,test="LRT") 

write.table(dff.S.bg.Kim.PCA, file="without_low_qual/HM_expression_phylom.csv", sep="\t",  quote=F, row.names=F)

#--> Analyse represented in Figure 2A and Figure S4