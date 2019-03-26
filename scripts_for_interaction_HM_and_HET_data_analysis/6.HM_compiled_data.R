####################################################################################################
#     Compile data for proportion HM study among singletons, WGD, SSD and double duplication      #
####################################################################################################
require(gitter)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggsignif)
library(grid)

###################################
##### All S. cerevisiae genes ##### 
###################################
genes <- read.table("data/genes-liste.csv", header = F)
colnames(genes) <- 'orf'


###################################
#####         BIOGRID         ##### 
###################################
bg <- read.table("data/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.166.tab2.txt", 
                 header = T,fill=T, quote="", sep="\t",stringsAsFactors=FALSE, comment.char = "")
#consider only physical interactions and high throughput
bg.p = bg %>% filter(Experimental.System.Type =="physical")#, Throughput=="High Throughput")

#look at the techniques to make sure we do not use things that do not reflect ppis
tech <- unique(bg.p$Experimental.System)

#remove techniques that are not ppis
bg.p <- filter(bg.p, Experimental.System!="Co-fractionation" & 
                 Experimental.System!="Affinity Capture-RNA" & 
                 Experimental.System!="Protein-RNA" & 
                 Experimental.System!="Co-purification" &
                 Experimental.System!="Proximity Label-MS" & 
                 Experimental.System!="Co-localization")

bg.p <- select(bg.p, Systematic.Name.Interactor.A,Systematic.Name.Interactor.B,
               Throughput, Experimental.System,Pubmed.ID)

#SELECT ONLY STUDIES WITH HM TESTED
bg.p = bg.p %>%  mutate(HM = ifelse(Systematic.Name.Interactor.A==Systematic.Name.Interactor.B, 1, 0))
Pubmed.ID <- bg.p %>% group_by(Pubmed.ID) %>% mutate(HM.in.study = ifelse(sum(HM) >= 1, 1, 0))
Pubmed.ID.with.HM <- filter(Pubmed.ID, HM.in.study==1)

#SELECT ONLY STUDIES WITH both prey and bait for a same orf
Pubmed.ID.bait <- unique(select(Pubmed.ID.with.HM, Pubmed.ID, Experimental.System, Systematic.Name.Interactor.A))
colnames(Pubmed.ID.bait) <- c("Pubmed.ID", "Experimental.System", "orf")
Pubmed.ID.prey <- unique(select(Pubmed.ID.with.HM, Pubmed.ID, Experimental.System, Systematic.Name.Interactor.B))
colnames(Pubmed.ID.prey) <- c("Pubmed.ID", "Experimental.System", "orf")

Pubmed.ID.bait$bait = 1
Pubmed.ID.prey$prey =1  

Pubmed.ID.bait.prey <- full_join(Pubmed.ID.bait, Pubmed.ID.prey, by=c("Pubmed.ID", "Experimental.System", "orf"))
Pubmed.ID.bait.prey <- filter(Pubmed.ID.bait.prey, !is.na(bait) & !is.na(prey))

HM <- filter(bg.p, Systematic.Name.Interactor.A==Systematic.Name.Interactor.B)
HM <- select(HM, Pubmed.ID, Experimental.System, Systematic.Name.Interactor.A, HM)
colnames(HM) <- c("Pubmed.ID", "Experimental.System", "orf", "HM")

Pubmed.ID.bait.prey.HM <- full_join(Pubmed.ID.bait.prey, HM, by=c("Pubmed.ID", "Experimental.System", "orf"))
Pubmed.ID.bait.prey.HM$HM[is.na(Pubmed.ID.bait.prey.HM$HM)] <- 0
Pubmed.ID.bait.prey.HM <- Pubmed.ID.bait.prey.HM %>% mutate(ExpRef = paste(Experimental.System,"(",Pubmed.ID, ")", sep="")) 
Pubmed.ID.bait.prey.HM.2 <- Pubmed.ID.bait.prey.HM %>% group_by(orf) %>% 
  summarise(comp_list = paste(sort(unique(c(as.character(ExpRef)))), collapse=";"),
            n.study = n(),
            n.HM = sum(HM))


#######################################
##### HM data from Kim et al 2019 ##### 
#######################################
Kim <- read.table("data/orf_tested_by_Kim_et_al_2019.txt", header = T, sep="\t")
Kim.falsePos <- read.table("data/fasle_positiv_Kim_et_al_2019.txt", header = T, sep="\t")
Kim.falsePos <- select(Kim.falsePos, ORF)
Kim.falsePos$false <- 'yes'
#remove potential false positive
Kim <- left_join(Kim, Kim.falsePos, by="ORF")
Kim <- filter(Kim, is.na(false))

HM.Kim <- read.table("data/HM_Kim_et_al_2019.txt", header = T, sep="\t")
HM.Kim <- select(HM.Kim, ORF)
HM.Kim$HM.Kim <- 1
Kim <- left_join(Kim, HM.Kim, by="ORF")
Kim$HM.Kim[is.na(Kim$HM.Kim)] <- 0
Kim <- select(Kim, ORF, HM.Kim)
colnames(Kim) <- c("orf", "HM.Kim")

Pubmed.ID.bait.prey.HM.Kim <- full_join(Pubmed.ID.bait.prey.HM.2, Kim, by="orf")

#annotate duplication statue
WGD <- read.table("data/WGD.csv", header = F, sep=";")
WGD$Duplication = "wgd"
SSD <- read.table("data/duplication_SDS_1paire.txt", header = F, sep="\t")
SSD$Duplication = "ssd"

P1.WGD <- select(WGD, V1, Duplication) 
colnames(P1.WGD) <- c("orf", "Duplication")
P2.WGD <- select(WGD, V2, Duplication)
colnames(P2.WGD) <- c("orf", "Duplication")

P1.SSD <- select(SSD, V1, Duplication) 
colnames(P1.SSD) <- c("orf", "Duplication")
P2.SSD <- select(SSD, V2, Duplication)
colnames(P2.SSD) <- c("orf", "Duplication")

df.2 <- rbind(P1.WGD, P2.WGD, P1.SSD, P2.SSD)
df.2 <- unique(df.2)
df.2 %<>% group_by(orf) %>% summarise(type_para = ifelse(n()==1, Duplication,
                                                         ifelse(n()>1, "ssd_wgd", 'pb')))

df.2 <- unique(df.2)

#add successive duplications
all_dupli <- read.table("data/all_duplicated_genes.txt", header = F)
all_dupli$dupli <- 'unknown'

P1.all_dupli <- select(all_dupli, V1, dupli) 
colnames(P1.all_dupli) <- c("orf", "dupli")
P2.all_dupli <- select(all_dupli, V2, dupli)
colnames(P2.all_dupli) <- c("orf", "dupli")
all_dupli2 <- unique(rbind(P1.all_dupli, P2.all_dupli))

all_dupli2 <- full_join(all_dupli2, df.2, by="orf")
all_dupli2 <- all_dupli2 %>% group_by(orf) %>% summarise(type_para= ifelse(is.na(type_para), "ssd-successiv", type_para))


Pubmed.ID.bait.prey.HM.Kim <- left_join(Pubmed.ID.bait.prey.HM.Kim, all_dupli2, by="orf")
Pubmed.ID.bait.prey.HM.Kim$type_para[is.na(Pubmed.ID.bait.prey.HM.Kim$type_para)] <- "S"
Pubmed.ID.bait.prey.HM.Kim <- Pubmed.ID.bait.prey.HM.Kim %>% group_by(orf) %>% 
  mutate(HM.bg.kim = ifelse(is.na(n.HM), HM.Kim,
                            ifelse(is.na(HM.Kim) & n.HM > 0, 1, 
                            ifelse(is.na(HM.Kim) & n.HM == 0, 0, 
                            ifelse(!is.na(n.HM) & !is.na(HM.Kim) & (n.HM + HM.Kim) > 0, 1, 0))))) %>% as.data.frame()


############################################
##### HM data from Stynen et al., 2008 ##### 
############################################

dff <- read.table("data/MedianValuesHomoMichnick_dupli_statue_2019_02_22_AM.tab", sep="\t", header=T)
#no any successiv ssd information in dff
dff <- select(dff, -type_para)
dff <- left_join(dff, all_dupli2, by='orf')
dff$type_para[is.na(dff$type_para)] <- 'S'

#remove proteins interacting with everythings that Tarassov et al., 2008
prot_remove <- read.table("data/Prot_removed.csv", header=T)
colnames(prot_remove) <- "orf"
prot_remove$remove_MatA <- "remove"
dff <- left_join(dff, prot_remove, by="orf")
dff <- filter(dff, is.na(remove_MatA))

dff.S.bg.Kim <- full_join(dff, Pubmed.ID.bait.prey.HM.Kim, by=c("orf", 'type_para'))
dff.S.bg.Kim <- dff.S.bg.Kim %>% group_by(orf) %>% mutate(HM.bg.kim.S = ifelse(is.na(interaction), HM.bg.kim,
                                                                               ifelse(is.na(HM.bg.kim), interaction,
                                                                               ifelse(!is.na(HM.bg.kim) & !is.na(interaction) & (HM.bg.kim + interaction) > 0, 1, 0)))) 

#adding PCA results
summary_table <- read.table("output/PCA_completed_by_Biog_Sty_2019_02.csv", sep="\t", header = T)
HM1 <- select(summary_table, P1, HM_P1)
HM2 <- select(summary_table, P2, HM_P2)
colnames(HM1) <- c("orf", "HM.PCA")
colnames(HM2) <- c("orf", "HM.PCA")
HM.PCA <- rbind(HM1, HM2)
HM.PCA <- left_join(HM.PCA, all_dupli2, by='orf')
HM.PCA$type_para[is.na(HM.PCA$type_para)] <- 'S'

dff.S.bg.Kim.PCA <- full_join(dff.S.bg.Kim, HM.PCA, by=c("orf", "type_para"))
dff.S.bg.Kim.PCA <- dff.S.bg.Kim.PCA %>% group_by(orf) %>% mutate(HM.bg.kim.S.PCA = ifelse(is.na(HM.bg.kim.S), HM.PCA,
                                                                                    ifelse(is.na(HM.PCA), HM.bg.kim.S,
                                                                                    ifelse(!is.na(HM.PCA) & !is.na(HM.bg.kim.S) & (HM.PCA + HM.bg.kim.S) > 0, 1, 0)))) 

write.table(dff.S.bg.Kim.PCA, file = "output/HM.data.csv", sep="\t",  quote=F, row.names=F)

