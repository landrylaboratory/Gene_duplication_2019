################################################################################################################
#            Script to evaluate age of Small Scall Duplication events with orthologs from phyloDB              #
###############################################################################################################

rm(list=ls())

setwd("/Users/axellemarchant/Documents/postdoc_Landry/AMarchant_2016-2019/papier_AMarchant_2019")

library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(circlize)
library(VennDiagram)
library(mixtools)
library(magrittr)
library(tidytree)
library(ape)

#We kept only SSDs duplicated one time, we removed successive duplications
SSD <- read.table("data/duplication_SDS_1paire.txt", sep="\t", header=F)
colnames(SSD) <- c("P1", "P2")
SSD$pair <- apply(cbind(as.character(SSD$P1), 
                             as.character(SSD$P2)), 
                       1, function(x) paste(sort(x), collapse="."))

####   Age estimation with more distant sp --> use of phylodb orthology###

####################################Phylomdb_Tidytree####################################
names_phylodb <- read.table("data/phylomDB/yeast_SN_genes_names.txt", sep="\t", header=F)
colnames(names_phylodb) <- c("phylodb_P1", "P1")
SSD <- left_join(SSD, names_phylodb, by="P1")
colnames(names_phylodb) <- c("phylodb_P2", "P2")
SSD <- left_join(SSD, names_phylodb, by="P2")
write.table(SSD, file="output/allSSD_phylodb_names_2019_02.txt", sep="\t",  quote=F, row.names=F)

Tree <- read.table("data/phylomDB/best_trees.txt", sep="\t", header=F)
colnames(Tree) <- c("phylodb_P1", "evolutionary_model_P1", "likelihood_value_P1", "tree_P1")
SSD_tree <- left_join(SSD, Tree, by="phylodb_P1")
colnames(Tree) <- c("phylodb_P2", "evolutionary_model_P2", "likelihood_value_P2", "tree_P2")
SSD_tree <- left_join(SSD_tree, Tree, by="phylodb_P2")

SSD_tree <- select(SSD_tree, P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2)
write.table(SSD_tree, file="output/allSSD_phylodb_tree_2019_02.txt", sep="\t",  quote=F, row.names=F)


older_orthP1 <- function(phylodb_P1, phylodb_P2, tree_P1)
{
  tree_P1 <- gsub(")0.[0-9]*:", "):", as.character(tree_P1))
  tree_P1 <- read.tree(text=tree_P1)
  
  x <- as_data_frame(as_tibble(tree_P1))
  node_P1 <- min(filter(x, label==phylodb_P1)$node)
  #some pairs without both paralogs in the same tree  
  if (nrow(filter(x, label==phylodb_P2)) > 0)
  {
    node_P2 <- min(filter(x, label==phylodb_P2)$node)
    commun_ancestors_P1P2 <-intersect(ancestor(x, node_P1)$parent, ancestor(x, node_P2)$parent)
    more_recent_commun_ancestors_P1P2 <- max(commun_ancestors_P1P2)
    
    #We want the more recent commun ancestor below 
    #the more recent commun ancestor between the 2 paralogs
    a <- ancestor(x, node_P1)$parent
    CA_P1 <- ifelse (length(a[a>more_recent_commun_ancestors_P1P2])>=1,
                     min(a[a > more_recent_commun_ancestors_P1P2]),
                     max(a))
  }
  
  if (nrow(filter(x, label==phylodb_P2)) == 0)
  {
    CA_P1 <- min(ancestor(x, node_P1)$parent)
  }
  
  #We next want the farest species of P1 in this undergroup
  Colnames <- c("parent", "node", "branch.length", "label")
  dfP1 = data.frame(matrix(vector(), 0, 4, dimnames=list(c(), Colnames)))
  head(dfP1)
  
  for (i in CA_P1:(parent(x, node_P1)$node))
  {
    dfP1 <- rbind(dfP1, child(x, i))
  }
  
  dfP1 <- filter(dfP1, !is.na(label))
  
  #We obtain the list of species in the undergroup
  #We want the more distant from P1 --> node the more distant
  #not work --> use of known tree and sp
  sp <- colsplit(dfP1$label, "[_]", names=c("ortholog_P1", "sp"))
  sp <- dplyr::select(sp, sp)
  dfP1 <- cbind(dfP1, sp)
  
  dfP1 <- dfP1 %>%
    mutate (group_phyloP1 = ifelse  (sp=="YEAST" | sp=="SACKU" | sp=="SACMI" | sp=="SACPA" |
                                       sp=="SACBA" | sp=="CANGA" | sp=="SACCA" | sp=="VANPO", 1,
                                     ifelse  (sp=="KLULA"| sp=="KLUWA"| sp=="SACKL"| sp=="ASHGO", 2,
                                              ifelse  (sp=="CANAL"| sp=="CANDU"| sp=="CANTR"| sp=="LODEL"| sp=="PICGU"| sp=="PICST"| sp=="DEBHA"| sp=="YARLI" | sp=="CLALS", 3,
                                                       ifelse  (sp=="ASPCL"| sp=="ASPFL"| sp=="ASPFU"| sp=="ASPNG"| sp=="ASPOR"| sp=="ASPTE"| sp=="EMENI"| sp=="NEOFI"| 
                                                                  sp=="COCIM"| sp=="UNCRE"| sp=="AJECA"| sp=="FUSOX"| sp=="GIBZE"| sp=="NECHA"| sp=="GIBMO"| sp=="TRIRE"|
                                                                  sp=="CHAGB"| sp=="PODAN"| sp=="NEUCR"| sp=="MAGGR"| sp=="SCLSC"| sp=="BOTFU"| sp=="PHANO"| sp=="383855", 4,
                                                                ifelse  (sp=="SCHJP"| sp=="SCHPO"| sp=="PNECA", 5,
                                                                         ifelse  (sp=="PHACH"| sp=="POSPM"| sp=="COPCI"| sp=="LACBI"| sp=="CRYNE"| sp=="PUCGR"| sp=="365493"| sp=="USTMA", 6, 
                                                                                  ifelse  (sp=="RHIOR"| sp=="PHYBL"| sp=="BATDE"| sp=="ENCCU", 7,
                                                                                           ifelse  (sp=="HUMAN", 8,
                                                                                                    ifelse  (sp=="ARATH", 9,
                                                                                                             ifelse  (is.na(sp), NA, "PB"))))))))))) %>% as.data.frame()
  
  dfP1$group_phyloP1=as.numeric(as.character(dfP1$group_phyloP1))
  age_dupliP1<-max(dfP1$group_phyloP1)
  age_dupliP1
}

older_orthP2 <- function(phylodb_P1, phylodb_P2, tree_P2)
{
  tree_P2 <- gsub(")0.[0-9]*:", "):", as.character(tree_P2))
  tree_P2 <- read.tree(text=tree_P2)
  
  x <- as_data_frame(as_tibble(tree_P2))
  node_P2 <- min(filter(x, label==phylodb_P2)$node)
  
  #some pairs without both paralogs in the same tree  
  if (nrow(filter(x, label==phylodb_P1)) > 0)
  {
    node_P1 <- min(filter(x, label==phylodb_P1)$node)
    commun_ancestors_P2P1 <-intersect(ancestor(x, node_P2)$parent, ancestor(x, node_P1)$parent)
    more_recent_commun_ancestors_P2P1 <- max(commun_ancestors_P2P1)
    
    #We want the more recent commun ancestor below 
    #the more recent commun ancestor between the 2 paralogs
    a <- ancestor(x, node_P2)$parent
    CA_P2 <- ifelse (length(a[a>more_recent_commun_ancestors_P2P1])>=1,
                     min(a[a > more_recent_commun_ancestors_P2P1]),
                     max(a))
  }
  
  if (nrow(filter(x, label==phylodb_P1)) == 0)
  {
    CA_P2 <- min(ancestor(x, node_P2)$parent)
  }
  
  
  #We next want the farest species of P2 in this undergroup
  Colnames <- c("parent", "node", "branch.length", "label")
  dfP2 = data.frame(matrix(vector(), 0, 4, dimnames=list(c(), Colnames)))
  head(dfP2)
  
  for (i in CA_P2:(parent(x, node_P2)$node))
  {
    dfP2 <- rbind(dfP2, child(x, i))
  }
  
  dfP2 <- filter(dfP2, !is.na(label))
  
  #We obtain the list of species in the undergroup
  #We want the more distant from P2 --> node the more distant
  #not work --> use of known tree and sp
  sp <- colsplit(dfP2$label, "[_]", names=c("ortholog_P2", "sp"))
  sp <- dplyr::select(sp, sp)
  dfP2 <- cbind(dfP2, sp)
  
  
  dfP2 <- dfP2 %>%
    mutate (group_phyloP2 = ifelse  (sp=="YEAST" | sp=="SACKU" | sp=="SACMI" | sp=="SACPA" |
                                       sp=="SACBA" | sp=="CANGA" | sp=="SACCA" | sp=="VANPO", 1,
                                     ifelse  (sp=="KLULA"| sp=="KLUWA"| sp=="SACKL"| sp=="ASHGO", 2,
                                              ifelse  (sp=="CANAL"| sp=="CANDU"| sp=="CANTR"| sp=="LODEL"| sp=="PICGU"| sp=="PICST"| sp=="DEBHA"| sp=="YARLI" | sp=="CLALS", 3,
                                                       ifelse  (sp=="ASPCL"| sp=="ASPFL"| sp=="ASPFU"| sp=="ASPNG"| sp=="ASPOR"| sp=="ASPTE"| sp=="EMENI"| sp=="NEOFI"| 
                                                                  sp=="COCIM"| sp=="UNCRE"| sp=="AJECA"| sp=="FUSOX"| sp=="GIBZE"| sp=="NECHA"| sp=="GIBMO"| sp=="TRIRE"|
                                                                  sp=="CHAGB"| sp=="PODAN"| sp=="NEUCR"| sp=="MAGGR"| sp=="SCLSC"| sp=="BOTFU"| sp=="PHANO"| sp=="383855", 4,
                                                                ifelse  (sp=="SCHJP"| sp=="SCHPO"| sp=="PNECA", 5,
                                                                         ifelse  (sp=="PHACH"| sp=="POSPM"| sp=="COPCI"| sp=="LACBI"| sp=="CRYNE"| sp=="PUCGR"| sp=="365493"| sp=="USTMA", 6, 
                                                                                  ifelse  (sp=="RHIOR"| sp=="PHYBL"| sp=="BATDE"| sp=="ENCCU", 7,
                                                                                           ifelse  (sp=="HUMAN", 8,
                                                                                                    ifelse  (sp=="ARATH", 9,
                                                                                                             ifelse  (is.na(sp), NA, "PB"))))))))))) %>% as.data.frame()
  
  dfP2$group_phyloP2=as.numeric(as.character(dfP2$group_phyloP2))
  age_dupliP2<-max(dfP2$group_phyloP2)
  age_dupliP2
}


SSD_tree <- dplyr::select(SSD_tree, P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2)
SSD_tree <- filter(SSD_tree, !is.na(tree_P1) & !is.na(tree_P2))
SSD_age_phylomedb <- SSD_tree %>% group_by(P1, P2) %>% summarise(age_dupliP1 = older_orthP1(phylodb_P1, phylodb_P2, tree_P1),
                                                                 age_dupliP2 = older_orthP2(phylodb_P1, phylodb_P2, tree_P2))

SSD_age_phylomedb %<>% group_by(P1, P2) %>% mutate(age_dupli = max(age_dupliP1, age_dupliP2)) %>% as.data.frame()

write.table(SSD_age_phylomedb, "output/allSSD_age_phylomedb_2019_02_AM.txt", sep="\t",  quote=F, row.names=F)

#Cluster all PCA SSD datas


##########Test of the script with WGD for which group ages should be 
########## 1 (genome duplication origin) or 2 (hybridization origin)
WGD <- read.table("data/WGD.csv", header = F, sep = ";")
colnames(WGD) <- c("P1", "P2")

names_phylodb <- read.table("data/phylomDB/yeast_SN_genes_names.txt", sep="\t", header=F)
colnames(names_phylodb) <- c("phylodb_P1", "P1")

WGD <- left_join(WGD, names_phylodb, by="P1")
colnames(names_phylodb) <- c("phylodb_P2", "P2")
WGD <- left_join(WGD, names_phylodb, by="P2")
write.table(WGD, file="output/WGD_phylodb_names.txt", sep="\t",  quote=F, row.names=F)

Tree <- read.table("data/phylomDB/best_trees.txt", sep="\t", header=F)
colnames(Tree) <- c("phylodb_P1", "evolutionary_model_P1", "likelihood_value_P1", "tree_P1")
WGD_tree <- left_join(WGD, Tree, by="phylodb_P1")
colnames(Tree) <- c("phylodb_P2", "evolutionary_model_P2", "likelihood_value_P2", "tree_P2")
WGD_tree <- left_join(WGD_tree, Tree, by="phylodb_P2")

WGD_tree <- dplyr::select(WGD_tree, P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2)
write.table(WGD_tree, file="output/WGD_phylodb_tree.txt", sep="\t",  quote=F, row.names=F)

WGD_tree <- filter(WGD_tree, !is.na(tree_P1) & !is.na(tree_P2))
WGD_age_phylomedb <- WGD_tree %>% group_by(P1, P2) %>% summarise(age_dupliP1 = older_orthP1(phylodb_P1, phylodb_P2, tree_P1),
                                                                 age_dupliP2 = older_orthP2(phylodb_P1, phylodb_P2, tree_P2))
WGD_age_phylomedb %<>% group_by(P1, P2) %>% mutate(age_dupli = max(age_dupliP1, age_dupliP2)) %>% as.data.frame()


pb_def_age <- filter(WGD_age_phylomedb, age_dupli > 2)
(nrow(pb_def_age)/nrow(WGD_age_phylomedb))*100 # --> 12.04% d'erreur

#idea: to see if in the age group defined, 
# we observe species with at least both paralogs
# -> 1) estimate of the duplication group 2) verification
WGD_age <- dplyr::select(WGD_age_phylomedb, P1, P2, age_dupli)
WGD_tree_age <- left_join(WGD_tree, WGD_age, by=c("P1", "P2"))


n_orth <- function(phylodb_P1, phylodb_P2, tree_P1, tree_P2, age_dupli)
{
  tree_P1 <- gsub(")0.[0-9]*:", "):", as.character(tree_P1))
  tree_P1 <- read.tree(text=tree_P1)
  
  y <- as_data_frame(as_tibble(tree_P1))
  
  tree_P2 <- gsub(")0.[0-9]*:", "):", as.character(tree_P2))
  tree_P2 <- read.tree(text=tree_P2)
  
  z <- as_data_frame(as_tibble(tree_P2))
  
  x <- rbind(y, z)
  x <- filter(x, !is.na(label))
  
  #We obtain the list of species in the undergroup
  #We want the more distant from P1 --> node the more distant
  #not work --> use of known tree and sp
  sp <- colsplit(x$label, "[_]", names=c("ortholog", "sp"))
  sp <- dplyr::select(sp, sp)
  x <- cbind(x, sp)
  
  n_orth_per_sp <- x %>% group_by(sp) %>% summarise(n=length(unique(label)))
  
  
  max_n_orth_for_grp <- ifelse  (age_dupli==1, max(n_orth_per_sp$n[n_orth_per_sp$sp=="YEAST"], n_orth_per_sp$n[n_orth_per_sp$sp=="SACKU"], 
                                                   n_orth_per_sp$n[n_orth_per_sp$sp=="SACMI"], n_orth_per_sp$n[n_orth_per_sp$sp=="SACPA"],
                                                   n_orth_per_sp$n[n_orth_per_sp$sp=="SACBA"], n_orth_per_sp$n[n_orth_per_sp$sp=="CANGA"], 
                                                   n_orth_per_sp$n[n_orth_per_sp$sp=="SACCA"], n_orth_per_sp$n[n_orth_per_sp$sp=="VANPO"]), 
                                 ifelse  (age_dupli==2, max(n_orth_per_sp$n[n_orth_per_sp$sp=="KLULA"], n_orth_per_sp$n[n_orth_per_sp$sp=="KLUWA"], 
                                                            n_orth_per_sp$n[n_orth_per_sp$sp=="SACKL"], n_orth_per_sp$n[n_orth_per_sp$sp=="ASHGO"]),
                                          ifelse  (age_dupli==3, max(n_orth_per_sp$n[n_orth_per_sp$sp=="CANAL"], n_orth_per_sp$n[n_orth_per_sp$sp=="CANDU"], 
                                                                     n_orth_per_sp$n[n_orth_per_sp$sp=="CANTR"], n_orth_per_sp$n[n_orth_per_sp$sp=="LODEL"], 
                                                                     n_orth_per_sp$n[n_orth_per_sp$sp=="PICGU"], n_orth_per_sp$n[n_orth_per_sp$sp=="PICST"], 
                                                                     n_orth_per_sp$n[n_orth_per_sp$sp=="DEBHA"], n_orth_per_sp$n[n_orth_per_sp$sp=="YARLI"] , 
                                                                     n_orth_per_sp$n[n_orth_per_sp$sp=="CLALS"]), 
                                                   ifelse  (age_dupli==4, max(n_orth_per_sp$n[n_orth_per_sp$sp=="ASPCL"], n_orth_per_sp$n[n_orth_per_sp$sp=="ASPFL"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="ASPFU"], n_orth_per_sp$n[n_orth_per_sp$sp=="ASPNG"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="ASPOR"], n_orth_per_sp$n[n_orth_per_sp$sp=="ASPTE"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="EMENI"], n_orth_per_sp$n[n_orth_per_sp$sp=="NEOFI"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="COCIM"], n_orth_per_sp$n[n_orth_per_sp$sp=="UNCRE"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="AJECA"], n_orth_per_sp$n[n_orth_per_sp$sp=="FUSOX"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="GIBZE"], n_orth_per_sp$n[n_orth_per_sp$sp=="NECHA"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="GIBMO"], n_orth_per_sp$n[n_orth_per_sp$sp=="TRIRE"],
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="CHAGB"], n_orth_per_sp$n[n_orth_per_sp$sp=="PODAN"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="NEUCR"], n_orth_per_sp$n[n_orth_per_sp$sp=="MAGGR"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="SCLSC"], n_orth_per_sp$n[n_orth_per_sp$sp=="BOTFU"], 
                                                                              n_orth_per_sp$n[n_orth_per_sp$sp=="PHANO"], n_orth_per_sp$n[n_orth_per_sp$sp=="383855"]),
                                                            ifelse  (age_dupli==5, max(n_orth_per_sp$n[n_orth_per_sp$sp=="SCHJP"], n_orth_per_sp$n[n_orth_per_sp$sp=="SCHPO"], 
                                                                                       n_orth_per_sp$n[n_orth_per_sp$sp=="PNECA"]),
                                                            ifelse  (age_dupli==6, max(n_orth_per_sp$n[n_orth_per_sp$sp=="PHACH"], n_orth_per_sp$n[n_orth_per_sp$sp=="POSPM"], 
                                                                                       n_orth_per_sp$n[n_orth_per_sp$sp=="COPCI"], n_orth_per_sp$n[n_orth_per_sp$sp=="LACBI"], 
                                                                                       n_orth_per_sp$n[n_orth_per_sp$sp=="CRYNE"], n_orth_per_sp$n[n_orth_per_sp$sp=="PUCGR"], 
                                                                                       n_orth_per_sp$n[n_orth_per_sp$sp=="365493"], n_orth_per_sp$n[n_orth_per_sp$sp=="USTMA"]),
                                                            ifelse  (age_dupli==7, max(n_orth_per_sp$n[n_orth_per_sp$sp=="RHIOR"], n_orth_per_sp$n[n_orth_per_sp$sp=="PHYBL"], 
                                                                                       n_orth_per_sp$n[n_orth_per_sp$sp=="BATDE"], n_orth_per_sp$n[n_orth_per_sp$sp=="ENCCU"]), 
                                                            ifelse  (age_dupli==8, n_orth_per_sp$n[n_orth_per_sp$sp=="HUMAN"],
                                                            ifelse  (age_dupli==9, n_orth_per_sp$n[n_orth_per_sp$sp=="ARATH"],
                                                           ifelse  (is.na(age_dupli), NA, "PB")))))))))) %>% as.data.frame()
  max_n_orth_for_grp
}

WGD_n_orth <- WGD_tree_age %>% group_by(P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2, age_dupli) %>% 
  summarise(max_n_orth_for_grp = n_orth(phylodb_P1, phylodb_P2, tree_P1, tree_P2, age_dupli)$.)

#We loop back until we have at least 2 duplicates for one species of the group or we arrive at age group 0
#Maximum 10 bearings

i=0
WGD_n_orth %<>% group_by(P1, P2) %>% mutate(age_dupli2= age_dupli) %>% as.data.frame()
WGD_n_orth$max_n_orth_for_grp <- as.numeric(as.character(WGD_n_orth$max_n_orth_for_grp))

while (i < 9)
{
  WGD_n_orth %<>% group_by(P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2, age_dupli) %>% 
                  summarise(age_dupli2= ifelse((max_n_orth_for_grp == 1 & age_dupli2 > 0), age_dupli2-1, age_dupli2))
  WGD_n_orth %<>% group_by(P1, P2) %>% 
    mutate(max_n_orth_for_grp = n_orth(phylodb_P1, phylodb_P2, tree_P1, tree_P2, age_dupli2)$.) %>% as.data.frame()
  i = i+1
} 


WGD_n_orth_simpl <- dplyr::select(WGD_n_orth, P1, P2, age_dupli, age_dupli2)
#Adding missing grouping variables: `phylodb_P1`, `tree_P1`, `phylodb_P2`, `tree_P2`
pb_def_age <- filter(WGD_n_orth_simpl, age_dupli2 > 2)
(nrow(pb_def_age)/nrow(filter(WGD_n_orth_simpl,!is.na(age_dupli2))))*10 #only 0.88% of error

###We do the correction for SSDs
SSD_tree<- read.table("output/allSSD_phylodb_tree_2019_02.txt", sep="\t", header=T)
SSD_age_phylomedb <- read.table("output/allSSD_age_phylomedb_2019_02_AM.txt", sep="\t", header=T)

SSD_age <- dplyr::select(SSD_age_phylomedb, P1, P2, age_dupli)
SSD_tree_age <- left_join(SSD_tree, SSD_age, by=c("P1", "P2"))
SSD_tree_age$age_dupli <- as.numeric(as.character(SSD_tree_age$age_dupli))
SSD_tree_age <- filter(SSD_tree_age, !is.na(age_dupli))
SSD_tree_age$phylodb_P1 <- as.character(SSD_tree_age$phylodb_P1)
SSD_tree_age$phylodb_P2 <- as.character(SSD_tree_age$phylodb_P2)

#1 case with the tree absent for P2, we removed it
SSD_tree_age <- filter(SSD_tree_age, !is.na(tree_P2))

SSD_n_orth <- SSD_tree_age %>% group_by(P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2, age_dupli) %>% 
  summarise(max_n_orth_for_grp = n_orth(phylodb_P1, phylodb_P2, tree_P1, tree_P2, age_dupli))


i=0
SSD_n_orth %<>% group_by(P1, P2) %>% mutate(age_dupli2= age_dupli) %>% as.data.frame()
SSD_n_orth$max_n_orth_for_grp <- as.numeric(as.character(SSD_n_orth$max_n_orth_for_grp))

while (i < 9)
{
  SSD_n_orth %<>% group_by(P1, phylodb_P1, tree_P1, P2, phylodb_P2, tree_P2, age_dupli) %>% 
    summarise(age_dupli2= ifelse((max_n_orth_for_grp == 1 & age_dupli2 > 0), age_dupli2-1, age_dupli2))
  SSD_n_orth %<>% group_by(P1, P2) %>% 
    mutate(max_n_orth_for_grp = n_orth(phylodb_P1, phylodb_P2, tree_P1, tree_P2, age_dupli2)$.) %>% as.data.frame()
  i = i+1
} 


SSD_n_orth_simpl <- select(SSD_n_orth, P1, P2, age_dupli, age_dupli2, max_n_orth_for_grp)
nrow(filter(SSD_n_orth_simpl,age_dupli!=age_dupli2))
#72 cases

SSD_n_orth_simpl $pair <- apply(cbind(as.character(SSD_n_orth_simpl $P1), 
                        as.character(SSD_n_orth_simpl $P2)), 
                  1, function(x) paste(sort(x), collapse="."))

write.table(SSD_n_orth_simpl ,file = "output/SSD_corrected_age_2019_02_AM.csv", sep="\t",  quote=F, row.names=F)
