###########################################################################################
#                             Script to filter PCA data                               #
###########################################################################################
rm(list=ls())

require(gitter)
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


df <- read.table("output/PCA_data_2019_02.tab", sep="\t", header=T)

colstocharacter = c("Duplication", "Mat_A", "Mat_alpha", 
                    "replicat_plate", "Condition", 
                    "PCA", "taille")
df[,colstocharacter] = apply(df[,colstocharacter], 2, function(x) as.character(x))

colstonumeric = c("Plate", "col", "row", "bloc", "taille")
df[,colstonumeric] = apply(df[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

###########################################################################################
#                               Sort the dataframe                                        #
###########################################################################################
##Empty flag with Gitter correspond to a correct colony
df$flag = as.character(df$flag)
df$flag[is.na(df$flag)] <- "OK"
df = df %>% mutate(flag = ifelse(flag == "", "OK",flag))
df$flag = as.factor(df$flag)

#look at distribution sizes flag
ggplot(df, aes(x=log2(taille+1), colour=flag)) + geom_density()

#after looking at the images, flag C appears to be promiscuous and also concerns rather small colonies, 
#which means more missing 
#data for non interacting pairs so should not eliminate them
#remove only S, and S,C

df %<>% mutate(taille = ifelse(flag == "S", NA,taille))
df %<>% mutate(taille = ifelse(flag == "S,C", NA,taille))

#recheck
ggplot(df, aes(x=log2(taille+1), colour=flag)) + geom_density()

df %<>% mutate(adjust = taille+1)
df %<>% mutate(log2 = log2(adjust))


###########################################################################################
#                           Evaluation of distribution size                               #
###########################################################################################
#Boxplot in log2(size) in function of the different replicate plates
#pdf("distributuon_size.pdf")
#df %>% ggplot(aes(x=as.factor(Plate), y=log2, stat="identity")) + 
#  geom_jitter(shape=16, position=position_jitter(0.2)) +
#  geom_boxplot(notch=TRUE, alpha=0.5)
#dev.off()

df %>% ggplot(aes(x=as.factor(Condition), y=log2, file=Plate, stat="identity")) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(notch=TRUE, alpha=0.5)

#remove borders and empty (df = 107520 rows)
df.noB <- df %>% filter(row!=1, row!=2, row!=32,row!=31, col!=1, col!=2, col!=47, col!=48, Mat_A!="empty")
nrow(df.noB) #85440

pair_df.noB <- filter(df.noB, Mat_A!=Mat_alpha) %>%
  group_by(Mat_A, Mat_alpha, Plate, bloc,replicat_plate, Condition, PCA) %>%
  mutate(pair=apply(cbind(as.character(Mat_A), 
                          as.character(Mat_alpha)), 
                    1, function(x) paste(sort(x), collapse=".")))
pair_df.noB <- select(pair_df.noB , pair, Plate, bloc, replicat_plate, Condition, PCA) %>% unique(.) %>% as.data.frame()
pair_df.noB <- select(pair_df.noB , pair, Plate, bloc, replicat_plate, Condition, PCA) %>% unique(.) %>% as.data.frame()

df.noB <- left_join(df.noB, pair_df.noB, by=c("Plate", "bloc", "replicat_plate", "Condition", "PCA"))

#remove strains not good in PCA1 and PCA2
reconstructed <- read.table("data/coordonnees_good_strains_young_SSD3.csv", header=T)
reconstructed$pair = apply(cbind(as.character(reconstructed$P1), 
                                     as.character(reconstructed$P2)), 
                               1, function(x) paste(sort(x), collapse="."))
reconstructed %<>% select(., project, pair) %>% filter(., project=="reconstructed")
df.noB <- left_join(df.noB, reconstructed, by="pair")
df.noB %<>% mutate (delete = ifelse((project=="reconstructed" & PCA!="PCA3"), "yes", "no"))
df.noB %<>% filter(., (is.na(delete) | delete=="no")) %>% select(., -project, -delete)
nrow(df.noB) #nrow=78640


#Size distribution
colplot=c("#00A44A", "#28F040", "#0083B3", "#00F5FF", "#003B63", "#7D0A90", "#CB10EA",  "#C98EF6", 
          "#8C0011", "#FF011B", "#FF7605", "#FDEF33")

df.noB %>% filter(!is.na(log2)) %$% plot(density(log2[Condition=="S2"]), col=colplot[1], 
                                         main = "Size of colonies : density depending steps")
df.noB %>% filter(!is.na(log2)) %$% lines(density(log2[Condition=="MTX1"]), col=colplot[3])
df.noB %>% filter(!is.na(log2)) %$% lines(density(log2[Condition=="MTX2"]), col=colplot[8])
legend("topleft", c("S2", "MTX1", "MTX2"), 
       lty=c(1,1), lwd=c(2.0,2.0), 
       col=c(colplot[1], colplot[3],colplot[8]))


###########################################################################################
#     Sort MTX result in function of sources diploids plates and DMSO                     #
###########################################################################################

#PCA1 and PCA2
PCA1_2 <- filter(df.noB, PCA=="PCA1" | PCA=="PCA2")
#Check if there are colonies which grow in some S2 or S4 replic plate and no other
a <- PCA1_2 %>% filter(log2 > 8.5 & (Condition=="S2" | Condition=="S4")) %>% group_by(Plate, col, row) %>%
  summarise(grow = length(log2)) %>% as.data.frame()
b <- PCA1_2 %>% filter(log2 < 8.5 & (Condition=="S2" | Condition=="S4")) %>% group_by(Plate, col, row) %>%
  summarise(no_grow = length(log2)) %>% as.data.frame()
comp_replicS = full_join(a, b, by=c("Plate", "col", "row")) %>% group_by(Plate, col, row) %>%
  mutate(comp = if(is.na(grow) | is.na(no_grow)) "ok" else "no") %>% as.data.frame()
nrow(filter(comp_replicS, comp=="no")) #6 --> (6/3376)*100=0.18%


df.batch1 <- PCA1_2 %>% filter((Condition=="S2" | Condition =="MTX2") & (Plate==1 | Plate==2) & (replicat_plate <= 3))

#same replic number between S2 and MTX for batch 1
a = df.batch1 %>% filter(Condition=="S2") %>% select(log2, replicat_plate, row, col)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
b = df.batch1 %>% filter(Condition=="MTX2")
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0.23%
df.batch1 = inner_join(b, a, by=c("replicat_plate", "row", "col")) %>%
  rename(log2 = log2.x) %>% rename(log2S = log2.y) %>%
  mutate(log2 = ifelse(log2S < 8.5, NA, log2)) %>% as.data.frame()
(nrow(filter(df.batch1 , is.na(log2)))/nrow(df.batch1 ))*100 #0.11

filterS <- function(a,b)
{
  inner_join(b, a, by=c("PCA", "Plate", "row", "col")) %>%
    dplyr::rename(., log2 = log2.x) %>% dplyr::rename(., log2S = log2.y) %>%
    mutate(log2 = ifelse(log2S < 8.5, NA, log2)) %>% as.data.frame()
}

df.batch2 <- df.noB %>% filter((Plate==1 | Plate==2) &  (PCA=="PCA1" | PCA=="PCA2") &
                                 (Condition=="S4" | Condition =="DMSO2" | (Condition=="MTX2" & replicat_plate > 3)))

#For batch 2 and PCA2, more complexe : see figure
#DMSO.R1 --> from S4R2
a = df.batch2 %>% filter(Condition=="S4" & replicat_plate==2) %>% select(log2, PCA, Plate, row, col)
b = df.batch2 %>% filter(Condition=="DMSO2" & replicat_plate==1)
df.batch2.DMSO.R1 = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0
(nrow(filter(df.batch2.DMSO.R1 , is.na(log2)))/nrow(df.batch2.DMSO.R1))*100 #3.77%

#DMSO.R2--> from S4R3
a = df.batch2 %>% filter(Condition=="S4" & replicat_plate==3) %>% select(log2, PCA, Plate, row, col)
b = df.batch2 %>% filter(Condition=="DMSO2" & replicat_plate==2)
df.batch2.DMSO.R2 = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0
(nrow(filter(df.batch2.DMSO.R2 , is.na(log2)))/nrow(df.batch2.DMSO.R2 ))*100 #3.77%
#DMSO.R3 --> from S4R5
a = df.batch2 %>% filter(Condition=="S4" & replicat_plate==5) %>% select(log2, PCA, Plate, row, col)
b = df.batch2 %>% filter(Condition=="DMSO2" & replicat_plate==3)
df.batch2.DMSO.R3 = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0
(nrow(filter(df.batch2.DMSO.R3 , is.na(log2)))/nrow(df.batch2.DMSO.R3))*100 #3.81%
#MTX.R4 --> from S4R6
a = df.batch2 %>% filter(Condition=="S4" & replicat_plate==6) %>% select(log2, PCA, Plate, row, col)
b = df.batch2 %>% filter(Condition=="MTX2" & replicat_plate==4)
df.batch2.MTX.R4 = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0.18
(nrow(filter(df.batch2.MTX.R4, is.na(log2)))/nrow(df.batch2.MTX.R4))*100 #3.99%
#MTX.R5 --> from S4R8
a = df.batch2 %>% filter(Condition=="S4" & replicat_plate==8) %>% select(log2, PCA, Plate, row, col)
b = df.batch2 %>% filter(Condition=="MTX2" & replicat_plate==5)
df.batch2.MTX.R5 = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0.09
(nrow(filter(df.batch2.MTX.R5, is.na(log2)))/nrow(df.batch2.MTX.R5))*100 #3.90%
#MTX.R6 --> from S4R9
a = df.batch2 %>% filter(Condition=="S4" & replicat_plate==9) %>% select(log2, PCA, Plate, row, col)
b = df.batch2 %>% filter(Condition=="MTX2" & replicat_plate==6)
df.batch2.MTX.R6 = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0.04
(nrow(filter(df.batch2.MTX.R6, is.na(log2)))/nrow(df.batch2.MTX.R6))*100 #3.77%

df.PCA2 <- df.noB %>% filter(Plate==3)
#DMSO --> from S2R2
a = df.PCA2 %>% filter(Condition=="S2" & replicat_plate==2) %>% select(log2, PCA, Plate, row, col)
b = df.PCA2 %>% filter(Condition=="DMSO2")
df.PCA2.DMSO = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0
(nrow(filter(df.PCA2.DMSO, is.na(log2)))/nrow(df.PCA2.DMSO))*100 #0%

#MTX --> from S2R3
a = df.PCA2 %>% filter(Condition=="S2" & replicat_plate==3) %>% select(log2, PCA, Plate, row, col)
b = df.PCA2 %>% filter(Condition=="MTX2")
df.PCA2.MTX = filterS(a,b)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0.09
(nrow(filter(df.PCA2.MTX, is.na(log2)))/nrow(df.PCA2.MTX))*100 #0.09%


# PCA3 
PCA3 <- filter(df.noB, PCA=="PCA3")
#same number of replicates between S2 and MTX for batch3
a = PCA3 %>% filter(Condition=="S2") %>% select(log2, row, col, Plate)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
b = PCA3 %>% filter(Condition=="MTX2")
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0
df.batch3 = inner_join(b, a, by=c("row", "col", "Plate")) %>%
  dplyr::rename(., log2 = log2.x) %>% dplyr::rename(., log2S = log2.y) %>%
  mutate(log2 = ifelse(log2S < 8.5, NA, log2)) %>% as.data.frame()
(nrow(filter(df.batch3 , is.na(log2)))/nrow(df.batch3 ))*100 #10.48

df.filtrS <- rbind(df.batch1, df.batch2.DMSO.R1, df.batch2.DMSO.R2, df.batch2.DMSO.R3,
                   df.batch2.MTX.R4, df.batch2.MTX.R5, df.batch2.MTX.R6,
                   df.PCA2.DMSO, df.PCA2.MTX, df.batch3)
(nrow(filter(df.filtrS, is.na(log2)))/nrow(df.filtrS))*100 #3.95%

#summarize per PPI
df.filtrS %<>% mutate(ppi= interaction(Mat_A,Mat_alpha))


###########################################################################################
#                                     Remove outliers                                     #
###########################################################################################

#I remove the na before calulating the total otherwise we
#a not the true number of replicats without NA
MTX <- filter(df.filtrS, Condition=="MTX2") %>% 
  select(Duplication, pair, ppi, Mat_A, Mat_alpha, Plate, bloc, col, row, log2)
(nrow(filter(MTX, is.na(log2)))/nrow(MTX))*100 #4.21

df.MTXnoOut = MTX  %>% 
  filter(!is.na(log2)) %>% group_by(ppi) %>%
  mutate(qt25 = quantile(log2, 0.25),
         qt75 = quantile(log2, 0.75),
         iqr = IQR(log2)) %>%
  ungroup()

#df.MTXnoOut %<>% filter(log2 < qt75 + iqr) %>% filter( log2 > qt25 - iqr)
#removes only small outliers
df.MTXnoOut %<>%  filter( log2 > qt25 - iqr)

#keep if at least 4 replicates after removing outliers
df.MTXnoOut %<>% group_by(ppi) %>% mutate (total = n()) %>% filter(total>=4)

df.med <- df.MTXnoOut %>% group_by(Duplication, pair, ppi, Mat_A, Mat_alpha, Plate) %>% 
  summarise(med_MTX=median(log2), replicats = n()) %>% as.data.frame()

df.med %<>% mutate(type = ifelse(Mat_A==Mat_alpha, "HM",
                          ifelse(Mat_A!=Mat_alpha, "HET", "PB")))
ggplot(df.med, aes(x=log2(med_MTX), colour=type)) + geom_density()

write.table(df.med, file="output/PCA_med_data_2019_02.tab", sep="\t",  quote=F, row.names=F)
