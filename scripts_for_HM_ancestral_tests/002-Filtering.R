# This script sorts and filters the data
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
library(gdata)

df <- read.table("df_2018_12_03_DIAS.tab", sep="\t", header=T,stringsAsFactors = F)

###########################################################################################
#                               Sort the dataframe                                        #
###########################################################################################

df$flag = as.character(df$flag)
df$flag[is.na(df$flag)] <- "OK"
df = df %>% mutate(flag = ifelse(flag == "", "OK",flag))
df$flag = as.factor(df$flag)

#look at distribution sizes flag
ggplot(df, aes(x=log2(taille), colour=flag)) + geom_density()

#after looking at the images, flag C appears to be promiscuous and also concerns rather small colonies, 
#which means more missing 
#data for non interacting pairs so should not eliminate them
#remove only S, and S,C

df %<>% mutate(taille = ifelse(flag == "S", NA,taille))
df %<>% mutate(taille = ifelse(flag == "S,C", NA,taille))

#recheck
ggplot(df, aes(x=log2(taille), colour=flag)) + geom_density()

df %<>% mutate(adjust = taille+1)
df %<>% mutate(log2 = log2(adjust))


###########################################################################################
#                           Evaluation of distribution size                               #
###########################################################################################


pdf("../figures/distribution_size.pdf")
df %>% ggplot(aes(x=as.factor(Etape), y=log2, file=Plate, stat="identity")) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(notch=TRUE, alpha=0.5)
dev.off()

#remove borders and empty
df.noB <- df %>% filter(row!=1, row!=2, row!=32,row!=31,  col!=1, col!=2, col!=47, col!=48, Mat_A!="empty")

#Size distribution
colplot=c("#00A44A", "#28F040", "#0083B3", "#00F5FF", "#003B63", "#7D0A90", "#CB10EA",  "#C98EF6", 
          "#8C0011", "#FF011B", "#FF7605", "#FDEF33")

pdf("../figures/size_density.pdf")
df.noB %>% filter(!is.na(log2)) %$% plot(density(log2[Etape=="MTX1"]), col=colplot[1], 
                                         main = "Size of colonies : density depending steps")
df.noB %>% filter(!is.na(log2)) %$% lines(density(log2[Etape=="MTX2"]), col=colplot[8])

legend("topleft", c("mtx1", "mtx2"), 
       lty=c(1,1), lwd=c(2.0,2.0), 
       col=c(colplot[1], colplot[8]))
dev.off()


###########################################################################################
#                                  filter data                                            #
###########################################################################################

#condense both steps in a single data.frame
a = df.noB %>% filter(Etape=="S2") %>% select(row, col, Plate,log2)
(nrow(filter(a, is.na(log2)))/nrow(a))*100 #0
b = df.noB %>% filter(Etape=="MTX2")
(nrow(filter(b, is.na(log2)))/nrow(b))*100 #0.011
df.batch1 = inner_join(b, a, by=c("row", "col", "Plate")) %>%
  dplyr::rename(log2 = log2.x) %>% dplyr::rename(log2S = log2.y) %>%
  mutate(log2 = ifelse(log2S < 8.5, NA, log2)) %>% as.data.frame()
(nrow(filter(df.batch1 , is.na(log2)))/nrow(df.batch1 ))*100 #16.266


write.table(df.batch1, file="df.med_2018_10_09_DIAS.tab", sep="\t",  quote=F, row.names=F)

