# This script calculate Zscores and generates final figures

source("theme_Publication.R")
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
library(ggExtra)
library(xlsx)
rm(list = ls())
df.med <- read.table("df.med_2018_10_09_DIAS.tab", sep="\t", header=T,stringsAsFactors = F)

###########################################################################################
#                              Calculate Zscores                                          #
###########################################################################################
data <- df.med
DMSO <- c(7,11,15,19)
MTXA <-  c(8,12,16,20)
MTXB <-  c(9,13,17,21)
MTXC <-   c(10,14,18,22)
data %<>% mutate(Condition = "normal")%>%
            mutate(Condition = ifelse(Plate %in% DMSO, "DMSO",Condition))%>%
            mutate(Condition = ifelse(Plate %in% MTXA, "mtxA",Condition))%>%
            mutate(Condition = ifelse(Plate %in% MTXB, "mtxB",Condition))%>%
            mutate(Condition = ifelse(Plate %in% MTXC, "mtxC",Condition))
graphics.off()           

#calculate mean and sd
dat <- filter(data, log2 > 6.5 & Condition =="normal")
test <- normalmixEM(dat$log2)
mu <- test$mu[1]
sd <- test$sigma[1]
th <- mu+2.5*sd
mu
sd
th
data %<>% mutate (Zscore_mtx2 = (log2 - mu) / sd)
data %<>% filter(!is.nan(Zscore_mtx2))

# Generate files with the data --------------------------------------------
cdt <- c("normal","DMSO", "mtxA","mtxB","mtxC")
co <- c("normal","DMSO", "mtxA(150ug/ml)","mtxB(175ug/ml)","mtxC(200ug/ml)")
alld <- NULL
for(i in 1:5){
  cond <- cdt[i]
  dta1 <- data %>% filter(Condition == cond) 
  dta1 %<>% select(Mat_A, Mat_alpha,log2) %>%
    mutate(interaction = paste(Mat_A,Mat_alpha, sep = "/"))  %>%
    group_by(interaction) %>%
    summarise(colony_size = median(log2), n = n(),colony_size_sd = sd(log2)) 
   
  dta2 <- data %>% filter(Condition == cond) 
  dta2 %<>% select(Mat_A, Mat_alpha,Zscore_mtx2) %>%
    mutate(interaction = paste(Mat_A,Mat_alpha, sep = "/"))  %>%
    group_by(interaction) %>%
    summarise(Zscore = median(Zscore_mtx2), n = n(),Zscore_sd = sd(Zscore_mtx2)) 
  dta <- left_join(dta1,dta2, by = c("interaction", "n"))
  
  dta %<>% filter(!is.na(colony_size)) %>% 
    separate(interaction, c("Mat_A", "Mat_alpha"), "/") %>%
    mutate(Condition = co[i]) 
    alld <- rbind(alld,dta)
 
}

alld %<>% select(Mat_A,Mat_alpha,Condition,"n", colony_size,colony_size_sd,Zscore, Zscore_sd)

dias_data <- alld %>% filter(Condition == "normal")
# Order genes by orthologous group ----------------------------------------------
f <- "paralogspairs_orthologs_bg.csv"
temp <- read.csv(file = f,header = T,stringsAsFactors=FALSE, fileEncoding = "UTF-8-BOM")
Sce1 <- temp %>% select(pair_no,bg, motif,P1) %>% mutate(Species = "Sce1") %>%
  rename("Gene" = P1 )  
Sce2 <- temp %>% select(pair_no, bg, motif,P2) %>% mutate(Species = "Sce2")%>%
  rename("Gene" = P2 )  
Lkluy <-temp %>% select(pair_no,bg, motif,Lkluy) %>% mutate(Species = "Lkluy")%>%
  rename("Gene" = Lkluy )  
Zrou <- temp %>% select(pair_no, bg, motif,Zrou) %>% mutate(Species = "Zrou")%>%
  rename("Gene" = Zrou )  
badg <- rbind(Sce1,Sce2,Lkluy,Zrou) %>% arrange(pair_no)



###########################################################################################
#                             Generate final figures                                    #
###########################################################################################
data <- df.med
d <- filter(badg, pair_no == 13)
glimpse(d)
dta <- filter(dias_data,  Mat_A %in% d$Gene &  Mat_alpha %in% d$Gene)
dta %<>% select(Mat_A, Mat_alpha,Zscore) %>%
  mutate(interaction = paste(Mat_A,Mat_alpha, sep = "/"))  %>%
  group_by(interaction) %>%
  summarise(mean = mean(Zscore), n = n(),sd = sd(Zscore)) 
dta %<>% filter(!is.na(mean))
dta %<>% separate(interaction, c("Mat_A", "Mat_alpha"), "/")
a <- dta %>% select(Mat_A,Mat_alpha,mean)
a %<>% spread(Mat_alpha, mean)
aa <- as.matrix(apply(a[,-1],2,as.numeric))
row.names(aa) <- a$Mat_A

row.names(aa) <- c("Lkluy","Scer (P2)","Scer (P1)","Zrou")
colnames(aa) <- c("Lkluy","Scer (P2)","Scer (P1)","Zrou")
# Select range of the colorbar
colors = c(seq(0,12,length = 15)) 
# Select color palette 
my_palette <- colorRampPalette(c("grey88","grey66", "grey40", "black"))(n = 14)
# heatmap plot
file <- paste("../figures/TOM70-TOM71",".pdf",sep = "")
pdf(file,width=7.5, height=6.5)
heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
          symm=F,symkey=F,symbreaks=T, scale="none", 
          offsetRow = -0.3, offsetCol = -0.3, key.xlab="Z-score", key.title = "",
          cexRow=1.6,cexCol=1.6,margins=c(6.5,7.5),srtCol=45, dendrogram = "none",
          lhei=c(1.6,7), keysize=1.2, cex.main = 1,
          main="TOM70 (P1) / TOM71 (P2)")
dev.off()


###### second figure
graphics.off()
d <- filter(badg, pair_no ==3)
glimpse(d)
cond <- "normal"
dta <- filter(data, Condition == cond & Mat_A %in% d$Gene &  Mat_alpha %in% d$Gene)
dta %<>% select(Mat_A, Mat_alpha,Zscore_mtx2) %>%
  mutate(interaction = paste(Mat_A,Mat_alpha, sep = "/"))  %>%
  group_by(interaction) %>%
  summarise(mean = mean(Zscore_mtx2), n = n(),sd = sd(Zscore_mtx2)) 
dta %<>% filter(!is.na(mean))
dta %<>% separate(interaction, c("Mat_A", "Mat_alpha"), "/")
a <- dta %>% select(Mat_A,Mat_alpha,mean)
a %<>% spread(Mat_alpha, mean)
aa <- as.matrix(apply(a[,-1],2,as.numeric))
row.names(aa) <- a$Mat_A

row.names(aa) <- c("Lkluy","Scer (P2)","Scer (P1)","Zrou")
colnames(aa) <- c("Lkluy","Scer (P2)","Scer (P1)","Zrou")

# Select range of the colorbar
colors = c(seq(0,12,length = 15)) 

# heatmap plot

file <- paste("../figures/TAL1-NQM1",".pdf",sep = "")
pdf(file,width=7.5, height=6.5)
heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
          symm=F,symkey=F,symbreaks=T, scale="none", 
          offsetRow = -0.3, offsetCol = -0.3, key.xlab="Z-score", key.title = "",
          cexRow=1.6,cexCol=1.6,margins=c(6.5,7.5),srtCol=45, dendrogram = "none",
          lhei=c(1.6,7), keysize=1.2, cex.main = 1,
          main="TAL1 (P1) / NQM1 (P2)")
dev.off()
