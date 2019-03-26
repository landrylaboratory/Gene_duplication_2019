rm(list=ls())
# setwd("D:/Users/dias/Dropbox/Project_HeteroHomodimers/pca/results/")
setwd("C:/Users/Diana Ascencio/Dropbox/Project_HeteroHomodimers/pca/results/")

source("../scripts/theme_Publication.R")
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

df.med <- read.table("df.med_2018_10_09_DEY_DIAS.tab", sep="\t", header=T,stringsAsFactors = F)
df <- read.table("df_2018_12_03_DIAS.tab", sep="\t", header=T, stringsAsFactors = F)

# Define the different conditions used in the experiment
data <- df.med
DMSO <- c(7,11,15,19)
MTXA <-  c(8,12,16,20)
MTXB <-  c(9,13,17,21)
MTXC <-   c(10,14,18,22)
# data%<>% mutate(log2 = log2(taille)) %>% filter(Etape == "MTX2") %>%
         # mutate(type = ifelse(Mat_A == "border","border","interior" ))
data %<>% mutate(Condition = "normal")%>%
            mutate(Condition = ifelse(Plate %in% DMSO, "DMSO",Condition))%>%
            mutate(Condition = ifelse(Plate %in% MTXA, "mtxA",Condition))%>%
            mutate(Condition = ifelse(Plate %in% MTXB, "mtxB",Condition))%>%
            mutate(Condition = ifelse(Plate %in% MTXC, "mtxC",Condition))
graphics.off()           
# Plot data distributions using different tresholds 
theme_set(theme_Publication())
dt <-filter(df, log2(taille) > 0 & Plate %in% c(2,3,4,5,6,7))

g <- ggplot(dt,aes(log2(taille))) + 
  geom_density(aes(fill=factor(Etape), alpha=0.8)) + 
  labs(title="Density plot", 
       x="Colony size, log2",
       fill="Stage")+
  scale_alpha_continuous(guide = FALSE) 
ggsave("../figures/Colony_size_density_etapes.pdf",g)

g+ theme_Publication()


g <- ggplot(data, aes(log2))
g + geom_density(aes(fill=factor(Condition), alpha=0.8)) + 
  labs(title="Density plot", 
              x="log2(size)",
       fill="Etape")

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

d <- data
g <- ggplot(d, aes(log2))
g + geom_density(aes(fill=Condition), size=1, alpha = 0.6)  + 
  geom_vline(xintercept = c(test$mu[1],test$mu[1]+2.5*test$sigma[1]),color="red", linetype="dashed", size=1) +
  labs(title="Colony size", 
       x="log2(size)",
       fill="Etape")+
  facet_grid(Condition ~ ., scales = "free_y")

g <- ggplot(d, aes( (log2 - mu) / sd))
g +  geom_density(aes(fill = Condition), size=1, alpha = 0.6) + 
  labs(title="Z-scores", 
       x="Zscore",
       fill="Etape")+
  geom_vline(xintercept =2,color="red", linetype="dashed", size=1)+
  facet_grid(Condition ~ ., scales = "free_y")


file <- "../figures/Colony_size_density_plots.pdf"
dt <-filter(df, log2(taille) > 0 & Plate %in% c(2,3,4,5,6,7))
theme_set(theme_classic())
g <- ggplot(dt,aes(log2(taille))) + geom_density(aes(fill=factor(Etape), alpha=0.8)) + 
  geom_vline(xintercept = c(test$mu[1],test$mu[1]+2.5*test$sigma[1]))+
  labs(title="Density plot", 
       x="Colony size, log2",
       fill="Stage")+
  scale_alpha_continuous(guide = FALSE) 
ggsave(file,plot = g)
g+theme_Publication()
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

dan_data <- alld %>% filter(Condition!= "normal")

write.csv(dan_data, file= "DEY_PCA_data_11122018.csv", quote=F,row.names=FALSE)

dan_data <- alld %>% filter(Condition == "normal")

write.csv(dan_data, file= "DIAS_PCA_data_11122018.csv", quote=F,row.names=FALSE)


# clustering by MTX concentration--------------------------------------------------------------
  
graphics.off()
for (i in 1:5) {
        d <- alld
        cond <- co[i]
        a<- d  %>% filter(Condition == cond ) %>% 
        select(Mat_A, Mat_alpha,Zscore) 
        a %<>% spread(Mat_alpha, Zscore)
        aa <- as.matrix(apply(a[,-1],2,as.numeric))
        row.names(aa) <- a$Mat_A
        
        colors = c(seq(-10,10,length = 10))
        my_palette <- colorRampPalette(c("grey", "white", "red"))(n = 9)
        heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
                  symm=F,symkey=F,symbreaks=T, scale="none",
                  offsetRow = -0.3, offsetCol = -0.4,
                  cexRow=1,cexCol=1,margins=c(11,11.5),srtCol=60,
                  main=co[i])
        }
  
  row.order <- hclust(dist(aa))$order # clustering rows
  col.order <- hclust(dist(t(aa)))$order # clustering cols
  dat_new <- aa[row.order, col.order] # re-order matrix accoring to clustering
  df_molten_dat <- melt(as.matrix(dat_new)) # reshape into dataframe
  names(df_molten_dat)[c(1:2)] <- c("Mat_A", "Mat_alpha")
  df_molten_dat %<>% mutate(score = ifelse(value<(-10), -10,value))
  uplim <- ceiling(max(df_molten_dat$score,na.rm = T))
  
  ggplot(data =df_molten_dat,
         aes(x = Mat_alpha, y = Mat_A, fill =score)) + 
    geom_raster() +
    scale_fill_distiller(type = "div",palette = 6, direction = -1, limits = c(-uplim,uplim), na.value = "white") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    ggtitle(cond)  


# Compare with ppi reported_biogrid --------------------------------------------------------
graphics.off()
library(xlsx)
biog <- read.xlsx(file = "../DataBases/ScPRS_BioGRID_integrated_11222018DEY.xlsx",sheetIndex = 1,startRow = 5,header = T)

dtb <- select(biog,Bait,Prey,sum) 
a <- dtb %>% select(Bait,Prey,sum) 
a %<>% spread(Prey, sum)
aa <- as.matrix(apply(a[,-1],2,as.numeric))
row.names(aa) <- a$Bait

colors = c(seq(0,13,length = 10))

my_palette <- colorRampPalette(c("white", "lightblue", "purple"))(n = 9)
heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
            symm=F,symkey=F,symbreaks=T, scale="none",cellnote=aa,  notecex=1.0,
            notecol="gray", na.color=par("bg"))

row.order <- hclust(dist(aa))$order # clustering rows
col.order <- hclust(dist(t(aa)))$order # clustering cols
dat_new <- aa[row.order, col.order] # re-order matrix accoring to clustering
df_molten_dat <- melt(as.matrix(dat_new)) # reshape into dataframe
names(df_molten_dat)[c(1:2)] <- c("Bait", "Prey")
# df_molten_dat %<>% mutate(score = ifelse(value>1, 1,value))
uplim <- ceiling(max(df_molten_dat$value,na.rm = T))

ggplot(data =df_molten_dat,
       aes(x = Bait, y = Prey, fill =value)) + 
  geom_raster() +
  scale_fill_distiller(direction = 1, limits = c(0,uplim), na.value = "white") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle("Biogrid")  


# Separate problematic genes ----------------------------------------------
f <- "../GeneLists/15paralogspairs_orthologs_bg.xlsx"
temp <- read.xlsx(file = f,sheetIndex = 1,header = T,stringsAsFactors=FALSE)
Sce1 <- temp %>% select(pair_no, Badgenes, motif,P1) %>% mutate(Species = "Sce1") %>%
        rename("Gene" = P1 )  
Sce2 <- temp %>% select(pair_no, Badgenes, motif,P2) %>% mutate(Species = "Sce2")%>%
        rename("Gene" = P2 )  
Lkluy <-temp %>% select(pair_no, Badgenes, motif,Lkluy) %>% mutate(Species = "Lkluy")%>%
        rename("Gene" = Lkluy )  
Zrou <- temp %>% select(pair_no, Badgenes, motif,Zrou) %>% mutate(Species = "Zrou")%>%
        rename("Gene" = Zrou )  
badg <- rbind(Sce1,Sce2,Lkluy,Zrou) %>% arrange(pair_no)

# Select pair of genes to plot
graphics.off()
d <- filter(badg, pair_no == 13)
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

# Select range of the colorbar
colors = c(seq(-15,15,length = 20)) 
# Select color palette 
my_palette <- colorRampPalette(c("grey", "white", "purple"))(n = 19)
# heatmap plot
heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
          symm=F,symkey=F,symbreaks=T, scale="none",
          offsetRow = -0.3, offsetCol = -0.3, key.xlab="Z-score",
          cexRow=1,cexCol=1,margins=c(6.5,7.5),srtCol=45,dendrogram = "none", Rowv = FALSE,
          lhei=c(1,6), keysize=0.9, key.par = list(cex=0.6))


# Manual clustering to use ggplot
row.order <- hclust(dist(aa))$order # clustering rows
col.order <- hclust(dist(t(aa)))$order # clustering cols
dat_new <- aa[row.order, col.order] # re-order matrix accoring to clustering
df_molten_dat <- melt(as.matrix(dat_new)) # reshape into dataframe
names(df_molten_dat)[c(1:2)] <- c("Mat_A", "Mat_alpha")
df_molten_dat %<>% mutate(score = ifelse(value<(-10), -10,value))
uplim <- ceiling(max(df_molten_dat$score,na.rm = T))

ggplot(data =df_molten_dat,
       aes(x = Mat_alpha, y = Mat_A, fill =score)) + 
  geom_raster() +
  scale_fill_distiller(type = "div",palette = 6, direction = -1, limits = c(-uplim,uplim), na.value = "white") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle(cond)  


# Create a data file with all useful information --------------------------



graphics.off()
finals <- NULL
for(i in 1:15){
    d <- filter(badg, pair_no == i)
    dt <- filter(dan_data, Mat_A %in% d$Gene &  Mat_alpha %in% d$Gene)
    if (nrow(dt)>0){
      dt$pair_no <- d$pair_no[1]
      dt$motif <- d$motif[1]
      dt$problematic <- d$Badgenes[1]
      finals <- rbind(finals,dt)
    }
    
}

write.csv(finals, file= "DIAS_PCA_data_21122018.csv", quote=F,row.names=FALSE)

# figure_paper ------------------------------------------------------------

graphics.off()
d <- filter(badg, pair_no == 13)
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
# Select color palette 
my_palette <- colorRampPalette(c("gray90", "#aaddff",  "#1ea5ff", "#14a1ff"))(n = 14)
# heatmap plot
file <- paste("../figures/blue_TOM70-TOM71",".pdf",sep = "")
pdf(file,width=7.5, height=6.5)
h1 <- heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
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
# Select color palette 
# my_palette <- colorRampPalette(c("grey", "white", "purple"))(n = 14)
# heatmap plot

file <- paste("../figures/blue_TAL1-NQM1",".pdf",sep = "")
pdf(file,width=7.5, height=6.5)
heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
          symm=F,symkey=F,symbreaks=T, scale="none", 
          offsetRow = -0.3, offsetCol = -0.3, key.xlab="Z-score", key.title = "",
          cexRow=1.6,cexCol=1.6,margins=c(6.5,7.5),srtCol=45, dendrogram = "none",
          lhei=c(1.6,7), keysize=1.2, cex.main = 1,
          main="TAL1 (P1) / NQM1 (P2)")
dev.off()
