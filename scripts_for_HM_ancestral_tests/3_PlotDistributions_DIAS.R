
rm(list=ls())
setwd("C:/Users/Diana Ascencio/Dropbox/Project_HeteroHomodimers/pca/results/")
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

df.med <- read.table("df.med_2018_10_09_DEY_DIAS.tab.tab", sep="\t", header=T,stringsAsFactors = F)
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

          
# Plot data distributions using different tresholds 
theme_set(theme_classic())
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
data %<>% mutate (Zscore_mtx2 = (log2 - mu) / sd)
data %<>% filter(!is.nan(Zscore_mtx2))

d <- data
pdf("../figures/Condition_density_DEYDIAS_11122018.pdf")
g <- ggplot(d, aes(log2))
g + geom_density(aes(fill=Condition), size=1, alpha = 0.6)  + 
  geom_vline(xintercept = c(test$mu[1],test$mu[1]+2.5*test$sigma[1]),color="red", linetype="dashed", size=1) +
  labs(title="Colony size", 
       x="log2(size)",
       fill="Etape")+
  facet_grid(Condition ~ ., scales = "free_y")
dev.off()

pdf("../figures/Zscore_density_DEYDIAS_11122018.pdf")
g <- ggplot(d, aes( (log2 - mu) / sd))
g +  geom_density(aes(fill = Condition), size=1, alpha = 0.6) + 
  labs(title="Z-scores", 
       x="Zscore",
       fill="Etape")+
  geom_vline(xintercept =2,color="red", linetype="dashed", size=1)+
  facet_grid(Condition ~ ., scales = "free_y")
dev.off()

# Generate files with the data --------------------------------------------
cdt <- c("normal","DMSO", "mtxA","mtxB","mtxC")
co <- c("normal","DMSO", "mtxA_150ug","mtxB_175ug","mtxC_200ug")
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

write.csv(dan_data, file= "DIAS_PCA_data_11122018.csv", quote=F,row.names=FALSE)


# clustering --------------------------------------------------------------
graphics.off()
colors = c(seq(-10,10,length = 10))
my_palette <- colorRampPalette(c("grey", "white", "red"))(n = 9)

for (i in 1:5) {
        d <- alld
        cond <- co[i]
        a<- d  %>% filter(Condition == cond ) %>% 
        select(Mat_A, Mat_alpha,Zscore) 
        a %<>% spread(Mat_alpha, Zscore)
        aa <- as.matrix(apply(a[,-1],2,as.numeric))
        row.names(aa) <- a$Mat_A
        
        # pdf(paste("../figures/",co[i],"_Heatmap_DEYDIAS_11122018.pdf",sep = ""))
        heatmap.2(aa, col=my_palette, breaks = colors, density.info="none", trace="none", 
                  symm=F,symkey=F,symbreaks=T, scale="none",
                  offsetRow = -0.3, offsetCol = -0.4,
                  cexRow=1,cexCol=1,margins=c(11,11.5),srtCol=60,
                  main=co[i])
        # dev.off()
        }
  
 



