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
library(reshape2)
library(magrittr)
library(gdata)
pat <- "C:/Users/Diana Ascencio/Dropbox/Project_HeteroHomodimers/pca/results/"
# 1) Process all photos with gitter.batch()
setwd(pat)

gitter.batch(image.files = "../raw_data/MTX1_j4/", plate.format = 1536, grid.save = "gridded/", verbose = "p", remove.noise = T)
gitter.batch(image.files = "../raw_data/MTX2_j4/", plate.format = 1536, grid.save = "gridded/", verbose = "p",remove.noise = T)
gitter.batch(image.files = "../raw_data/S2/", plate.format = 1536, grid.save = "gridded/", verbose = "p",remove.noise = T)

# 2) Create an info_phot an index of array files 

fl <- list.files(pattern = ".JPG.dat",full.names = F)
fl <- paste(pat,fl,sep = "")
fs2 <- fl[45:53]
fs2 <- c(fs2[4:9],rep(fs2[1],5),rep(fs2[2],5),rep(fs2[3],6))
fls <- c(fs2,fl[1:44])

pat2 <- "C:/Users/Diana Ascencio/Dropbox/Project_HeteroHomodimers/pca/array_files/"
gpath <- list.files(path = "../array_files/",pattern = "dhfr")
gpath <- paste(pat2,gpath,sep = "")
d12 <- gpath[c(1,3)]
d3 <- gpath[-c(1,3)]
dhfr3 <- c(rep(d3[2:4],2),rep(d3[1],16))
dhfr3 <- rep(dhfr3,3)
dhfr12 <- c(rep(d12[2],6),rep(d12[1],16))
dhfr12 <- rep(dhfr12,3)
plt <- rep(1:22,3)
etp <- c(rep("S2",22),rep("MTX1",22),rep("MTX2",22))
indx <- data.frame(fls,dhfr12,dhfr3,plt,etp)
colnames(indx) <- c("File","MATa","MATalpha","Plate","Etape")

write.csv(indx,file = "info_photo_2018_11_06_DIAS.csv", quote=F,row.names=FALSE)

## Importation du fichier contenant les informations sur les photos (info_photo.csv). 
##  Analyse directement en J4 / derniere photo

ip=read.csv(file="info_photo_2018_11_06_DIAS.csv", sep = ",",header=T, stringsAsFactor=F, strip.white = T)
head(ip)
tail(ip)
nrow(ip) 


# 3) Create a dataframe that contains all the information of the screening
## To create a dataframe that will be use in the next steps of the analysis
## This structure will contain the information: "Duplication", "Mat_A", "Mat_alpha", "Plate", 
## "col", "row", "bloc", "replicat_plate", "Etape", "Condition", "taille", "flag"
## Importation of a file that contains the information files(info_photo.csv). 

Colnames = c("Mat_A", "Mat_alpha", "Plate", "col", "row", "bloc", "Etape", "taille", "flag")

df = data.frame(matrix(vector(), 0, 9, dimnames=list(c(), Colnames)))
head(df)

for (i in 1:nrow(ip))
{
  data = read.table(as.character(ip[i,1]), sep="\t", header=F, stringsAsFactor=F, strip.white = T)
  head(data)
  nrow(data)
  data<-data[order(data[,2], data[,1]),]
  DHFR12 = read.table(as.character(ip[i,2]),sep=",", header=T, stringsAsFactor=F, strip.white = T)
  head(DHFR12)
  #need to format that as well, time, od
  DHFR3 = read.table(as.character(ip[i,3]),sep=",", header=T,  strip.white = T)
  head(DHFR3)
  Plate = ip[i,4]
  head(Plate)
  Etape=ip[i,5]
  head(Etape)
  
  subdf <- cbind(DHFR12[,6], DHFR3[,6], rep(Plate, nrow(data)), DHFR3[,c(2,3,1)], rep(Etape, nrow(data)), data[,3], data[,5])
  colnames(subdf) <- Colnames
  
  df = rbind(df, subdf)
}

head(df)
tail(df)
nrow(df) # 101376
1536* nrow(ip) # 101376 --> OK
write.table(df, file="df_2018_12_03_DIAS.tab", sep="\t",  quote=F, row.names=F)
