# This script analize 1532 yeast colony plate pictures using gitter_batch() 
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

# 1) Process all photos with gitter.batch()
gitter.batch(image.files = "../raw_data/MTX1_j4/", plate.format = 1536, grid.save = "gridded/", verbose = "p", remove.noise = T)
gitter.batch(image.files = "../raw_data/MTX2_j4/", plate.format = 1536, grid.save = "gridded/", verbose = "p",remove.noise = T)
gitter.batch(image.files = "../raw_data/S2/", plate.format = 1536, grid.save = "gridded/", verbose = "p",remove.noise = T)

# 2) load info_photo an index of array files 

ip=read.csv(file="info_photo_2018_11_06_DIAS.csv", sep = ",",header=T, stringsAsFactor=F, strip.white = T)
head(ip)
tail(ip)
nrow(ip) 
ip <- indx

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
