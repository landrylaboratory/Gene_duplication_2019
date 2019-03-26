#Create a dataframe of raw data
rm(list=ls())

setwd("/Users/axellemarchant/Documents/postdoc_Landry/AMarchant_2016-2019/papier_AMarchant_2019")

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

## Importing files that contains the plate plan 1536
## PCA1_pl1_A.csv : PCA1 - plate 1- DHFR[1,2]
## PCA1_pl2_A.csv : PCA1 - plate 2- DHFR[1,2]
## PCA2_A.csv : PCA2- DHFR[1,2]
## PCA1_pl1_alpha.csv : PCA1 - plate 1- DHFR[3]
## PCA1_pl2_alpha.csv : PCA1 - plate 2- DHFR[3]
## PCA2_alpha.csv : PCA2- DHFR[3]

#   Colonne 1: Origin of duplication (WGD ? SSD)
#   Colonne 2: Prot.Name (ex: YDL048C)
#   Colonne 3: Plate_1536 number (1 ? 4) 
#   Colonne 4: Column_1536 number (1 ? 48)
#   colonne 5: Row_1536 number (1 ? 32)
#   colonne 6: Numbloc number (NA, 1 ? 308)


PCA1_pl1_A=read.csv("data/PCA1_pl1_A.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA1_pl1_alpha=read.csv("data/PCA1_pl1_alpha.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA1_pl2_A=read.csv("data/PCA1_pl2_A.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA1_pl2_alpha=read.csv("data/PCA1_pl2_alpha.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA2_A=read.csv("data/PCA2_A.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA2_alpha=read.csv("data/PCA2_alpha.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA3_pl1_A=read.csv("data/PCA3_pl1_A.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)
PCA3_pl1_alpha=read.csv("data/PCA3_pl1_alpha.csv", sep="\t", header=T, stringsAsFactor=F, strip.white = T)

## Importing the file containing the information about the photos
#   Colonne 1: File (MTX2_J1_0001.txt)
#   Colonne 2: Mat_A (PCA1_pl1_A ? PCA2_A)
#   Colonne 3: Mat_alpha (PCA1_pl1_a ? PCA2_a)
#   Colonne 4: Plate (1 ? 3)
#   Colonne 6: replicat_plate (1 ? 3)
#   Colonne 7: Etape (S2, DMSO1, MTX1, DMSO2, MTX2)
#   Colonne 8: Condition (gal, eth, NaCl, 30C, 37C)

ip=read.table(file="data/info_photo_ttes_les_batch_30C_no_random.csv", header=T, sep=",")
head(ip)
tail(ip)
nrow(ip) #70

###################################################################################
###################################################################################
###### 2. Creation of the dataframe which will contain all the basic information 
######    of the screen
###################################################################################
###################################################################################


## Creation of an empty dataframe that will be filled in the next steps.
## It will contain the following information: "Duplication", "Mat_A", "Mat_alpha", "Plate",
## "col", "row", "block", "replicat_plate", "step", "condition", "size", "flag"


Colnames = c("Duplication", "Mat_A", "Mat_alpha", "Plate", "col", "row", "bloc", 
             "replicat_plate", "Condition", "PCA", "taille", "flag")

df = data.frame(matrix(vector(), 0, 12, dimnames=list(c(), Colnames)))
head(df)

for (i in 1:nrow(ip))
{
  data = read.table(as.character(ip[i,1]), sep="\t", header=F, stringsAsFactor=F, strip.white = T)
  head(data)
  nrow(data)
  data<-data[order(data[,2], data[,1]),]
  DHFR12 = read.table(as.character(ip[i,2]),sep="\t", header=T, stringsAsFactor=F, strip.white = T)
  head(DHFR12)
  DHFR3 = read.table(as.character(ip[i,3]),sep="\t", header=T, stringsAsFactor=F, strip.white = T)
  head(DHFR3)
  replicat_plate = ip[i,5]
  head(replicat_plate)
  Cond=ip[i,6]
  head(Cond)
  PCA=ip[i,7]
  head(PCA)
  subdf <- cbind(DHFR12[,1], DHFR12[,2], DHFR3[,2], DHFR3[,3:6], rep(replicat_plate, nrow(data)), 
                 rep(Cond, nrow(data)), rep(PCA, nrow(data)), data[,3], data[,5])
  colnames(subdf) <- Colnames
  
  df = rbind(df, subdf)
}

head(df)
tail(df)
nrow(df) # 107520
1536* nrow(ip) #107520 --> OK

write.table(df, file="output/PCA_data_2019_02.tab", sep="\t",  quote=F, row.names=F)
