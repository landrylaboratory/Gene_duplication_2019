#install.packages("tidyr")
#install.packages("reshape")
library(tidyr)
library(reshape)
library(dplyr)
library(xlsx)
rm(list=ls())
setwd("C:/Users/Diana Ascencio/Dropbox/Project_HeteroHomodimers/pca/array_files/")
dhfr12 <- read.xlsx("DIAS_dest_array_12_20180920.xls",sheetIndex = 1,startRow = 2,stringsAsFactors=FALSE)

# rearrangenment matrix
plt<- matrix(1:1536, nrow = 32,ncol = 48)
colm <- seq(1,48,2)
rowm <- seq(1,32,2)

A <- plt[rowm,colm] 
B <- plt[rowm,colm+1]
C <- plt[rowm+1,colm]
D <- plt[rowm+1,colm+1]
 
 temp <- NULL 
 for (c in 1:48) {
  for (r in 1:32) {
  
   temp <-  rbind(temp, c(c,r))
   }
 }
  
  pl1536 <- temp
  a <- as.vector(A)
  indA <- cbind(a,pl1536[a,])
  colnames(indA) <- c("no1536","c1536","r1536")
  a <- as.vector(B)
  indB <- cbind(a,pl1536[a,])
  colnames(indB) <- c("no1536","c1536","r1536")
  a <- as.vector(C)
  indC <- cbind(a,pl1536[a,])
  colnames(indC) <- c("no1536","c1536","r1536")
  a <- as.vector(D)
  indD <- cbind(a,pl1536[a,])
  colnames(indD) <- c("no1536","c1536","r1536")
  
  

# Generate collection index -----------------------------------------------

  
  # generate the plate for my DHFR[1,2] strains
  pl <- rep("A",384)
  pltA <- cbind(indA,dhfr12,pl)
  pl <- rep("B",384)
  pltB <- cbind(indB,dhfr12,pl)
  pl <- rep("C",384)
  pltC <- cbind(indC,dhfr12,pl)
  pl <- rep("D",384)
  pltD <- cbind(indD,dhfr12,pl)
  
  tp <- rbind(pltA,pltB,pltC,pltD)
  plt_dhfr12 <- arrange(tp,no1536)
  write.csv(plt_dhfr12,file = "plate_dhfr12.csv",  quote=F, row.names=F)
  
  index1536 <- tp %>% select(no1536,c1536,r1536,pl)
  
  # generate the plate for my DHFR[3] strains
  p1 <-  read.xlsx("DIAS_dest_array1_3_20180920.xls",sheetIndex = 1,startRow = 5,stringsAsFactors=FALSE)
  p2 <- read.xlsx("DIAS_dest_array2_3_20180920.xls",sheetIndex = 1,startRow = 5,stringsAsFactors=FALSE)
  p3 <- read.xlsx("DIAS_dest_array3_3_20180920.xls",sheetIndex = 1,startRow = 5,stringsAsFactors=FALSE)
  p4 <-  read.xlsx("DIAS_dest_array4_3_20180920.xls",sheetIndex = 1,startRow = 5,stringsAsFactors=FALSE)
  p5 <- read.xlsx("DIAS_dest_array5_3_20180920.xls",sheetIndex = 1,startRow = 5,stringsAsFactors=FALSE)
  p6 <- read.xlsx("DIAS_dest_array6_3_20180920.xls",sheetIndex = 1,startRow = 5,stringsAsFactors=FALSE)
  
  
  #plates 1-3 dhfr[3]
  
  pl1 <- rbind(p1,p2,p3,p4)
  pl2 <- rbind(p5,p6,p1,p2)
  pl3 <- rbind(p3,p4,p5,p6)
  
  
  plt1_dhfr3 <- cbind(index1536,pl1) %>% arrange(no1536) %>%
    select(no1536,c1536,r1536,c,r,Gene,pl) %>%
    mutate(Gene = ifelse(Gene == "Boder", "border",Gene)) %>%
    mutate(Gene = ifelse(Gene == "NA", "blank",Gene)) 
  
  plt2_dhfr3 <- cbind(index1536,pl2) %>% arrange(no1536)%>%
    select(no1536,c1536,r1536,c,r,Gene,pl) %>%
    mutate(Gene = ifelse(Gene == "Boder", "border",Gene)) %>%
    mutate(Gene = ifelse(Gene == "NA", "blank",Gene)) 
  
  plt3_dhfr3 <- cbind(index1536,pl3) %>% arrange(no1536)%>%
    select(no1536,c1536,r1536,c,r,Gene,pl) %>%
    mutate(Gene = ifelse(Gene == "Boder", "border",Gene)) %>%
    mutate(Gene = ifelse(Gene == "NA", "blank",Gene)) 

 
  
  
  
  # plates after mating
  
  plate1 <- left_join(plt_dhfr12,plt1_dhfr3, by = c("no1536","c1536","r1536","c","r", "pl")) %>%
            select(c(1:5,7,6,8)) %>%  mutate(plate = "plate_1")   
            colnames(plate1)[7:8] <- c("dhfr12","dhfr3")
  
  plate2 <- left_join(plt_dhfr12,plt2_dhfr3, by = c("no1536","c1536","r1536","c","r", "pl")) %>%
            select(c(1:5,7,6,8)) %>%  mutate(plate = "plate_2")       
            colnames(plate2)[7:8] <- c("dhfr12","dhfr3")
  
  plate3 <- left_join(plt_dhfr12,plt3_dhfr3, by = c("no1536","c1536","r1536","c","r", "pl")) %>%
            select(c(1:5,7,6,8))  %>%  mutate(plate = "plate_2")      
            colnames(plate3)[7:8] <- c("dhfr12","dhfr3")
 
  
  
  all_plates <- rbind(plate1,plate2,plate3)
  
  write.csv(plate1,file = "plate1.csv",  quote=F, row.names=F)
  write.csv(plate2,file = "plate2.csv",  quote=F, row.names=F)
  write.csv(plate2,file = "plate3.csv",   quote=F, row.names=F)
  write.csv(plt1_dhfr3,file = "plate1_dhfr3.csv", quote=F, row.names=F)
  write.csv(plt2_dhfr3,file = "plate2_dhfr3.csv", quote=F, row.names=F)
  write.csv(plt3_dhfr3,file = "plate3_dhfr3.csv",  quote=F, row.names=F)
  
  write.csv(all_plates, file = "all_plates.csv",  quote=F, row.names=F)
  

# Check if the array is correct -------------------------------------------

  tim12 <-  grep("YER062C", x = all_plates$dhfr12)
  tim3 <-  grep("YER062C", x = all_plates$dhfr3)
  tim <- intersect(tim12,tim3)
  all_plates[tim,]
  
  
  tp <- grep("23_5_D", all_plates$tag)
  crim <- intersect(rm,cm)
    