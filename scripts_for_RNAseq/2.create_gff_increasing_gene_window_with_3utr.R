library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(magrittr)
rm(list=ls())

setwd("/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Genomes/S288C_reference_genome_R64-2-1_20150113")
S288c <- read.table
S288c <- read.table("saccharomyces_cerevisiae_R64-2-1_20150113_AM.gff", header = F, sep="\t")
colnames(S288c) <- c("chr", "strain", "feature", "start", "end", ".1", "fram", ".2", "ID", "Name")
S288c$Name <- str_replace_all(S288c$Name,"ID=","")
S288c.gene <- filter(S288c, feature=="gene")
S288c.chr <- filter(S288c, feature=="chromosome")

UTR <- read.table("../Nagalakshmi_2008/TableS4.csv", header = T, sep="\t")
UTR <- select(UTR, Name, SGD_Start, SGD_End, X5UTR_Start, X3UTR_End)
S288c.gene.utr <- left_join(S288c.gene, UTR, by="Name")

#check if SGD start / end are similar between both databases
filter(S288c.gene.utr, start!=min(SGD_Start, SGD_End) | end!=max(SGD_Start, SGD_End)) #ok

#increase window of gene with UTR --> +/- = to exons
S288c.gene.utr <- S288c.gene.utr %>% mutate(start_review = ifelse(!is.na(X5UTR_Start) & (X5UTR_Start < start), X5UTR_Start, 
                                                           ifelse(!is.na(X3UTR_End) & (X3UTR_End < start), X3UTR_End, start)))
S288c.gene.utr <- S288c.gene.utr %>% mutate(end_review = ifelse(!is.na(X5UTR_Start) & (X5UTR_Start > end), X5UTR_Start, 
                                                         ifelse(!is.na(X3UTR_End) & (X3UTR_End > end), X3UTR_End, end)))

S288c.gene.utr <- select(S288c.gene.utr, "chr", "strain", "feature", "start_review", "end_review", ".1", "fram", ".2", "Name")
colnames(S288c.gene.utr) <- c("chr", "strain", "feature", "start", "end", ".1", "fram", ".2", "Name")

S288c_without_gene <- filter(S288c, feature!="gene")
S288c_without_gene <- select(S288c_without_gene, "chr", "strain", "feature", "start", "end", ".1", "fram", ".2", "Name")

S288c.utr <- rbind(S288c_without_gene, S288c.gene.utr)
S288c.utr <- S288c.utr[with(S288c.gene.utr, order(chr, start)),]

write.table(S288c.utr, file="S288c_with_extended_with_utr.gff", sep="\t",  quote=F, row.names=F)

S288c.chr$group <- "Sequence"
S288c.chr <- select(S288c.chr, "chr", "strain", "feature", "start", "end", ".1", "fram", ".2", "group", "Name")
S288c.gene.utr$group <- "systematic"
S288c.gene.utr <- select(S288c.gene.utr, "chr", "strain", "feature", "start", "end", ".1", "fram", ".2", "group", "Name")
S288c.chr.gene.utr <- rbind(S288c.chr, S288c.gene.utr)

S288c.chr.gene.utr <- S288c.chr.gene.utr[with(S288c.chr.gene.utr, order(chr, start)),]
S288c.chr.gene.utr <- select(S288c.chr.gene.utr, "chr", "strain", "feature", "start", "end", ".1", "fram", ".2", "group", "Name")
write.table(S288c.chr.gene.utr, file="S288c_genes_extended_with_utr.gff", sep="\t",  quote=F, row.names=F)


setwd("/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/RNAseq/alignment_bwa_S288C_chr_renamed_SGD")
