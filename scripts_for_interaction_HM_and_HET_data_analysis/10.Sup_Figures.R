rm(list=ls())

setwd("/Users/axellemarchant/Documents/postdoc_Landry/AMarchant_2016-2019/papier_AMarchant_2019")

require(gitter)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggsignif)
library(grid)
library("cowplot")
library(ggsignif)
library(ggpubr)
library("tidyverse")
library(UpSetR)
cl <- colors()[]

summary_table <- read.delim('output/TableS1.csv', header=T, sep="\t")


#FigS1
Interact <- read.table("output/comp_Biogrid.csv", sep="\t",  header=T)

Biogrid <-
  ggplot(Interact, aes(x=as.factor(N.report.in.biogrid), y=median.size.PCA)) + 
  geom_jitter(width = 0.2) + 
  geom_violin(alpha=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Number of times reported in BioGRID") + 
  ylab("Median colony size (pixels, log2)") +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')


HM_Stynen_our_PCA <- read.table("output/comp_Stynen.csv", sep="\t",  header=T)

Michnick <- 
  ggplot(HM_Stynen_our_PCA, aes(x=med_MTX_our_PCA, y=med_MTX_Stynen)) + 
  geom_point() + geom_smooth(method=lm) +
  ylab("Median colony size (pixels, log2) from Stynen") + 
  xlab("Median colony size (pixels, log2) from this study") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  annotation_custom(grobTree(textGrob("Spearman: r = 0.58 p-value < 2.2e-16", 
                                      x=0.05,  y=0.98, hjust=0, gp=gpar(fontsize=12))))

HM_Tarassov_our_PCA <- read.table("output/comp_Tarassov.csv", sep="\t",  header=T)

Tarassov <-
  ggplot(HM_Tarassov_our_PCA, aes(x=med_MTX_our_PCA, y=log2(intensity_Tarassov))) + 
  geom_point() + geom_smooth(method=lm) +
  ylab("Integrated colony pixel intensity (log2) from Tarassov") + 
  xlab("Median colony size (pixels, log2) from this study") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  annotation_custom(grobTree(textGrob("Spearman: r = 0.55 p-value = 1.86e-08", 
                                      x=0.05,  y=0.98, hjust=0, gp=gpar(fontsize=12))))


ggsave(file="sup_fig/FigS1.pdf", width=7, height=15, dpi=500)
plot_grid(Biogrid, Michnick, Tarassov, labels = c("A", "B", "C"), nrow=3, ncol=1)
dev.off()


######FigS2
HM1 <- select(summary_table, P1, HM1.PCA, HM1.PDB, HM1.Kim, HM1.Stynen, list.bg.HM1.Pubmed.id)
HM1$ppi <- apply(cbind(as.character(HM1$P1), as.character(HM1$P1)), 
                 1, function(x) paste(x, collapse="."))
HM1 <- select(HM1, ppi, HM1.PCA, HM1.PDB, HM1.Kim, HM1.Stynen, list.bg.HM1.Pubmed.id)
colnames(HM1) <- c('ppi', 'PCA_this_study', 'PDB', 'BiFC', 'PCA.S', 'list.bg')
HM2 <- select(summary_table, P2, HM2.PCA, HM2.PDB, HM2.Kim, HM2.Stynen, list.bg.HM2.Pubmed.id)
HM2$ppi <- apply(cbind(as.character(HM2$P2), as.character(HM2$P2)), 
                 1, function(x) paste(x, collapse="."))
HM2 <- select(HM2, ppi, HM2.PCA, HM2.PDB, HM2.Kim, HM2.Stynen, list.bg.HM2.Pubmed.id)
colnames(HM2) <- c('ppi', 'PCA_this_study', 'PDB', 'BiFC', 'PCA.S', 'list.bg')
HM <- rbind(HM1, HM2)
HM$PDB[is.na(HM$PDB)] <- 0
HM$BiFC[is.na(HM$BiFC)] <- 0
HM$PCA_this_study[is.na(HM$PCA_this_study)] <- 0

HM <- HM %>% mutate(Affinity_Capture_MS=ifelse(is.na(list.bg), 0,
                                        ifelse(str_detect(list.bg, "Affinity Capture-MS"), 1, 0)))
HM <- HM %>% mutate(Reconstituted_Complex=ifelse(is.na(list.bg), 0,
                                          ifelse(str_detect(list.bg, "Reconstituted Complex"), 1, 0)))
HM <- HM %>% mutate(Biochemical_Activity=ifelse(is.na(list.bg), 0,
                                                ifelse(str_detect(list.bg, "Biochemical Activity"), 1, 0)))
HM <- HM %>% mutate(Far_Western=ifelse(is.na(list.bg), 0,
                                       ifelse(str_detect(list.bg, "Far Western"), 1, 0)))
HM <- HM %>% mutate(Protein_peptide=ifelse(is.na(list.bg), 0,
                                           ifelse(str_detect(list.bg, "Protein-peptide"), 1, 0)))
HM <- HM %>% mutate(Affinity_Capture_Lumi=ifelse(is.na(list.bg), 0,
                                                 ifelse(str_detect(list.bg, "Affinity Capture-Luminescence"), 1, 0)))
HM <- HM %>% mutate(Affinity_Capture_Western=ifelse(is.na(list.bg), 0,
                                                    ifelse(str_detect(list.bg, "Affinity Capture-Western"), 1, 0)))
HM <- HM %>% mutate(Two_hybrid=ifelse(is.na(list.bg), 0,
                                      ifelse(str_detect(list.bg, "Two-hybrid"), 1, 0)))
HM <- HM %>% mutate(Co_crystal_Structure=ifelse(is.na(list.bg), 0,
                                                ifelse(str_detect(list.bg, "Co-crystal Structure"), 1, 0)))
HM <- HM %>% mutate(FRET=ifelse(is.na(list.bg), 0,
                                ifelse(str_detect(list.bg, "FRET"), 1, 0)))
HM <- HM %>% mutate(PCA.bg=ifelse(is.na(list.bg), 0,
                                  ifelse(str_detect(list.bg, "PCA"), 1, 0)))
HM <- HM %>% mutate(PCA = ifelse(is.na(PCA.S), PCA.bg,
                                 ifelse(PCA.S==1 | PCA.bg==1, 1, 0)))

HM <- select(HM, -PCA.S, -PCA.bg, -list.bg)

#FigS2A
upset(HM,  sets = colnames(HM)[2:15], sets.bar.color = "pink", point.size = 1.5, order.by = "freq",
      mainbar.y.label = "Method intersections for HM", sets.x.label = "HMs per method")


HET <- select(summary_table, pair, HET.ourPCA, list.bg.HET.Pubmed.id, HET.PDB)
colnames(HET) <- c('ppi', 'PCA_this_study', 'list.bg', 'PDB')
HET$PCA_this_study[is.na(HET$PCA_this_study)] <- 0


HET <- HET %>% mutate(Affinity_Capture_MS=ifelse(is.na(list.bg), 0,
                                                 ifelse(str_detect(list.bg, "Affinity Capture-MS"), 1, 0)))
HET <- HET %>% mutate(Reconstituted_Complex=ifelse(is.na(list.bg), 0,
                                                   ifelse(str_detect(list.bg, "Reconstituted Complex"), 1, 0)))
HET <- HET %>% mutate(Biochemical_Activity=ifelse(is.na(list.bg), 0,
                                                  ifelse(str_detect(list.bg, "Biochemical Activity"), 1, 0)))
HET <- HET %>% mutate(Far_Western=ifelse(is.na(list.bg), 0,
                                         ifelse(str_detect(list.bg, "Far Western"), 1, 0)))
HET <- HET %>% mutate(Protein_peptide=ifelse(is.na(list.bg), 0,
                                             ifelse(str_detect(list.bg, "Protein-peptide"), 1, 0)))
HET <- HET %>% mutate(Affinity_Capture_Lumi=ifelse(is.na(list.bg), 0,
                                                   ifelse(str_detect(list.bg, "Affinity Capture-Luminescence"), 1, 0)))
HET <- HET %>% mutate(Affinity_Capture_Western=ifelse(is.na(list.bg), 0,
                                                      ifelse(str_detect(list.bg, "Affinity Capture-Western"), 1, 0)))
HET <- HET %>% mutate(Two_hybrid=ifelse(is.na(list.bg), 0,
                                        ifelse(str_detect(list.bg, "Two-hybrid"), 1, 0)))
HET <- HET %>% mutate(Co_crystal_Structure=ifelse(is.na(list.bg), 0,
                                                  ifelse(str_detect(list.bg, "Co-crystal Structure"), 1, 0)))
HET <- HET %>% mutate(FRET=ifelse(is.na(list.bg), 0,
                                  ifelse(str_detect(list.bg, "FRET"), 1, 0)))
HET <- HET %>% mutate(PCA=ifelse(is.na(list.bg), 0,
                                 ifelse(str_detect(list.bg, "PCA"), 1, 0)))
HET <- select(HET, -list.bg)

#FigS2B
  upset(HET,  sets = colnames(HET)[2:14], sets.bar.color = "purple", point.size = 1.5, order.by = "freq",
        mainbar.y.label = "Method intersections for HET", sets.x.label = "HETs per method")

#####FigS3
#FigS3A

HM1.PCA <- select(summary_table, P1, HM1.PCA, log10.RNAseq.MTX.P1, HM1bg.S.PDB.Kim)
HM2.PCA <- select(summary_table, P2, HM2.PCA, log10.RNAseq.MTX.P2, HM2bg.S.PDB.Kim)
colnames(HM1.PCA) <- c("paralog", "HM", "log10_RNAseq_MTX", "HMbg.S.PDB.Kim")
colnames(HM2.PCA) <- c("paralog", "HM", "log10_RNAseq_MTX", "HMbg.S.PDB.Kim")
HM <- rbind(HM1.PCA, HM2.PCA)
HM <- HM %>% .[complete.cases(.), ]
HM <- filter(HM, HMbg.S.PDB.Kim==1)
lines <- (ksmooth(HM$log10_RNAseq_MTX, HM$HM, "normal", 
                  bandwidth = 0.8, x.points = seq(0, 5, by=0.001))) %>% as.data.frame()

FigS3A <-
  ggplot(HM, aes(x=log10_RNAseq_MTX, y=HM)) + geom_point() + geom_line(data=lines, aes(x=x, y=y))+
  xlab("log10(mapped RNAseq read MTX) of HM previously reported") +  ylab("P(HM)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, lineheight=0.4, hjust=0.5), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.1, 0.95), 
        legend.title = element_text(size=15), legend.text = element_text(size=12))

#FigS3B
df <- filter(summary_table, !is.na(med.MTX.HM.P1) & !is.na(med.MTX.HM.P2)) %>% 
  mutate(diff_expP2P1 = log10.RNAseq.MTX.P2-log10.RNAseq.MTX.P1,
         tot_exp_P1P2 = log10.RNAseq.MTX.P2+log10.RNAseq.MTX.P1)

FigS3B <- ggplot(df, aes(x=med.MTX.HM.P1, y=med.MTX.HM.P2, color=diff_expP2P1))+
  geom_point(aes(size = 2^tot_exp_P1P2))+
  scale_color_gradient2(midpoint=0, low="blue", mid="gray",
                        high="red", space ="Lab" ) +
  xlab("PCA score of paralog 1 (P1)") + 
  ylab("PCA score of paralog 2 (P2)") +
  guides( colour = guide_colorbar(title = "Difference of expression P2-P1"),
          size = guide_legend("Total expression P1+P2")) +
  scale_size_continuous(breaks = c(1, 4, 16, 64, 256), 
                        labels = c(log2(1), log2(4), log2(16), log2(64), log2(256)))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0, 12),
        legend.title = element_text(size=15), legend.text = element_text(size=12))


#FigS3C
HM1_rep_undec <- filter(summary_table, HM1.PCA==0 & HM1bg.S.PDB.Kim==1)
HM1_rep_undec <- select(HM1_rep_undec, P1, HM1.PCA, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
colnames(HM1_rep_undec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")
HM2_rep_undec <- filter(summary_table, HM2.PCA==0 & HM2bg.S.PDB.Kim==1)
HM2_rep_undec <- select(HM2_rep_undec, P2, HM2.PCA, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM2_rep_undec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")

HM_rep_undec <- rbind(HM1_rep_undec, HM2_rep_undec)

HM1_rep_dec <- filter(summary_table, HM1.PCA==1 & HM1bg.S.PDB.Kim==1)
HM1_rep_dec <- select(HM1_rep_dec, P1, HM1.PCA, log10.RNAseq.MTX.P1, log10.abundance.paxDB.P1)
colnames(HM1_rep_dec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")
HM2_rep_dec <- filter(summary_table, HM2.PCA==1 & HM2bg.S.PDB.Kim==1)
HM2_rep_dec <- select(HM2_rep_dec, P2, HM2.PCA, log10.RNAseq.MTX.P2, log10.abundance.paxDB.P2)
colnames(HM2_rep_dec) <- c("orf", "HM.PCA", "log10_RNAseq_MTX", "log10_abundance_paxDB")

HM_rep_dec <- rbind(HM1_rep_dec, HM2_rep_dec)
wilcox.test(HM_rep_undec$log10_RNAseq_MTX, HM_rep_dec$log10_RNAseq_MTX)
#p-value = 1.42e-05
wilcox.test(HM_rep_undec$log10_abundance_paxDB, HM_rep_dec$log10_abundance_paxDB)
#p-value = 1.651e-05

HM_rep_tested <- rbind(HM_rep_dec, HM_rep_undec)

my_comparisons <- list(c("0", "1"))


library(ggpubr)

FigS3C <-
  ggboxplot(HM_rep_tested,  x="HM.PCA", y="log10_RNAseq_MTX")+ 
  stat_compare_means(comparisons = my_comparisons, group.by = "HM.PCA", method = "wilcox.test") +
  ylab("log10(number of mapped RNAseq reads in MTX)") + xlab("Detection by PCA of HM previously reported") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')


library("cowplot")
ggsave(file="sup_fig/FigS3.pdf", width=14, height=7, dpi=500)
plot_grid(FigS3A, FigS3B, FigS3C, labels = c("A", "B", "C"), nrow=3, ncol=1)
dev.off()


#FigS4
dff.S.bg.Kim.PCA <- read.table('output/HM_expression.csv', header=T, sep='\t')
my_comparisons <- list( c("S", "ssd"), c("S", "wgd"), c("ssd", "wgd"))

FigS4A <-
  ggplot(filter(dff.S.bg.Kim.PCA,type_para!="ssd_wgd") , aes(x=type_para, y=log10_RNA_med_MTX, fill=type_para)) + 
  geom_boxplot()  +
  scale_fill_manual(guide=FALSE,values = c(cl[16],cl[144],cl[129],cl[310])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  scale_x_discrete(labels=c(c("S","SSD","WGD"))) +
  ylab("mRNA level in MTX") + xlab("Groups of genes")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0.08, 0.9), legend.title = element_text(size=15), 
        legend.text = element_text(size=12))

FigS4B <-
  ggplot(filter(dff.S.bg.Kim.PCA, type_para!="ssd_wgd" & type_para!="ssd-successiv") , 
         aes(x=type_para, y=log(paxdb), fill=type_para)) + 
  geom_boxplot()  +
  scale_fill_manual(guide=FALSE,values = c(cl[16],cl[144],cl[129],cl[310])) +
  scale_x_discrete(labels=c(c("S","SSD","WGD"))) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  ylab("Protein abundance (log)") + xlab("Groups of genes")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0.08, 0.9), legend.title = element_text(size=15), 
        legend.text = element_text(size=12))


ggsave(file="sup_fig/FigureS4.pdf", width=14, height=7, dpi=500)
plot_grid(FigS4A, FigS4B,
          labels = c("A", "B"), nrow=1, ncol = 2, label_size = 12)
dev.off()

#######FigS5#####
#Fig5A
my_comparisons <- list(c("Homeologs","True_ohnologs"))
summary_table$Origin.of.WGDs <- factor(summary_table$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

FigS5A <-
  ggplot(filter(summary_table, !is.na(Origin.of.WGDs)), aes(x=Origin.of.WGDs, y=pid, fill=Origin.of.WGDs)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  size = 4) +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Origin of WGDs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')


#Fig5B
HM <- filter(summary_table, motif.categories=="HM" & !is.na(age.group)) %>% 
  group_by(age.group, Duplication) %>%  summarise(HM=length(motif.categories)) %>% as.data.frame()
HM_HET <-  filter(summary_table, motif.categories=="HM&HET" & !is.na(age.group)) %>% 
  group_by(age.group, Duplication) %>%  summarise("HM_HET"=length(motif.categories)) %>% as.data.frame()
HM$age.group <- as.character(HM$age.group)
HM_HET$age.group <- as.character(HM_HET$age.group)
HM$HM <- as.numeric(as.character(HM$HM))
HM_HET$HM_HET <- as.numeric(as.character(HM_HET$HM_HET))
ratio_HM_HET <- left_join(HM, HM_HET, by=c("age.group", "Duplication")) 
ratio_HM_HET$HM_HET[is.na(ratio_HM_HET$HM_HET)] <- 0
ratio_HM_HET <- ratio_HM_HET %>% mutate(ratioHETsurHM=(HM_HET/HM))
ratio_HM_HET <- ratio_HM_HET %>% group_by(age.group, Duplication) %>%
  mutate(effectif=paste(c(HM_HET,HM), collapse = "/"))

FigS5B <-
  ggplot(ratio_HM_HET, aes(x=as.factor(age.group), y=as.numeric(ratioHETsurHM),  fill=Duplication)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  ylab("HM&HET / HM ratio") + xlab("Age group") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  geom_text(aes(label=effectif), vjust=-0.3, size=3.5)


#Fig5C
my_comparisons <- list( c("HM", "HET"), c("HET", "HM&HET"), c("HET", "NI"),
                        c("HM", "HM&HET"), c("HM", "NI"),
                        c("HM&HET", "NI"))
FigS5C <- ggboxplot(filter(summary_table, Duplication=="WGD", !is.na(Origin.of.WGDs)),  x="motif.categories", y="pid", fill="Origin.of.WGDs", facet.by = "Origin.of.WGDs")+ 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test") +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')


#Fig5D
FigS5D <- 
  ggplot(filter(summary_table, Duplication=="SSD" & !is.na(age.group)), 
         aes(x=as.character(age.group), y=pid)) + 
  geom_boxplot() + 
  ylab("Pairwise amino acid sequence identity (%)") + xlab("SSD age group") + 
  scale_fill_manual(values = cl[129]) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')



ggsave(file="sup_fig/Fig5S.pdf", width=14, height=14, dpi=500)
plot_grid(FigS5A, FigS5B, FigS5C, FigS5D, labels = c("A", "B", "C", "D"), ncol=2, nrow = 2)
dev.off()


#FigS7
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}

WGD <- filter(summary_table, Duplication=="WGD", !is.na(Origin.of.WGDs))
WGD$int <- interaction(WGD$motif.categories, WGD$Origin.of.WGDs)

my_comparisons <- list(c("HM.Homeologs", "HM&HET.Homeologs"), c("HM.True_ohnologs", "HM&HET.True_ohnologs"))

FigS7A <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.mol.fct.P1P2, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) +
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("GO molecular function similarity") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=1.1)



FigS7B <-WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.bio.proc.P1P2, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("GO biological process similarity") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=1.1)


FigS7C <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.pheno.P1P2, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("Similarity of phenotypes") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=0.55)


FigS7D  <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=med.gi.cor, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("Correlation of genetic interaction profile") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=0.3)


ggsave(file="sup_fig/FigureS7.pdf", width=10, height=10, dpi=500)
plot_grid(FigS7A, FigS7B, FigS7C, FigS7D, labels=c("A", "B", "C", "D"), nrow=2, align=c("h"))
dev.off()



#FigS8
tot <- summary_table %>% group_by(Duplication) %>% filter(., !is.na(motif.number)) %>%  
  summarise(tot = length(pair))
n <- summary_table %>% group_by(Duplication, motif.number) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.number)) %>% as.data.frame()

freq <- full_join(n, tot, by="Duplication") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n

#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))
diff_SSD_WGD_per_motif.simple <- freq %>% group_by(motif.number) %>% 
  #summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()


FigS8A<-
  ggplot(freq, aes(as.factor(motif.number), freq)) + 
  geom_bar(aes(fill = Duplication), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  labs(x=("Interaction motifs"), y=("Percentage (%)"), title = ("Duplications")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75),xend=c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), 
                              y=c(16, 40, 15, 10, 17.5, 25), 
                              annotation=c("0.29", "3.32e-03", "2.97e-02", "0.38", "9.99e-3", "4.79e-4")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.9, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.background = element_blank())

tot <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs) %>% filter(., !is.na(motif.number)) %>%  
  summarise(tot = length(pair))
n <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs, motif.number) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.number)) %>% as.data.frame()

freq <- full_join(n, tot, by="Origin.of.WGDs") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n


#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

diff_Inffered.topo_per_motif <- freq %>% group_by(motif.number) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
#motif.number p.value_fisherTest
#1            1         1.00000000
#2            2         0.04222921
#3            3         1.00000000
#4            4         0.40946337
#5            5         0.22081352
#6            6         0.10468644

freq$Origin.of.WGDs <- factor(freq$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

FigS8B <-
  ggplot(freq, aes(as.factor(motif.number), freq)) + 
  geom_bar(aes(fill = Origin.of.WGDs), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[430],cl[490]), labels = c("Homeologs", "True ohnologs")) + 
  labs(x=("Interaction motifs"), title = ("Origin of WGDs")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75),xend=c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), 
                              y=c(15.5, 32, 11, 15.5, 20, 33.5), 
                              annotation=c("1.00", "4.22e-02", " 1.00 ", "0.41", "0.22", "0.10")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.title.x = element_text(size=15), 
        axis.text.y= element_text(size=12), axis.title.y = element_text(size=15), 
        #axis.title.y = element_blank(), axis.ticks.y =  element_blank(), axis.text.y= element_blank(), 
        plot.title = element_text(size=15, hjust = 0.5), 
        legend.position = c(0.3, 0.94), legend.title = element_blank(), legend.text = element_text(size=12), 
        legend.background = element_blank())

ggsave(file="sup_fig/FigS8.pdf", width=7, height=7, dpi=500)
plot_grid(SupMotigNubr, SupMotigNubrWGD, nrow=1, ncol=2)
dev.off()



###FigS13

my_comparisons <- list(c("Homeologs","True_ohnologs"))
summary_table$Origin.of.WGDs <- factor(summary_table$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

FigS13A <-
  ggplot(filter(summary_table, !is.na(Origin.of.WGDs)), aes(x=Origin.of.WGDs, y=expression.correl.coeff, fill=Origin.of.WGDs)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 4) +
  ylab("Correlation coefficient (r)") + xlab("Origin of WGDs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')

my_comparisons <- list( c("HM", "HET"), c("HET", "HM&HET"), c("HET", "NI"),
                        c("HM", "HM&HET"), c("HM", "NI"),
                        c("HM&HET", "NI"))

FigS13B <- ggboxplot(filter(summary_table, Duplication=="WGD", !is.na(Origin.of.WGDs)),  x="motif.categories", 
                     y="expression.correl.coeff", fill="Origin.of.WGDs", facet.by = "Origin.of.WGDs")+ 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test") +
  ylab("Correlation coefficient (r)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')

WGD_coexp_pi_het <- summary_table %>% 
  filter(!is.na(Origin.of.WGDs) & (motif.categories=="HM" | motif.categories=="HM&HET")) %>%
  filter(!is.na(expression.correl.coeff)) %>%
  filter(!is.na(pid)) %>%
  arrange(pid) %>%
  mutate(window_pi = cut_interval(pid,6)) %>%
  group_by(window_pi, motif.categories, Duplication, Origin.of.WGDs) %>%
  summarise(mean.coexp = mean(expression.correl.coeff, na.rm=T),
            median.coexp = median(expression.correl.coeff, na.rm=T),
            ci = 1.96*sd(expression.correl.coeff, na.rm=T)/sqrt(n()),
            mean_interval = mean(pid),
            npoints=n())

FigS13C <- ggplot(WGD_coexp_pi_het, aes(x=mean_interval,y=mean.coexp, col=as.factor(motif.categories)))+
  geom_point(shape=15, size=4, position=position_dodge(.3))+ facet_grid(.~Origin.of.WGDs) +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position = c(0,0.9)) +
  ylab("Correlation coefficient (r)") +
  xlab("Pairwise amino acid sequence identity (%)") +
  geom_errorbar(aes(ymin=mean.coexp-ci, ymax=mean.coexp+ci), width=0.2, size=0.5,
                position=position_dodge(.3))+
  guides(color=guide_legend(""))+
  geom_abline(intercept = 0.0, slope = 0, color="grey", linetype="dotted")


#FigS13D
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}

WGD <- filter(summary_table, Duplication=="WGD", !is.na(Origin.of.WGDs))
WGD$int <- interaction(WGD$motif.categories, WGD$Origin.of.WGDs)

my_comparisons <- list(c("HM.Homeologs", "HM&HET.Homeologs"), c("HM.True_ohnologs", "HM&HET.True_ohnologs"))

FigS13D <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.TF, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("Similarity of transcription factors") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=0.7)


FigS13E<-WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.cell.comp.P1P2, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("GO cellular component similarity") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=1.1)

FigS13F  <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.loc, fill=as.factor(Origin.of.WGDs)))+
  scale_fill_manual(values = c(cl[490],cl[430])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("Similarity of localization") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.Homeologs" = "HM", "HM&HET.Homeologs" = "HM&HET", "HM.True_ohnologs" = "HM", "HM&HET.True_ohnologs"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test", label.y=1.1, size=4.5)


ggsave(file="sup_fig/FigS13.pdf", width=14, height=21, dpi=500)
plot_grid(FigS13A, FigS13B, FigS13C, FigS13D, FigS13E, FigS13F, 
          labels=c("A", 'B', 'C', "D", "E", "F"), nrow=3, align=c("h"))
dev.off()


#FigS17
df.med.seuil <- read.table("output/PCA_med_seuil_data_2019_02.tab", sep='\t', header=T)
FigureS17 <- ggplot(df.med.seuil, aes(x=Zscore), color="blue") + 
  geom_density() + geom_vline(aes(xintercept=2.5),
                              color="red", linetype="dashed", size=1) +
  labs(x=("z-score"), y=("Frequency")) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))

ggsave(file="sup_fig/FigS17.pdf", width=7, height=7, dpi=500)
FigureS17
dev.off()
