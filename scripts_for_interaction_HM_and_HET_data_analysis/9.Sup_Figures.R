#############################################################################
#                    Supplementary figures script                           #
#############################################################################

#Script to generate supplementary figures from 6.HM_compiled_data.R and 
#8.HM_HET_expression_duplication.r scripts output

rm(list=ls())

setwd("dir")

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

summary_table <- read.delim('without_low_qual/TableS1.csv', header=T, sep="\t")

#FigureS1
#Comparison of PCA data generated in this study with previously published data.

Interact <- read.table("without_low_qual/comp_Biogrid.csv", sep="\t",  header=T)

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


HM_Stynen_our_PCA <- read.table("without_low_qual/comp_Stynen.csv", sep="\t",  header=T)

Michnick <- 
  ggplot(HM_Stynen_our_PCA, aes(x=med_MTX_our_PCA, y=med_MTX_Stynen)) + 
  geom_point() + geom_smooth(method=lm) +
  ylab("Median colony size \n (pixels, log2) from Stynen") + 
  xlab("Median colony size (pixels, log2) from this study") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  annotation_custom(grobTree(textGrob("Spearman: r = 0.59 p-value < 2.2e-16", 
                                      x=0.05,  y=0.98, hjust=0, gp=gpar(fontsize=12))))

HM_Tarassov_our_PCA <- read.table("without_low_qual/comp_Tarassov.csv", sep="\t",  header=T)

Tarassov <-
  ggplot(HM_Tarassov_our_PCA, aes(x=med_MTX_our_PCA, y=log2(intensity_Tarassov))) + 
  geom_point() + geom_smooth(method=lm) +
  ylab("Integrated colony pixel \n intensity (log2) from Tarassov") + 
  xlab("Median colony size (pixels, log2) from this study") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  annotation_custom(grobTree(textGrob("Spearman: r = 0.49 p-value = 1.58e-05", 
                                      x=0.05,  y=0.98, hjust=0, gp=gpar(fontsize=12))))


ggsave(file="without_low_qual/FigureS1.pdf", width=7, height=15, dpi=500)
plot_grid(Biogrid, Michnick, Tarassov, labels = c("A", "B", "C"), nrow=3, ncol=1)
dev.off()


######FigureS2
#Intersections of detected HMs and HETs from this study and previously reported HMs and HETs.

#FigureS2A (HM)
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


upset(HM,  sets = colnames(HM)[2:15], sets.bar.color = "pink", point.size = 1.5, order.by = "freq",
      mainbar.y.label = "Number of HM interactions", sets.x.label = "HMs per method")


#FigureS2B (HET)
HET <- select(summary_table, pair, HET.ourPCA, list.bg.HET.Pubmed.id, HET.PDB)
colnames(HET) <- c('ppi', 'PCA_this_study', 'list.bg', 'PDB')
HET$PCA_this_study[is.na(HET$PCA_this_study)] <- 0
HET$PDB[is.na(HET$PDB)] <- 0


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

  upset(HET,  sets = colnames(HET)[2:13], sets.bar.color = "purple", point.size = 1.5, order.by = "freq",
        mainbar.y.label = "Number of HET interactions", sets.x.label = "HETs per method")
  
  
#####FigureS3
#Association between expression and the proability of HM detection by PCA in this study.

#FigureS3A
#The plot shows the detection probability of HMs  as a function of mRNA abundance for previously reported HMs.
HM1.PCA <- select(summary_table, P1, HM1.PCA, log10.RNAseq.MTX.P1, HM1bg.S.PDB.Kim)
HM2.PCA <- select(summary_table, P2, HM2.PCA, log10.RNAseq.MTX.P2, HM2bg.S.PDB.Kim)
colnames(HM1.PCA) <- c("paralog", "HM", "log10_RNAseq_MTX", "HMbg.S.PDB.Kim")
colnames(HM2.PCA) <- c("paralog", "HM", "log10_RNAseq_MTX", "HMbg.S.PDB.Kim")
HM <- rbind(HM1.PCA, HM2.PCA)
HM <- HM %>% .[complete.cases(.), ]
HM <- filter(HM, HMbg.S.PDB.Kim==1)
lines <- (ksmooth(HM$log10_RNAseq_MTX, HM$HM, "normal", 
                  bandwidth = 0.8, x.points = seq(0, 5, by=0.001))) %>% as.data.frame()

FigureS3A <-
  ggplot(HM, aes(x=log10_RNAseq_MTX, y=HM)) + geom_point() + geom_line(data=lines, aes(x=x, y=y))+
  xlab("mRNA abundance of HM previously reported \n log10 (number of mapped RNAseq reads in MTX)") +  ylab("HM detection probability") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, lineheight=0.4, hjust=0.5), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.1, 0.95), 
        legend.title = element_text(size=15), legend.text = element_text(size=12))

#FigureS3B
#To show that difference in HM formation between paralogs results in part from their differential mRNA abundance,
#the PCA score of paralog 1 (P1) is compared to the PCA score of paralog 2 (P2)

df <- filter(summary_table, !is.na(med.MTX.HM.P1) & !is.na(med.MTX.HM.P2)) %>% 
  mutate(diff_expP2P1 = log10.RNAseq.MTX.P2-log10.RNAseq.MTX.P1,
         tot_exp_P1P2 = log10.RNAseq.MTX.P2+log10.RNAseq.MTX.P1)

FigureS3B <- ggplot(df, aes(x=med.MTX.HM.P1, y=med.MTX.HM.P2, color=diff_expP2P1))+
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


#FigureS3C
#Comparison of expression levels of previously reported HMs for 
#HMs undetected and detected in the PCA experiments performed in this study.

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
#p-value =  4.029e-06
wilcox.test(HM_rep_undec$log10_abundance_paxDB, HM_rep_dec$log10_abundance_paxDB)
#p-value = 9.823e-05

HM_rep_tested <- rbind(HM_rep_dec, HM_rep_undec)

my_comparisons <- list(c("0", "1"))


library(ggpubr)

FigureS3C <-
  ggboxplot(HM_rep_tested,  x="HM.PCA", y="log10_RNAseq_MTX")+ 
  stat_compare_means(comparisons = my_comparisons, group.by = "HM.PCA", method = "wilcox.test") +
  ylab("log10 (number of mapped \n RNAseq reads in MTX)") + xlab("Detection of known HM") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')


library("cowplot")
ggsave(file="without_low_qual/FigureS3.pdf", width=5, height=15, dpi=500)
plot_grid(FigureS3A, FigureS3B, FigureS3C, labels = c("A", "B", "C"), nrow=3, ncol=1)
dev.off()


#FigureS4
#Comparison of expression between singletons and duplicates
dff.S.bg.Kim.PCA <- read.table('without_low_qual/HM_expression.csv', header=T, sep='\t')
dff.S.bg.Kim.PCA$pair <- apply(cbind(as.character(dff.S.bg.Kim.PCA$orf), 
                                     as.character(dff.S.bg.Kim.PCA$prey)), 
                               1, function(x) paste(sort(x), collapse="."))
dff.S.bg.Kim.PCA <- left_join(dff.S.bg.Kim.PCA, wgd_filter, by='pair')
dff.S.bg.Kim.PCA <- filter(dff.S.bg.Kim.PCA, is.na(wgd.filter))

my_comparisons <- list( c("S", "ssd"), c("S", "wgd"), c("ssd", "wgd"))

#FigureS4A : mRNA level
FigureS4A <-
  ggplot(filter(dff.S.bg.Kim.PCA,type_para!="ssd_wgd") , aes(x=type_para, y=log10_RNA_med_MTX, fill=type_para)) + 
  geom_boxplot()  +
  scale_fill_manual(guide=FALSE,values = c(cl[16],cl[144],cl[129])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  scale_x_discrete(labels=c(c("S","SSD","WGD"))) +
  ylab("mRNA level in MTX") + xlab("Groups of genes")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0.08, 0.9), legend.title = element_text(size=15), 
        legend.text = element_text(size=12))

#FigureS4B : protein level
FigureS4B <-
  ggplot(filter(dff.S.bg.Kim.PCA, type_para!="ssd_wgd" & type_para!="ssd-successiv") , 
         aes(x=type_para, y=log(paxdb), fill=type_para)) + 
  geom_boxplot()  +
  scale_fill_manual(guide=FALSE,values = c(cl[16],cl[144],cl[129])) +
  scale_x_discrete(labels=c(c("S","SSD","WGD"))) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  ylab("Protein abundance (log10)") + xlab("Groups of genes")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0.08, 0.9), legend.title = element_text(size=15), 
        legend.text = element_text(size=12))


ggsave(file="without_low_qual/FigureS4.pdf", width=10, height=5, dpi=500)
plot_grid(FigureS4A, FigureS4B,
          labels = c("A", "B"), nrow=1, ncol = 2, label_size = 12)
dev.off()

#######FigureS5#####
#Linl interaction motifs and percentage of pairwise amino acid sequence identity between paralogs

#FigureS5A
#Pairs of paralogs were clustered in 6 pairwise amino acid sequence identity groups and 
#their distribution (in percentage) of these groups were compared between SSDs and WGDs
proportion_pid <- summary_table %>%
  filter(!is.na(phylomDB.similarity)) %>%
  arrange(phylomDB.similarity) %>%
  mutate(window_pi = cut_interval(phylomDB.similarity,6)) %>%
  group_by(window_pi, Duplication) 


tot <- proportion_pid %>% group_by(Duplication) %>% filter(., !is.na(window_pi)) %>%  
  summarise(tot = length(pair))
n <- proportion_pid %>% group_by(Duplication, window_pi) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(window_pi)) %>% as.data.frame()

freq <- full_join(n, tot, by="Duplication") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n

#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

diff_SSD_WGD_per_window_pi <- freq %>% group_by(window_pi) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
diff_SSD_WGD_per_window_pi

#result
#window_pi p.value_fisherTest
#1 [5.09,20.9]       3.400104e-15
#2 (20.9,36.7]       2.555369e-02
#3 (36.7,52.5]       8.106963e-05
#4 (52.5,68.4]       1.857072e-03
#5 (68.4,84.2]       1.251342e-05
#6  (84.2,100]       8.871017e-02

FigureS5A <-
  ggplot(freq, aes(window_pi, freq)) + 
  geom_bar(aes(fill = Duplication), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  labs(x=("Pairwise amino acid sequence identity (%)"), y=("Percentage (%)"), title = ("Duplication")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75),
                              xend=c(1.25, 2.25, 3.25, 4.25, 5.25, 6.25), 
                              y=c(41, 34.5, 27.5, 17, 15, 12), 
                              annotation=c("3.40e-15", "2.55e-02", '8.11e-05', '1.86e-03', '1.25e-05', '8.87e-02')),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.15, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.background = element_blank())



#FigureS5B
#Percentage of the number of paralog pairs forming HMs and HETs among the total number of paralog pairs
#forming HM (HM and HM&HET) is shown as a function of the percentage of pairwise amino acid sequence identity
HM.HET_pi <- summary_table %>% filter(motif.categories=="HM" | motif.categories=="HM&HET") %>%
  filter(!is.na(phylomDB.similarity)) %>%
  arrange(phylomDB.similarity) %>%
  mutate(window_pi = cut_interval(phylomDB.similarity,6),
         HM = ifelse(motif.categories=="HM", 1, 0),
         HM.HET = ifelse(motif.categories=="HM&HET", 1, 0),
         dupli = ifelse(Duplication=='SSD', 'SSD',
                 ifelse(Duplication=='WGD' & Origin.of.WGDs=='Homeologs', 'Homeologs',
                 ifelse(Duplication=='WGD'& Origin.of.WGDs=='True_ohnologs', 'True_ohnologs', 'NA')))) %>%
  group_by(window_pi, Duplication) %>%
  summarise(percentHM.HET = sum(HM.HET)/(sum(HM)+sum(HM.HET))*100,
            HM=sum(HM),
            HM.HET=sum(HM.HET),
            percent.HM.HET=100*sum(HM.HET)/(sum(HM)+sum(HM.HET)),
            effectif=paste(c(sum(HM.HET),(sum(HM)+sum(HM.HET))), collapse = "/")) %>% as.data.frame()


FigureS5B <-
  ggplot(HM.HET_pi, aes(x=as.factor(window_pi), y=as.numeric(percent.HM.HET),  fill=Duplication)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144], cl[129])) +
  ylab("100 * HM&HET / HM+HM&HET") + xlab("Pairwise amino acid sequence identity (%)") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  geom_text(aes(label=effectif), vjust=-0.3, size=3.5)

#FigureS5C
#Percentage of pairwise amino acid sequence identity between paralogs for each motif number 
#1HM: shows one homomer only
#2HM: shows both homomers
#1HM.HET: shows one homomer and the heteromer
#2HM.HET: shows both homomers and the heteromer

my_comparisons <- list( c("1HM", "2HM"), c("1HM", "1HM.HET"), c("1HM", "2HM.HET"),
                        c("2HM", "1HM.HET"), c("2HM", "2HM.HET"),
                        c("1HM.HET", "2HM.HET"))

summary_table$motif.number<-factor(summary_table$motif.number, levels=c("1HM", "2HM", "1HM.HET", "2HM.HET"))

FigureS5C <- ggboxplot(filter(summary_table, motif.number=="1HM" | motif.number=="2HM" | motif.number=="1HM.HET" | motif.number=="2HM.HET"),  
                       x="motif.number", y="phylomDB.similarity", fill="Duplication", facet.by = "Duplication")+ 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')


#FigureS5D
#The percentage of pairwise amino acid sequence identity among true onhologs and homeologs
my_comparisons <- list(c("Homeologs","True_ohnologs"))
summary_table$Origin.of.WGDs <- factor(summary_table$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))


FigureS5D <-
  ggplot(filter(summary_table, !is.na(Origin.of.WGDs)), aes(x=Origin.of.WGDs, y=phylomDB.similarity, fill=Origin.of.WGDs)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",  size = 4) +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Origin of WGDs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')

#FigureS5E
#Percentage of pairwise amino acid  pairwise sequence  identity between paralogs for HM and HM&HET each motifs for true ohnologs and homeologs
FigureS5E <- ggboxplot(filter(summary_table, Duplication=="WGD" & 
                               !is.na(Origin.of.WGDs) & 
                               (motif.categories=='HM' | motif.categories=='HM&HET')),  
                        x="motif.categories", y="phylomDB.similarity", fill="Origin.of.WGDs", facet.by = "Origin.of.WGDs")+ 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = c("HM", "HM&HET"), group.by = "Origin.of.WGDs", method = "wilcox.test") +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')
wilcox.test(summary_table$phylomDB.similarity[summary_table$motif.categories=='HM' & summary_table$Origin.of.WGDs=='Homeologs'],
            summary_table$phylomDB.similarity[summary_table$motif.categories=='HM&HET' & summary_table$Origin.of.WGDs=='Homeologs']) 
#p-value = 0.0002377
wilcox.test(summary_table$phylomDB.similarity[summary_table$motif.categories=='HM' & summary_table$Origin.of.WGDs=='True_ohnologs'],
            summary_table$phylomDB.similarity[summary_table$motif.categories=='HM&HET' & summary_table$Origin.of.WGDs=='True_ohnologs']) 
#p-value = 0.04385

FigureS5DE <- plot_grid(FigureS5D, FigureS5E, labels = c("D", "E"), ncol=2, nrow = 1)


ggsave(file="without_low_qual/FigureS5.test.pdf", width=15, height=10, dpi=500)
plot_grid(FigureS5A, FigureS5B, FigureS5C, FigureS5DE, labels = c("A", "B", "C", ''), ncol=2, nrow = 2)
dev.off()

########FigureS8########
#Comparison of functional similarity between HM and HM&HET pairs with more details than Figure 3
#(all pairs represented)
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}


summary_table$int <- interaction(summary_table$motif.categories, summary_table$Duplication)
my_comparisons <- list(c("HM.SSD", "HM&HET.SSD"), c("HM.WGD", "HM&HET.WGD"))


#GO molecular functions
FigureS8A <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.mol.fct.P1P2*100, fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=110, size=4.5) +
  ylab("GO molecular function similarity (%)") + xlab("Interaction motifs")+
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

#GO biological processes
FigureS8B <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET")  %>% 
  ggplot(., aes(x=as.factor(int), y=sim.bio.proc.P1P2*100, fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+ 
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("GO biological process similarity (%)") +  xlab("Interaction motifs") +
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=110, size=4.5)

#growth phenotypes
FigureS8C <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.pheno.P1P2*100, fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("Phenotype similarity (%)") +  xlab("Interaction motifs") + 
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(strip.text.x = element_text(size = 16,face='bold'))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=110, size=4.5)#+

#correlation of genetic interaction profiles
FigureS8D <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
  ggplot(., aes(x=as.factor(int), y=as.numeric(med.gi.cor), fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=0.55, size=4.5) +
  ylab("Correlation of genetic interaction profile") + xlab("Interaction motifs") +
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")


ggsave(file="without_low_qual/FigureS8.pdf", width=10, height=10, dpi=500)
plot_grid(FigureS8A, FigureS8B, FigureS8C, FigureS8D, labels=c("A", "B", "C", "D"), nrow=2, align=c("h"))
dev.off()



#FigureS9
#Same than figure S9 among WGDs for homeologs and true onhologs
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}

WGD <- filter(summary_table, Duplication=="WGD", !is.na(Origin.of.WGDs))
WGD$int <- interaction(WGD$motif.categories, WGD$Origin.of.WGDs)

my_comparisons <- list(c("HM.Homeologs", "HM&HET.Homeologs"), c("HM.True_ohnologs", "HM&HET.True_ohnologs"))

#GO molecular functions
FigureS9A <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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

#GO biological processes 
FigureS9B <-WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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

#growth phenotypes
FigureS9C <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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

#correlation of genetic interaction profiles
FigureS9D  <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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


ggsave(file="without_low_qual/FigureS9.pdf", width=10, height=10, dpi=500)
plot_grid(FigureS9A, FigureS9B, FigureS9C, FigureS9D, labels=c("A", "B", "C", "D"), nrow=2, align=c("h"))
dev.off()

#FigureS10
##Function similarity could be due to sequence similarity
#We checked functional similarity between paralogs as a function of their pairwise amino acid sequence identity
#to control if the functional effect observed on motif interaction in figure S8 wasn't due just for sequence similarity

funct <- filter(summary_table, motif.categories !="NI" & motif.categories !="HET") %>%
  select(., pair, Duplication, motif.categories, sim.mol.fct.P1P2, sim.bio.proc.P1P2, sim.pheno.P1P2, med.gi.cor, sim.cell.comp.P1P2, sim.loc, sim.TF, phylomDB.similarity)
colnames(funct) <- c('pair', 'Duplication', 'motif.categories', 'Molecular.function', 
                     'Biological.process', 'Phenotype', 'med.gi.cor', "cellular.comp", "localization", "transcription.Factor", 'phylomDB.similarity')

funct %<>% gather('Molecular.function', 'Biological.process', 'Phenotype',  'med.gi.cor', "cellular.comp", "localization", "transcription.Factor", key=GO, value= similarity)

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), ...)
}

#GO molecular functions
plot1 <- ggplot(filter(funct, GO=='Biological.process', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories)) + facet_grid(.~Duplication) +
  geom_point() + ylab("GO biological process similarity (%)") + xlab("Pairwise amino acid sequence identity (%)") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) + binomial_smooth()

#GO biological processes
plot2 <- ggplot(filter(funct, GO=='Molecular.function', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories))  + facet_grid(.~Duplication) +
  geom_point() +   ylab("GO molecular function similarity (%)") + xlab("Pairwise amino acid sequence identity (%)") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) +binomial_smooth()

#growth phenotypes
plot3 <- ggplot(filter(funct, GO=='Phenotype', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories))  + facet_grid(.~Duplication) +
  geom_point() +   ylab("Phenotype similarity (%)") + xlab("Pairwise amino acid sequence identity (%)") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) +binomial_smooth()

#correlation of genetic interaction profiles
plot4 <- ggplot(filter(funct, GO=='med.gi.cor', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories)) +
  geom_point() +   ylab("Correlation of genetic interaction profile") + xlab("Pairwise amino acid sequence identity (%)")  + facet_grid(.~Duplication) +
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") + 
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) + geom_smooth(method=lm)

ggsave(file="without_low_qual/FigureS10.pdf", width=10, height=10, dpi=500)
plot_grid(plot1, plot2, plot3, plot4, labels=c("A", "B", "C", "D"), nrow=2, align=c("h"))
dev.off()

#GO cellular component
plot5 <- ggplot(filter(funct, GO=='cellular.comp', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories)) + facet_grid(.~Duplication) +
  geom_point() + ylab("GO cellular component similarity (%)") + xlab("Pairwise amino acid sequence identity (%)2") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) + binomial_smooth()

#Localization
plot6 <- ggplot(filter(funct, GO=='localization', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories))  + facet_grid(.~Duplication) +
  geom_point() +   ylab("Localization (%)") + xlab("Pairwise amino acid sequence identity (%)") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) +binomial_smooth()

#Transcription factor
plot7 <- ggplot(filter(funct, GO=='transcription.Factor', !is.na(similarity)), aes(x=phylomDB.similarity, y=similarity, color=motif.categories))  + facet_grid(.~Duplication) +
  geom_point() +   ylab("Transcription Factor (%)") + xlab("Pairwise amino acid sequence identity (%)") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) +binomial_smooth()

ggsave(file="without_low_qual/FigureS20.pdf", width=15, height=5, dpi=500)
plot_grid(plot5, plot6, plot7, labels=c("A", "B", "C", "D"), nrow=1, align=c("h"))
dev.off()


#To test if they are an effect of function independently to sequence similarity
#we processed to GLM tests:

#For SSD
SSD.BP <- filter(funct, Duplication=='SSD' & GO=='Biological.process')
glm1.SSD.BP<-glm(SSD.BP$motif.categories ~ 
                   SSD.BP$similarity+
                   SSD.BP$phylomDB.similarity, family=binomial)

summary(glm1.SSD.BP)
export_summs(glm1.SSD.BP)
anova(glm1.SSD.BP,test="LRT") 


SSD.MF <- filter(funct, Duplication=='SSD' & GO=='Molecular.function')
glm1.SSD.MF<-glm(SSD.MF$motif.categories ~ 
                   SSD.MF$similarity+
                   SSD.MF$phylomDB.similarity, family=binomial)

summary(glm1.SSD.MF)
export_summs(glm1.SSD.MF)
anova(glm1.SSD.MF,test="LRT") 

SSD.Phe <- filter(funct, Duplication=='SSD' & GO=='Phenotype')
glm1.SSD.Phe<-glm(SSD.Phe$motif.categories ~ 
                    SSD.Phe$similarity+
                    SSD.Phe$phylomDB.similarity, family=binomial)

summary(glm1.SSD.Phe)
export_summs(glm1.SSD.Phe)
anova(glm1.SSD.Phe,test="LRT") 

SSD.gi.cor <- filter(funct, Duplication=='SSD' & GO=='med.gi.cor')
glm1.SSD.gi.cor<-glm(SSD.gi.cor$motif.categories ~ 
                       SSD.gi.cor$similarity+
                       SSD.gi.cor$phylomDB.similarity, family=binomial)

summary(glm1.SSD.gi.cor)
export_summs(glm1.SSD.gi.cor)
anova(glm1.SSD.gi.cor,test="LRT") 

SSD.cel.comp <- filter(funct, Duplication=='SSD' & GO=='cellular.comp')
glm1.SSD.cel.comp<-glm(SSD.cel.comp$motif.categories ~ 
                       SSD.cel.comp$similarity+
                       SSD.cel.comp$phylomDB.similarity, family=binomial)

summary(glm1.SSD.cel.comp)
export_summs(glm1.SSD.cel.comp)
anova(glm1.SSD.cel.comp,test="LRT") 


SSD.loc <- filter(funct, Duplication=='SSD' & GO=='localization')
glm1.SSD.loc<-glm(SSD.loc$motif.categories ~ 
                         SSD.loc$similarity+
                         SSD.loc$phylomDB.similarity, family=binomial)

summary(glm1.SSD.loc)
export_summs(glm1.SSD.loc)
anova(glm1.SSD.loc,test="LRT") 

SSD.TF <- filter(funct, Duplication=='SSD' & GO=="transcription.Factor")
glm1.SSD.TF<-glm(SSD.TF$motif.categories ~ 
                    SSD.TF$similarity+
                    SSD.TF$phylomDB.similarity, family=binomial)

summary(glm1.SSD.TF)
export_summs(glm1.SSD.TF)
anova(glm1.SSD.TF,test="LRT") 

#For WGD
WGD.BP <- filter(funct, Duplication=='WGD' & GO=='Biological.process')
glm1.WGD.BP<-glm(WGD.BP$motif.categories ~ 
                   WGD.BP$similarity+
                   WGD.BP$phylomDB.similarity, family=binomial)

summary(glm1.WGD.BP)
export_summs(glm1.WGD.BP)
anova(glm1.WGD.BP,test="LRT") 

WGD.MF <- filter(funct, Duplication=='WGD' & GO=='Molecular.function')
glm1.WGD.MF<-glm(WGD.MF$motif.categories ~ 
                   WGD.MF$similarity+
                   WGD.MF$phylomDB.similarity, family=binomial)

summary(glm1.WGD.MF)
export_summs(glm1.WGD.MF)
anova(glm1.WGD.MF,test="LRT") 

WGD.cel.comp <- filter(funct, Duplication=='WGD' & GO=='cellular.comp')
glm1.WGD.cel.comp<-glm(WGD.cel.comp$motif.categories ~ 
                         WGD.cel.comp$similarity+
                         WGD.cel.comp$phylomDB.similarity, family=binomial)

summary(glm1.WGD.cel.comp)
export_summs(glm1.WGD.cel.comp)
anova(glm1.WGD.cel.comp,test="LRT") 


WGD.Phe <- filter(funct, Duplication=='WGD' & GO=='Phenotype')
glm1.WGD.Phe<-glm(WGD.Phe$motif.categories ~ 
                    WGD.Phe$similarity+
                    WGD.Phe$phylomDB.similarity, family=binomial)

summary(glm1.WGD.Phe)
export_summs(glm1.WGD.Phe)
anova(glm1.WGD.Phe,test="LRT") 

WGD.gi.cor <- filter(funct, Duplication=='WGD' & GO=='med.gi.cor')
glm1.WGD.gi.cor<-glm(WGD.gi.cor$motif.categories ~ 
                       WGD.gi.cor$similarity+
                       WGD.gi.cor$phylomDB.similarity, family=binomial)

summary(glm1.WGD.gi.cor)
export_summs(glm1.WGD.gi.cor)
anova(glm1.WGD.gi.cor,test="LRT") 

WGD.loc <- filter(funct, Duplication=='WGD' & GO=='localization')
glm1.WGD.loc<-glm(WGD.loc$motif.categories ~ 
                    WGD.loc$similarity+
                    WGD.loc$phylomDB.similarity, family=binomial)

summary(glm1.WGD.loc)
export_summs(glm1.WGD.loc)
anova(glm1.WGD.loc,test="LRT") 

WGD.TF <- filter(funct, Duplication=='WGD' & GO=="transcription.Factor")
glm1.WGD.TF<-glm(WGD.TF$motif.categories ~ 
                   WGD.TF$similarity+
                   WGD.TF$phylomDB.similarity, family=binomial)

summary(glm1.WGD.TF)
export_summs(glm1.WGD.TF)
anova(glm1.WGD.TF,test="LRT") 

#Try with Guan et al 2007 data
Guan <- read.table('data/Guan_data')
colnames(Guan) <- c('P1', 'P2', 'Guan.similarity')
Guan$pair <- apply(cbind(as.character(Guan$P1), 
                         as.character(Guan$P2)), 
                   1, function(x) paste(sort(x), collapse="."))
summary_table <- left_join(summary_table, Guan, by='pair')

funct <- filter(summary_table, motif.categories !="NI" & motif.categories !="HET") %>%
  select(., pair, Duplication, motif.categories, sim.mol.fct.P1P2, sim.bio.proc.P1P2, sim.pheno.P1P2, med.gi.cor, Guan.similarity)
colnames(funct) <- c('pair', 'Duplication', 'motif.categories', 'Molecular.function', 
                     'Biological.process', 'Phenotype', 'med.gi.cor', 'Guan.similarity')

funct %<>% gather('Molecular.function', 'Biological.process', 'Phenotype',  'med.gi.cor', key=GO, value= similarity)

plot1 <- ggplot(filter(funct, GO=='Biological.process'), aes(x=Guan.similarity, y=similarity, color=motif.categories)) + facet_grid(.~Duplication) +
  geom_point() + ylab("GO biological process similarity (%)") + xlab("Guan similarity P1P2") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) + binomial_smooth()

plot2 <- ggplot(filter(funct, GO=='Molecular.function'), aes(x=Guan.similarity, y=similarity, color=motif.categories))  + facet_grid(.~Duplication) +
  geom_point() +   ylab("GO molecular function similarity (%)") + xlab("Guan similarity P1P2") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) +binomial_smooth()

plot3 <- ggplot(filter(funct, GO=='Phenotype'), aes(x=Guan.similarity, y=similarity, color=motif.categories))  + facet_grid(.~Duplication) +
  geom_point() +   ylab("Phenotype similarity (%)") + xlab("Guan similarity P1P2") + 
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) +binomial_smooth()

plot4 <- ggplot(filter(funct, GO=='med.gi.cor'), aes(x=Guan.similarity, y=similarity, color=motif.categories)) +
  geom_point() +   ylab("Correlation of genetic interaction profile") + xlab("Guan similarity P1P2")  + facet_grid(.~Duplication) +
  theme(axis.text.x= element_text(size=9), axis.text.y= element_text(size=9), 
        axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
        plot.title = element_text(size=12, hjust = 0.5), legend.position = "none") + 
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple")) + geom_smooth(method=lm)

plot_grid(plot1, plot2, plot3, plot4, labels=c("A", "B", "C", "D"), nrow=2, align=c("h"))

#--> similar results (figure not showed)

#FigureS11
#Distribution of interaction motifs for SSDs, WGDs and the two types of WGDs
#separating 1 homolog observation to 2 homologs

tot <- summary_table %>% group_by(Duplication) %>% 
  filter(., motif.number=='1HM' | motif.number=='2HM' | motif.number=='1HM.HET' | motif.number=='2HM.HET' ) %>%  
  summarise(tot = length(pair))
n <- summary_table %>% filter(., motif.number=='1HM' | motif.number=='2HM' | motif.number=='1HM.HET' | motif.number=='2HM.HET') %>%  
  group_by(Duplication, motif.number) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.number)) %>% as.data.frame()

freq <- full_join(n, tot, by="Duplication") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n

#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))
diff_SSD_WGD_per_motif.simple <- freq %>% group_by(motif.number) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
diff_SSD_WGD_per_motif.simple
#motif.number p.value_fisherTest
#1          1HM         0.14889317
#2      1HM.HET         0.07308352
#3          2HM         0.28018431
#4      2HM.HET         0.05410324



FigureS11A<-
  ggplot(freq, aes(as.factor(motif.number), freq)) + 
  geom_bar(aes(fill = Duplication), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  labs(x=("Interaction motifs"), y=("Percentage (%)"), title = ("Duplications")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75),xend=c(1.25, 2.25, 3.25, 4.25), 
                              y=c(45, 17.5, 23, 33.5), 
                              annotation=c("0.14", " 0.07 ", "0.28", "0.05")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.9, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.background = element_blank())

tot <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs) %>% 
  filter(., motif.number=='1HM' | motif.number=='2HM' | motif.number=='1HM.HET' | motif.number=='2HM.HET' ) %>%  
  summarise(tot = length(pair))
n <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs, motif.number) %>% 
  summarise(n = length(pair)) %>% 
  filter(., motif.number=='1HM' | motif.number=='2HM' | motif.number=='1HM.HET' | motif.number=='2HM.HET' )  %>% as.data.frame()

freq <- full_join(n, tot, by="Origin.of.WGDs") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n


#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

diff_Inffered.topo_per_motif <- freq %>% group_by(motif.number) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
#motif.number p.value_fisherTest
#1          1HM         0.03208741
#2          2HM         0.74629934
#3      1HM.HET         0.30564231
#4      2HM.HET         0.18334996

freq$Origin.of.WGDs <- factor(freq$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

FigureS11B <-
  ggplot(freq, aes(as.factor(motif.number), freq)) + 
  geom_bar(aes(fill = Origin.of.WGDs), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[430],cl[490]), labels = c("Homeologs", "True ohnologs")) + 
  labs(x=("Interaction motifs"), title = ("Origin of WGDs")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75),xend=c(1.25, 2.25, 3.25, 4.25), 
                              y=c(44, 16, 27, 44), 
                              annotation=c("0.03", "0.75", "0.31", "0.18")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.title.x = element_text(size=15), 
        axis.text.y= element_text(size=12), axis.title.y = element_text(size=15), 
        plot.title = element_text(size=15, hjust = 0.5), 
        legend.position = c(0.3, 0.94), legend.title = element_blank(), legend.text = element_text(size=12), 
        legend.background = element_blank())

ggsave(file="without_low_qual/FigureS11.pdf", width=5, height=5, dpi=500)
plot_grid(FigureS11A, FigureS11B, nrow=1, ncol=2)
dev.off()


###FigureS17


#FigureS17A
#The correlation coefficients between the expression profile of paralog pairs 
#are compared among the different interaction motifs for SSDs and WGDs
#using Ihmel et al., 2004 data
FigureS18A <-
  ggboxplot(filter(summary_table, motif.categories=='HM' | motif.categories=='HM&HET'), x="motif.categories", y="exp.correl.coeff.microarray", 
            fill="Duplication", facet.by = "Duplication", panel.labs.font = list(size = 11))+ 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  stat_compare_means(comparisons = c("HM", "HM&HET"), group.by = "Duplication", method = "t.test") +
  ylab("Correlation coefficient (r) from \n Ihmels et al.'s microarray data") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

SSD <- filter(summary_table, Duplication=='SSD')
t.test(SSD$exp.correl.coeff.microarray[SSD$motif.categories=='HM'], 
       SSD$exp.correl.coeff.microarray[SSD$motif.categories=='HM&HET']) #p-value = 0.006532

WGD <- filter(summary_table, Duplication=='WGD')
t.test(WGD$exp.correl.coeff.microarray[WGD$motif.categories=='HM'], 
       WGD$exp.correl.coeff.microarray[WGD$motif.categories=='HM&HET']) #p-value = 0.00613

#FigureS17B
#Correlation of expression based on Ihmels et al. 2004 data between paralogs forming 
#HM&HET or only HM as a function of their amino acid sequence identity
coexp_pi_het <- summary_table %>% filter(motif.categories=="HM" | motif.categories=="HM&HET") %>%
  filter(!is.na(exp.correl.coeff.microarray)) %>%
  filter(!is.na(phylomDB.similarity)) %>%
  arrange(phylomDB.similarity) %>%
  mutate(window_pi = cut_interval(phylomDB.similarity,6)) %>%
  group_by(window_pi, motif.categories, Duplication) %>%
  summarise(mean.coexp = mean(exp.correl.coeff.microarray, na.rm=T),
            median.coexp = median(exp.correl.coeff.microarray, na.rm=T),
            ci = 1.96*sd(exp.correl.coeff.microarray, na.rm=T)/sqrt(n()),
            mean_interval = mean(phylomDB.similarity),
            npoints=n())

FigureS17B <- 
  ggplot(coexp_pi_het,aes(x=mean_interval,y=mean.coexp, col=as.factor(motif.categories)))+
  geom_point(shape=15, size=4, position=position_dodge(.3))+ facet_grid(.~Duplication) +
  scale_color_manual(labels = c("HM", "HM&HET"), values = c("pink","purple"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.position = c(0,0.9)) +
  ylab("Correlation coefficient (r) from \n Ihmels et al.'s microarray data") +
  xlab("Pairwise amino acid sequence identity (%)") +
  geom_errorbar(aes(ymin=mean.coexp-ci, ymax=mean.coexp+ci), width=0.2, size=0.5,
                position=position_dodge(.3))+
  guides(color=guide_legend(""))+
  geom_abline(intercept = 0.0, slope = 0, color="grey", linetype="dotted")

ggsave(file="without_low_qual/FigureS17.pdf", width=10, height=5, dpi=500)
plot_grid(FigureS17A, FigureS17B, labels = c("A", 'B'), nrow=1, ncol = 2)
dev.off()

#FigureS18
# Interaction motifs and similarity of functions for SSDs and WGDs
# With more details than Figure 6 (all pairs shown)

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}


summary_table$int <- interaction(summary_table$motif.categories, summary_table$Duplication)
my_comparisons <- list(c("HM.SSD", "HM&HET.SSD"), c("HM.WGD", "HM&HET.WGD"))

#Similarity of GO cellular components
FigureS18A <- summary_table %>% filter(motif.categories =="HM" | motif.categories =="HM&HET") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.TF*100, fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=75, size=4.5) +
  ylab("Transcription factor similarity (%)") +  xlab("Interaction motifs") + 
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")


#Similarity of GO cellular components.
FigureS18B <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
  ggplot(., aes(x=as.factor(int), y=sim.cell.comp.P1P2*100, fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+ 
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("GO cellular component similarity (%)") +  xlab("Interaction motifs") +
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=110, size=4.5)

#Similarity of localization.
FigureS18C <- summary_table %>% filter((motif.categories =="HM" | motif.categories =="HM&HET") & !is.na(sim.loc)) %>% 
  ggplot(., aes(x=as.factor(int), y=as.numeric(sim.loc*100), fill=as.factor(Duplication)))+
  scale_fill_manual(values = c(cl[144],cl[129])) +
  geom_violin(alpha=0.5)+
  geom_jitter(width = 0.2, color="grey30", alpha=0.5) + 
  stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7, size=4.5, fontface="bold") +
  ylab("Similarity of localization (%)") +  xlab("Interaction motifs") + 
  scale_x_discrete(labels=c("HM.SSD" = "HM", "HM&HET.SSD" = "HM&HET", "HM.WGD" = "HM", "HM&HET.WGD"= "HM&HET")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test", label.y=110, size=4.5)



ggsave(file="without_low_qual/FigureS18.pdf", width=15, height=5, dpi=500)
plot_grid(FigureS18A, FigureS18B, FigureS18C, labels = c("A", 'B', "C"), nrow=1, ncol = 3)
dev.off()


###FigureS19
# Expression of whole-genome duplicates (WGDs) and consequences on interaction motifs.

my_comparisons <- list(c("Homeologs","True_ohnologs"))
summary_table$Origin.of.WGDs <- factor(summary_table$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

#Correlation coefficients between expression profiles across conditions are compared between homeologs and true ohnologs. 
FigureS19A <-
  ggplot(filter(summary_table, !is.na(Origin.of.WGDs)), aes(x=Origin.of.WGDs, y=exp.correl.coeff.RNAseq, fill=Origin.of.WGDs)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 4) +
  ylab("Correlation coefficient (r)") + xlab("Origin of WGDs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')

my_comparisons <- list(c("HM", "HM&HET"))

#Correlation coefficients (r) between the expression profiles of paralogs pair across growth conditions 
#are compared among the different interaction motifs for homeologs and ohnologs.
FigureS19B <- ggboxplot(filter(summary_table, Duplication=="WGD" & !is.na(Origin.of.WGDs) & (motif.categories=="HM" | motif.categories=="HM&HET")),  x="motif.categories", 
                     y="exp.correl.coeff.RNAseq", fill="Origin.of.WGDs", facet.by = "Origin.of.WGDs")+ 
  scale_fill_manual(values = c(cl[490],cl[430])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Origin.of.WGDs", method = "wilcox.test") +
  ylab("Correlation coefficient (r)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')

WGD_coexp_pi_het <- summary_table %>% 
  filter(!is.na(Origin.of.WGDs) & (motif.categories=="HM" | motif.categories=="HM&HET")) %>%
  filter(!is.na(exp.correl.coeff.RNAseq)) %>%
  filter(!is.na(phylomDB.similarity)) %>%
  arrange(phylomDB.similarity) %>%
  mutate(window_pi = cut_interval(phylomDB.similarity,6)) %>%
  group_by(window_pi, motif.categories, Duplication, Origin.of.WGDs) %>%
  summarise(mean.coexp = mean(exp.correl.coeff.RNAseq, na.rm=T),
            median.coexp = median(exp.correl.coeff.RNAseq, na.rm=T),
            ci = 1.96*sd(exp.correl.coeff.RNAseq, na.rm=T)/sqrt(n()),
            mean_interval = mean(phylomDB.similarity),
            npoints=n())
#Percentage of amino acid sequence identity between paralogs forming HM&HET 
#motif or only HM as a function of their correlation of expression for homeologs and ohnologs.
FigureS19C <- ggplot(WGD_coexp_pi_het, aes(x=mean_interval,y=mean.coexp, col=as.factor(motif.categories)))+
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


#FigureS19D
# Similarity of transcription factor binding sites
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}

WGD <- filter(summary_table, Duplication=="WGD", !is.na(Origin.of.WGDs))
WGD$int <- interaction(WGD$motif.categories, WGD$Origin.of.WGDs)

my_comparisons <- list(c("HM.Homeologs", "HM&HET.Homeologs"), c("HM.True_ohnologs", "HM&HET.True_ohnologs"))

FigureS19D <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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

#Similarity of GO cellular components
FigureS19E<-WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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


# Similarity of localization
FigureS19F  <- WGD %>% filter(motif.categories =="HM&HET" | motif.categories =="HM") %>% 
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


ggsave(file="without_low_qual/FigureS19.pdf", width=14, height=21, dpi=500)
plot_grid(FigureS19A, FigureS19B, FigureS19C, FigureS19D, FigureS19E, FigureS19F, 
          labels=c("A", 'B', 'C', "D", "E", "F"), nrow=3, align=c("h"))
dev.off()


#FigureS22
#Density of colony size converted to z-score
df.med.seuil <- read.table("without_low_qual/PCA_med_seuil_data_2019_02.tab", sep='\t', header=T)
FigureS22 <- ggplot(df.med.seuil, aes(x=Zscore), color="blue") + 
  geom_density() + geom_vline(aes(xintercept=2.5),
                              color="red", linetype="dashed", size=1) +
  labs(x=("z-score"), y=("Frequency")) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15))

ggsave(file="without_low_qual/FigureS22.pdf", width=7, height=7, dpi=500)
FigureS22
dev.off()


