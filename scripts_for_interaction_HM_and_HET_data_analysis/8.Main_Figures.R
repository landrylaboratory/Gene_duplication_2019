#############################################################################
#                             Main figures script                           #
#############################################################################

#Script to generate main figures from 6.HM_compiled_data.R and 
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
library("tidyverse")
library(ggpubr)
cl <- colors()[]


summary_table <- read.delim('TableS3.csv', header=T, sep="\t")

######## Figure2 #########

# Figure2A
# This figure compares percentage of homomeric proteins among singleton, SSDs, WGDs and double duplication
dff.S.bg.Kim.PCA <- read.table("TableS1.csv", sep="\t", header=T)
dff.S.bg.Kim.PCA$pair <- apply(cbind(as.character(dff.S.bg.Kim.PCA$orf), 
                             as.character(dff.S.bg.Kim.PCA$prey)), 
                       1, function(x) paste(sort(x), collapse="."))
dff.S.bg.Kim.PCA <- left_join(dff.S.bg.Kim.PCA, wgd_filter, by='pair')

#keeping successive SSDs
dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para== "ssd-successiv"] <- "ssd"
ct_hom <- table(droplevels(dff.S.bg.Kim.PCA)$type_para, dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA)

dbg <- data.frame(class=c("No HM S", "HM S", 
                          "No HM SSD", "HM SSD",
                          "No HM 2D", "HM 2D",
                          "No HM WGD", "HM WGD"),
                  count=c(ct_hom[1,1], ct_hom[1,2],
                          ct_hom[2,1],ct_hom[2,2],
                          ct_hom[3,1],ct_hom[3,2],
                          ct_hom[4,1],ct_hom[4,2]))
ct_hom <-t(ct_hom)
chisq.test(ct_hom) #p-value  2.2e-16

fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA, as.factor(dff.S.bg.Kim.PCA$type_para), workspace=2e8)
#p-value < 2.2e-16
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"]), workspace=2e8)
#p-value < 2.2e-16
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"]), workspace=2e8)
#p-value =   1.618e-05
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"]), workspace=2e8)
#p-value = 0.001695

tot <- sum(dbg$count)
dbg <- dbg %>% separate(class, c("HMstatue", "dupli"), "HM ")
dbg$HMstatue[(dbg$HMstatue)==""] <- "Yes"
dbg$HMstatue[(dbg$HMstatue)=="No "] <- "No"
dbg <- dbg %>% group_by(dupli) %>% mutate(sum_per_dupli=sum(count)) %>% as.data.frame()
dbg_percent <- filter(dbg, HMstatue=="Yes") %>% group_by(dupli) %>% 
  summarise(percentHM = (count/sum_per_dupli)*100) %>% as.data.frame()

genome <- filter(dbg, HMstatue=="Yes") %>% group_by(HMstatue) %>% 
  summarise(sum_per_genome = (sum(count)/tot)*100) %>% as.data.frame()
genome$HMstatue[(genome$HMstatue)=="Yes"] <- "Total"
colnames(genome) <- c("dupli", "percentHM")
dbg_percent <- rbind(genome, dbg_percent)

dbg_percent$n <- c(1,5,2,3,4)
dbg_percent$dupli <- factor(dbg_percent$dupli, levels = dbg_percent$dupli[order(dbg_percent$n)]) 


Figure2A <-
  ggplot(data=dbg_percent, aes(x=dupli, y=percentHM, fill=dupli)) + 
  geom_bar(stat="identity") + geom_text(aes(label=paste(round(percentHM,2))), vjust=-0.3, size=4) +
  scale_fill_manual(guide=FALSE,values = c(cl[285],cl[16],cl[144],cl[129],cl[310])) +
  geom_signif(xmin=c(2, 2, 2),
              xmax=c(3, 4, 5), 
              y_position=c(42, 46, 50), 
              annotation=c("<2.0e-16","1.6e-05","1.7e-03"), textsize=4) +
  ylab("Percentage of homomers (%)") + xlab("Groups of genes")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0.08, 0.9), legend.title = element_text(size=15), 
        legend.text = element_text(size=12))


# Figure2B, C, D
# Those figures came from analyse of orthologs HM formation 
# (see folder scripts_for_HM_ancestral_tests)
Figure2BC <- ggdraw() + draw_image("data/Figure2BC_Diana.png")
Figure2D <- ggdraw() + draw_image("data/Figure2D.png")

# Figure2E
#This figure compare the percentage of HM and HM&HET 
#pannel 1: comparition SSDs / WGDs
tot <- summary_table %>% group_by(Duplication) %>% filter(., !is.na(motif.categories)) %>%  
  summarise(tot = length(pair))
n <- summary_table %>% group_by(Duplication, motif.categories) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.categories)) %>% as.data.frame()

freq <- full_join(n, tot, by="Duplication") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n

#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))

diff_SSD_WGD_per_motif.simple <- freq %>% group_by(motif.categories) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
diff_SSD_WGD_per_motif.simple

#result
#motif.categories p.value_fisherTest
#1              HET        1.000000000
#2               HM        0.084312003
#3           HM&HET        0.002504755
#4               NI        0.107998862

freq <- filter(freq, motif.categories=='HM&HET' | motif.categories=='HM')

Figure2E1 <-
  ggplot(freq, aes(motif.categories, freq)) + 
  geom_bar(aes(fill = Duplication), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  labs(x=("Interaction motifs"), y=("Percentage (%)"), title = ("Duplication")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75),xend=c(1.25, 2.25), 
                              y=c(45, 45), 
                              annotation=c("0.08", "2.50e-03")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.15, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.background = element_blank())

#pannel 2: comparition homeologs / true ohnologs

tot <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs) %>% filter(., !is.na(motif.categories)) %>%  
  summarise(tot = length(pair))
n <- filter(summary_table, !is.na(Origin.of.WGDs)) %>% group_by(Origin.of.WGDs, motif.categories) %>% 
  summarise(n = length(pair)) %>% filter(., !is.na(motif.categories)) %>% as.data.frame()

freq <- full_join(n, tot, by="Origin.of.WGDs") %>% as.data.frame()
freq$freq <- (freq$n/freq$tot)*100
freq$tot_n <- freq$tot-freq$n


#stat
colstonumeric = c(3:6)
freq[,colstonumeric] = apply(freq[,colstonumeric], 2, function(x) as.numeric(as.character(x)))
diff_Inffered.topo_per_motif <- freq %>% group_by(motif.categories) %>% 
  summarise(p.value_fisherTest = (fisher.test(matrix(c(n, tot_n), nrow = 2)))$p.value) %>% as.data.frame()
diff_Inffered.topo_per_motif

#result:
#motif.categories p.value_fisherTest
#1              HET         0.40318878
#2               HM         0.07316155
#3           HM&HET         0.01652867
#4               NI         0.79815288


freq$Origin.of.WGDs <- factor(freq$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

freq <- filter(freq, motif.categories=='HM&HET' | motif.categories=='HM')

Figure2E2 <-
  ggplot(freq, aes(motif.categories, freq)) + 
  geom_bar(aes(fill = Origin.of.WGDs), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[430],cl[490]), labels = c("Homeologs", "True ohnologs")) + 
  labs(x=("Interaction motifs"), title = ("Origin of WGDs")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75),xend=c(1.25, 2.25), 
                              y=c(43, 54), 
                              annotation=c("0.07", "1.65e-02")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), #axis.text.y= element_blank(), 
        axis.title.x = element_text(size=15), #axis.title.y = element_blank(), axis.ticks.y =  element_blank(),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.3, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12), 
        legend.background = element_blank())


Figure2E <- plot_grid(Figure2E1, Figure2E2, labels = c("", ""), align = 'h')


# Figure 2F
# Comparison of pairwise amino acid sequence identity between paralogs for HM and HM&HET
# For SSD and WGD

my_comparisons <- list( c("HM", "HM&HET"))

Figure2F <-
ggplot(filter(summary_table, motif.categories=="HM"| motif.categories=="HM&HET"), 
          aes(x=motif.categories, y=phylomDB.similarity, fill=Duplication))+ 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none') +
  geom_violin() +  facet_grid(. ~ Duplication) + 
  geom_jitter(shape=16, color="#999999", position=position_jitter(0.2)) + geom_boxplot(width=0.2) 

SSD <- filter(summary_table, Duplication=='SSD')
wilcox.test(SSD$phylomDB.similarity[SSD$motif.categories=='HM'] , SSD$phylomDB.similarity[SSD$motif.categories=='HM&HET'])
WGD <- filter(summary_table, Duplication=='WGD')
wilcox.test(WGD$phylomDB.similarity[WGD$motif.categories=='HM'] , WGD$phylomDB.similarity[WGD$motif.categories=='HM&HET'])
#p-value = 0.0001324
Figure2_top <- plot_grid(Figure2A, Figure2BC, Figure2D, labels = c("A", "", "D"), ncol=3)
Figure2_bottom <- plot_grid(Figure2E, Figure2F, labels=c('E', 'F'), ncol=2)

ggsave(file="Figure2.pdf", width=14, height=14, dpi=500)
plot_grid(Figure2_top, Figure2_bottom, ncol=1)
dev.off()


#Figure 3
#This figure linkes loss of heteromerization between paralogs with an increase of functional divergence
funct <- filter(summary_table, motif.categories !="NI" & motif.categories !="HET") %>%
  select(., Duplication, motif.categories, sim.mol.fct.P1P2, sim.bio.proc.P1P2, sim.pheno.P1P2, med.gi.cor)
colnames(funct) <- c('Duplication', 'motif.categories', 'Molecular.function', 'Biological.process', 'Phenotype', 'Cor.genetic.interaction')
funct$Cor.genetic.interaction <- funct$Cor.genetic.interaction*10

funct %<>% gather('Molecular.function', 'Biological.process', 'Phenotype', 'Cor.genetic.interaction', key=GO, value= similarity)

funct_resum <- funct %>% group_by(Duplication, GO, motif.categories) %>% 
  summarise(mean.similarity = mean(similarity, na.rm=TRUE)*100) %>% as.data.frame()

funct_resum %<>% mutate(p.value = ifelse(Duplication=='SSD' & GO=='Molecular.function', 
               wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='Molecular.function' & funct$motif.categories=='HM'],
                           funct$similarity[funct$Duplication=='SSD' & funct$GO=='Molecular.function' & funct$motif.categories=='HM&HET'])$p.value,
               ifelse(Duplication=='SSD' & GO=='Biological.process', 
                wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='Biological.process' & funct$motif.categories=='HM'],
                            funct$similarity[funct$Duplication=='SSD' & funct$GO=='Biological.process' & funct$motif.categories=='HM&HET'])$p.value,
               ifelse(Duplication=='SSD' & GO=='Phenotype', 
                wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='Phenotype' & funct$motif.categories=='HM'],
                            funct$similarity[funct$Duplication=='SSD' & funct$GO=='Phenotype' & funct$motif.categories=='HM&HET'])$p.value,
                ifelse(Duplication=='SSD' & GO=='Cor.genetic.interaction', 
                       wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='Cor.genetic.interaction' & funct$motif.categories=='HM'],
                                   funct$similarity[funct$Duplication=='SSD' & funct$GO=='Cor.genetic.interaction' & funct$motif.categories=='HM&HET'])$p.value,
              ifelse(Duplication=='WGD' & GO=='Molecular.function', 
               wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='Molecular.function' & funct$motif.categories=='HM'],
                           funct$similarity[funct$Duplication=='WGD' & funct$GO=='Molecular.function' & funct$motif.categories=='HM&HET'])$p.value,
              ifelse(Duplication=='WGD' & GO=='Biological.process', 
               wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='Biological.process' & funct$motif.categories=='HM'],
                           funct$similarity[funct$Duplication=='WGD' & funct$GO=='Biological.process' & funct$motif.categories=='HM&HET'])$p.value,
              ifelse(Duplication=='WGD' & GO=='Phenotype', 
               wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='Phenotype' & funct$motif.categories=='HM'],
                          funct$similarity[funct$Duplication=='WGD' & funct$GO=='Phenotype' & funct$motif.categories=='HM&HET'])$p.value,
               ifelse(Duplication=='WGD' & GO=='Cor.genetic.interaction', 
                      wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='Cor.genetic.interaction' & funct$motif.categories=='HM'],
                                  funct$similarity[funct$Duplication=='WGD' & funct$GO=='Cor.genetic.interaction' & funct$motif.categories=='HM&HET'])$p.value, NA)))))))))

funct_resum %<>% mutate(p.value.range = ifelse(p.value > 0.05, 'NS',
                                        ifelse(p.value < 0.05 & p.value > 0.01, '< 0.05',
                                        ifelse(p.value < 0.01 & p.value > 0.001, '< 0.01',
                                        ifelse(p.value < 0.001, '< 0.001', 'pb')))))

ggsave(file="Figure3.pdf", width=10, height=5, dpi=500)
ggplot(funct_resum, 
       aes(motif.categories, mean.similarity, group = GO, color = as.factor(GO))) +
  geom_point(size=4) +
  geom_line(color = 'grey20', alpha = 4/10, aes(linetype=p.value.range)) + facet_grid(. ~ Duplication) +
  labs(x = 'Interaction motifs', y ='Similarity (%)', color='') 
dev.off()




#######Figure6#######

#Figure6A
#This figure shows that The loss of heteromerization between paralogs may result from regulatory divergence

#Figure6A
#The correlation coefficients between the expression profile of paralog pairs 
#are compared among the different interaction motifs for SSDs and WGDs
Figure6A <-
  ggboxplot(filter(summary_table, motif.categories=='HM' | motif.categories=='HM&HET'), x="motif.categories", y="exp.correl.coeff.microarray", 
            fill="Duplication", facet.by = "Duplication", panel.labs.font = list(size = 11))+ 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  stat_compare_means(comparisons = c("HM", "HM&HET"), group.by = "Duplication", method = "t.test") +
  ylab("Correlation coefficient (r)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

SSD <- filter(summary_table, Duplication=='SSD')
t.test(SSD$exp.correl.coeff.microarray[SSD$motif.categories=='HM'], 
            SSD$exp.correl.coeff.microarray[SSD$motif.categories=='HM&HET']) # 0.0159

WGD <- filter(summary_table, Duplication=='WGD')
t.test(WGD$exp.correl.coeff.microarray[WGD$motif.categories=='HM'], 
            WGD$exp.correl.coeff.microarray[WGD$motif.categories=='HM&HET']) # 0.860


#Figure6B
#The correlation of expression between paralogs forming HM&HET or only HM
#as a function of their amino acid sequence identity. 
#The correlation of expression data was binned into six equal categories. 

funct <- filter(summary_table, motif.categories !="NI" & motif.categories !="HET") %>%
  select(., Duplication, motif.categories, sim.TF, sim.cell.comp.P1P2, sim.loc)
colnames(funct) <- c('Duplication', 'motif.categories', 'Transcription factor', 'GO cellular component', 'Localization')

funct %<>% gather('Transcription factor', 'GO cellular component', 'Localization', key=GO, value= similarity)

funct_resum <- funct %>% group_by(Duplication, GO, motif.categories) %>% 
  summarise(mean.similarity = mean(similarity, na.rm=TRUE)*100) %>% as.data.frame()

funct_resum %<>% mutate(p.value = ifelse(Duplication=='SSD' & GO=='Transcription factor', 
                                         wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='Transcription factor' & funct$motif.categories=='HM'],
                                                     funct$similarity[funct$Duplication=='SSD' & funct$GO=='Transcription factor' & funct$motif.categories=='HM&HET'])$p.value,
                                  ifelse(Duplication=='SSD' & GO=='GO cellular component', 
                                         wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='GO cellular component' & funct$motif.categories=='HM'],
                                                     funct$similarity[funct$Duplication=='SSD' & funct$GO=='GO cellular component' & funct$motif.categories=='HM&HET'])$p.value,
                                  ifelse(Duplication=='SSD' & GO=='Localization', 
                                         wilcox.test(funct$similarity[funct$Duplication=='SSD' & funct$GO=='Localization' & funct$motif.categories=='HM'],
                                                     funct$similarity[funct$Duplication=='SSD' & funct$GO=='Localization' & funct$motif.categories=='HM&HET'])$p.value,
                                  ifelse(Duplication=='WGD' & GO=='Transcription factor', 
                                         wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='Transcription factor' & funct$motif.categories=='HM'],
                                                     funct$similarity[funct$Duplication=='WGD' & funct$GO=='Transcription factor' & funct$motif.categories=='HM&HET'])$p.value,
                                  ifelse(Duplication=='WGD' & GO=='GO cellular component', 
                                         wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='GO cellular component' & funct$motif.categories=='HM'],
                                                     funct$similarity[funct$Duplication=='WGD' & funct$GO=='GO cellular component' & funct$motif.categories=='HM&HET'])$p.value,
                                  ifelse(Duplication=='WGD' & GO=='Localization', 
                                         wilcox.test(funct$similarity[funct$Duplication=='WGD' & funct$GO=='Localization' & funct$motif.categories=='HM'],
                                                     funct$similarity[funct$Duplication=='WGD' & funct$GO=='Localization' & funct$motif.categories=='HM&HET'])$p.value, NA)))))))

funct_resum %<>% mutate(p.value.range = ifelse(p.value > 0.05, 'NS',
                                        ifelse(p.value < 0.05 & p.value > 0.01, '< 0.05',
                                        ifelse(p.value < 0.01 & p.value > 0.001, '< 0.01',
                                        ifelse(p.value < 0.001, '< 0.001', 'pb')))))

Figure6B <-
ggplot(funct_resum, 
       aes(motif.categories, mean.similarity, group = GO, color = as.factor(GO))) +
  geom_point(size=4) +
  geom_line(color = 'grey20', alpha = 4/10, aes(linetype=p.value.range)) + facet_grid(. ~ Duplication) +
  labs(x = 'Interaction motifs', y ='Similarity (%)', color='') 



#Figure6C
#The similarity of transcription factor binding sites, 
#GO cellular component and GFP-based localization 
#are compared between HM and HM&HET for each SSDs and WGDs

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


Figure6C <- 
  ggplot(coexp_pi_het,aes(x=mean_interval,y=mean.coexp, col=as.factor(motif.categories)))+
  geom_point(shape=15, size=4, position=position_dodge(.3))+ facet_grid(.~Duplication) +
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

Figure6AC <- plot_grid(Figure6A, Figure6C, labels = c("A", 'B'), nrow=1, ncol = 2)


ggsave(file="Figure6.pdf", width=10, height=10, dpi=500)
plot_grid(Figure6AC, Figure6B, labels = c("", "C"), nrow=2, ncol = 1)
dev.off()
