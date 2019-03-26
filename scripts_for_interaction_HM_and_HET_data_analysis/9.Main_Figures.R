rm(list=ls())

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


summary_table <- read.delim('output/TableS1.csv', header=T, sep="\t")

######## Fig2 #########

# Fig2A
dff.S.bg.Kim.PCA <- read.table("output/HM.data.csv", sep="\t", header=T)

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
#2.497e-15
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd"]), workspace=2e8)
#p-value < 2.2e-16
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="wgd"]), workspace=2e8)
#p-value =  9.875e-06
fisher.test(dff.S.bg.Kim.PCA$HM.bg.kim.S.PCA[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"], 
            as.factor(dff.S.bg.Kim.PCA$type_para[dff.S.bg.Kim.PCA$type_para=="S" | dff.S.bg.Kim.PCA$type_para=="ssd_wgd"]), workspace=2e8)
#p-value = 0.03717

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


Fig2A <-
  ggplot(data=dbg_percent, aes(x=dupli, y=percentHM, fill=dupli)) + 
  geom_bar(stat="identity") + geom_text(aes(label=paste(round(percentHM,2))), vjust=-0.3, size=4) +
  scale_fill_manual(guide=FALSE,values = c(cl[285],cl[16],cl[144],cl[129],cl[310])) +
  geom_signif(xmin=c(2, 2, 2),
              xmax=c(3, 4, 5), 
              y_position=c(44, 48, 52), 
              annotation=c("<2.0e-16","9.875e-06","2.438e-13"), textsize=4) +
  ylab("Percentage of homomers (%)") + xlab("Groups of genes")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        legend.position = c(0.08, 0.9), legend.title = element_text(size=15), 
        legend.text = element_text(size=12))


# Fig2B, C, D
Fig2BC <- ggdraw() + draw_image("data/Fig2BC_Diana.png")
Fig2D <- ggdraw() + draw_image("data/Fig2D.png")

# Fig2E
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


#result
#motif.categories p.value_fisherTest
#1              HET       3.836343e-01
#2               HM       2.529560e-05
#3           HM&HET       1.630493e-06
#4               NI       2.860273e-01


Fig2E1 <-
  ggplot(freq, aes(motif.categories, freq)) + 
  geom_bar(aes(fill = Duplication), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  labs(x=("Interaction motifs"), y=("Percentage (%)"), title = ("Duplication")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75),xend=c(1.25, 2.25, 3.25, 4.25), 
                              y=c(12, 55, 43, 18), 
                              annotation=c("3.84e-01", "2.53e-05", "1.63e-06", "0.29")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.15, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.background = element_blank())

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
#result:
#motif.categories p.value_fisherTest
#1              HET         0.40946337
#2               HM         0.04951537
#3           HM&HET         0.01095495
#4               NI         1.00000000


freq$Origin.of.WGDs <- factor(freq$Origin.of.WGDs, levels = c("Homeologs", "True_ohnologs"))

Fig2E2 <-
  ggplot(freq, aes(motif.categories, freq)) + 
  geom_bar(aes(fill = Origin.of.WGDs), stat="identity", position="dodge") + 
  scale_fill_manual(values = c(cl[430],cl[490]), labels = c("Homeologs", "True ohnologs")) + 
  labs(x=("Interaction motifs"), title = ("Origin of WGDs")) +
  geom_signif(textsize=16, stat="identity",
              data=data.frame(x=c(0.75, 1.75, 2.75, 3.75),xend=c(1.25, 2.25, 3.25, 4.25), 
                              y=c(16, 43, 54, 16), 
                              annotation=c("0.41", "4.95e-02", "1.09e-02", "1.00")),
              aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), #axis.text.y= element_blank(), 
        axis.title.x = element_text(size=15), #axis.title.y = element_blank(), axis.ticks.y =  element_blank(),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = c(0.3, 0.94), 
        legend.title = element_blank(), legend.text = element_text(size=12), 
        legend.background = element_blank())


Fig2E <- plot_grid(Fig2E1, Fig2E2, labels = c("", ""), align = 'h')


# Fig2F
my_comparisons <- list( c("HM", "HET"), c("HET", "HM&HET"), c("HET", "NI"),
                        c("HM", "HM&HET"), c("HM", "NI"),
                        c("HM&HET", "NI"))
Fig2F <-
  ggboxplot(summary_table,  x="motif.categories", y="pid", fill="Duplication", facet.by = "Duplication")+ 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  ylab("Pairwise amino acid sequence identity (%)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position='none')

ggsave(file="main_fig/Figure2_previous.pdf", width=14, height=21, dpi=500)
plot_grid(Fig2A, Fig2BC, Fig2D, Fig2E, Fig2F, labels = c("A", "", "D", "E", "F"), ncol=3, nrow = 2)
dev.off()


########FigS3########

fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}


summary_table$int <- interaction(summary_table$motif.categories, summary_table$Duplication)
my_comparisons <- list(c("HM.SSD", "HM&HET.SSD"), c("HM.WGD", "HM&HET.WGD"))



Fig3A <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
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

Fig3B <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET")  %>% 
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

Fig3C <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
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

Fig3D <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
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


ggsave(file="main_fig/Figure3.pdf", width=10, height=10, dpi=500)
plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, labels=c("A", "B", "C", "D"), nrow=2, align=c("h"))
dev.off()


#######Fig6#######

#Fig6A

my_comparisons <- list( c("HM", "HET"), c("HET", "HM&HET"), c("HET", "NI"),
                        c("HM", "HM&HET"), c("HM", "NI"),
                        c("HM&HET", "NI"))

Fig6A <-
  ggboxplot(filter(summary_table, !is.na(motif.categories)), x="motif.categories", y="expression.correl.coeff", 
            fill="Duplication", facet.by = "Duplication", panel.labs.font = list(size = 11))+ 
  scale_fill_manual(values = c(cl[144],cl[129])) +
  stat_compare_means(comparisons = my_comparisons, group.by = "Duplication", method = "wilcox.test") +
  ylab("Correlation coefficient (r)") + xlab("Interaction motifs") + 
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15), axis.title.y = element_text(size=15),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")



#Fig6B
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=round(mean(x,na.rm=T), 2)))}


summary_table$int <- interaction(summary_table$motif.categories, summary_table$Duplication)
my_comparisons <- list(c("HM.SSD", "HM&HET.SSD"), c("HM.WGD", "HM&HET.WGD"))

Fig6B <- summary_table %>% filter(motif.categories =="HM" | motif.categories =="HM&HET") %>% 
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


#Fig6C
coexp_pi_het <- summary_table %>% filter(motif.categories=="HM" | motif.categories=="HM&HET") %>%
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


Fig6C<- 
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

Fig6D <- summary_table %>% filter(motif.categories !="NI" & motif.categories !="HET") %>% 
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


Fig6E <- summary_table %>% filter((motif.categories =="HM" | motif.categories =="HM&HET") & !is.na(sim.loc)) %>% 
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



ggsave(file="main_fig/Fi6.pdf", width=10, height=15, dpi=500)
plot_grid(Fig6A, Fig6B, Fig6C, Fig6D, Fig6E, labels = c("A", 'B', "C", 'D', 'E'), nrow=3, ncol = 2)
dev.off()
