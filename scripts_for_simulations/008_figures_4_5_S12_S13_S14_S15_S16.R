#########################################################################################
####                          Figures simulations                                    ####
#### This script plots the figures on the results of the simulations, that is:       ####
#### - Figure 4                                                                      ####
#### - Figure 5                                                                      ####
#### - Figure S12                                                                    ####
#### - Figure S13                                                                    ####
#### - Figure S14                                                                    ####
#### - Figure S15                                                                    ####
#### - Figure S16                                                                    ####
#########################################################################################

# Libraries 
library(ggplot2)
library(tidyverse)
library(magrittr)
library(cowplot)
library(grid)
library(gridExtra)
library(ggpubr)
library(Cairo)

#### Set the working directory to the folder with the script

setwd('/path/to/scripts_for_simulations')

#### Figure 4 ####

# Load the data
data_fig_4_s12_1A82 <- read.table('Data/data_fig_4_s12_1A82.tsv',
                                 sep = '\t', header = T)
data_fig_4_s12_1M38 <- read.table('Data/data_fig_4_s12_1M38.tsv',
                                 sep = '\t', header = T)
data_fig_4_s12_2JKY <- read.table('Data/data_fig_4_s12_2JKY.tsv',
                                 sep = '\t', header = T)
data_fig_4_s12_2O1V <- read.table('Data/data_fig_4_s12_2O1V.tsv',
                                 sep = '\t', header = T)
data_fig_4_s12_3D8X <- read.table('Data/data_fig_4_s12_3D8X.tsv',
                                 sep = '\t', header = T)
data_fig_4_s12_4FGW <- read.table('Data/data_fig_4_s12_4FGW.tsv',
                                 sep = '\t', header = T)

data_fig_4_s12 <- rbind(data_fig_4_s12_1A82, data_fig_4_s12_1M38, data_fig_4_s12_2JKY,
                       data_fig_4_s12_2O1V, data_fig_4_s12_3D8X, data_fig_4_s12_4FGW)

# Load the diagram for the simulations
fig_4A <- ggdraw() + draw_image('Data/Figure4-A.png')

# Relevel to adjust order in legend
data_fig_4_s12$ComplexType = factor(data_fig_4_s12$ComplexType, levels = c('HM AA', 'HET AB', 'HM BB'))

# Plot each of the scenarios for 1M38
# No selection #

p_1M38_nosel <- data_fig_4_s12 %>% 
  filter(PDB == '1M38', Scenario == 'Neutral evolution') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) + 
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on both HMs #

p_1M38_sel_AA_BB <- data_fig_4_s12 %>% 
  filter(PDB == '1M38', Scenario == 'Selection on both HMs') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on HET #

p_1M38_sel_AB <- data_fig_4_s12 %>% 
  filter(PDB == '1M38', Scenario == 'Selection on HET AB') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) + 
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on HM AA #

p_1M38_sel_AA <- data_fig_4_s12 %>% 
  filter(PDB == '1M38', Scenario == 'Selection on HM AA') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on HM BB #

p_1M38_sel_BB <- data_fig_4_s12 %>% 
  filter(PDB == '1M38', Scenario == 'Selection on HM BB') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Save the panels for figure 4 #
final_fig4 <- plot_grid(fig_4A, p_1M38_nosel, p_1M38_sel_AA_BB,
                        p_1M38_sel_AB, p_1M38_sel_AA, p_1M38_sel_BB,
                        nrow = 3, 
                        labels = c('A', 'B', 'C', 'D', 'E', 'F'))

ggsave(filename = 'Figure4.pdf',
       width = 14, height = 21, dpi = 500, plot = final_fig4, device = cairo_pdf)

#### Figure 5 ####

# Load the data
data_fig_5_s14 <- read.table('Data/data_fig_5_s14.tsv',
                             h = T, sep = '\t')

# Save helper variables
y_min <- -5
y_max <- 10
x_min <- -5
x_max <- 15

# Load the PDB figure
fig_1m38 <- ggdraw() + draw_image('Data/1M38_new_2.png')

fig_1M38_sel_AA_BB <- data_fig_5_s14 %>% 
  filter(Scenario == 'Selection on both HMs', PDB == '1M38') %>%
  ggplot(aes(y = het_ddg_be, x = homo_ddg_be, color = Verdict, shape = out_of_bounds)) + 
  facet_grid(~Scenario) +
  scale_color_manual(values = c('red', 'black')) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  geom_point(aes(alpha = Verdict), size = 2) +
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = "bold"), axis.title.y = element_text(size=15, face = "bold"))+
  theme(plot.title = element_text(size = 20, face = "bold",family="Arial", hjust = 0.5),
        legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=15), legend.position = "top",
        legend.justification = "center")+
  guides(color = guide_legend(override.aes = list(size=3))) +
  ylab("ΔΔG HET (Kcal/mol)") + xlab("ΔΔG HM (Kcal/mol)") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  stat_cor(method = 'pearson', label.x.npc = 0.35,
           label.y.npc = 0.1, show.legend = FALSE,
           inherit.aes = FALSE,
           aes(x = het_ddg_be, y = homo_ddg_be, label = gsub(x = ..label.., replacement = 'r', pattern = 'R')),
           size = 6
  ) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white'))

fig_1M38_sel_AB <- data_fig_5_s14 %>% 
  filter(Scenario == 'Selection on HET AB', PDB == '1M38') %>%
  ggplot(aes(y = het_ddg_be, x = homo_ddg_be, color = Verdict, shape = out_of_bounds)) + 
  facet_grid(~Scenario) +
  scale_color_manual(values = c('red', 'black')) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  geom_point(aes(alpha = Verdict), size = 2) +
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = "bold"), axis.title.y = element_text(size=15, face = "bold"))+
  theme(plot.title = element_text(size = 20, face = "bold",family="Arial", hjust = 0.5),
        legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=15), legend.position = "top",
        legend.justification = "center")+
  guides(color = guide_legend(override.aes = list(size=3))) +
  ylab("ΔΔG HET (Kcal/mol)") + xlab("ΔΔG HM (Kcal/mol)") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  stat_cor(method = 'pearson', label.x.npc = 0.35,
           label.y.npc = 0.1, show.legend = FALSE,
           inherit.aes = FALSE,
           aes(x = het_ddg_be, y = homo_ddg_be, label = gsub(x = ..label.., replacement = 'r', pattern = 'R')),
           size = 6
  ) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white'))

all_plots <- plot_grid(
 fig_1m38, fig_1M38_sel_AA_BB, fig_1M38_sel_AB,
 ncol = 3, nrow = 1, labels = c('', 'A', 'B'))

ggsave(filename = 'Figure5.pdf',
       width = 21, height = 7, dpi = 500, device = cairo_pdf, plot = all_plots)


#### Figure S12 ####

# Load the figures of the PDB structures
fig_1a82 <- ggdraw() + draw_image('Data/1A82_new.png')
fig_2o1v <- ggdraw() + draw_image('Data/2O1V_new.png')
fig_1m38 <- ggdraw() + draw_image('Data/1M38_new.png')
fig_4fgw <- ggdraw() + draw_image('Data/4FGW_new.png')
fig_3d8x <- ggdraw() + draw_image('Data/3D8X_new.png')
fig_2jky <- ggdraw() + draw_image('Data/4FGW_new.png')

# Relevel to adjust order in legend
data_fig_4_s12$ComplexType = factor(data_fig_4_s12$ComplexType, levels = c('HM AA', 'HET AB', 'HM BB'))

#### Save figure S12 ####

p_all <- data_fig_4_s12 %>% ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(Scenario ~ PDB, labeller = label_wrap_gen(width = 15)) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=20), axis.text.y= element_text(size=20), 
        axis.title.x = element_text(size=25, face = 'bold'), axis.title.y = element_text(size=25, face = 'bold'))+
  theme(
    legend.key = element_rect(fill = "white"),
    legend.text=element_text(size=25), legend.position = "none",
    legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-60,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 25, face = 'bold'),
        strip.text.y = element_text(size = 25, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  ) 


legend <- get_legend(p_all + theme(legend.position = 'top', 
                                   legend.text=element_text(size=25),
                                   legend.justification = "center"
) +
  guides(color = guide_legend(override.aes = list(size=7))))

header_pdb <- plot_grid(NULL, fig_1a82, fig_1m38, fig_2jky, fig_2o1v, fig_3d8x, fig_4fgw, NULL, nrow = 1,
                        rel_widths = c(0.2, 1, 1, 1, 1, 1, 1, 0.1))

fig_s12 <- plot_grid(header_pdb, legend, p_all,
                    rel_heights = c(0.5, 0.05, 3), ncol = 1)

ggsave(filename = 'FigureS12.pdf',
       width = 28, height = 28, dpi = 500, device = cairo_pdf, plot = fig_s12)

#### Figure S13 ####

# Load the data
data_fig_s13_1 <- read_delim(file = 'Data/data_fig_s13_1.tsv', delim = '\t')
data_fig_s13_2 <- read_delim(file = 'Data/data_fig_s13_2.tsv', delim = '\t')
data_fig_s13_3 <- read_delim(file = 'Data/data_fig_s13_3.tsv', delim = '\t')
data_fig_s13_4 <- read_delim(file = 'Data/data_fig_s13_4.tsv', delim = '\t')
data_fig_s13_5 <- read_delim(file = 'Data/data_fig_s13_5.tsv', delim = '\t')

# Put them together
data_fig_s13 <- rbind(data_fig_s13_1, data_fig_s13_2, data_fig_s13_3, data_fig_s13_4, data_fig_s13_5) %>%
  mutate(Scenario = gsub(pattern = 'Selection on both homodimers', x = Scenario, replacement = 'Selection on both HMs'),
         Scenario = gsub(pattern = 'Selection on heterodimer AB', x = Scenario, replacement = 'Selection on HET AB'),
         Scenario = gsub(pattern = 'Selection on homodimer AA', x = Scenario, replacement = 'Selection on HM AA'),
         Scenario = gsub(pattern = 'Selection on homodimer BB', x = Scenario, replacement = 'Selection on HM BB'),
         ComplexType = gsub(pattern = 'Homodimer AA', x = ComplexType, replacement = 'HM AA'),
         ComplexType = gsub(pattern = 'Homodimer BB', x = ComplexType, replacement = 'HM BB'),
         ComplexType = gsub(pattern = 'Heterodimer AB', x = ComplexType, replacement = 'HET AB'))

data_fig_s13$ComplexType = factor(data_fig_s13$ComplexType, levels = c('HM AA', 'HET AB', 'HM BB'))

# Plot
p_all <- data_fig_s13 %>% ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(Scenario ~ params, labeller = label_wrap_gen(width = 30)) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "violet"),
                      name = "")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), 
        axis.title.x = element_text(size=20, face = 'bold'), axis.title.y = element_text(size=20, face = 'bold'))+
  theme(
    legend.key = element_rect(fill = "white"),
    legend.text=element_text(size=20), legend.position = "none",
    legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-55,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=4),
        strip.text.x = element_text(size = 20, face = 'bold'),
        strip.text.y = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  ) 


legend <- get_legend(p_all + theme(legend.position = 'top', 
                                   legend.text=element_text(size=20),
                                   legend.justification = "center"
) +
  guides(color = guide_legend(override.aes = list(size=5))))

fig_s13 <- plot_grid(legend, p_all,
                    rel_heights = c(0.05, 3), ncol = 1)

ggsave(filename = 'FigureS13.pdf',
       width = 21, height = 21, dpi = 500, device = cairo_pdf, plot = fig_s13)

#### Figure S14 ####

# Load the data
data_fig_5_s14 <- read.table('Data/data_fig_5_s14.tsv',
                             h = T, sep = '\t')

# Load the figures of the PDB structures
fig_1a82 <- ggdraw() + draw_image('Data/1A82_new.png')
fig_2o1v <- ggdraw() + draw_image('Data/2O1V_new.png')
fig_1m38 <- ggdraw() + draw_image('Data/1M38_new.png')
fig_4fgw <- ggdraw() + draw_image('Data/4FGW_new.png')
fig_3d8x <- ggdraw() + draw_image('Data/3D8X_new.png')
fig_2jky <- ggdraw() + draw_image('Data/4FGW_new.png')

p_all_facet <- data_fig_5_s14 %>% ggplot(aes(y = het_ddg_be, x = homo_ddg_be, color = Verdict, shape = out_of_bounds)) + 
  facet_grid(Scenario ~ PDB, labeller = label_wrap_gen(width = 15)) +
  scale_color_manual(values = c('red', 'black')) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  geom_point(aes(alpha = Verdict), size = 2) +
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), 
        axis.title.x = element_text(size=20, face = "bold"), axis.title.y = element_text(size=20, face = "bold"))+
  theme(plot.title = element_text(size = 20, face = "bold",family="Arial", hjust = 0.5),
        legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "none",
        legend.justification = "center")+
  guides(color = guide_legend(override.aes = list(size=3))) +
  ylab("ΔΔG HET (Kcal/mol)") + xlab("ΔΔG HM (Kcal/mol)") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  stat_cor(method = 'pearson', label.x.npc = 0.35,
           label.y.npc = 0.1, show.legend = FALSE,
           inherit.aes = FALSE,
           aes(x = het_ddg_be, y = homo_ddg_be, label = gsub(x = ..label.., replacement = 'r', pattern = 'R')),
           size = 6
  ) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size = 1),
        strip.text.x = element_text(size = 20, face = 'bold'),
        strip.text.y = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white'))

header_pdb <- plot_grid(NULL, fig_1a82, fig_1m38, fig_2jky, fig_2o1v, fig_3d8x, fig_4fgw, NULL, nrow = 1,
                        rel_widths = c(0.2, 1, 1, 1, 1, 1, 1, 0.1))

legend <- get_legend(p_all_facet + theme(legend.position = 'top', 
                                         legend.text=element_text(size=20)) +
                       guides(color = guide_legend(override.aes = list(size=4))))

fig_s14 <- plot_grid(header_pdb, legend, p_all_facet, rel_heights = c(1, 0.2, 4), ncol = 1)

ggsave(filename = 'FigureS14.pdf',
       width = 28, height = 14, dpi = 500, device = cairo_pdf, plot = fig_s14)

#### Figure S15 ####

# Load the data for panel A
data_fig_s15A <- read.table('Data/data_fig_s15A.tsv',
                            h = T, sep = '\t')

# Relevel for order in legend
data_fig_s15A$Complex <- factor(data_fig_s15A$Complex, levels = c('HM', 'HET'))

p15_A <- data_fig_s15A %>% 
  ggplot(aes(x = quadrant, y = ddg, fill = Complex)) +
  geom_point(aes(colour = Complex),
             position = position_jitterdodge(jitter.width = 0.25)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_fill_manual(values = c('#737373', '#d9d9d9')) +
  scale_colour_manual(values = c('#737373', '#d9d9d9'), guide = 'none') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold', size = 15),
        axis.title.y = element_text(face = 'bold', size = 15),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
  xlab('Class of mutation') + ylab('ΔΔG (kcal/mol)') + labs(fill = '') +
  ylim(-5,10)
p15_A

# Load the data for panel B
data_fig_s15B <- read.table('Data/data_fig_s15B.tsv',
                            h = T, sep = '\t')

# Relevel for order in legend
data_fig_s15B$Complex <- factor(data_fig_s15B$Complex, levels = c('HM', 'HET'))

p15_B <- data_fig_s15B %>% 
  ggplot(aes(x = quadrant, y = ddg, fill = Complex)) +
  geom_point(aes(colour = Complex),
             position = position_jitterdodge(jitter.width = 0.25)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_fill_manual(values = c('#737373', '#d9d9d9')) +
  scale_colour_manual(values = c('#737373', '#d9d9d9'), guide = 'none') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold', size = 15),
        axis.title.y = element_text(face = 'bold', size = 15),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
  xlab('Class of mutation') + ylab('ΔΔG (kcal/mol)') + labs(fill = '') +
  ylim(-5,10)
p15_B

p_effect_sizes <- plot_grid(p15_A, NULL, p15_B, labels = c('A', NULL, 'B'), nrow = 3, rel_heights = c(1, 0.1, 1))

ggsave(filename = 'FigureS15.pdf',
       device = cairo_pdf, width = 7, height = 14, plot = p_effect_sizes, dpi = 500)

#### Figure S16 ####

# Load the data
data_fig_s16 <- read.table('Data/data_fig_s16.tsv',
                           h = T, sep = '\t')

fig_s16 <- data_fig_s16 %>% 
  ggplot(aes(x = quadrant, y = accepted_percentage_of_total, fill = selection)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = min, ymax = max), 
                width = 0.5, size = 0.3,
                position = position_dodge(0.9)) +
  scale_fill_manual(values = c('#737373', '#d9d9d9')) + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5)) +
  xlab("Effect on HMs") + ylab("Fixed mutations (%)") + labs(fill = 'Selection on')
fig_s16

ggsave(filename = 'FigureS16.pdf',
       width = 10,height = 7, plot = fig_s16, device = cairo_pdf, dpi = 500)


