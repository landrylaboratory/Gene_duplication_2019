#############################################
####            Figure S6                ####
#### This script plots the conservation  ####
#### of sequences within interfaces      ####
#### based on structural data from the   ####
#### PDB and motifs from our PCA data.   ####
#############################################

library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(Cairo)
library(grid)
library(gridExtra)

#### Set the working directory ####

setwd('/path/to/scripts_for_interface_conservation')

#### Panel A ####

# Load the data
pdb_data_p1 <- read.table('Data/data_fig_s6A.tsv', h = T, sep = '\t')

# Plot
p1 <- pdb_data_p1 %>%
  ggplot(aes(x = region, y = seq_ident, fill = motif_categories, colour = motif_categories, group = region)) +
  geom_point(aes(colour = motif_categories),
             position = position_jitterdodge(jitter.width = 0.05)
  ) +
  geom_line(aes(group = interaction(pair_final)), alpha = 0.7) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  facet_wrap(~motif_categories) + 
  scale_fill_manual(values = c('#eea2ac', 'purple', '#eea2ac', 'purple')) + 
  scale_colour_manual(values = c('#eea2ac', 'purple', '#eea2ac', 'purple')) + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold', size = 15), axis.title.y = element_text(face = 'bold', size = 15),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')) +
  stat_compare_means(comparisons = comparison_list, method = "wilcox.test", 
                     paired = TRUE) +
  xlab("Protein region") + ylab("Pairwise amino acid sequence identity (%)") + labs(fill = '')
p1

#### Panel B ####

# Load the data
pdb_data_p2 <- read.table('Data/data_fig_s6B.tsv', h = T, sep = '\t')

# Plot
p2 <- pdb_data_p2 %>% 
  filter(motif_categories %in% c('HM', 'HM&HET')) %>%
  ggplot(aes(x = motif_categories, y = cons_score, fill = motif_categories, colour = motif_categories)) +
  geom_jitter(width = 0.15) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_colour_manual(values = c('#eea2ac', 'purple', '#eea2ac', 'purple')) + 
  scale_fill_manual(values = c('#eea2ac', 'purple', '#eea2ac', 'purple')) + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold', size = 15), axis.title.y = element_text(face = 'bold', size = 15),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  stat_compare_means(comparisons = comparison_list, method = "wilcox.test", 
                     paired = FALSE) +
  xlab("Interaction motif") + ylab("Relative conservation score") + labs(fill = '')
p2

#### Save figure S6 ####

p_int_conservation <- plot_grid(p1, NULL, p2, nrow = 3, rel_heights = c(1, 0.05, 1), labels = c('A', NULL, 'B'))
p_int_conservation
ggsave(filename = 'FigureS6.pdf',
      device = cairo_pdf, width = 7, height = 14, plot = p_int_conservation, dpi = 500)




