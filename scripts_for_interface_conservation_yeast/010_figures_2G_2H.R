#############################################
####            Figures 2G, 2H           ####
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

setwd(<path_to_scripts_for_interface_conservation_yeast>)

#### Load the data ####
### These are the structures we removed because of the following reasons:
# 1.- The crystallized domain is only present in one of the paralogs:
# - 1zy4: The kinase domain is present in YDR283C but not in YPR033C
# - 1ygh: The HAT domain is present in YGR252W but not in YLR399C
# - 1m4z: The BAH domain is present in YML065W but not in YJL194W
# 2.- The crystallized fragment is too small
# - 2rkl (53 residues)
# 3.- The subunits do not come into contact in the crystallized structures
# - 4q0w
# - 5jea

removed <- c('4q0w_1', '2rkl_1', '1zy4_1', '1ygh_1', '1m4z_1', '5jea_1')

pdb_data <- read.table('Data/Paralogs_sequence_identities_dist_completeness_2019-06-19.txt',
                      h = T, sep = '\t')

pdb_data %<>% filter(!(HET %in% removed),
                     !(Strict_HM_P1 %in% removed),
                     !(Strict_HM_P2 %in% removed))

axelles_data <- read.table('Data/summary_PCA_results_per_pairs_2019-01-07_clean.csv',
                           h = T, sep = ',')

#### Add a pair column to the PDB data ####

pdb_data$pair <- paste(pdb_data$P1_ID, pdb_data$P2_ID, sep = '.')

pdb_data %<>% 
  mutate(pair_final = ifelse(pair %in% axelles_data$pair,
                     pair,
                     paste(P2_ID, P1_ID, sep = '.')),
         completeness = PDB_length*100/Full_protein_length) %>%
  filter(completeness >= 50, completeness <= 100)

# Now, I should add the category to the pdb data
pdb_data_final <- left_join(x = pdb_data, y = axelles_data %>% select(pair, motif_categories_PCA, motif_categories),
                            by = c("pair" = "pair"))

table(pdb_data_final$motif_categories)
# HET     HM HM&HET     NI 
#   0     33     28      1 

# Apply the filter
# There is a homomer (YGR203W) that has a crystal structure but is not reported in our data set.
pdb_data_final %<>% filter(!(is.na(motif_categories))) %>%
  filter(or(Full_sequence_identity >= 20, Same_phylome == 1)) %>%
  mutate(motif_categories = gsub(x = motif_categories, pattern = "NI", replacement = "HM")) 


table(pdb_data_final$motif_categories)
# HM HM&HET      
# 30     28      

pdb_data_final %<>%
  gather(key = region, value = seq_ident, Full_sequence_identity, PDB_sequence_identity,
         Interface_sequence_identity, Non_interface_sequence_identity)

comparison_list = list(c('Full sequence', 'Interface'))

pdb_data_p1 <- pdb_data_final %>% 
  mutate(region = gsub(x = region, pattern = "Full_sequence_identity", replacement = 'Full sequence'),
         region = gsub(x = region, pattern = "PDB_sequence_identity", replacement = 'PDB sequence'),
         region = gsub(x = region, pattern = "Non_interface_sequence_identity", replacement = 'Non-interfaces'),
         region = gsub(x = region, pattern = "Interface_sequence_identity", replacement = 'Interface'))

jitter_b <- runif(nrow(pdb_data_p1), -0.1, 0.1)

pdb_data_p1 %<>%
  filter(motif_categories != 'Only crystal', region %in% c('Full sequence', 'Interface'))

#### Plot figure 2G ####

p1 <- pdb_data_p1 %>%
  filter(motif_categories != 'Only crystal', region %in% c('Full sequence', 'Interface')) %>%
  ggplot(aes(x = region, y = seq_ident, fill = motif_categories, colour = motif_categories, group = region)) +
  geom_point(aes(colour = motif_categories),
             position = position_jitterdodge(jitter.width = 0.0)
  ) +
  geom_line(aes(group = interaction(pair_final, Compared_interface)), alpha = 0.7) +
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

pdb_data_final <- left_join(x = pdb_data %>% filter(or(Full_sequence_identity >= 20, Same_phylome == 1)), 
                            y = axelles_data %>% select(pair, motif_categories_PCA, motif_categories),
                            by = c("pair" = "pair"))

#### Plot Figure 2H ####

pdb_data_final %<>% 
  mutate(cons_score = Interface_sequence_identity / Non_interface_sequence_identity)

comparison_list = list(c('HM', 'HM&HET'))

p2 <- pdb_data_final %>% 
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

p_int_conservation <- plot_grid(p1, NULL, p2, nrow = 3, rel_heights = c(1, 0.05, 1), labels = c('A', NULL, 'B'))

ggsave(filename = 'Figures_2G_2H.pdf',
       device = cairo_pdf, width = 7, height = 14, plot = p_int_conservation, dpi = 500)

