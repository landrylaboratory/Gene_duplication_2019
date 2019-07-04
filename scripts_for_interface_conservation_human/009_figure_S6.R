#############################################
####    Interface conservation human     ####
#### This script plots the conservation  ####
#### of sequences within interfaces      ####
#### based on structural data from the   ####
#### PDB and motifs from BioGRID data.   ####
#############################################

library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(Cairo)

# Set directory
setwd(<path_to_scripts_for_interface_conservation_human>)

#### Load the data ####

pdb_data <- read_delim('Data/Paralogs_sequence_identities_dist_best_isoforms_2019-06-19.txt', delim = '\t')

rohans_data_biogrid_intact <- read_delim('Data/dparalog_uniprotpdbids_gene_ids_biogrid_intact.tsv', delim = '\t')
rohans_data_biogrid_intact %<>% select(`gene1 id`, `gene2 id`, `gene1 name`, `gene2 name`,
                                      `heteromer or not (direct biogrid intact)`, 
                                      `homomer or not (direct biogrid intact) gene1`, 
                                      `homomer or not (direct biogrid intact) gene2`)

colnames(rohans_data_biogrid_intact) <- c("Gene_1", "Gene_2", "Gene_1_name", "Gene_2_name", "HET", "HM_P1", "HM_P2")

rohans_data_biogrid_intact %<>%
mutate(motif = ifelse(and(HET == 0, and(HM_P1 == 0, HM_P2 == 0)), 'NI',
                     ifelse(and(HET == 0, or(HM_P1 == 1, HM_P2 == 1)), "HM",
                            ifelse(and(HET == 1, and(HM_P1 == 0, HM_P2 == 0)), "HET", "HM&HET"))))

rohans_data_condensed <- unique(rohans_data_biogrid_intact)

#### Put the data together ####

# Calculate the completeness
pdb_data %<>% mutate(completeness = round(PDB_length*100/Full_protein_length, 2))

# Name P1 and P2 the same in both tables, with P1 and P2 sorted by gene ID by alphabetical order
pdb_data %<>% rowwise() %>% filter(completeness >= 50) %>%
  mutate(Gene_1_new = ifelse(Gene_1 < Gene_2, Gene_1, Gene_2),
         Gene_2_new = ifelse(Gene_1 > Gene_2, Gene_1, Gene_2),
         
         P1_ID_new = ifelse(Gene_1 < Gene_2, P1_ID, P2_ID),
         P2_ID_new = ifelse(Gene_1 > Gene_2, P1_ID, P2_ID),
         
         Strict_HM_P1_new = ifelse(Gene_1 < Gene_2, Strict_HM_P1, Strict_HM_P2),
         Strict_HM_P2_new = ifelse(Gene_1 > Gene_2, Strict_HM_P1, Strict_HM_P2),
         
         Compared_interface_new = ifelse(Gene_1 < Gene_2, Compared_interface,
                                         ifelse(Compared_interface == "P1_HM", "P2_HM",
                                                ifelse(Compared_interface == "P2_HM", "P1_HM",
                                                       ifelse(Compared_interface == "HET_P1", "HET_P2",
                                                              ifelse(Compared_interface == "HET_P2", "HET_P1", "NA")))))
         )

rohans_data_condensed %<>% rowwise() %>%
  mutate(Gene_1_new = ifelse(Gene_1 < Gene_2, Gene_1, Gene_2),
         Gene_2_new = ifelse(Gene_1 > Gene_2, Gene_1, Gene_2),
         
         Gene_1_name_new = ifelse(Gene_1 < Gene_2, Gene_1_name, Gene_2_name),
         Gene_2_name_new = ifelse(Gene_1 > Gene_2, Gene_1_name, Gene_2_name),
         
         HM_P1_new = ifelse(Gene_1 < Gene_2, HM_P1, HM_P2),
         HM_P2_new = ifelse(Gene_1 > Gene_2, HM_P1, HM_P2)
         )


gathered_data <- inner_join(x = pdb_data %>% 
                              mutate(PDB = ifelse(Compared_interface_new == 'P1_HM', Strict_HM_P1_new, 
                                                  ifelse(Compared_interface_new == 'P2_HM', Strict_HM_P2_new,
                                                         HET))) %>%
                              select(Gene_1_new, Gene_2_new, P1_ID_new, P2_ID_new, 
                                                 Full_sequence_identity, Interface_sequence_identity, 
                                                 Non_interface_sequence_identity, PDB_sequence_identity,
                                                 Compared_interface_new, Same_phylome, PDB),
                         y = rohans_data_condensed %>% select(-Gene_1, -Gene_2, -Gene_1_name, -Gene_2_name, 
                                                              -HM_P1, -HM_P2),
                         by = c("Gene_1_new" = "Gene_1_new", "Gene_2_new" = "Gene_2_new")) %>%
  unique() %>% # There are some duplicate rows
  filter(!(is.na(motif))) %>%
  filter(or(Full_sequence_identity >= 20, Same_phylome == 1))

# There are still some pairs of paralogs that are repeated because different isoform IDs that
# map to identical sequences are compared.
final_data <- c()
already_seen <- c()
for(line_num in 1:nrow(gathered_data)){
  pair = paste(gathered_data$Gene_1_new[line_num], gathered_data$Gene_2_new[line_num], gathered_data$Compared_interface_new[line_num], sep = '.')
  if(!(pair %in% already_seen)){
    # Add it to the new table
    final_data <- rbind(final_data, gathered_data[line_num,])
    already_seen <- c(already_seen, pair)
  }
}
  
# Include the data from the PDB structures in the motifs if they are still zeros
# This adds the data from one structure whose HM had not been seen in BioGRID but 
# is crystallized
table(final_data$Compared_interface_new, final_data$motif)
final_data2 <- final_data %>% 
  mutate(HM_P1_new = ifelse(Compared_interface_new == 'P1_HM', 1, HM_P1_new),
         HM_P2_new = ifelse(Compared_interface_new == 'P2_HM', 1, HM_P2_new),
         HET = ifelse(Compared_interface_new == 'HET_P1' || Compared_interface_new == 'HET_P2', 1, HET)) %>%
  mutate(motif = ifelse(and(HET == 0, and(HM_P1_new == 0, HM_P2_new == 0)), 'NI',
                        ifelse(and(HET == 0, or(HM_P1_new == 1, HM_P2_new == 1)), "HM",
                               ifelse(and(HET == 1, and(HM_P1_new == 0, HM_P2_new == 0)), "HET", "HM&HET"))))
           
table(final_data2$motif)
# HM HM&HET 
# 40     25 

pdb_data_p1 <- final_data2 %>%
  gather(key = region, value = seq_ident, Full_sequence_identity, PDB_sequence_identity,
         Interface_sequence_identity, Non_interface_sequence_identity) %>%
  mutate(region = gsub(x = region, pattern = "Full_sequence_identity", replacement = 'Full sequence'),
         region = gsub(x = region, pattern = "PDB_sequence_identity", replacement = 'PDB sequence'),
         region = gsub(x = region, pattern = "Non_interface_sequence_identity", replacement = 'Non-interfaces'),
         region = gsub(x = region, pattern = "Interface_sequence_identity", replacement = 'Interface'))

comparison_list = list(c('Full sequence', 'Interface'))

pdb_data_p1 %<>% 
  unite(col = pair, sep = '.', Gene_1_new, Gene_2_new, remove = FALSE) %>%
  filter(!(PDB %in% c('1cki_1', '2oan_1', '2c6q_1', '2jcs_1'))) # Incorrectly assigned structures
  
#### The plot for sequence identities (Panel A) ####

p1 <- pdb_data_p1 %>% rowwise() %>%
  filter(motif != 'HET', region %in% c('Full sequence', 'Interface')) %>%
  ggplot(aes(x = region, y = seq_ident, fill = motif, colour = motif, group = region)) +
  geom_point(aes(colour = motif),
             position = position_jitterdodge(jitter.width = 0.0)
  ) +
  geom_line(aes(group = interaction(pair, Compared_interface_new)), alpha = 0.7) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  facet_wrap(~motif) + 
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

#### The plot for the interface conservation score (panel B) ####

final_data2 %<>% filter(or(Full_sequence_identity >= 20, Same_phylome == 1)) %>%
  mutate(cons_score = Interface_sequence_identity / Non_interface_sequence_identity) %>%
  filter(!(PDB %in% c('1cki_1', '2oan_1', '2c6q_1', '2jcs_1'))) # Incorrectly assigned structures

comparison_list = list(c('HM', 'HM&HET'))

p2 <- final_data2 %>% 
  filter(motif %in% c('HM', 'HM&HET')) %>%
  ggplot(aes(x = motif, y = cons_score, fill = motif, colour = motif)) +
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

ggsave(filename = 'FigureS6.pdf',
       device = cairo_pdf, width = 7, height = 14, plot = p_int_conservation, dpi = 500)
