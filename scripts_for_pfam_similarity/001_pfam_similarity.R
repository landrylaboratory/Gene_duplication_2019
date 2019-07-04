###############################################
####      PFAM domain assignments          ####
#### I will use this script to derive      ####
#### Jaccard indices for our pairs of      ####
#### yeast paralogs. This script also      ####
#### figure S7                             ####
###############################################

# Load libraries

library(tidyverse)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)

# Set working directory (path to the folder that contains this script)
setwd(<path_to_scripts_for_pfam_similarity>)

# Load the data on Pfam assignments
pfam_yeast <- read_delim('Data/yeast_proteome_pfam.tab', delim = '\t')

# Get a column with the systematic yeast IDs
pfam_yeast %<>% rowwise() %>%
  mutate(systematic_ID_1 = str_extract(`Gene names`, 'Y.[LR]...[CW](-.)?'),
         systematic_ID_2 = str_extract(`Gene names`, 'Q[0-9]+')) %>%
  mutate(systematic_ID_final = ifelse(is.na(systematic_ID_1) && !(is.na(systematic_ID_2)), systematic_ID_2,
                                      systematic_ID_1))

# Load the lists of paralogs
ssd_pairs <- read_delim('Data/duplication_SDS_1paire.txt', delim = '\t', col_names = c('P1', 'P2'))
wgd_pairs <- read_delim('Data/WGD.csv', delim = ';', col_names = c('P1', 'P2'))

# Concatenate the lists
all_duplicates <- bind_rows(ssd_pairs %>% mutate(Dup_type = 'SSD'), wgd_pairs %>% mutate(Dup_type = 'WGD'))

# Add the information about Pfam domains to each of the duplicates
all_duplicates_pfam <- left_join(x = all_duplicates, y = pfam_yeast %>% select(`Cross-reference (Pfam)`, systematic_ID_final),
                                  by = c("P1" = "systematic_ID_final"))
all_duplicates_pfam <- left_join(x = all_duplicates_pfam, y = pfam_yeast %>%  select(`Cross-reference (Pfam)`, systematic_ID_final),
                                 by = c("P2" = "systematic_ID_final"))
colnames(all_duplicates_pfam) <- c('P1', 'P2', 'Dup_type', 'Pfam_P1', 'Pfam_P2')

# Function for similarity
similarity <- function(col1,col2){
  
  if (!is.na(col1) & !is.na(col2) & col1 !="" & col2 !="") {
    a <- unlist(strsplit(col1, split=";"))
    b <- unlist(strsplit(col2, split=";"))
    
    res = intersect(a,b)
    resl = length(res)
    tot = length(union(a,b))
    
    if (resl>0) {
      return(resl/tot)
    } 
    else {
      return(0)
    }
  }
  
  else{
    return(NA)
  }
}
similarity <- Vectorize(similarity)

# Removing NAs
all_duplicates_pfam_jaccard <- all_duplicates_pfam %>% filter(!(is.na(Pfam_P1)) && !(is.na(Pfam_P2))) %>%
  mutate(jaccard_pfam = similarity(Pfam_P1, Pfam_P2))
all_duplicates_pfam_jaccard %<>% filter(!(is.na(jaccard_pfam)))

#### Add data on interaction motifs ####

# Load summary table
interaction_data <- read_delim('Data/summary_PCA_results_per_pairs_2019-01-07_clean.csv',
                                                      delim = ',')

# Add data on interactions 
all_duplicates_pfam_interactions <- inner_join(x = all_duplicates_pfam_jaccard, y = interaction_data %>% select(P1, P2, motif_categories),
                                               by = c("P1" = "P1", "P2" = "P2"))


comparison_list <- list(c('HM', 'HM&HET'),
                        c('HM&HET', 'HM'))

cl <- colors()

p1 <- all_duplicates_pfam_interactions %>% 
  filter(motif_categories %in% c('HM', 'HM&HET')) %>%
  group_by(motif_categories, Dup_type) %>%
  mutate(group_mean = round(mean(jaccard_pfam), 2)) %>%
  ggplot(aes(x = motif_categories, y = jaccard_pfam, fill = Dup_type)) +
  # scale_fill_manual(values = c("yellow", "lightblue")) +
  scale_fill_manual(values = c(cl[144], cl[129])) +
  facet_wrap(~Dup_type) +
  geom_violin(alpha = 0.5) +
  geom_jitter(width = 0.25, alpha = 0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        # strip.text.x = element_text(size = 15, face = 'bold'),
        # strip.background = element_rect(fill = 'white'),
        panel.spacing.x = unit(0, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  stat_compare_means(comparisons = comparison_list, method = "wilcox.test", 
                     paired = FALSE, position = 'top', label.y.npc = 1,
                     label.y = 1.05) +
  xlab("Motif") + ylab("Pfam term similarity") + labs(fill = '') +
  stat_summary(geom = 'point', fun.y = mean, colour = 'red') +
  geom_text(aes(y = group_mean - 0.02, label = group_mean), fontface = 'bold') +
  scale_x_discrete(labels = c('HM', 'HM&HET', 'HM', 'HM&HET'))
p1

#### Add info on sequence identity ####

seq_ident <- read_delim('Data/Paralogs_new_sequence_identities_2019-06-11_phylome_0003.txt', 
                       delim = '\t')

all_duplicates_pfam_seq_ident <- inner_join(x = all_duplicates_pfam_interactions, y = seq_ident %>% select(-Duplication),
                                                         by = c("P1" = "P1", "P2" = "P2"))

quasibinomial_smooth <- function(...){
  geom_smooth(method = "glm", method.args = list(family = "quasibinomial"), ...)
}

p2 <- all_duplicates_pfam_seq_ident %>%
  mutate(jaccard_pfam = jaccard_pfam) %>%
  filter(or(Sequence_identity > 20, Same_phylome == 1)) %>%
  filter(motif_categories %in% c('HM', 'HM&HET')) %>%
  ggplot(aes(x = Sequence_identity, y = jaccard_pfam, group = motif_categories, colour = motif_categories)) +
  geom_point() +
  facet_wrap(~Dup_type) +
  quasibinomial_smooth() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.line = element_line(size = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'grey'),
        panel.spacing = unit(2, "lines")) +
  xlab("Pairwise amino acid sequence identity (%)") + ylab("Pfam term similarity") + labs(fill = '') +
  scale_colour_manual(values = c("pink", "purple"))
p2

p_pfam_similarity <- plot_grid(p1, p2, nrow = 2, labels = c('A', 'B'))
ggsave(plot = p_pfam_similarity, filename = 'Figure_S7.pdf', dpi = 500, width = 7, height = 14, device = cairo_pdf)
