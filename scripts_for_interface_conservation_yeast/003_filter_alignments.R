###########################################################
####          003_filter_alignments                    ####
#### This script will be used to identify the matching ####
#### structures from the PDB for each of the paralogs. ####
###########################################################

library(tidyverse)
library(magrittr)

# Set the working directory
setwd('/path/to/scripts_for_interface_conservation')

# Load the data for the new alignments with all SSD and WGD proteins
datos_WGD <- read.table('Data/WGD_alignment_to_PDB.aln',
                        h = F)
datos_SSD <- read.table('Data/SSD_alignment_to_PDB.aln',
                        h = F)

colnames(datos_WGD) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                     'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
colnames(datos_SSD) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                         'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

# Get all the lines for which the pDB complex has 100% identity to the paralog's sequence,
# that is, the PDB files that can certainly be assigned to a paralog.
PDB_matches_WGD <- datos_WGD %>% filter(pident == 100)
PDB_matches_SSD <- datos_SSD %>% filter(pident == 100)

#### Save this formatting part to a function ####
format_matches <- function(PDB_matches){

  # Separate the PDB_ID from the chain within the PDB file.
  sseqid_PDB <- rep(0, nrow(PDB_matches))
  sseqid_chain <- rep(0, nrow(PDB_matches))
  
  for(pos in 1:length(PDB_matches$sseqid)){
    id <- as.character(PDB_matches$sseqid[pos])
    split_list <- strsplit(id,'_')
    
    sseqid_PDB[pos] <- split_list[[1]][1]
    sseqid_chain[pos] <- split_list[[1]][2]
  }
  
  PDB_matches$sseqid_PDB <- sseqid_PDB
  PDB_matches$sseqid_chain <- sseqid_chain
  
  # Rearrange the columns
  PDB_matches <- PDB_matches[,c(1,13,14, 3:12)]
  return(PDB_matches)
}

PDB_matches_SSD <- format_matches(PDB_matches_SSD)
PDB_matches_WGD <- format_matches(PDB_matches_WGD)

# Read the tables for paralogs from SSD and WGD
pairs_SSD <- read.table('Data/duplication_SDS_1paire.txt', h = F)
pairs_WGD <- read.csv2('Data/WGD.csv', h = F)

# Assign column names to the tables of paralogs
colnames(pairs_SSD) <- c('P1', 'P2')
colnames(pairs_WGD) <- c('P1', 'P2')

# I will join the hits into pairs of paralogs by saving the following lines as a function
join_best_hits_pairs <- function(best_hits, paralog_table){
  
  # Retrieve the second paralogs when the one matched to a PDB structure was P1
  best_hits_P1 <- left_join(x = best_hits, y = paralog_table, by = c("qseqid" = "P1"))
  colnames(best_hits_P1)[1] <- 'P1'
  colnames(best_hits_P1)[15] <- 'P2'

  
  # Retrieve the second paralogs when the one matched to a PDB structure was P2
  best_hits_P2 <- left_join(x = best_hits, y = paralog_table, by = c("qseqid" = "P2"))
  colnames(best_hits_P2)[1] <- 'P1'
  colnames(best_hits_P2)[15] <- 'P2'
  
  # Remove the rows in which we got NAs from each of the two joins and merge the data frames
  best_hits_pairs_P1 <- best_hits_P1 %>% ungroup() %>% filter(!(is.na(P2)))
  best_hits_pairs_P2 <- best_hits_P2 %>% ungroup() %>% filter(!(is.na(P2)))
  
  all_best_hits_pairs <- rbind(best_hits_pairs_P1, best_hits_pairs_P2)
  return(all_best_hits_pairs)
}

#### Get a table with all the alignment results ####

PDB_matches_SSD %<>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE))
PDB_matches_WGD %<>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE))

all_hits_pairs_SSD <- join_best_hits_pairs(PDB_matches_SSD, pairs_SSD)
all_hits_pairs_WGD <- join_best_hits_pairs(PDB_matches_WGD, pairs_WGD)

# Check if both paralogs map well to the same structure
check_repeat_candidates <- function(all_best_hits_pairs){

  check_P1 <- table(all_best_hits_pairs$P1)
  check_P2 <- table(all_best_hits_pairs$P2)
  
  # Look for paralogs that appear on both sides of the table
  names_P1 <- names(check_P1[check_P1 != 0])
  names_P2 <- names(check_P2[check_P2 != 0])
  check_repeat <- names_P1[names_P1 %in% names_P2]
  
  # Filter for entries that include the paralogs that could be repeated
  potential_repeats <- all_best_hits_pairs %>% filter(P1 %in% check_repeat)
  
  table_pot_repeats <- table(potential_repeats$sseqid_PDB)
  
  repeats <- table_pot_repeats[table_pot_repeats > 1]
  check_repeats <- all_best_hits_pairs %>% filter(sseqid_PDB %in% names(repeats)) %>% arrange(sseqid_PDB, sseqid_chain)
  
  # Check which chains from specific PDB IDs map to more than one protein
  undistinguishable <- check_repeats %>% group_by(sseqid_PDB, sseqid_chain) %>% summarise(uniq = n_distinct(P1)) %>% filter(uniq > 1)
  
  # Get all the information for such entries
  undistinguishable2 <- check_repeats %>% filter(and(sseqid_PDB %in% undistinguishable$sseqid_PDB, sseqid_chain %in% undistinguishable$sseqid_chain))
  
  return(undistinguishable2)
}

# Look for undistinguishable hits
undistinguishable_all_hits_SSD <- check_repeat_candidates(all_hits_pairs_SSD)
undistinguishable_all_hits_WGD <- check_repeat_candidates(all_hits_pairs_WGD)

#### Check which proteins are in the list of inspecific interactions ####
removed_prot <- read.table('Data/Prot_removed.csv',
                           h = T)

# Write the tables of matches for SSD and WGD after removing the undistinguishable structures, 
# add columns indicating if the chains are in the sets of unspecific interactions
final_set_SSD <- all_best_hits_pairs_SSD %>%
  filter(!(P1 %in% undistinguishable_SSD$P1)) %>% # Remove undistinguishable pairs
  arrange(sseqid_PDB, sseqid_chain, P1) %>% # Organize the data
  mutate(P1_unspecific = ifelse(P1 %in% removed_prot$prot_removed, 1, 0), P2_unspecific = ifelse(P2 %in% removed_prot$prot_removed, 1, 0)) # Add columns about unspecific interactions

final_set_WGD <- all_best_hits_pairs_WGD %>%
  filter(!(P1 %in% undistinguishable_WGD$P1)) %>%
  arrange(sseqid_PDB, sseqid_chain, P1) %>%
  mutate(P1_unspecific = ifelse(P1 %in% removed_prot$prot_removed, 1, 0), P2_unspecific = ifelse(P2 %in% removed_prot$prot_removed, 1, 0))

final_set_all_SSD <- all_hits_pairs_SSD %>% 
  filter(!(P1 %in% undistinguishable_all_hits_SSD$P1)) %>%
  arrange(sseqid_PDB, sseqid_chain, P1) %>%
  mutate(P1_unspecific = ifelse(P1 %in% removed_prot$prot_removed, 1, 0), P2_unspecific = ifelse(P2 %in% removed_prot$prot_removed, 1, 0))

final_set_all_WGD <- all_hits_pairs_WGD %>% 
  filter(!(P1 %in% undistinguishable_all_hits_WGD$P1)) %>%
  arrange(sseqid_PDB, sseqid_chain, P1) %>%
  mutate(P1_unspecific = ifelse(P1 %in% removed_prot$prot_removed, 1, 0), P2_unspecific = ifelse(P2 %in% removed_prot$prot_removed, 1, 0))

# Concatenate the two tables by adding an extra column with the type of duplication
final_set_SSD$Dup_type <- rep('SSD', nrow(final_set_SSD))
final_set_WGD$Dup_type <- rep('WGD', nrow(final_set_WGD))

complete_set <- rbind(final_set_SSD, final_set_WGD)

#### Write the tables to files ####

# All the data
write.table(x = complete_set,
            file = 'Data/PDB_matches_SSD_WGD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Lists of undistinguishable paralogs
write.table(x = undistinguishable_SSD,
            file = 'Data/undistinguishable_SSD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

write.table(x = undistinguishable_WGD,
            file = 'Data/undistinguishable_WGD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Save the whole table of matches without undistinguishable paralogs
# Add the type of duplication
final_set_all_SSD$Dup_type <- rep('SSD', nrow(final_set_all_SSD))
final_set_all_WGD$Dup_type <- rep('WGD', nrow(final_set_all_WGD))

complete_set_all <- rbind(final_set_all_SSD, final_set_all_WGD)
write.table(x = complete_set_all,
            file = 'Data/PDB_matches_all_chains_SSD_WGD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)
