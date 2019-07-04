###########################################################
#### This script will be used to identify the matching ####
#### structures from the PDB for each of the paralogs. ####
###########################################################

library(tidyverse)
library(magrittr)

# Load the data for the original alignment
# datos <- read.table('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/alignment_to_PDB.aln',
# h = F)

# Load the data for the new alignments with all SSD and WGD proteins
# datos_WGD <- read.table('/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/WGD_alignment_to_PDB.aln',
#                         h = F)
data_human <- read.table('/home/axelle/Dropbox/Hiver2019/Paralog_interference/Human_paralogs/001_PDB_alignment/PDB_matches_human_paralogs.aln')

# colnames(datos) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
#                          'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
# colnames(datos_WGD) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
#                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
colnames(data_human) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                         'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

# Get all the lines for which the pDB complex has 100% identity to the paralog's sequence,
# that is, the PDB files that can certainly be assigned to a paralog.
# PDB_matches <- datos %>% filter(pident == 100)
# PDB_matches_WGD <- datos_WGD %>% filter(pident == 100)
PDB_matches_human <- data_human %>% filter(pident == 100)

#### Save this formatting part to a function ####
format_matches <- function(PDB_matches){

  # I will separate the PDB_ID from the chain within the PDB file.
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

# PDB_matches <- format_matches(PDB_matches)
# PDB_matches_WGD <- format_matches(PDB_matches_WGD)
PDB_matches_human <- format_matches(PDB_matches_human)

# There are 5 structures for which the best alignment is under 50 bp long
# Previously, I had foud 9 of them. This could be due to the fact that I did not use my reference proteome
# test <- PDB_matches %>% group_by(qseqid) %>% summarise(best_hit = max(length))
# bad <- test %>% filter(best_hit < 50)

# There are only 8 proteins with an alignment length < 50 aa for SSD paralogs
# and 7 for WGD paralogs.
# bad_WGD <- PDB_matches_WGD %>% group_by(qseqid) %>% summarise(best_hit = max(length)) %>% filter(best_hit < 50)
bad_human <- PDB_matches_human %>% group_by(qseqid) %>% summarise(best_hit = max(length)) %>% filter(best_hit < 50)

# Read the table with the pairs of paralogs
# paralog_table <- read.table('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/Data/pairs_seq_2018_08_09_AM.csv',
#                            h = T)

# Read the tables for paralogs from SSD and WGD
# pairs_WGD <- read.csv2('/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/000_data/WGD.csv', h = F)
pairs_human <- read.table('/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/000_data/duplication_SDS_1paire.txt', h = F)

# If I restrict the matches to the QSbio dataset, I get only 26 structures. If I use all the matches from
# the alignment, I find structures for 117 IDs. However, this would require me to check manually the structures to make sure
# they don't have only one chain.
# test <- PDB_matches %>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE))
# best_hits <- test %>% filter(rank == 1)

best_hits_SSD <- PDB_matches_SSD %>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE)) %>% filter(rank == 1)
best_hits_WGD <- PDB_matches_WGD %>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE)) %>% filter(rank == 1)

# Let's try to work with all 117 best hits. I can extract the list of structures and download the biological assemblies from
# the PDB. Afterwards, I can just count the number of chains.
# Now, it is just a matter of organizing the data and going back to the table that had the pairs of paralogs together
# and identifying the best structure for each paralog pair.

# Assign column names to the tables of paralogs
colnames(pairs_SSD) <- c('P1', 'P2')
colnames(pairs_WGD) <- c('P1', 'P2')

# I will join the best hits into pairs of paralogs by saving the followig lines as a function
join_best_hits_pairs <- function(best_hits, paralog_table){
  # Rename qseqid as P1
  # colnames(best_hits)[1] = 'P1'
  
  # Retrieve the second paralogs when the one I matched to a PDB structure was bait
  best_hits_P1 <- left_join(x = best_hits, y = paralog_table, by = c("qseqid" = "P1"))
  colnames(best_hits_P1)[1] <- 'P1'
  colnames(best_hits_P1)[15] <- 'P2'
  # best_hits_P1_clean <- best_hits_P1[,c(1:5, 15:21)]
  # colnames(best_hits_P1_clean)[1] = 'P1'
  # colnames(best_hits_P1_clean)[12] = 'P2'
  
  # Retrieve the second paralogs when the one I matched to a PDB structure was prey
  best_hits_P2 <- left_join(x = best_hits, y = paralog_table, by = c("qseqid" = "P2"))
  colnames(best_hits_P2)[1] <- 'P1'
  colnames(best_hits_P2)[15] <- 'P2'
  # best_hits_P2_clean <- best_hits_P2[,c(1:5, 15:21)]
  # colnames(best_hits_P2_clean)[1] = 'P1'
  # colnames(best_hits_P2_clean)[12] = 'P2'
  
  # Now, we can just remove the rows in which we got NAs from each of the two joins and merge the data frames
  best_hits_pairs_P1 <- best_hits_P1 %>% ungroup() %>% filter(!(is.na(P2)))
  best_hits_pairs_P2 <- best_hits_P2 %>% ungroup() %>% filter(!(is.na(P2)))
  
  all_best_hits_pairs <- rbind(best_hits_pairs_P1, best_hits_pairs_P2)
  return(all_best_hits_pairs)
}

all_best_hits_pairs_SSD <- join_best_hits_pairs(best_hits_SSD, pairs_SSD)
all_best_hits_pairs_WGD <- join_best_hits_pairs(best_hits_WGD, pairs_WGD)

#### Let's also try to get a table with all the alignment results, not only the best hits ####

PDB_matches_SSD %<>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE))
PDB_matches_WGD %<>% group_by(qseqid) %>% mutate(rank = order(bitscore, decreasing = TRUE))

all_hits_pairs_SSD <- join_best_hits_pairs(PDB_matches_SSD, pairs_SSD)
all_hits_pairs_WGD <- join_best_hits_pairs(PDB_matches_WGD, pairs_WGD)

# I will write a function that will help me get the candidates for repetition
# Could it be possible that both paralogs align very well to the
# sequence in the PDB structure?
check_repeat_candidates <- function(all_best_hits_pairs){
  # Let's see if I get IDs on both columns
  check_P1 <- table(all_best_hits_pairs$P1)
  check_P2 <- table(all_best_hits_pairs$P2)
  
  # Let's see if I get pairs for which the same protein appears on both sides. This could mean that I found a structure for the
  # two paralogs of that pair, or that a gene has been duplicated twice
  names_P1 <- names(check_P1[check_P1 != 0])
  names_P2 <- names(check_P2[check_P2 != 0])
  check_repeat <- names_P1[names_P1 %in% names_P2]
  
  # Filter for entries that include the paralogs that could be repeated
  potential_repeats <- all_best_hits_pairs %>% filter(P1 %in% check_repeat)
  
  table_pot_repeats <- table(potential_repeats$sseqid_PDB)
  
  # Cases for which I could potentially not distinguish the two paralogs. This is looking for PDB structures that appear more than 
  repeats <- table_pot_repeats[table_pot_repeats > 1]
  check_repeats <- all_best_hits_pairs %>% filter(sseqid_PDB %in% names(repeats)) %>% arrange(sseqid_PDB, sseqid_chain)
  
  # When looking at the repeats after arranging them, we can see that some paralogs are mapped to the same chain in the same
  # structure. We can get them by seeing which combinations of sseqid_PDB and sseqid_chain map to more than one distinct
  # P1 ID
  
  # Check which chains from specific PDB IDs map to more than one protein
  undistinguishable <- check_repeats %>% group_by(sseqid_PDB, sseqid_chain) %>% summarise(uniq = n_distinct(P1)) %>% filter(uniq > 1)
  
  # Get all the information for such entries
  undistinguishable2 <- check_repeats %>% filter(and(sseqid_PDB %in% undistinguishable$sseqid_PDB, sseqid_chain %in% undistinguishable$sseqid_chain))
  
  return(undistinguishable2)
}

# For the best hits, this filter removes a pair of paralogs if they are undistinguishable in their best hits.
undistinguishable_SSD <- check_repeat_candidates(all_best_hits_pairs_SSD)
undistinguishable_WGD <- check_repeat_candidates(all_best_hits_pairs_WGD)

# For all the hits, this filter removes a pair of paralogs if they are undistinguishable in at least one of their hits.
undistinguishable_all_hits_SSD <- check_repeat_candidates(all_hits_pairs_SSD)
undistinguishable_all_hits_WGD <- check_repeat_candidates(all_hits_pairs_WGD)

#### Check which proteins are in the list of unspecific interactions ####
removed_prot <- read.table('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/Final_table/Prot_removed.csv',
                           h = T)
# removed_SSD <- all_best_hits_pairs_SSD %>% filter(or(P1 %in% removed_prot$prot_removed, P2 %in% removed_prot$prot_removed))
# removed_WGD <- all_best_hits_pairs_WGD %>% filter(or(P1 %in% removed_prot$prot_removed, P2 %in% removed_prot$prot_removed))


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

#### Some checks before writing the tables ####

# Check for paralogs that are mapped twice to structures
ugh <- table(final_set_SSD$P1)
ugh[ugh > 1]
# named integer(0), so there are no repeated PDB IDs

ugh <- table(final_set_WGD$P1)
ugh[ugh > 1]
# named integer(0)

# Let's check if there are PDB structures whose chains map to more than one protein
ugh <- table(final_set_SSD$sseqid_PDB)
ugh[ugh > 1]
# 3jaq 3jct 3o8o 4arz 4bzi 4xr7 5jcs 5mc6 5mrf 5nrl 5sva 5u5q 5w66 5wvk 5wyk 
# 2    2    2    2    2    2    2    5   12    2    3    3    3    5   13 

ugh <- table(final_set_WGD$sseqid_PDB)
ugh[ugh > 1]
# 4v7h 5m1j 5mc6 5t6r 5tgm 5wyk 
# 2    2    6   19    2    6 

# We can check how many of the entries we have are matched to proteins with unspecific interactions
# The most abundant groups are always the ones without unspecific interactions (0,0)
table(final_set_SSD$P1_unspecific, final_set_SSD$P2_unspecific)
#     0   1
# 0 152  14
# 1  20  30

table(final_set_WGD$P1_unspecific, final_set_WGD$P2_unspecific)
#    0  1
# 0 99  6
# 1 13 45

# I will concatenate the two tables by adding an extra column with the type of duplication
final_set_SSD$Dup_type <- rep('SSD', nrow(final_set_SSD))
final_set_WGD$Dup_type <- rep('WGD', nrow(final_set_WGD))

complete_set <- rbind(final_set_SSD, final_set_WGD)

#### Write the tables to files ####
# I can write them separately
# write.table(x = final_set_SSD,
#            file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/PDB_matches_SSD.txt',
#            sep = '\t', quote = F, row.names = F, col.names = T)

# write.table(x = final_set_WGD,
#             file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/PDB_matches_WGD.txt',
#             sep = '\t', quote = F, row.names = F, col.names = T)

# Or together
write.table(x = complete_set,
            file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/PDB_matches_SSD_WGD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Let's save the lists of undistinguishable paralogs, too
write.table(x = undistinguishable_SSD,
            file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/undistinguishable_SSD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

write.table(x = undistinguishable_WGD,
            file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/undistinguishable_WGD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)

# Save the whole table of matches without undistinguishable paralogs
# I will concatenate the two tables by adding an extra column with the type of duplication
final_set_all_SSD$Dup_type <- rep('SSD', nrow(final_set_all_SSD))
final_set_all_WGD$Dup_type <- rep('WGD', nrow(final_set_all_WGD))

complete_set_all <- rbind(final_set_all_SSD, final_set_all_WGD)
write.table(x = complete_set_all,
            file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/All_paralogs/001_PDB_alignment/PDB_matches_all_chains_SSD_WGD.txt',
            sep = '\t', quote = F, row.names = F, col.names = T)






#### Code I am not using anymore ####

# I ended up thinking that this function was really not a good idea when I can just use the merge function more efficiently

# Define my function to match tables
match_columns <- function(query, target_lookup, target_result, unique_bool){
  dummy <- data.frame(lookup = target_lookup, result = target_result)
  dummy$lookup <- as.character(dummy$lookup)
  dummy$result <- as.character(dummy$result)
  values <- dummy %>% filter(lookup == query)
  if(is.na(values) || nrow(values) == 0){
    final_value <- NA
  }else{
    if(unique_bool){
      final_value <- str_c(unique(values$result), collapse = ';')
    }else{
      final_value <- str_c(values$result, collapse = ';')
    }
    return(final_value)
  }
}

qsbio_matches$Gene_id <- lapply(X = as.character(qsbio_matches$V2),
                                FUN = match_columns,
                                PDB_matches$sseqid_PDB, 
                                PDB_matches$qseqid, TRUE)
qsbio_matches <- qsbio_matches[,c(10, 1:9)]
colnames(qsbio_matches) <- c('Gene_id', 'Bio_assembly', 'Pdb_id', 'PiQSi_score', 'QSalign_score',
                             'Symmetry', 'Number_of_subunits', 'h_90', 'Error_prob', 'Confidence')

#### Some lines used to check which files were not downloaded by the PDB server tool ####

list_downloaded <- read.table('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/list_downloaded.txt', h = F)
list_wanted <- read.table('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/final_set_noQSbio_list.txt', h = F)
list_wanted$V1[which(!(list_wanted$V1 %in% list_downloaded$V1))]

#### Filter for the list of structures I downloaded ####

available_structures <- final_set %>% filter(sseqid_PDB %in% list_downloaded$V1) %>% arrange(sseqid_PDB)
pairs <- table(available_structures$paire)
pairs[pairs > 0]
length(pairs[pairs > 0])

# Get IDs that are repeated. However, when checking the structures both of them mapped to monomers.
repeated_ids <- final_set %>% filter(P1 %in% names(table(available_structures$P1)[table(available_structures$P1) > 1]))

# Write available structures
write.table(x = available_structures, file = '/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/available_structures.csv',
            quote = F, sep = ',', row.names = F, col.names = T)


#### Stuff I am not using ####

# I used to check which PDB structures were present in QSbio based on filtering done by another script. I was only looking for PDB structures
# matching proteins from the PCA tests, though.
qsbio_matches <- read.csv2('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/QSbio_matches_rerun.csv',
                           h = F)


colnames(qsbio_matches) <- c('Bio_assembly', 'Pdb_id', 'PiQSi_score', 'QSalign_score',
                             'Symmetry', 'Number_of_subunits', 'h_90', 'Error_prob', 'Confidence')

# I can filter out the PDB IDs that have only one subunit
qsbio_matches %<>% filter(Number_of_subunits > 1)

qsbio_matches$Gene_id <- lapply(X = as.character(qsbio_matches$V2),
                                FUN = match_columns,
                                PDB_matches$sseqid_PDB, 
                                PDB_matches$qseqid, TRUE)

# Indicate all the proteins that have structures in the PDB
merged_set <- inner_join(x = qsbio_matches, y = PDB_matches, by = c("Pdb_id" = "sseqid_PDB"))

# Some features of the data set
table(merged_set$Number_of_subunits)
# 2    3    4    6 
# 1078  164   76   18 

table(merged_set$Symmetry)
# C2   C3   C4   C6   C7   D2   D3   D4   D5   D6   D8  NPS   NS Octa Tetr 
# 1070  146    0    0    0   76   18    0    0    0    0    0   26    0    0 

# Check for how many proteins I have a PDB file with more than one subunit annotated in the QSbio survey
check_ids <- table(merged_set$qseqid)
check_ids[check_ids != 0]
length(check_ids[check_ids != 0])

# Let's do the same for PDB matches
check_pdb_ids <- table(PDB_matches$qseqid)
check_pdb_ids <- check_pdb_ids[check_pdb_ids != 0]
length(check_pdb_ids[check_pdb_ids != 0])

#### As I had to correct the sort in this part, I should check the differences between the two tables ####

final_table_old <- read.table('/Users/angelcisneros/Documents/Ete2018/Paralog_interference/Axelles_paralogs/PDB_search/rerun_2018-11-12/final_set_noQSbio.tab', 
                              sep = '\t', h =T) %>% arrange(sseqid_PDB, sseqid_chain, P1)

common_entries <- intersect(final_table_old$P1, final_set$P1)
