#### This script is useful for formatting the tables of PDB structures that matched ####
#### yeast paralogs from either WGD or SSD.                                         ####

library(tidyverse)
library(magrittr)
library(ggplot2)

#### Set the working directory ####

setwd('/path/to/scripts_for_interface_conservation')

#### Format the table of PDB matches ####

final_PDB_matches <- read.table('Data/paralogs_PDB_structures_complete2.txt',
                                h = T, sep = '\t')

final_PDB_matches$Monomer_P1 <- ifelse(is.na(final_PDB_matches$Monomer_P1), 
                                      0,
                                      1)
final_PDB_matches$Monomer_P2 <- ifelse(is.na(final_PDB_matches$Monomer_P2), 
                                       0,
                                       1)

final_PDB_matches$HM_P1 <- ifelse(is.na(final_PDB_matches$HM_P1), 
                                       0,
                                       1)

final_PDB_matches$HM_P2 <- ifelse(is.na(final_PDB_matches$HM_P2), 
                                  0,
                                  1)

final_PDB_matches$HET <- ifelse(is.na(final_PDB_matches$HET), 
                                  0,
                                  1)

complex_with_others <- final_PDB_matches %>% filter(Monomer_P1 == 0,
                             Monomer_P2 == 0,
                             HM_P1 == 0,
                             HM_P2 == 0,
                             HET == 0)

relevant_complexes <- final_PDB_matches %>% 
  filter(!(P1_ID %in% complex_with_others$P1_ID))

complexes_checked <- relevant_complexes %>% 
  filter(or(HM_P1 == 1, or(HM_P2 == 1, HET == 1)))

# Let's just write the complete table with this format
write.table(x = final_PDB_matches, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/paralogs_PDB_structures_reformatted.txt')

#### Reformat the matrix that has strict HMs ####
final_PDB_matches <- read.table('Data/paralogs_PDB_structures_with_strict.txt',
                                h = T, sep = '\t')

struct_one_monomer_oneHM <- final_PDB_matches %>% filter(
  or(and(!(is.na(Monomer_P1)), !(is.na(HM_P1))),
     and((!is.na(Monomer_P2)), !(is.na(HM_P2)))
  )
)

# Write them
write.table(x = struct_one_monomer_oneHM, file = 'Data/paralogs_PDB_structures_one_monomer_one_HM.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')

final_PDB_matches$Monomer_P1 <- ifelse(is.na(final_PDB_matches$Monomer_P1), 
                                       0,
                                       1)
final_PDB_matches$Monomer_P2 <- ifelse(is.na(final_PDB_matches$Monomer_P2), 
                                       0,
                                       1)

final_PDB_matches$HM_P1 <- ifelse(is.na(final_PDB_matches$HM_P1), 
                                  0,
                                  1)

final_PDB_matches$HM_P2 <- ifelse(is.na(final_PDB_matches$HM_P2), 
                                  0,
                                  1)

final_PDB_matches$Strict_HM_P1 <- ifelse(is.na(final_PDB_matches$Strict_HM_P1), 
                                  0,
                                  1)

final_PDB_matches$Strict_HM_P2 <- ifelse(is.na(final_PDB_matches$Strict_HM_P2), 
                                  0,
                                  1)

final_PDB_matches$HET <- ifelse(is.na(final_PDB_matches$HET), 
                                0,
                                1)

final_PDB_matches$other_HET_P1 <- ifelse(is.na(final_PDB_matches$other_HET_P1),
                                         0,
                                         1)
final_PDB_matches$other_HET_P2 <- ifelse(is.na(final_PDB_matches$other_HET_P2),
                                         0,
                                         1)


# Write the complete table with this format
write.table(x = final_PDB_matches, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/paralogs_PDB_structures_reformatted_with_strict.txt')

complex_with_others <- final_PDB_matches %>% filter(Monomer_P1 == 0,
                                                    Monomer_P2 == 0,
                                                    HM_P1 == 0,
                                                    HM_P2 == 0,
                                                    HET == 0)

relevant_complexes <- final_PDB_matches %>% 
  filter(!(P1_ID %in% complex_with_others$P1_ID))

complexes_checked <- relevant_complexes %>% 
  filter(or(HM_P1 == 1, or(HM_P2 == 1, HET == 1)))

two_HMs <- relevant_complexes %>% filter(and(HM_P1 == 1, HM_P2 == 1))
