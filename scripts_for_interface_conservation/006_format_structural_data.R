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

# There are 10 pairs for which I have HET!
sum(!(is.na(final_PDB_matches$HET)))

final_PDB_matches$HET <- ifelse(is.na(final_PDB_matches$HET), 
                                  0,
                                  1)
table(final_PDB_matches$Monomer_P1, final_PDB_matches$HM_P1)
#     0   1
# 0 128  45
# 1  77  14

table(final_PDB_matches$Monomer_P2, final_PDB_matches$HM_P2)
#     0   1
# 0 137  48
# 1  63  16

# There are 36 pairs for which I could only find complexes with a third protein
complex_with_others <- final_PDB_matches %>% filter(Monomer_P1 == 0,
                             Monomer_P2 == 0,
                             HM_P1 == 0,
                             HM_P2 == 0,
                             HET == 0)

table(final_PDB_matches$Duplication_type)

relevant_complexes <- final_PDB_matches %>% 
  filter(!(P1_ID %in% complex_with_others$P1_ID))
table(relevant_complexes$Duplication_type)
# This many pairs have at least a monomeric structure or a structure for
# one of the complexes we are looking for
# SSD WGD 
# 142  86

complexes_checked <- relevant_complexes %>% 
  filter(or(HM_P1 == 1, or(HM_P2 == 1, HET == 1)))
table(complexes_checked$Duplication_type)
# This many pairs have at least an HM or a HET of paralogs
# SSD WGD 
# 81  40 

# Let's just write the complete table with this format
write.table(x = final_PDB_matches, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/paralogs_PDB_structures_reformatted.txt')

#### Let's do the same reformatting for the matrix that has strict HMs ####
final_PDB_matches <- read.table('Data/paralogs_PDB_structures_with_strict.txt',
                                h = T, sep = '\t')

struct_one_monomer_oneHM <- final_PDB_matches %>% filter(
  or(and(!(is.na(Monomer_P1)), !(is.na(HM_P1))),
     and((!is.na(Monomer_P2)), !(is.na(HM_P2)))
  )
)

# Write them to check them later
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


# Let's just write the complete table with this format
write.table(x = final_PDB_matches, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/paralogs_PDB_structures_reformatted_with_strict.txt')

#### Let's check the counts of structures we found ####

table(final_PDB_matches$Duplication_type)
# SSD WGD 
# 156  91

# There are 10 pairs for which I have HET!
sum(final_PDB_matches$HET)

table(final_PDB_matches$Monomer_P1, final_PDB_matches$HM_P1)
#     0   1
# 0 127  47
# 1  61  12



table(final_PDB_matches$Monomer_P2, final_PDB_matches$HM_P2)
#     0   1
# 0 133  53
# 1  50  11

# There are 41 pairs for which I could only find complexes with a third protein
complex_with_others <- final_PDB_matches %>% filter(Monomer_P1 == 0,
                                                    Monomer_P2 == 0,
                                                    HM_P1 == 0,
                                                    HM_P2 == 0,
                                                    HET == 0)

table(final_PDB_matches$Duplication_type)

relevant_complexes <- final_PDB_matches %>% 
  filter(!(P1_ID %in% complex_with_others$P1_ID))
table(relevant_complexes$Duplication_type)
# This many pairs have at least a monomeric structure or a structure for
# one of the complexes we are looking for
# SSD WGD 
# 130  76

complexes_checked <- relevant_complexes %>% 
  filter(or(HM_P1 == 1, or(HM_P2 == 1, HET == 1)))
table(complexes_checked$Duplication_type)
# This many pairs have at least an HM or a HET of paralogs
# SSD WGD 
# 81  40 

# Let's see for how many proteins we found two HMs
two_HMs <- relevant_complexes %>% filter(and(HM_P1 == 1, HM_P2 == 1))

table(final_PDB_matches$Strict_HM_P1)
#   0   1 
# 205  42 

table(final_PDB_matches$Strict_HM_P2)
#   0   1 
# 208  39 

# How many of them have HM?
table(relevant_complexes$HM_P1)
#   0   1 
# 147  59

table(relevant_complexes$HM_P2)
#   0   1 
# 142  64

table(relevant_complexes$HM_P1, relevant_complexes$HM_P2)
#    0  1
# 0 91 56
# 1 51  8

