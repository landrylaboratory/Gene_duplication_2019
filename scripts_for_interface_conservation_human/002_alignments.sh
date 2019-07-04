#!/bin/bash

# This script performs the alignments of the longest isoform file to the
# sequences in the PDB to assign PDB structures to the pairs of paralogs.

blastp -outfmt 6 -query <path_to_longest_isoform_file> -db <path_to_pdb_seqres> -num_threads 4 -evalue 0.000001 -out Data/PDB_matches_human_paralogs.aln

