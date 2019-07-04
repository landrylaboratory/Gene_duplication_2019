#!/bin/bash

# This script performs the alignments of the SSD and WGD pairs to the
# sequences in the PDB to assign PDB structures to the pairs of paralogs.

blastp -outfmt 6 -query Data/WGD_fasta.fa -db Data/pdb_seqres.txt -num_threads 4 -evalue 0.000001 -out Data/WGD_alignment_to_PDB.aln
blastp -outfmt 6 -query Data/SSD_fasta.fa -db Data/pdb_seqres.txt -num_threads 4 -evalue 0.000001 -out Data/SSD_alignment_to_PDB.aln
