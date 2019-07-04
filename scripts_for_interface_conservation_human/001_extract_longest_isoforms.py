#### This script will process the human reference proteome to obtain the longest isoforms
#### for duplicated genes as to BLAST them to the set of PDB sequences.

# Load libraries
import csv
import re
import os
import glob
from collections import OrderedDict
import math
from Bio.PDB import *
from Bio import SeqIO
from Bio.Seq import *
from Bio.SeqRecord import *
from Bio.Align.Applications import MuscleCommandline
from shutil import copyfile
import copy

human_proteome = SeqIO.parse(<path_to_reference_proteome>, 'fasta')

# Initialize the dictionary
proteome_dict = OrderedDict()

for record in human_proteome:
    for entry in record.description.split(' '):
        # Extract the entry with the gene ID
        if entry.startswith('gene:'):
            # Retrieve the gene ID without the '.X' at the end
            gene_id = entry.split(':')[1].split('.')[0]
            
            if type(proteome_dict.get(gene_id, -1)) == SeqRecord:
                # Make sure I save the longest isoform
                if len(proteome_dict[gene_id].seq) < len(record.seq):
                    proteome_dict[gene_id] = record
            else:
                # Save the first sequence I see for this gene ID
                proteome_dict[gene_id] = record

# Save the list of paralogs as a dictionary so that I can access them quickly
paralog_dict = OrderedDict()
paralog_table_handle = open('Data/dparalog_uniprotpdbids_gene_ids.tsv', 'r')
paralog_table = csv.reader(paralog_table_handle, delimiter = '\t')

# Skip header
header = paralog_table.next()

# Loop through all the entries and add the IDs of the paralogs to the dictionary
for entry in paralog_table:
    p1 = entry[1]
    p2 = entry[2]
    
    if paralog_dict.get(p1, -1) == -1:
        paralog_dict[p1] = 1
    if paralog_dict.get(p2, -1) == -1:
        paralog_dict[p2] = 1
    
records = []
for gene_id, record in proteome_dict.items():
    if paralog_dict.get(gene_id, -1) == 1:
        records.append(record)
    
SeqIO.write(records, <path_to_longest_isoform_file>, 'fasta')



