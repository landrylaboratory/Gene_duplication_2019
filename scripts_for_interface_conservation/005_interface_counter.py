
# coding: utf-8

# # Interface counter
# 
# This script will receive a FASTA file with the interaces as lowercase residues and will return the counts of residues from each chain that are a part of interfaces and the total counts of residues.
# 
# 

# In[12]:


# Load libraries
from Bio import SeqIO
import glob
import csv


# In[13]:


# Parse the file
handle = SeqIO.parse('/home/axelle/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/013_FASTA_interfaces/dist_regions_1a3w.fasta', 'fasta')

for sequence in handle:
    print sequence.description
    description_list = sequence.description.split(' | ')
    chain = description_list[1].split(' ')[1]
    print chain
    print sequence.seq


# In[14]:


# A function to count lowercase letters
# Taken from https://stackoverflow.com/questions/10953189/count-lower-case-characters-in-a-string
def n_lower_chars(string):
    return sum([c.islower() for c in string])


# In[15]:


handle = SeqIO.parse('/home/axelle/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/013_FASTA_interfaces/dist_regions_1a3w.fasta', 'fasta')

# Prepare a list to contain all the information
interface_list = []

for sequence in handle:
    description_list = sequence.description.split(' | ')
    chain = description_list[1].split(' ')[1]
    interface_count = n_lower_chars(sequence.seq)
    all_seq_count = len(sequence.seq)
    interface_list.append(', '.join([chain, str(interface_count), str(all_seq_count)]))

'; '.join(interface_list)


# Now that I have the code to count the interface residues in a FASTA file, I can just extend this to all of the files in the folder.

# In[16]:


outfile_handle = open('/home/axelle/Documents/Hiver2019/Paralog_interference/Interface_counts/interface_counts.txt', 'w')
writer = csv.writer(outfile_handle, delimiter = '\t')

file_list = glob.glob('/home/axelle/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/013_FASTA_interfaces/*fasta')

for fasta_file in file_list:
    prot = fasta_file.split('/')[-1].split('_')[2][0:4]
    handle = SeqIO.parse(fasta_file, 'fasta')
    
    interface_list = []

    for sequence in handle:
        description_list = sequence.description.split(' | ')
        chain = description_list[1].split(' ')[1]
        interface_count = n_lower_chars(sequence.seq)
        all_seq_count = len(sequence.seq)
        interface_list.append(', '.join([chain, str(interface_count), str(all_seq_count)]))

    counts = '; '.join(interface_list)
    
    new_row = [prot, counts]
    writer.writerow(new_row)
    
outfile_handle.close()
    

