
# coding: utf-8

# # Gather PDB info
# 
# This script will create a table summarizing the information I can get from the PDB files. Such a final table will have the following columns:
# 
# - P1_ID
# - P2_ID
# - Duplication type
# - Monomer_P1 (the ID of a PDB structure that shows the P1 monomer or NA if there is none in my set)
# - Monomer_P2 (the ID of a PDB structure that shows the P2 monomer or NA if there is none in my set)
# - HM_P1 (the ID of a PDB structure whose biological assembly is not a monomer AND has at least two copies of P1)
# - HM_P2 (the ID of a PDB structure whose biological assembly is not a monomer AND has at least two copies of P2)
# - HET (the ID of a PDB structure whose biological assembly is not a monomer AND has at least one copy of each of P1 and P2)
# - P1_unspecific (1 if P1 is in the list of structures with unspecific interactions in Tarassov's data, 0 otherwise)
# - P2_unspecific (1 if P2 is in the list of structures with unspecific interactions in Tarassov's data, 0 otherwise)

# ## Overall plan
# 
# - Write a script to extract the following from each of the PDB files I could download
#     - REMARK 350: Number of subunits in assembly 1 (monomeric, dimeric, etc)
#     - Resolution (I can optionally use this as a filter, say make sure all the structures have a resolution of 3 Å or better)
# - Load table with data about the files I could not download
# - Load table with alignment data for all SSD, WGD hits
# - Apply filters
# - Match the tables

# In[1]:


# Load libraries
import csv
import glob
import re
from collections import OrderedDict


# ### Working with the information from the PDB files

# In[2]:


dict_text_2_numbers = {
     'MONOMERIC': 1, 
     'DIMERIC': 2,
     'TRIMERIC': 3,
     'TETRAMERIC': 4,
     'PENTAMERIC': 5,
     'HEXAMERIC': 6,
     'HEPTAMERIC': 7,
     'OCTAMERIC': 8,
     'NONAMERIC': 9,
     'DECAMERIC': 10,
     'UNDECAMERIC': 11,
     'DODECAMERIC': 12,
     'TRIDECAMERIC': 13,
     'TETRADECAMERIC': 14,
     'PENTADECAMERIC': 15,
     'HEXADECAMERIC': 16,
     'HEPTADECAMERIC': 17,
     'OCTADECAMERIC': 18,
     'NONADECAMERIC': 19,
     'EICOSAMERIC': 20
}


# In[3]:


def extract_pdb_data(pdb_file, dict_text_2_numbers):
    '''This function will receive the path to a PDB file and extract its information as a list with:
    - The PDB ID
    - The biological assembly assigned by the authors
    - The total number of subunits in that assembly
    - The structure's resolution, if applicable
    - A dictionary containing the IDs of those chains and how many times they appear in the assembly
    '''
    handle = open(pdb_file, 'r')
    pdb_id = pdb_file.split('/')[-1][0:4]
    resolution = 'NA'
    quit_bool = False
    
    # This dictionary will tell me how many times each chain is found in the selected biological assembly
    chains_dict = OrderedDict()
    
    # Loop through the lines to look for REMARK 350
    for line in handle:
        if line.startswith('EXPDTA'):
            # I will split the line on the experimental data with at least two spaces
            expdata = re.split(' [ ]+', line)[1]    
        if line.startswith('REMARK   2 RESOLUTION.'):
            # Extract the resolution
            res_match = re.search('([0-9\.]+)[ ]+ANGSTROM', line)
            if res_match:
                resolution = res_match.group(1)
        if line.startswith('REMARK 350'):
            # Then, I want to start reading.
            # I want to read what the biological assembly is, and then check what the one determined by the authors is
            # I will stop reading once I find one that was determined by the authors
            match_assembly = re.search('BIOMOLECULE:[ ]+([0-9]+)', line)
            if match_assembly:   
                if quit_bool:
                    # If I have all the data for the assembly selected by the authors and I find a new assembly, 
                    # I will stop and just quit
                    subunit_number = dict_text_2_numbers.get(subunit_number, subunit_number)
                    return [pdb_id, bio_assembly, subunit_number, expdata, resolution, chains_dict]
                else:
                    # Otherwise, I read the next one because I still have not found the one assigned by the authors
                    bio_assembly = int(match_assembly.group(1))
            
            # Now that we know with which biological assembly we are working, we can check if it is the one determined
            # by the authors.
            match_author = re.search('AUTHOR DETERMINED BIOLOGICAL UNIT: ([a-zA-Z0-9]+)', line)
            if match_author:
                subunit_number = match_author.group(1)
                
                # This means that we have now found the assembly determined by the authors
                quit_bool = True
            
            # I will need to check which chains are in this biological assembly
            match_chains_1 = re.search('APPLY THE FOLLOWING TO CHAINS: ([a-zA-Z0-9, ]+)', line)
            if match_chains_1:
                chains_assembly = match_chains_1.group(1).strip()
            
            # Sometimes the chains don't fit in a single line (example: 2ja7)
            match_chains_2 = re.search('AND CHAINS: ([a-zA-Z0-9, ]+)', line)
            if match_chains_2:
                chains_assembly = chains_assembly + ' ' + match_chains_2.group(1).strip()
            
            # Sometimes a single chain is used to obtain two chains (example: 3qps)
            # I need to look for the BIOMT line
            match_biomt = re.search('BIOMT\d   (\d)', line)
            if match_biomt:
                all_chains = chains_assembly.split(', ')
                for chain in all_chains:
                    chains_dict[chain] = int(match_biomt.group(1))
            
            
    # If I reach the end of the file without having an author determined unit, I will consider it to be monomeric
    # This is something I can check easily as to see if there are exceptions. 
    # An example of a PDB file that does not state the authors' preference is 1a1d, but it
    # was obtained with solution NMR
    if quit_bool:
        # Make sure I convert the assembly types from text to numbers
        subunit_number = dict_text_2_numbers.get(subunit_number, subunit_number)
        return [pdb_id, bio_assembly, subunit_number, expdata, resolution, chains_dict]
    else:
        return [pdb_id, 1, 1, expdata, resolution, chains_dict]


# In[4]:


# extract_pdb_data('/Users/intermilan1102/Dropbox/All_paralogs/002_PDB_structures/1qso.pdb')
print extract_pdb_data('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/002_PDB_structures/1a1d.pdb', dict_text_2_numbers)
print extract_pdb_data('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/002_PDB_structures/2ja7.pdb', dict_text_2_numbers)


# In[5]:


# Now let's just loop through all the structures to extract these data
file_list = glob.glob('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/002_PDB_structures/*pdb')
handle_out = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/PDB_structures.txt', 'w')
structure_dict = OrderedDict()

## I decided I would not keep writing to a file. Instead, I will save all these data in a big dictionary
# writer = csv.writer(handle_out, delimiter = '\t')
# header = ['PDB_ID', 'Biological_assembly', 'Number_of_subunits', 'Technique', 'Resolution']
# writer.writerow(header)
for pdb_file in file_list:
    # print pdb_file
    new_line = extract_pdb_data(pdb_file, dict_text_2_numbers)
    
    # The key will be the PDB ID and everything else will be in the values
    structure_dict[new_line[0]] = new_line[1:6]
#     writer.writerow(new_line)

# handle_out.close()


# In[6]:


structure_dict['2ja7']


# I will just need to remember that empty dictionaries will mean that there is no information for transformations on chains in the PDB file. As such, I will just consider every chain to appear only once.

# ### Loading the data for the structures I could not download

# In[7]:


# Load the table with all the info about the structures
handle_in = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/not_downloaded_reformatted.txt', 'r')
reader = csv.reader(handle_in, delimiter = '\t')

for line in reader:
    # I could also filter here based on the technique and the resolution if needed
    # I can ask Rong what he thinks about this
    
    pdb_id = line[0]
    bio_assembly = line[1]
    subunit_num = line[2]
    technique = line[3]
    resolution = line[4]
    comments = line[5]
    
    # Skip the first line
    if pdb_id == 'PDB_ID':
        continue
    
    # I will add all of this to the dictionary with the info on the downloaded structures
    structure_dict[pdb_id] = [int(bio_assembly), int(subunit_num), technique, resolution, OrderedDict()]
 


# In[8]:


structure_dict['1a3w']


# In[9]:


# I wanted to test this structure because YOR224C maps to chains H and T. These are on separate copies of the same
# complex, so only one of them should count when looking at the biological assembly. I got this right because
# assembly one only contains chain H, and so only chain H appears in this dictionary.
structure_dict['2ja7']


# For the alignment data, I would like to have the following structure:
# - Main Dictionary
#     - First level: Dictionary with pairs of paralogs sorted alphabetically as keys
#         - Second level: Dictionary with P1 and P2 as keys
#             - Third level: Dictionaries with PDB IDs that were matched by P1 and P2, respectively, as keys
#                 - Values: List of chains from each PDB ID that were matched by the respective paralog
#        

# In[10]:


# Load the alignment data as a dictionary
alignment_dict = OrderedDict()
handle_in = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/PDB_matches_all_chains_SSD_WGD.txt', 'r')
reader = csv.reader(handle_in, delimiter = '\t')

unspecific_dict = OrderedDict()
dup_type_dict = OrderedDict()

for line in reader:
    P1 = line[0]
    P2 = line[14]
    pair_list = [P1, P2]
    pair_list.sort()
    pair = tuple(pair_list)
    
    P1_match_PDB = line[1]
    P1_match_chain = line[2]
    
    P1_unspecific = line[15]
    P2_unspecific = line[16]
    
    dup_type = line[17]
    
    # Save the data about unspecific interactions and duplication types
    unspecific_dict[P1] = P1_unspecific
    unspecific_dict[P2] = P2_unspecific
    
    dup_type_dict[pair] = dup_type
    
    # Skip the first line
    if P1 == 'P1':
        continue
    
    # Now that I extracted all the info, I can start filling the dictionary
    if alignment_dict.get(pair, -1) == -1:
        # First level
        alignment_dict[pair] = OrderedDict()
        
        # Second level
        alignment_dict[pair][P1] = OrderedDict()
        alignment_dict[pair][P2] = OrderedDict()
        
        # Third level
        alignment_dict[pair][P1][P1_match_PDB] = [P1_match_chain]
    
    # Now, we can consider the case that the whole hierarchy is there but this is a new PDB structure
    elif alignment_dict[pair][P1].get(P1_match_PDB, -1) == -1:
        alignment_dict[pair][P1][P1_match_PDB] = [P1_match_chain]
        
    # Otherwise, everything is there and we just need to append to the list of matching chains
    else:
        alignment_dict[pair][P1][P1_match_PDB].append(P1_match_chain)
        


# Two examples of how the dictionary works. The first one shows the hierarchy and the second one shows how I have a list of four chains within the same structure that map to the same paralog.

# In[11]:


print alignment_dict.keys()[0]
print alignment_dict[('YNL172W', 'YOR224C')].keys()[0]
print alignment_dict[('YNL172W', 'YOR224C')]['YOR224C']
print alignment_dict[('YNL172W', 'YOR224C')]['YOR224C']['1a1d']

print alignment_dict[('YNL172W', 'YOR224C')]['YNL172W'].keys()


# In[12]:


print alignment_dict[('YDR256C', 'YGR088W')]['YDR256C']['1a4e']
print alignment_dict[('YDR256C', 'YGR088W')]['YGR088W'].keys()


# In[13]:


# This is just an example of how the unspecific dict works
unspecific_dict['YGR088W']


# In[14]:


# An example of the dictionary of the duplication types
dup_type_dict[('YDR256C', 'YGR088W')]


# ### Getting the final table

# This part should be a loop through the dictionaries that should allow me to look at each of the paralog pairs, then the structures that matched each of them, and the number of chains within each structure that mapped to them.
# 
# I will also use the dictionaries on the structures' data and the unspecific interactions to complete the table.

# The code for this should be similar to what I used to have but checking if the chains I am considering
# are a part of the biological assembly the authors selected and how many times they appear.

# In[15]:


len(structure_dict.keys())


# In[19]:


# handle_out = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures.txt', 'w')
handle_out = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures_complete2.txt', 'w')
writer = csv.writer(handle_out, delimiter = '\t')
header = ['P1_ID', 'P2_ID', 'Duplication_type','Monomer_P1', 'Monomer_P2', 'HM_P1', 'HM_P2', 'HET', 'other_HET_P1', 'other_HET_P2', 'P1_unspecific', 'P2_unspecific']
writer.writerow(header)

# Loop through the pairs of paralogs
for pair in alignment_dict.keys():
    P1 = pair[0]
    P2 = pair[1]
    
    dup_type = dup_type_dict[pair]
    P1_unspecific = unspecific_dict[P1]
    P2_unspecific = unspecific_dict[P2]

    monomers_P1_list = []
    monomers_P2_list = []
    HM_P1_list = []
    HM_P2_list = []
    other_HET_P1_list = []
    other_HET_P2_list = []
    HET_list = []
    
    # Let's work with P1 and the list of structures whose chains match it
    for structure, chains in alignment_dict[pair][P1].items():
        # I can check how many chains from that structure mapped to it
        # To know the number of matches, I need to look at each of the chains that this paralog matched
        # and how many times they appear in the biological assembly according to the structure dict.
        
        # Skip structures that were solved with solution NMR
        if structure_dict[structure][2] == 'SOLUTION NMR':
            continue
        # However, I must remember that if the dictionary is empty, every chain is considered a part
        # of the biological assembly. This is the case for some monomers and the big complexes.
        elif len(structure_dict[structure][4].keys()) == 0:
            matches_in_structure = len(chains)
        else:
            matches_in_structure = 0
            for chain in chains:
                # Let's see how many chains in the assembly derive from each of these matches
                # If this chain is not in the assembly, count it as a zero
                chain_matches = structure_dict[structure][4].get(chain, 0)
                matches_in_structure = matches_in_structure + chain_matches
        
        # I can now check if this is a monomer (the assembly has only one chain AND there is one match)
        if structure_dict[structure][1] == 1 and matches_in_structure == 1:
            # I will save the monomer for P1 as the PDB_ID followed by an underscore and the number of the assembly
            monomer_P1 = structure + '_' + str(structure_dict[structure][0])
            monomers_P1_list.append(monomer_P1)
        # or an HM (the assembly has more than one chain AND this paralog matches more than one chain)
        elif structure_dict[structure][1] > 1 and matches_in_structure > 1:
            HM_P1 = structure + '_' + str(structure_dict[structure][0])
            HM_P1_list.append(HM_P1)
        # or a HET but not of paralogs (other_HET)
        elif structure_dict[structure][1] > 1 and matches_in_structure == 1:
            other_HET_P1 = structure + '_' + str(structure_dict[structure][0])
            other_HET_P1_list.append(other_HET_P1)
    
    # Let's work with P2 and the list of structures whose chains match it
    for structure, chains in alignment_dict[pair][P2].items():
        # I can check how many chains from that structure mapped to it
        # To know the number of matches, I need to look at each of the chains that this paralog matched
        # and how many times they appear in the biological assembly according to the structure dict.
        
        # Skip structures that were solved with solution NMR
        if structure_dict[structure][2] == 'SOLUTION NMR':
            continue
        # However, I must remember that if the dictionary is empty, every chain is considered a part
        # of the biological assembly. This is the case for some monomers and the big complexes.
        elif len(structure_dict[structure][4].keys()) == 0:
            matches_in_structure = len(chains)
        else:
            matches_in_structure = 0
            for chain in chains:
                # Let's see how many chains in the assembly derive from each of these matches
                # If this chain is not in the assembly, count it as a zero
                chain_matches = structure_dict[structure][4].get(chain, 0)
                matches_in_structure = matches_in_structure + chain_matches
            
        # I can now check if this is a monomer (the assembly has only one chain AND there is only one match)
        if structure_dict[structure][1] == 1 and matches_in_structure == 1:
            # I will save the monomer for P2 as the PDB_ID followed by an underscore and the number of the assembly
            monomer_P2 = structure + '_' + str(structure_dict[structure][0])
            monomers_P2_list.append(monomer_P2)
        # or an HM (the assembly has more than one chain AND this paralog matches more than one chain)
        elif structure_dict[structure][1] > 1 and matches_in_structure > 1:
            HM_P2 = structure + '_' + str(structure_dict[structure][0])
            HM_P2_list.append(HM_P2)
        # or a HET but not of paralogs (other_HET)
        elif structure_dict[structure][1] > 1 and matches_in_structure == 1:
            other_HET_P2 = structure + '_' + str(structure_dict[structure][0])
            other_HET_P2_list.append(other_HET_P2)
    
    # Get the list of HET
    # Start with the list of structures that have a match of P1 and P2
    P1_matches = alignment_dict[pair][P1].keys()
    P2_matches = alignment_dict[pair][P2].keys()
    
    # I can now loop through each of the structures in P1_matches and check if it also has matches for P2
    for candidate in P1_matches:
        if candidate in P2_matches:
            chains_P1 = alignment_dict[pair][P1][candidate]
            chains_P2 = alignment_dict[pair][P2][candidate]
            total_chains = structure_dict[candidate][1]
            assembly = structure_dict[candidate][0]
            
            chains_P1_match = 0
            if len(structure_dict[candidate][4].keys()) == 0:
                chains_P1_match = len(chains_P1)
            else:
                for chain in chains_P1:
                    P1_chain_matches = structure_dict[candidate][4].get(chain, 0)
                    chains_P1_match = chains_P1_match + P1_chain_matches 

            chains_P2_match = 0
            if len(structure_dict[candidate][4].keys()) == 0:
                chains_P2_match = len(chains_P2)
            else:
                for chain in chains_P2:
                    P2_chain_matches = structure_dict[candidate][4].get(chain, 0)
                    chains_P2_match = chains_P2_match + P2_chain_matches             
            
            # Then, I just have to check if:
            # There is at least one match for P1 AND 
            # There is at least one match for P2 AND
            # The assembly has at least as many subunits as the sum of matches of P1 and P2
            if chains_P1_match >= 1 and chains_P2_match >= 1 and total_chains >= (chains_P1_match + chains_P2_match):
                HET = candidate + '_' + str(assembly)
                HET_list.append(HET)
                # I will remove them from the lists of other HET if they are there
                if HET in other_HET_P1_list:
                    other_HET_P1_list.remove(HET)
                if HET in other_HET_P2_list:
                    other_HET_P2_list.remove(HET)
    
    # Once I have all of this, I can just join all the elements from each of the lists and put everything together to
    # write the line. If the lists are empty, I will write NA instead
    save_bool = False
    if len(monomers_P1_list) == 0:
        monomers_P1 = 'NA'
    else:
        monomers_P1 = ','.join(monomers_P1_list)
        save_bool = True
        
    if len(monomers_P2_list) == 0:
        monomers_P2 = 'NA'
    else:    
        monomers_P2 = ','.join(monomers_P2_list)
        save_bool = True
        
    if len(HM_P1_list) == 0:
        HM_P1 = 'NA'
    else:
        HM_P1 = ','.join(HM_P1_list)
        save_bool = True
        
    if len(HM_P2_list) == 0:
        HM_P2 = 'NA'
    else:
        HM_P2 = ','.join(HM_P2_list)
        save_bool = True
        
    if len(HET_list) == 0:
        HET = 'NA'
    else:
        HET = ','.join(HET_list)
        save_bool = True
        
    if len(other_HET_P1_list) == 0:
        other_HET_P1 = 'NA'
    else:
        other_HET_P1 = ','.join(other_HET_P1_list)
        save_bool = True
        
    if len(other_HET_P2_list) == 0:
        other_HET_P2 = 'NA'
    else:
        other_HET_P2 = ','.join(other_HET_P2_list)
        save_bool = True
        
    # Now, I can put everything together and write the new row
    if save_bool == True:
        new_row = [P1, P2, dup_type, monomers_P1, monomers_P2, HM_P1, HM_P2, HET, other_HET_P1, other_HET_P2, P1_unspecific, P2_unspecific]
        writer.writerow(new_row)

handle_out.close()


# In[20]:


new_row


# In[21]:


# I will do some tests to see if this is working properly.
# The final pair of paralogs has no monomers, HM, or HET. The only match I could find was for YKR029C as part of
# a dimeric complex with a protein that is not YJL105W.
print alignment_dict[('YJL105W', 'YKR029C')]['YJL105W']
print alignment_dict[('YJL105W', 'YKR029C')]['YKR029C']
print structure_dict['5tdr']
print structure_dict['5tdw']


# In[22]:


# Let's test: YDR066C YER139C WGD     NA      5c2y_1  NA      NA      NA      0       1
# This is correctly assigned as a monomer for YER139C
print alignment_dict[('YDR066C', 'YER139C')]['YDR066C']
print alignment_dict[('YDR066C', 'YER139C')]['YER139C']
print structure_dict['5c2y']


# In[23]:


# Let's test: YBR014C YDL010W WGD     NA      3l4n_1  NA      5j3r_1  NA      0       0
# These are correctly assigned as a monomeric structure and a homodimer
print alignment_dict[('YBR014C', 'YDL010W')]['YBR014C']
print alignment_dict[('YBR014C', 'YDL010W')]['YDL010W']
print structure_dict['3l4n']
print structure_dict['5j3r']


# ### Let's write some code to check which of the HMs are strict homomers

# To do this, I will just need to look at the entries in the final table and check their subunits with the structure and alignment dictionaries I loaded previously.

# In[24]:


# Load the final data table
handle = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures_complete2.txt', 'r')
table_reader = csv.reader(handle, delimiter = '\t')

# Skip headers
header = table_reader.next()


# This is the block I originally used to look for strict homomers in the pairs that had no HET

# In[ ]:


for line in table_reader:
    
    # I am especially interested in pairs that have two HM but no HET. Are they strict HMs or occurrences of more
    # than one copy in bigger complexes?
    paralog_1 = line[0]
    paralog_2 = line[1]
    HM_1 = line[5]
    HM_2 = line[6]
    HET = line[7]
    
    if HM_1 != 'NA' and HM_2 != 'NA' and HET == 'NA':
        # print line
        
        HM_1 = HM_1.split(',')
        HM_2 = HM_2.split(',')
        
        print 'Checking', paralog_1, 'and', paralog_2
        
        # Then, I should look at each of the structures in HM_1 and HM_2 and their chains based on the
        # alignment and structure dictionaries.
        for structure in HM_1:
            # Check which chains in that structure correspond to paralog 1
            chains = alignment_dict[(paralog_1, paralog_2)][paralog_1][structure[0:4]]
            # I can get the total number of chains that come from those chains
            total = 0
            
            print structure
            
            for chain in chains:
                # Get the total number of times this chain appears in the biological assembly
                # Some might not appear because they could be present in the file but in a different assembly
                chain_appears = structure_dict[structure[0:4]][4].get(chain,0)
                
                # Count the total number of chains
                total = total + chain_appears
            
            # Check if it is a strict HM. This would be the case if all the chains that form the structure were matches
            if total == structure_dict[structure[0:4]][1]:
                print structure, 'is a strict HM for', paralog_1
            
        # Repeat for HM_2
        for structure in HM_2:
            # Check which chains in that structure correspond to paralog 1
            chains = alignment_dict[(paralog_1, paralog_2)][paralog_2][structure[0:4]]
            # I can get the total number of chains that come from those chains
            total = 0
            
            print structure
            
            for chain in chains:
                # Get the total number of times this chain appears in the biological assembly
                # Some might not appear because they could be present in the file but in a different assembly
                chain_appears = structure_dict[structure[0:4]][4].get(chain, 0)
                
                # Count the total number of chains
                total = total + chain_appears 
                
            # Check if it is a strict HM. This would be the case if all the chains that form the structure were matches
            if total == structure_dict[structure[0:4]][1]:
                print structure, 'is a strict HM for', paralog_2
        
        print '------'


# These are the blocks I used to fill the table with info on the strict homomers

# In[25]:


# Load the final data table
handle = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures_complete2.txt', 'r')
table_reader = csv.reader(handle, delimiter = '\t')

# Skip headers
header = table_reader.next()

# Prepare the file to write
handle_writer = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures_with_strict2.txt', 'w')
writer = csv.writer(handle_writer, delimiter = '\t')

header = ['P1_ID', 'P2_ID', 'Duplication_type','Monomer_P1', 'Monomer_P2', 'HM_P1', 'HM_P2', 'HET', 'other_HET_P1', 'other_HET_P2', 'P1_unspecific', 'P2_unspecific', 'Strict_HM_P1', 'Strict_HM_P2']   
writer.writerow(header)


# In[26]:


for line in table_reader:
    
    # Are they strict HMs or occurrences of more
    # than one copy in bigger complexes?
    paralog_1 = line[0]
    paralog_2 = line[1]
    HM_1 = line[5]
    HM_2 = line[6]
    HET = line[7]
    strict_HM_1 = []
    strict_HM_2 = []
    
    HM_1 = HM_1.split(',')
    HM_2 = HM_2.split(',')

    # I just need to make sure that these paralogs are not matched to NAs
    if HM_1[0] == 'NA':
        # No HMs means there are no strict HMs
        line.append('NA')
    else:
        # Then, I should look at each of the structures in HM_1 and HM_2 and their chains based on the
        # alignment and structure dictionaries.
        for structure in HM_1:
            # Check which chains in that structure correspond to paralog 1
            chains = alignment_dict[(paralog_1, paralog_2)][paralog_1][structure[0:4]]
            # I can get the total number of chains that come from those chains
            total = 0

            for chain in chains:
                # Get the total number of times this chain appears in the biological assembly
                # Some might not appear because they could be present in the file but in a different assembly
                chain_appears = structure_dict[structure[0:4]][4].get(chain,0)

                # Count the total number of chains
                total = total + chain_appears

            # Check if it is a strict HM. This would be the case if all the chains that form the structure were matches
            if total == structure_dict[structure[0:4]][1]:
                # Then we have a strict HM for paralog 1
                strict_HM_1.append(structure)
        
        # Now, we should just add the data on strict homomers to the line and write to a file
        # If there are strict homomers we write them to the column. Otherwise, we write NA.
        if len(strict_HM_1) > 0:
            strict_HM_1_final = ','.join(strict_HM_1)
            line.append(strict_HM_1_final)
        else:
            line.append('NA')

    # Repeat for HM_2
    # I just need to make sure that these paralogs are not matched to NAs
    if HM_2[0] == 'NA':
        # No HMs means there are no strict HMs
        line.append('NA')
    else:
        for structure in HM_2:
            # Check which chains in that structure correspond to paralog 1
            chains = alignment_dict[(paralog_1, paralog_2)][paralog_2][structure[0:4]]
            # I can get the total number of chains that come from those chains
            total = 0

            for chain in chains:
                # Get the total number of times this chain appears in the biological assembly
                # Some might not appear because they could be present in the file but in a different assembly
                chain_appears = structure_dict[structure[0:4]][4].get(chain, 0)

                # Count the total number of chains
                total = total + chain_appears 

            # Check if it is a strict HM. This would be the case if all the chains that form the structure were matches
            if total == structure_dict[structure[0:4]][1]:
                # Then we have a strict HM for paralog 2
                strict_HM_2.append(structure)

        if len(strict_HM_2) > 0:
            strict_HM_2_final = ','.join(strict_HM_2)
            line.append(strict_HM_2_final)
        else:
            line.append('NA')
    
    writer.writerow(line)

handle_writer.close()


# In[27]:


a = '1'
b = a.split('a')
print b


# In[28]:


print structure_dict['3w4y']
print alignment_dict['YGR029W','YPR037C']['YGR029W']['3w4y']


# In[29]:


print alignment_dict[('YDR423C', 'YML007W')]
print structure_dict['1sse']


# In[30]:


print alignment_dict[('YDR409W', 'YOR156C')]
print structure_dict['2rnn']


# ## Let's use a loop to select the complexes with the best resolution for interface analyses

# In[31]:


handle = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures_with_strict2.txt', 'r')
reader = csv.reader(handle, delimiter = '\t')

handle_writer = open('~/Documents/Hiver2019/Paper_duplication/Submission_eLife_AC/Data/Interface_conservation/003_data_tables/paralogs_PDB_structures_best_structures2.txt', 'w')
writer = csv.writer(handle_writer, delimiter = '\t')

header = ['P1_ID', 'P2_ID', 'Duplication_type', 'HM_P1', 'HM_P2', 'HET', 'P1_unspecific', 'P2_unspecific', 'Strict_HM_P1', 'Strict_HM_P2']   
writer.writerow(header)


# In[32]:


# Skip the first line
header = reader.next()

for line in reader:
    P1_ID = line[0]
    P2_ID = line[1]
    dup_type = line[2]
    
    HM_P1 = line[5]
    best_HM_P1 = ['NA', 1000]
    
    HM_P2 = line[6]
    best_HM_P2 = ['NA', 1000]
    
    HET = line[7]
    best_HET = ['NA', 1000]
    
    P1_unspecific = line[10]
    P2_unspecific = line[11]
    
    strict_HM_P1 = line[12]
    best_strict_HM_P1 = ['NA', 1000]
    
    strict_HM_P2 = line[13]
    best_strict_HM_P2 = ['NA', 1000]
    
    bool_save = False
    
    # We have to look at the list of structures in each column, and select the crystal with the best resolution
    if HM_P1 != 'NA':
        for structure in HM_P1.split(','):
            # Get the resolution of such structure if it is a crystal
            structure_info = structure_dict[structure[0:4]]
            technique = structure_info[2]
            resolution = float(structure_info[3])
            if technique == 'X-RAY DIFFRACTION' and resolution < best_HM_P1[1]:
                best_HM_P1 = [structure, resolution]
                bool_save = True
    
    
    if HM_P2 != 'NA':
        for structure in HM_P2.split(','):
            # Get the resolution of such structure if it is a crystal
            structure_info = structure_dict[structure[0:4]]
            technique = structure_info[2]
            resolution = float(structure_info[3])
            if technique == 'X-RAY DIFFRACTION' and resolution < best_HM_P2[1]:
                best_HM_P2 = [structure, resolution]
                bool_save = True
                
    if HET != 'NA':
        for structure in HET.split(','):
            # Get the resolution of such structure if it is a crystal
            structure_info = structure_dict[structure[0:4]]
            technique = structure_info[2]
            resolution = float(structure_info[3])
            if technique == 'X-RAY DIFFRACTION' and resolution < best_HET[1]:
                best_HET = [structure, resolution]
                bool_save = True
            
    if strict_HM_P1 != 'NA':
        for structure in strict_HM_P1.split(','):
            # Get the resolution of such structure if it is a crystal
            structure_info = structure_dict[structure[0:4]]
            technique = structure_info[2]
            resolution = float(structure_info[3])
            if technique == 'X-RAY DIFFRACTION' and resolution < best_strict_HM_P1[1]:
                best_strict_HM_P1 = [structure, resolution]
                bool_save = True
            
    if strict_HM_P2 != 'NA':
        for structure in strict_HM_P2.split(','):
            # Get the resolution of such structure if it is a crystal
            structure_info = structure_dict[structure[0:4]]
            technique = structure_info[2]
            resolution = float(structure_info[3])
            if technique == 'X-RAY DIFFRACTION' and resolution < best_strict_HM_P2[1]:
                best_strict_HM_P2 = [structure, resolution]
                bool_save = True
                
    # Now that we have all of them, we can just write the info to the table
    if bool_save:
        new_row = [P1_ID, P2_ID, dup_type, best_HM_P1[0], best_HM_P2[0], best_HET[0], P1_unspecific, P2_unspecific, best_strict_HM_P1[0], best_strict_HM_P2[0]]
        writer.writerow(new_row)
        
handle_writer.close()

