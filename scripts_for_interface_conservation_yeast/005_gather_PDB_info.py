#!/usr/bin/env python
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


# Load libraries
import csv
import glob
import re
from collections import OrderedDict



# A folder to write the output
output_path = <path/to/output>

# A folder with the PDB files downloaded from the alignments
pdb_folder = <path/to/pdb/folder>


# ### Working with the information from the PDB files


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
    
    # This dictionary will save how many times each chain is found in the selected biological assembly
    chains_dict = OrderedDict()
    
    # Loop through the lines to look for REMARK 350
    for line in handle:
        if line.startswith('EXPDTA'):
            # Split the line on the experimental data with at least two spaces
            expdata = re.split(' [ ]+', line)[1]    
        if line.startswith('REMARK   2 RESOLUTION.'):
            # Extract the resolution
            res_match = re.search('([0-9\.]+)[ ]+ANGSTROM', line)
            if res_match:
                resolution = res_match.group(1)
        if line.startswith('REMARK 350'):

            # Stop reading on the one that was determined by the authors
            match_assembly = re.search('BIOMOLECULE:[ ]+([0-9]+)', line)
            if match_assembly:   
                if quit_bool:

                    subunit_number = dict_text_2_numbers.get(subunit_number, subunit_number)
                    return [pdb_id, bio_assembly, subunit_number, expdata, resolution, chains_dict]
                else:
                    
                    bio_assembly = int(match_assembly.group(1))
            
.
            match_author = re.search('AUTHOR DETERMINED BIOLOGICAL UNIT: ([a-zA-Z0-9]+)', line)
            if match_author:
                subunit_number = match_author.group(1)
                
                quit_bool = True
            
            # Check which chains are in this biological assembly
            match_chains_1 = re.search('APPLY THE FOLLOWING TO CHAINS: ([a-zA-Z0-9, ]+)', line)
            if match_chains_1:
                chains_assembly = match_chains_1.group(1).strip()
            
            # Sometimes the chains don't fit in a single line (example: 2ja7)
            match_chains_2 = re.search('AND CHAINS: ([a-zA-Z0-9, ]+)', line)
            if match_chains_2:
                chains_assembly = chains_assembly + ' ' + match_chains_2.group(1).strip()
            
            # Sometimes a single chain is used to obtain two chains (example: 3qps)
            # Look for the BIOMT line
            match_biomt = re.search('BIOMT\d   (\d)', line)
            if match_biomt:
                all_chains = chains_assembly.split(', ')
                for chain in all_chains:
                    chains_dict[chain] = int(match_biomt.group(1))
            
            
    # Structures without biological assemblies will be considered monomeric
    if quit_bool:
        # Make sure I convert the assembly types from text to numbers
        subunit_number = dict_text_2_numbers.get(subunit_number, subunit_number)
        return [pdb_id, bio_assembly, subunit_number, expdata, resolution, chains_dict]
    else:
        return [pdb_id, 1, 1, expdata, resolution, chains_dict]



# Loop through all the structures to extract these data
file_list = glob.glob(pdb_folder + '/*pdb')
handle_out = open(output_path + '/PDB_structures.txt', 'w')
structure_dict = OrderedDict()


for pdb_file in file_list:
    new_line = extract_pdb_data(pdb_file, dict_text_2_numbers)
    
    # The key will be the PDB ID and everything else will be in the values
    structure_dict[new_line[0]] = new_line[1:6]





# Load the alignment data as a dictionary
alignment_dict = OrderedDict()
handle_in = open('Data/PDB_matches_all_chains_SSD_WGD.txt', 'r')
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
    
    # Start filling the dictionary
    if alignment_dict.get(pair, -1) == -1:
        # First level
        alignment_dict[pair] = OrderedDict()
        
        # Second level
        alignment_dict[pair][P1] = OrderedDict()
        alignment_dict[pair][P2] = OrderedDict()
        
        # Third level
        alignment_dict[pair][P1][P1_match_PDB] = [P1_match_chain]
    
    # For new PDB structures for paralogs that were already added
    elif alignment_dict[pair][P1].get(P1_match_PDB, -1) == -1:
        alignment_dict[pair][P1][P1_match_PDB] = [P1_match_chain]
        
    # Append the chain to the list of chains
    else:
        alignment_dict[pair][P1][P1_match_PDB].append(P1_match_chain)
        



handle_out = open(output_path + '/paralogs_PDB_structures_complete2.txt', 'w')
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
    
    for structure, chains in alignment_dict[pair][P1].items():
        
        # Skip structures that were solved with solution NMR
        if structure_dict[structure][2] == 'SOLUTION NMR':
            continue
        elif len(structure_dict[structure][4].keys()) == 0:
            matches_in_structure = len(chains)
        else:
            matches_in_structure = 0
            for chain in chains:
                # Check the number of chains in the assembly that derive from this chain
                chain_matches = structure_dict[structure][4].get(chain, 0)
                matches_in_structure = matches_in_structure + chain_matches
        
        # Check if this is a monomer (the assembly has only one chain AND there is one match)
        if structure_dict[structure][1] == 1 and matches_in_structure == 1:
            
            monomer_P1 = structure + '_' + str(structure_dict[structure][0])
            monomers_P1_list.append(monomer_P1)
        # Check if this is a HM (the assembly has more than one chain AND this paralog matches more than one chain)
        elif structure_dict[structure][1] > 1 and matches_in_structure > 1:
            HM_P1 = structure + '_' + str(structure_dict[structure][0])
            HM_P1_list.append(HM_P1)
        # Check if this is a HET but not of paralogs (other_HET)
        elif structure_dict[structure][1] > 1 and matches_in_structure == 1:
            other_HET_P1 = structure + '_' + str(structure_dict[structure][0])
            other_HET_P1_list.append(other_HET_P1)
    
    # Work with P2 and the list of structures whose chains match it
    for structure, chains in alignment_dict[pair][P2].items():
        
        # Skip structures that were solved with solution NMR
        if structure_dict[structure][2] == 'SOLUTION NMR':
            continue
        elif len(structure_dict[structure][4].keys()) == 0:
            matches_in_structure = len(chains)
        else:
            matches_in_structure = 0
            for chain in chains:
                # Check the number of chains in the assembly that derive from this chain
                chain_matches = structure_dict[structure][4].get(chain, 0)
                matches_in_structure = matches_in_structure + chain_matches
            
        # Check if this is a monomer (the assembly has only one chain AND there is only one match)
        if structure_dict[structure][1] == 1 and matches_in_structure == 1:
            
            monomer_P2 = structure + '_' + str(structure_dict[structure][0])
            monomers_P2_list.append(monomer_P2)
        # Check if this is a HM (the assembly has more than one chain AND this paralog matches more than one chain)
        elif structure_dict[structure][1] > 1 and matches_in_structure > 1:
            HM_P2 = structure + '_' + str(structure_dict[structure][0])
            HM_P2_list.append(HM_P2)
        # Check if this is a HET but not of paralogs (other_HET)
        elif structure_dict[structure][1] > 1 and matches_in_structure == 1:
            other_HET_P2 = structure + '_' + str(structure_dict[structure][0])
            other_HET_P2_list.append(other_HET_P2)
    
    # Get the list of HET
    # Start with the list of structures that have a match of P1 and P2
    P1_matches = alignment_dict[pair][P1].keys()
    P2_matches = alignment_dict[pair][P2].keys()
    
    # Loop through each of the structures in P1_matches and check if it also has matches for P2
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
            
            # Check if:
            # There is at least one match for P1 AND 
            # There is at least one match for P2 AND
            # The assembly has at least as many subunits as the sum of matches of P1 and P2
            if chains_P1_match >= 1 and chains_P2_match >= 1 and total_chains >= (chains_P1_match + chains_P2_match):
                HET = candidate + '_' + str(assembly)
                HET_list.append(HET)
                # Remove these HET from the other_HET category if they are there
                if HET in other_HET_P1_list:
                    other_HET_P1_list.remove(HET)
                if HET in other_HET_P2_list:
                    other_HET_P2_list.remove(HET)
    
    # Assemble the table
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
        
    # NPut everything together and write the new row
    if save_bool == True:
        new_row = [P1, P2, dup_type, monomers_P1, monomers_P2, HM_P1, HM_P2, HET, other_HET_P1, other_HET_P2, P1_unspecific, P2_unspecific]
        writer.writerow(new_row)

handle_out.close()


# ### Look for strict homomers


# Load the final data table
handle = open(output_path + '/paralogs_PDB_structures_complete2.txt', 'r')
table_reader = csv.reader(handle, delimiter = '\t')

# Skip headers
header = table_reader.next()




# Load the final data table
handle = open(output_path +'/paralogs_PDB_structures_complete2.txt', 'r')
table_reader = csv.reader(handle, delimiter = '\t')

# Skip headers
header = table_reader.next()

# Prepare the file to write
handle_writer = open(output_path + '/paralogs_PDB_structures_with_strict2.txt', 'w')
writer = csv.writer(handle_writer, delimiter = '\t')

header = ['P1_ID', 'P2_ID', 'Duplication_type','Monomer_P1', 'Monomer_P2', 'HM_P1', 'HM_P2', 'HET', 'other_HET_P1', 'other_HET_P2', 'P1_unspecific', 'P2_unspecific', 'Strict_HM_P1', 'Strict_HM_P2']   
writer.writerow(header)




for line in table_reader:
    
    paralog_1 = line[0]
    paralog_2 = line[1]
    HM_1 = line[5]
    HM_2 = line[6]
    HET = line[7]
    strict_HM_1 = []
    strict_HM_2 = []
    
    HM_1 = HM_1.split(',')
    HM_2 = HM_2.split(',')

    if HM_1[0] == 'NA':
        line.append('NA')
    else:
        # Look at each of the structures in HM_1 and HM_2 and their chains based on the
        # alignment and structure dictionaries.
        for structure in HM_1:
            # Check which chains in that structure correspond to paralog 1
            chains = alignment_dict[(paralog_1, paralog_2)][paralog_1][structure[0:4]]
            # Get the total number of chains that come from those chains
            total = 0

            for chain in chains:
                # Get the total number of times this chain appears in the biological assembly
                # Some might not appear because they could be present in the file but in a different assembly
                chain_appears = structure_dict[structure[0:4]][4].get(chain,0)

                # Count the total number of chains
                total = total + chain_appears

            # Check if it is a strict HM. This would be the case if all the chains that form the structure were matches
            if total == structure_dict[structure[0:4]][1]:
                strict_HM_1.append(structure)
        
        # Now, add the data on strict homomers to the line and write to a file
        # If there are strict homomers write them to the column. Otherwise, write NA.
        if len(strict_HM_1) > 0:
            strict_HM_1_final = ','.join(strict_HM_1)
            line.append(strict_HM_1_final)
        else:
            line.append('NA')

    # Repeat for HM_2
    if HM_2[0] == 'NA':
        line.append('NA')
    else:
        for structure in HM_2:
            # Check which chains in that structure correspond to paralog 2
            chains = alignment_dict[(paralog_1, paralog_2)][paralog_2][structure[0:4]]
            # Get the total number of chains that come from those chains
            total = 0

            for chain in chains:
                # Get the total number of times this chain appears in the biological assembly
                # Some might not appear because they could be present in the file but in a different assembly
                chain_appears = structure_dict[structure[0:4]][4].get(chain, 0)

                # Count the total number of chains
                total = total + chain_appears 

            # Check if it is a strict HM. This would be the case if all the chains that form the structure were matches
            if total == structure_dict[structure[0:4]][1]:
                strict_HM_2.append(structure)

        if len(strict_HM_2) > 0:
            strict_HM_2_final = ','.join(strict_HM_2)
            line.append(strict_HM_2_final)
        else:
            line.append('NA')
    
    writer.writerow(line)

handle_writer.close()


# ## Select the complexes with the best resolution for interface analyses



handle = open(output_path + '/paralogs_PDB_structures_with_strict2.txt', 'r')
reader = csv.reader(handle, delimiter = '\t')

handle_writer = open(output_path + '/paralogs_PDB_structures_best_structures2.txt', 'w')
writer = csv.writer(handle_writer, delimiter = '\t')

header = ['P1_ID', 'P2_ID', 'Duplication_type', 'HM_P1', 'HM_P2', 'HET', 'P1_unspecific', 'P2_unspecific', 'Strict_HM_P1', 'Strict_HM_P2']   
writer.writerow(header)




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
    
    # Look at the list of structures in each column, and select the crystal with the best resolution
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
                
    # Now that we have all of them, write the info to the table
    if bool_save:
        new_row = [P1_ID, P2_ID, dup_type, best_HM_P1[0], best_HM_P2[0], best_HET[0], P1_unspecific, P2_unspecific, best_strict_HM_P1[0], best_strict_HM_P2[0]]
        writer.writerow(new_row)
        
handle_writer.close()

