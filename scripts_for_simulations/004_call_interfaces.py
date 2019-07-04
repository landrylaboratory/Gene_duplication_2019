import re
import os
from collections import OrderedDict
import math
import argparse
from call_interfaces_helper import *
import csv

# This script will receive an input PDB file whose interfaces it will call based on the distance to the other chain (Tsai, 1996)

parser = argparse.ArgumentParser(description = 'This script will receive an input PDB file whose interfaces it will call based on the distance to the other chain (Tsai, 1996)')
parser.add_argument('-i', type = str, help = 'The input PDB file', dest = 'infile')
parser.add_argument('-o', type = str, help = 'Name of the folder to save the output files.', dest = 'out_folder')

args = parser.parse_args()
pdb_file = args.infile
out_folder = args.out_folder
outfile_name = pdb_file.split('/')[-1]

out_pdb_dist = os.path.join(out_folder, 'dist_regions_'  + outfile_name)
out_dict = os.path.join(out_folder, 'dist_regions_' + outfile_name + '.txt')

if not os.path.exists(out_folder):
	os.makedirs(out_folder)

# Reference vdW radii are based on the output of FreeSASA for PDB entry 1aof. 
ref_vdw_dict = 'Data/1aof_sasa.pdb'

# Create the reference dictionary for van der Waals radii
vdw_dict = get_vdw_dict(ref_vdw_dict)

alpha_carbons = parse_coordinates(pdb_file, True)
chains = alpha_carbons.keys()

# Prepare an output file for the pairs of directly contacting residues
list_contacts = os.path.join(out_folder, 'dist_contacting_residues.txt')
handle_contacts = open(list_contacts, 'w')
writer_contacts = csv.writer(handle_contacts, delimiter = '\t')
all_contacts = []

distance_region_dict = OrderedDict()
for chain in chains:
    distance_region_dict[chain] = OrderedDict()

all_atoms = parse_coordinates(pdb_file, False)

alpha_carbon_residue_list = []

for chainA_num in range(len(chains) - 1):
    chainA = chains[chainA_num]

    for chainB_num in range(chainA_num + 1, len(chains)):
        chainB = chains[chainB_num]
        
        # Compare the distances for the alpha carbons
        for alphaA in alpha_carbons[chainA]:
            residueA = alphaA.residue[1:]
        
            for alphaB in alpha_carbons[chainB]:

                residueB = alphaB.residue[1:]

                # Use this code to distinguish regions
                # 0.00 = neither nearby nor interacting
                # 0.25 = interaction candidates (used to reduce the number of all-atom measurements)
                # 0.75 = nearby
                # 1.00 = interacting
                if distance_region_dict[chainA].get(residueA, -1) != 0.25:
                    distance_region_dict[chainA][residueA] = 0.00

                if distance_region_dict[chainB].get(residueB, -1) != 0.25:
                    distance_region_dict[chainB][residueB] = 0.00

                # Get candidates based on their alpha carbons
                if alphaA.measure_distance(alphaB) <= 10:

                    distance_region_dict[chainA][residueA] = 0.25
                    distance_region_dict[chainB][residueB] = 0.25

# Save a dictionary of the residues that make up each interface
interacting_dict = OrderedDict()
for chainA_num in range(len(chains) - 1):
    chainA = chains[chainA_num]
    for chainB_num in range(chainA_num + 1, len(chains)):
        chainB = chains[chainB_num]
        interacting_dict[(chainA, chainB)] = [[], []]            

# Look at all the atom pairwise comparisons between interacting candidates
for chainA_num in range(len(chains) - 1):
    chainA = chains[chainA_num]
    for atomA in all_atoms[chainA]:
        residueA = atomA.residue[1:] 
        # Only keep going if this is a candidate
        if distance_region_dict[chainA].get(residueA, 0) >= 0.25:
            for chainB_num in range(chainA_num + 1, len(chains)):           
                chainB = chains[chainB_num]
                for atomB in all_atoms[chainB]:
                    residueB = atomB.residue[1:]
                    # Only keep going if this is a candidate
                    if distance_region_dict[chainB].get(residueB, 0) >= 0.25:

                        # Skip hydrogen atoms
                        if not atomA.atomtype.startswith('H') and not atomB.atomtype.startswith('H'):
                            vdw_atomA = vdw_dict[residueA[-3:]][atomA.atomtype]
                            vdw_atomB = vdw_dict[residueB[-3:]][atomB.atomtype]
                            atom_distance = atomA.measure_distance(atomB)

                            # Check if the atoms are interacting
                            if atom_distance < vdw_atomA + vdw_atomB + 0.5:
                                distance_region_dict[chainA][residueA] = 1.00
                                distance_region_dict[chainB][residueB] = 1.00

                                # Use a tuple to save the interacting pair
                                new_res_A = residueA[0:-3] + chainA + residueA[-3:]
                                new_res_B = residueB[0:-3] + chainB + residueB[-3:]
                                contact_pair = [new_res_A, new_res_B]
                                contact_pair.sort()
                                
                                # Write them to the list if they have not been written before
                                if not contact_pair in all_contacts:
                                    writer_contacts.writerow(contact_pair)
                                    all_contacts.append(contact_pair)

                                # Save to the dictionary
                                if not residueA in interacting_dict[(chainA, chainB)][0]:
                                    interacting_dict[(chainA, chainB)][0].append(residueA)
                                if not residueB in interacting_dict[(chainA, chainB)][1]:
                                    interacting_dict[(chainA, chainB)][1].append(residueB)

# Use the interacting residues to find the nearby residues
all_alpha = []
for chain in chains:
    all_alpha = all_alpha + alpha_carbons[chain]

# Keep track of comparisons to avoid repeating
for pos1 in range(0,len(all_alpha)-1):

    alpha1 = all_alpha[pos1]
    residue1 = alpha1.residue[1:]
    chain1 = alpha1.residue[0]

    for pos2 in range(pos1+1, len(all_alpha)):
        alpha2 = all_alpha[pos2]
        residue2 = alpha2.residue[1:]
        chain2 = alpha2.residue[0]

        if alpha1.measure_distance(alpha2) <= 6:

            # If residueB is near an interacting residue, classify it as nearby
            if distance_region_dict[chain1][residue1] == 1.00 and distance_region_dict[chain2][residue2] < 0.75:
                distance_region_dict[chain2][residue2] = 0.75
                # Add the nearby residue to the interface dict
                for chain_pair, interface_lists in interacting_dict.items():
                    if residue1 in interface_lists[0] and chain2 == chain_pair[0]:
                        interacting_dict[chain_pair][0].append(residue2)
                    if residue1 in interface_lists[1] and chain2 == chain_pair[1]:
                        interacting_dict[chain_pair][1].append(residue2)

            # If residueA is near an interacting residue, classify it as nearby
            if distance_region_dict[chain2][residue2] == 1.00 and distance_region_dict[chain1][residue1] < 0.75:
                distance_region_dict[chain1][residue1] = 0.75
                # Add the nearby residue to the interface dict
                for chain_pair, interface_lists in interacting_dict.items():
                    if residue2 in interface_lists[0] and chain1 == chain_pair[0]:
                        interacting_dict[chain_pair][0].append(residue1)
                    if residue2 in interface_lists[1] and chain1 == chain_pair[1]:
                        interacting_dict[chain_pair][1].append(residue1)

# Set interaction candidates that do not interact to non interfaces
for chain in chains:
    for alpha_carbon in alpha_carbons[chain]:
        residue = alpha_carbon.residue[1:]
        chain = alpha_carbon.residue[0]
        if distance_region_dict[chain][residue] == 0.25:
            distance_region_dict[chain][residue] = 0.00

# Write the new PDB file with the regions
replace_b_factor(pdb_file, distance_region_dict, out_pdb_dist)

# Close the file for the contacts
handle_contacts.close()

int_dict_file = open(out_dict, 'w')
writer = csv.writer(int_dict_file, delimiter = '\t')
for key, value in interacting_dict.items():
    col1 = ','.join(key)
    col2_1 = ','.join(value[0])
    col2_2 = ','.join(value[1])
    col2 = col2_1 + ';' + col2_2
    writer.writerow([col1, col2])
int_dict_file.close()  

