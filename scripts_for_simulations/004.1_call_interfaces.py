import re
import os
from collections import OrderedDict
import math
import argparse
from interface_identifier_new_h import *
import csv

parser = argparse.ArgumentParser(description = 'This script will receive an input PDB file whose interfaces it will call based on two methods: Levy, 2010 (using solvent accessible surfaces) and Keskin, 2004 (using distances)')
parser.add_argument('-i', type = str, help = 'The input PDB file', dest = 'infile')
parser.add_argument('-o', type = str, help = 'Name of the folder to save the output files.', dest = 'out_folder')

args = parser.parse_args()
pdb_file = args.infile
out_folder = args.out_folder
outfile_name = pdb_file.split('/')[-1]

if not os.path.exists(out_folder):
	os.makedirs(out_folder)

# This dictionary is based on the output of FreeSASA for PDB entry 1aof. 
ref_vdw_dict = 'Data/1aof_sasa.pdb'

# Create the reference dictionary for van der Waals radii
vdw_dict = get_vdw_dict(ref_vdw_dict)

######################## Start finding interfaces with the distance definition #########################

alpha_carbons = parse_coordinates(pdb_file, True)
chains = alpha_carbons.keys()
chainA = chains[0]
chainB = chains[1]

distance_region_dict = {chainA : OrderedDict(),
           chainB : OrderedDict()}

all_atoms = parse_coordinates(pdb_file, False)

# Before measuring distances, I will make sure that all the residues have an alpha carbon, which I will
# take as a proxy of chain integrity

# I will compare the distances for the alpha carbons
for alphaA in alpha_carbons[chainA]:

    residueA = alphaA.residue[1:]

    for alphaB in alpha_carbons[chainB]:

        residueB = alphaB.residue[1:]

        # Check if alphaA has already been added to the dict
        # As I want to compare the "nearby" residues to the rim and the "interacting" residues
        # to the core, I will use this color code:
        # 0.00 = not nearby or interacting
        # 0.25 = interaction candidates (a placeholder to accelerate the search for interacting residues)
        # 0.75 = nearby
        # 1.00 = interacting
        if distance_region_dict[chainA].get(residueA, -1) != 0.25:
            distance_region_dict[chainA][residueA] = 0.00

        if distance_region_dict[chainB].get(residueB, -1) != 0.25:
            distance_region_dict[chainB][residueB] = 0.00

        # I decided to increase this candidate threshold to 10 because my new implementation was not
        # getting a pair of interacting residues
        if alphaA.measure_distance(alphaB) <= 10:
            # Then I add alphaB to the list of residues in the second chain that
            # are near alphaA in the first chain
            distance_region_dict[chainA][residueA] = 0.25
            distance_region_dict[chainB][residueB] = 0.25

# I will start looking at all the atom pairwise comparisons between interacting candidates
for atomA in all_atoms[chains[0]]:

    residueA = atomA.residue[1:]                        

    for atomB in all_atoms[chains[1]]:

        residueB = atomB.residue[1:]

        # To know if I must make the comparison, I need to check if the residues are nearby
        if distance_region_dict[chainA][residueA] == 0.25 and distance_region_dict[chainB][residueB] == 0.25:
            # I will need to check the distance to see if they are interacting or not
            # Remember that my convention for the PDB_coordinates objects' residue is
            # <position> + <residue_type>
            # Therefore, if I want to know the residue type, i need the last three letters
            vdw_atomA = vdw_dict[residueA[-3:]][atomA.atomtype]
            vdw_atomB = vdw_dict[residueB[-3:]][atomB.atomtype]
            atom_distance = atomA.measure_distance(atomB)

            # Now, I can check if the atoms are close enough to be considered as interacting
            if atom_distance < vdw_atomA + vdw_atomB + 0.5:
                # I can list them as interacting residues
                distance_region_dict[chainA][residueA] = 1.00
                distance_region_dict[chainB][residueB] = 1.00

# Now that I know which residues are interacting, I should go back and classify the nearby residues
# Nearby residues are those whose alpha carbons are within 6 A of the alpha carbon of an interacting atom,
# which could be the case for residues in the same chain 
# To do so, I will get a list of all alpha carbons
all_alpha = alpha_carbons[chainA] + alpha_carbons[chainB]

# I will define the cycles so that I only perform the lower triangle matrix of comparisons, that is,
# I will compare the first residue to all other residues but I will not return to it after that.
for pos1 in range(0,len(all_alpha)-1):

    alpha1 = all_alpha[pos1]
    residue1 = alpha1.residue[1:]
    chain1 = alpha1.residue[0]

    for pos2 in range(pos1+1, len(all_alpha)):
        alpha2 = all_alpha[pos2]
        residue2 = alpha2.residue[1:]
        chain2 = alpha2.residue[0]

        if alpha1.measure_distance(alpha2) <= 6:

            # If residueA is an interacting residue and residueB is not an interacting residue, I will
            # classify residueB as a nearby residue
            if distance_region_dict[chain1][residue1] == 1.00 and distance_region_dict[chain2][residue2] < 0.75:
                distance_region_dict[chain2][residue2] = 0.75

            # The same logic but the other way around
            if distance_region_dict[chain2][residue2] == 1.00 and distance_region_dict[chain1][residue1] < 0.75:
                distance_region_dict[chain1][residue1] = 0.75            

# I also need to remove the placeholder for interacting candidates (0.25)
for alpha_carbon in all_alpha:
    residue = alpha_carbon.residue[1:]
    chain = alpha_carbon.residue[0]
    if distance_region_dict[chain][residue] == 0.25:
        distance_region_dict[chain][residue] = 0.00

# Use the distance_region_dict to write the new PDB file with the regions
out_pdb_dist = os.path.join(out_folder, 'dist_regions_' + outfile_name)
replace_b_factor(pdb_file, distance_region_dict, out_pdb_dist)

####################### Finished finding interfaces with the distance definition #######################

################## Start finding interfaces with the surface accessibility definition ##################

# Get the sasa for each residue in the complex file and write a PDB with the relative surface areas as b-factors
complex_rsa_file = os.path.join(out_folder, outfile_name + '.rsa')
get_sasa(pdb_file, complex_rsa_file)

complex_rsa = parse_rsa(complex_rsa_file)
complex_rsa_pdb = os.path.join(out_folder, 'complex_rsa_' + outfile_name)
replace_b_factor(pdb_file, complex_rsa, complex_rsa_pdb)

## Isolate each chain and compute their residues sasa when found as monomers
change_dict = OrderedDict()
classifier_dict = OrderedDict()

for chain in complex_rsa.keys():
    # The following lines get the SASA for each of the individual chains and prepare the PDB file for visualization
    # chain_file will contain the coordinates of the chain in PDB format
    # chain_sasa_rsa will contain the chain's SASA in the rsa format
    # chain_sasa_pdb will be a PDB file with the coordinates and the SASA
    chain_file =  os.path.join(out_folder, 'chain_' + chain + '_' + outfile_name)
    os.system("grep -e 'ATOM' " + pdb_file + " | grep '^.\{21\}'" + chain + ' > ' + chain_file)
	
    chain_sasa_rsa = os.path.join(out_folder, 'chain_' + chain + '_' + outfile_name + '.rsa')
    get_sasa(chain_file, chain_sasa_rsa)
    parsed_chain = parse_rsa(chain_sasa_rsa)

    chain_sasa_pdb = os.path.join(out_folder, 'sasa_chain_' + chain + '_' + outfile_name)
    replace_b_factor(chain_file, parsed_chain, chain_sasa_pdb)
    
    # Now, I want to know the change in relative SASA    
    # Initialize a subdictionary for the change in relative SASA in this chain
    change_dict[chain] = OrderedDict()
    classifier_dict[chain] = OrderedDict()

    # Loop through the aminoacids and calculate the change
    for aminoacid in parsed_chain[chain].keys():
        # Just remember that, in this code, rsa is the dictionary that contains the SASA for the 
        # complex and parsed_chain is the dictionary that contains the SASA for the unbound chain
        change_dict[chain][aminoacid] = abs(parsed_chain[chain][aminoacid] - complex_rsa[chain][aminoacid])
    
        # Now, I need to identify the region of the protein in which this residue belongs (surface, interior,
        # interface rim, interface support, interface core)
        # I will use a color code to distinguish the 5 regions Levy proposed
        # The idea will be to use a "b-factor" of:
        # 0.0 for the surface (it is exposed and does not change when in the complex)
        if complex_rsa[chain][aminoacid] > 0.25 and change_dict[chain][aminoacid] == 0:
            classifier_dict[chain][aminoacid] = 0.00
            
        # 25.0 for the interior (it is not exposed and does not change when in the complex)
        elif complex_rsa[chain][aminoacid] <= 0.25 and change_dict[chain][aminoacid] == 0:
            classifier_dict[chain][aminoacid] = 0.25
        
        # 50.0 for the support (it is not exposed in the monomeric form but changes when in the complex)
        elif change_dict[chain][aminoacid] > 0 and parsed_chain[chain][aminoacid] <= 0.25:
            classifier_dict[chain][aminoacid] = 0.50
        
        # 75.0 for the rim (it is exposed in the complex but changes from the monomeric form)
        elif change_dict[chain][aminoacid] > 0 and complex_rsa[chain][aminoacid] > 0.25:
            classifier_dict[chain][aminoacid] = 0.75
        
        # 1.00 for the core (it changes from exposed in the monomer to not exposed in the complex)
        elif change_dict[chain][aminoacid] > 0 and parsed_chain[chain][aminoacid] > 0.25 and complex_rsa[chain][aminoacid] <= 0.25:
            classifier_dict[chain][aminoacid] = 1.00

# Write a new PDB file with the sasa_change
sasa_change_pdb = os.path.join(out_folder, 'sasa_change_' + outfile_name)
replace_b_factor(pdb_file, change_dict, sasa_change_pdb)

# Write a new PDB file with the final regions
region_file = os.path.join(out_folder, 'sasa_regions_' + outfile_name)
replace_b_factor(pdb_file, classifier_dict, region_file)

################ Finished finding interfaces with the surface accessibility definition #################
