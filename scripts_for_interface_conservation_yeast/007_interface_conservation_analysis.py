#!/usr/bin/env python
# coding: utf-8

# # Interface conservation analyses
# 
# This is the final script that will do the analyses of the sequence identity of paralogs at the interface. It will do the following steps:
# 
# - 1.- Go through the list of best structures as to find the PDB structures I need (strict HMs and HETs)
# - 2.- Get the needed biological assembly for each of them. 
# - 3.- Use such files to create dictionaries that remember which chains represent my paralogs of interest.
# - 4.- Identify the interfaces that include that chain. 
# - 5.- Retrieve the sequence of each paralog pair and use the phylomeDB phylogenies to get alignments to calculate sequence identity at the interface.


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


# Define paths to input folders
input_PDB = <path> # Folder with all the PDB files
phylomeDB_path = <path> # Path to the alignments from PhylomeDB
path_reference_proteome = <path> # Path to the reference proteome

# Define paths to output_folders
output_path = <path> # A parent folder for all the outputs
output_004_bio_assemblies = <path> # Biological assemblies for each input PDB
output_005_chain_match_tables = <path> # Chain match tables explaining which chains in the biological assembly come from which chain in the asymmetric unit
output_006_dist_interfaces = <path> # Files with the calculated interfaces
output_007_interface_dict = <path> # Files specifying which residues from each chain belong in interfaces
output_008_phylo_PDB_seq = <path> # FASTA files containing PhylomeDB alignments, the paralogs, and the PDB sequence
output_009_phylo_PDB_aln = <path> # FASTA files containing the aligned sequences from the output_008_phylo_seq folder
output_010_only_PDB_positions = <path> # FASTA files containing the parts of the alignment that correspond to the region observed in the PDB structure
output_011_only_interfaces = <path> # FASTA files containing the parts of the alignment that correspond to the interfaces observed in the PDB structure
output_012_only_non_interfaces = <path> # FASTA files containing the parts of the alignment that correspond to the non-interfaces observed in the PDB structure



def add2dict(dictionary, key, value):
    '''This is a small function I can use to add values to dictionaries. If the key is not
    already in the dictionary, then it initializes its value as a list with one element (value).
    Otherwise, it appends to that list.'''
    
    if dictionary.get(key, [-1]) == [-1]:
        dictionary[key] = value
    else:
        curr_val = list(dictionary[key])
        dictionary[key] = curr_val.append(value)
    
    return dictionary


# ## 1.- Map each of the paralog pairs to their corresponding structures


handle_in = open('Data/paralogs_PDB_structures_best_structures.txt', 'r')
reader = csv.reader(handle_in, delimiter = '\t')

paralogs2structures = OrderedDict()

# Skip header
header = reader.next()

for line in reader:
    P1 = line[0]
    P2 = line[1]
    HET = line[5]
    P1_HM = line[8]
    P2_HM = line[9]
    
    # Save the pair of paralogs as a key in the dictionary
    paralogs2structures = add2dict(paralogs2structures, (P1, P2), OrderedDict())
    
    # Add each of the complexes only if it is not NA
    if P1_HM != 'NA':
        paralogs2structures[(P1, P2)] = add2dict(paralogs2structures[(P1, P2)], 'P1_HM', P1_HM)
    if P2_HM != 'NA':
        paralogs2structures[(P1, P2)] = add2dict(paralogs2structures[(P1, P2)], 'P2_HM', P2_HM)
    if HET != 'NA':
        paralogs2structures[(P1, P2)] = add2dict(paralogs2structures[(P1, P2)], 'HET', HET)

handle_in.close()


# ## 2.- Get the needed biological assemblies

# Use a loop to call the script that generates the biological assemblies
for paralog_pair, struct_dict in paralogs2structures.items():
    for key, structure in struct_dict.items():
        assembly = structure[5:]
        pdb = structure[0:4]
        call_script = 'python ../scripts_for_simulations/001_generate_bio_assembly.py '
        arg1 = '-i ' + input_PDB + '/' + pdb + '.pdb '
        arg2 = '-o ' + output_004_bio_assemblies + '/' + structure + '.pdb '
        arg3 = '-n ' + assembly + ' '
        arg4 = '-t ' + output_005_chain_match_tables + '/' + pdb + '_chain_table.txt'
        command = call_script + arg1 + arg2 + arg3 + arg4
        os.system(command)


# ## 3.- Use the matching tables and the previous alignments to know which chains correspond to my proteins of interest

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
    
    # When adding a new PDB structure
    elif alignment_dict[pair][P1].get(P1_match_PDB, -1) == -1:
        alignment_dict[pair][P1][P1_match_PDB] = [P1_match_chain]
        
    # When adding chains to a PDB structure
    else:
        alignment_dict[pair][P1][P1_match_PDB].append(P1_match_chain)
        

# Load all the data
file_list = glob.glob(output_005_chain_match_tables + '/*')

chain_dict = OrderedDict()
for match_table in file_list:
    structure = match_table.split('/')[-1][0:4]
    handle = open(match_table, 'r')
    reader = csv.reader(handle, delimiter = '\t')
    chain_dict[structure] = OrderedDict()
    for line in reader:
        key = line[0]
        values = line[1].split(',')
        chain_dict[structure][key] = values


# ## 4.- Call the interfaces in each biological assembly

def parse_pdb_line(pdb_line):
    '''This function will receive a line from a PDB file and parse it as a list. It will do so based on the
    PDB format explanation from this site:

    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    '''
    atom = pdb_line[0:4].strip(' ')
    atom_num = pdb_line[6:11].strip(' ')
    atom_name = pdb_line[12:16].strip(' ')
    resname = pdb_line[17:20].strip(' ')
    chain = pdb_line[21]
    res_num = pdb_line[22:26].strip(' ')
    x = pdb_line[30:38].strip(' ')
    y = pdb_line[38:46].strip(' ')
    z = pdb_line[46:54].strip(' ')
    region = float(pdb_line[61:66].strip())

    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z, region]

##################################

def replace_b_factor(pdb_infile, in_dict, outfile):
    '''This function uses an input PDB file, selects the lines that correspond to the atoms, and replaces their
    b-factor values with the relative SASA from the dictionary. The output is written directly to a the outfile.'''
    in_handle = open(pdb_infile, 'r')
    out_handle = open(outfile, 'w')

    for line in in_handle:
        if line.startswith('ATOM'):

            parsed_line = parse_pdb_line(line)

            # Ignore DNA sequences
            if parsed_line[3] in ['DA', 'DG', 'DC', 'DT', 'A', 'C', 'G', 'T', 'U']:
                continue
        
            chain = parsed_line[4]
            resid = parsed_line[5] + parsed_line[3]

            dict_sasa = str(in_dict[chain].get(resid, 0.0)) +'0'
            final_sasa = (6-len(dict_sasa))*' ' + dict_sasa

            final_line = line.replace(line[60:66], final_sasa)
            out_handle.write(final_line)

    # Once the loop has finished, I can close the outfile
    out_handle.close()

##################################

# Define a class for PDB coordinates
class PDB_coordinates:
    def __init__(self, x_coord, y_coord, z_coord, residue, atomtype):
        '''The constructor for the PDB coordinates class'''
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
        self.residue = residue
        self.atomtype = atomtype

    def measure_distance(self, target):
        '''A function I can use to measure the distance between two PDB_coordinates objects.'''
        x_dist = (self.x - target.x)**2
        y_dist = (self.y - target.y)**2
        z_dist = (self.z - target.z)**2
        return math.sqrt(x_dist + y_dist + z_dist)

##################################

def parse_coordinates(infile, alpha_only):
    '''This function will parse a pdb file using the PDB_coordinates class.
    If the "alpha only" argument is set to True the function will only look at 
    alpha carbons. It will look at all atoms if it is set to False.
    '''
    handle = open(infile, 'r')
    coord_dict = OrderedDict()

    for line in handle:
        if line.startswith('ATOM'):

            atom_line = parse_pdb_line(line)

            # Ignore DNA sequences
            if atom_line[3] in ['DA', 'DG', 'DC', 'DT', 'A', 'C', 'G', 'T', 'U']:
                continue
            
            if alpha_only:

                if atom_line[2] == 'CA':
                    if coord_dict.get(atom_line[4], -1) == -1:
                        coord_dict[atom_line[4]] = []
                    # Add the atom to the list
                    x_coord = float(atom_line[6])
                    y_coord = float(atom_line[7])
                    z_coord = float(atom_line[8])
                    residue = atom_line[4] + atom_line[5] + atom_line[3]
                    coord_dict[atom_line[4]].append(PDB_coordinates(x_coord, y_coord, z_coord, residue, 'CA'))
            else:
                if coord_dict.get(atom_line[4], -1) == -1:
                        coord_dict[atom_line[4]] = []
                # Add the atom to the list
                x_coord = float(atom_line[6])
                y_coord = float(atom_line[7])
                z_coord = float(atom_line[8])
                atomtype = atom_line[2]
                residue = atom_line[4] + atom_line[5] + atom_line[3]
                coord_dict[atom_line[4]].append(PDB_coordinates(x_coord, y_coord, z_coord, residue, atomtype))

    return coord_dict

##################################

def get_vdw_dict(infile):
    '''This function parses the FreeSasa file I will use as a reference for the van der Waals radii.
    '''
    handle = open(infile, 'r')

    vdw_dict = {}

    for line in handle:
        if line.startswith('ATOM'):
            line_list = re.split('\s+', line)

            aminoacid = line_list[3]
            atomtype = line_list[2]
            vdw_radius = float(line_list[9])
            # Add the amino acid if it has not been registered
            if vdw_dict.get(aminoacid, -1) == -1:
                vdw_dict[aminoacid] = {'OXT' : 1.46}
            # Add the vdw radius for the current atomtype for this aminoacid if it has not been registered
            if vdw_dict[aminoacid].get(atomtype, -1) == -1:
                vdw_dict[aminoacid][atomtype] = vdw_radius

    # Now that I finished, I can return the dictionary
    return vdw_dict

ref_vdw_dict = '../scripts_for_simulations/Data/1aof_sasa.pdb'
vdw_dict = get_vdw_dict(ref_vdw_dict)


# Define the main function that will find interfaces
def distance_interfaces(pdb_file, out_pdb_dist):
    alpha_carbons = parse_coordinates(pdb_file, True)
    chains = alpha_carbons.keys()

    distance_region_dict = OrderedDict()
    for chain in chains:
        distance_region_dict[chain] = OrderedDict()

    all_atoms = parse_coordinates(pdb_file, False)

    alpha_carbon_residue_list = []
    
    for chainA_num in range(len(chains) - 1):
        chainA = chains[chainA_num]

        for chainB_num in range(chainA_num + 1, len(chains)):
            chainB = chains[chainB_num]
            
            # I will compare the distances for the alpha carbons
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

                            vdw_atomA = vdw_dict[residueA[-3:]][atomA.atomtype]
                            vdw_atomB = vdw_dict[residueB[-3:]][atomB.atomtype]
                            atom_distance = atomA.measure_distance(atomB)

                            # Check if the atoms are interacting
                            if atom_distance < vdw_atomA + vdw_atomB + 0.5:
                                # I can list them as interacting residues
                                distance_region_dict[chainA][residueA] = 1.00
                                distance_region_dict[chainB][residueB] = 1.00

                                # I can save them to the dictionary
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
    
    return(interacting_dict)   


# # The main block that calculates interfaces

# Get the list of files
file_list = glob.glob(output_004_bio_assemblies + '/*')

# Call interfaces
for pdb_file in file_list:
    print pdb_file
    structure = pdb_file.split('/')[-1][0:4]
    out_file = output_006_dist_interfaces + '/dist_regions_' + structure + '.pdb'
    interacting_dict = distance_interfaces(pdb_file, out_file)
    
    # Save the interacting dict to a file
    int_dict_file = open(output_007_interface_dict + '/' + structure + '_interface_dict.txt', 'w')
    writer = csv.writer(int_dict_file, delimiter = '\t')
    for key, value in interacting_dict.items():
        col1 = ','.join(key)
        col2_1 = ','.join(value[0])
        col2_2 = ','.join(value[1])
        col2 = col2_1 + ';' + col2_2
        writer.writerow([col1, col2])
    int_dict_file.close()


# ## 5.- Retrieve the interfaces and start looking at their conservation

# For all pairs of paralogs, the table that helps me match the pairs
integrated_table_file = 'Data/paralogs_PDB_structures_best_structures.txt'

P1_id = 'YDL022W'

phylomedb_path = phylomeDB_path
desired_regions = [0.75, 1.00]


# Load a dictionary that helps me map the protein IDs to PhylomeDB IDs
handle = open("Data/yeast_SN_genes_names.txt", 'r')
reader = csv.reader(handle, delimiter = '\t')

systematic2phylomeDB = OrderedDict()
phylomeDB2systematic = OrderedDict()
for line in reader:
    systematic2phylomeDB[line[1]] = line[0]
    phylomeDB2systematic[line[0]] = line[1]


# Define some functions
def align_seqs(in_fasta):
    '''This function will take the FASTA file with the PDB entry's sequences and
    align them. It will return a dictionary with the aligned sequences.
    '''
    muscle_cline = MuscleCommandline(input = in_fasta)

    stdout, stderr = muscle_cline()
    
    align_dict = {}
    for line in stdout.split('\n'):
        if line.startswith('>'):
            
            current_chain = line.split(' ')[0]
            
            if current_chain[1:] in align_dict.keys():
                current_chain = current_chain + '2'
            align_dict[current_chain[1:]] = ''
        else:
            align_dict[current_chain[1:]] = align_dict[current_chain[1:]] + line
        
    return align_dict

#########################################

def align2records(alignment):
    '''This function will receive a muscle alignment and change it into a list of records that can be saved.
    '''
    records = []
    for prot_id, sequence in alignment.items():
        new_record = SeqRecord(seq = Seq(sequence), id = prot_id, description = '')
        records.append(new_record)
    return records

#########################################

def percentage_identity(in_file):
    '''This function will receive a fasta file with two sequences and calculate their sequence identity.
    '''
    records = []
    for record in SeqIO.parse(in_file, 'fasta'):
        records.append(record)
    
    total_residues = len(record.seq)
    if total_residues == 0:
        return 0
    
    matches = 0.0
    dual_gaps = 0.0
    seq1 = records[0]
    seq2 = records[1]
    for pos in range(total_residues):
        if seq1[pos] == seq2[pos]:
            if seq1[pos] == '-':
                dual_gaps += 1.0
            else:
                matches = matches + 1.0
        
    return round(matches*100/(total_residues - dual_gaps), 2)

#########################################

def percentage_identity2(align_dict):
    '''This function will receive a dictionary with two sequences and calculate their sequence identity.
    '''
    seq1 = align_dict.values()[0]
    seq2 = align_dict.values()[1]
    
    total_residues = len(seq1)
    matches = 0.0
    dual_gaps = 0.0
    for pos in range(total_residues):
        if seq1[pos] == seq2[pos]:
            if seq1[pos] == '-':
                # This would mean that both sequences have gaps in this position. This happens when I extract
                # particular sequences from an MSA.
                dual_gaps += 1.0
            else:
                # Then they are the same and we have a match.
                matches = matches + 1.0
        
    return round(matches*100/(total_residues - dual_gaps), 2)

#########################################

def region_parser(interface_file, aa_dict):
    '''This function takes a PDB file in which regions have been called in order to identify the interfaces.
    It returns a dictionary in which residues map to a two-element list:
    - First position: True if they belong to the desired regions and False otherwise
    - Second position: The ID of the region to which they belong
    '''
    
    interface_file_handle = open(interface_file, 'r')

    # Create a dictionary that will store the chains, its residues, and the regions to which they belong
    pdb_dict = OrderedDict()
    start_dict = OrderedDict()
    end_dict = OrderedDict()
    
    sequence_dict = OrderedDict()

    # Save the previous position number to make sure I deal with alternative coordinates
    prev_pos = 'X'
    
    # Parse the file into the dictionary
    for line in interface_file_handle:
        pdb_line = parse_pdb_line(line)

        # Add the data for this residue to the dictionary
        res_id = pdb_line[3] + pdb_line[5]
        region = pdb_line[9]
        chain = pdb_line[4]
        atom_type = pdb_line[2]
        
        # Keep track of the start and the end positions
        curr_pos = int(pdb_line[5])

        # Look at alpha carbons to avoid repeating
        if atom_type == 'CA' and prev_pos != curr_pos:

            prev_pos = curr_pos
            
            # Keep track of the start and end positions
            if not chain in start_dict.keys():
                start_dict[chain] = curr_pos
            elif start_dict[chain] > curr_pos:
                start_dict[chain] = curr_pos

            if not chain in end_dict.keys():
                end_dict[chain] = curr_pos
            elif end_dict[chain] < curr_pos:
                end_dict[chain] = curr_pos

            # Save the regions to which each residue belongs
            if pdb_dict.get(chain, -1) == -1:
                pdb_dict[chain] = OrderedDict()
                pdb_dict[chain][res_id] = region
            else:
                pdb_dict[chain][res_id] = region

            # Save sequences
            if sequence_dict.get(chain, -1) == -1:
                sequence_dict[chain] = aa_dict[res_id[0:3]]
            else:
                sequence_dict[chain] += aa_dict[res_id[0:3]]

            
    return pdb_dict, start_dict, end_dict, sequence_dict

#########################################

def extract_PDB_region(alignment_dict, P1_id, output_path):
    ''' This function will receive the alignment parsed as a dictionary and return the region of the alignment
    that corresponds to residues in the PDB sequence. It will write the resulting sequences to a file in
    the 010_only_PDB_positions folder in the output path.
    '''
    # Using the MSA, extract the part of the sequences that match to the PDB.
    start = -1
    end = -1
    PDB_seq = alignment_dict[P1_id + "_pdb"].seq
    for pos in range(0,len(PDB_seq)):
        if PDB_seq[pos] != '-':
            # Assign this position to the end because it contains an amino acid
            end = pos
            # If this is the first one, save it as the start
            if start == -1:
                start = pos

    ### Knowing the start and end positions, save these sequences in the 007 folder.
    only_PDB_positions = []
    for key, record in alignment_dict.items():
        # Only save the records that are not the PDB sequence
        if key != P1_id + "_pdb":
            new_record = SeqRecord(record.seq[start:end], id = key, description = '')
            only_PDB_positions.append(new_record)

    SeqIO.write(only_PDB_positions, os.path.join(output_path, P1_id + "_only_PDB_pos.fasta"), "fasta")
    return(PDB_seq)

#########################################

def fasta2list(in_fasta):
    ''' This function receives a multi fasta file and saves its records to a list.
    '''
    phylo_records = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        phylo_records.append(record)
    return phylo_records

########################################

def fasta2dict(in_fasta):
    '''This function receives a multi fasta file and parses it into a dictionary.
    '''
    dict_out = OrderedDict()
    for record in SeqIO.parse(in_fasta, "fasta"):
        dict_out[record.id] = record
    return dict_out


# Add MSE (selenomethionine) as M and MLY (methylysine) as K
aa_dict = {
    'ALA':'A',
    'ARG':'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS':'C',
    'GLU':'E',
    'GLN':'Q',
    'GLY':'G',
    'HIS':'H',
    'ILE':'I',
    'LEU':'L',
    'LYS':'K',
    'MET':'M',
    'PHE':'F',
    'PRO':'P',
    'SER':'S',
    'THR':'T',
    'TRP':'W',
    'TYR':'Y',
    'VAL':'V',
    'MSE':'M',
    'MLY':'K'
}

aa2three_letter = OrderedDict()
for key, value in aa_dict.items():
    if key in ('MSE', 'MLY'):
        continue
    else:
        aa2three_letter[value] = key

desired_regions = [0.75, 1.00]


# Load the yeast reference proteome so that I can extract sequences easily
# Needed whenever the sequences are missing from the PhylomeDB alignments
proteome_dict = OrderedDict()
pdb_seqs = SeqIO.parse(path_reference_proteome, 'fasta')
for record in pdb_seqs:
    proteome_dict[record.id] = record



# I will define the functions to work with the PDB sequences here
def match_PDB_full_seq(align_dict, PDB_dict, full_seq_id, pdb_seq_id, aa2three_letter):
    '''This function will help me find assign PDB regions for the full sequence.
    '''
    # Retrieve the aligned sequences
    aligned_chain_pdb = align_dict[pdb_seq_id]
    aligned_chain_full = align_dict[full_seq_id]
    
    counter_pdb = 0
    counter_full = 0
    
    # Loop through the keys of the dictionary
    pdb_residues = PDB_dict.keys()
    
    match_list = []
    full_region_dict = OrderedDict()
    
    # Loop through the positions of the PDB chain
    for position in range(len(aligned_chain_pdb)):
        # Get the residue in this position 
        aligned_pdb = aligned_chain_pdb[position]
        
        if counter_pdb < len(pdb_residues):
            pdb_pos = pdb_residues[counter_pdb]
        else:
            pdb_pos = '---'
        
        # Get the residue in this position for the full sequence
        aligned_full = aligned_chain_full[position]
    
        # Look at the four possible scenarios with respect to gaps in the sequences
        if aligned_pdb == '-' and aligned_full == '-':
            continue
        elif aligned_pdb == '-' and aligned_full != '-':
            # This is a gap in the PDB sequence.
            full_pos = aa2three_letter[aligned_full] + str(counter_full)
            full_region_dict[full_pos] = -1
            match_list.append(['-', full_pos, -1])
            counter_full = counter_full + 1
        elif aligned_pdb != '-' and aligned_full == '-':
            # This is a gap in the full sequence.
            counter_pdb = counter_pdb + 1
        elif aligned_pdb != '-' and aligned_full != '-':
            # Make sure that the alignment position residue type matches the query
            if aa2three_letter.get(aligned_pdb, -1) != pdb_pos[0:3]:
                # This will identify unusual residues
                print 'Found something else', pdb_pos, aligned_pdb, structure
            else:
                region = PDB_dict.get(pdb_pos, '-')
                # Save a dictionary with the matches
                full_pos = aa2three_letter[aligned_full] + str(counter_full)
                match_list.append([pdb_pos, full_pos, region])

                # Add one to the PDB counter and to the full counter
                counter_pdb = counter_pdb + 1
                counter_full = counter_full + 1

                full_region_dict[full_pos] = region

        #####################################
    
    # Return the dictionary that maps each residue in the sequence to the available structural information
    return full_region_dict, match_list
    


# ## Use the PhylomeDB phylogenies to get sequence identities at the interface

def first_round_alignments2(phylo_dict, phylo_ID, paralog_ID, alignment_dict, chain_dict):
    '''This function takes the dictionaries and does the first alignment, that is, the sequence from the PDB to
    the full protein, which will help identify the position of the interfaces in that full protein.
    '''
    needed_chains = alignment_dict[(P1, P2)][paralog_ID][structure[0:4]]
    for needed_chain in needed_chains:
        new_chain = chain_dict[structure[0:4]][needed_chain][0]
        if new_chain != '':
            break
    pdb_sequence = pdb_sequence_dict[new_chain]
    
    # The pdb_dict already has all residues that belong to the interface of at least one of the subunits
    PDB_sequence_record = SeqRecord(Seq(pdb_sequence), id = paralog_ID + '_pdb', description = '')

    # Make records with these two sequences and align them
    full_seq_record = phylo_dict[phylo_ID]
    records = [full_seq_record, PDB_sequence_record]

    pdb_seq_id = PDB_sequence_record.id
    full_seq_id = full_seq_record.id

    # Write the PDB sequence and the full sequence for that protein
    outpath = output_008_phylo_PDB_seq + '/' + phylo_ID + '.fasta'
    SeqIO.write(records, outpath, 'fasta')

    # Align with MUSCLE
    aln_dict = align_seqs(outpath)
    
    # Use the alignment dict to know where the PDB chain matches the full sequence 
    full_region_dict, match_list = match_PDB_full_seq(aln_dict, pdb_dict[new_chain], full_seq_id, pdb_seq_id, aa2three_letter)
    
    return full_region_dict, match_list
    
    ####################################


def save_final_sequences2(phylo_dict, full_region_dict, phylo_ID_2, phylo_ID, aa2three_letter):
    '''This function transfers the annotation to the second paralog and writes the final files.'''

    full_region_dict_2, match_list_2 = match_PDB_full_seq(phylo_dict, full_region_dict, phylo_ID_2, phylo_ID, aa2three_letter)
    
    pdb_seqs = ['','']
    only_interfaces = ['','']
    only_non_interfaces = ['','']
    for position in match_list_2: 
        if position[2] >= 0.00:  
            # Add to the pdb sequences
            pdb_seqs[0] = pdb_seqs[0] + aa_dict[position[0][0:3]]
            pdb_seqs[1] = pdb_seqs[1] + aa_dict[position[1][0:3]]

            # Check region to decide where to add it
            if position[2] >= 0.75:
                # Add to interfaces
                only_interfaces[0] = only_interfaces[0] + aa_dict[position[0][0:3]]
                only_interfaces[1] = only_interfaces[1] + aa_dict[position[1][0:3]]
            else:
                only_non_interfaces[0] = only_non_interfaces[0] + aa_dict[position[0][0:3]]
                only_non_interfaces[1] = only_non_interfaces[1] + aa_dict[position[1][0:3]]

    # Save these sequences as records and write them
    # PDB seqs
    record_1 = SeqRecord(seq = Seq(pdb_seqs[0]), id = paralog_ID)
    record_2 = SeqRecord(seq = Seq(pdb_seqs[1]), id = paralog_ID_2)
    records = [record_1, record_2]
    out_pdb = output_010_only_PDB_positions + '/' + P1 + '_' + P2 + '_only_pdb_positions.fasta'
    SeqIO.write(records, out_pdb, 'fasta')

    # Only interfaces
    record_1 = SeqRecord(seq = Seq(only_interfaces[0]), id = paralog_ID)
    record_2 = SeqRecord(seq = Seq(only_interfaces[1]), id = paralog_ID_2)
    records = [record_1, record_2]
    out_interfaces = output_011_only_interfaces + '/' + P1 + '_' + P2 + '_only_interfaces.fasta'
    SeqIO.write(records, out_interfaces, 'fasta')

    # Only non-interfaces
    record_1 = SeqRecord(seq = Seq(only_non_interfaces[0]), id = paralog_ID)
    record_2 = SeqRecord(seq = Seq(only_non_interfaces[1]), id = paralog_ID_2)
    records = [record_1, record_2]
    out_non_interfaces = output_012_only_non_interfaces + '/' + P1 + '_' + P2 + '_only_non_interfaces.fasta'
    SeqIO.write(records, out_non_interfaces, 'fasta')
    
    return [out_pdb, out_interfaces, out_non_interfaces]


handle = open(integrated_table_file, 'r')
reader = csv.reader(handle, delimiter = '\t')

new_table = open(output_path + "/Paralogs_sequence_identities_dist.txt", 'w')
new_table_writer = csv.writer(new_table, delimiter = '\t')

headers = reader.next()
headers = headers + ["Compared_interface", "Full_sequence_identity", "PDB_sequence_identity", "Interface_sequence_identity", "Non_interface_sequence_identity"]
new_table_writer.writerow(headers)

distance_interface_path = output_006_dist_interfaces
interface_dict_path = output_007_interface_dict

# Loop through the lines
for line in reader:
    P1 = line[0]
    P2 = line[1]
    
    # Get the PDB files I need for that pair using the structure dictionary
    for complex_type, structure in paralogs2structures[(P1, P2)].items():
        
        # Get the PDB file that corresponds to this complex
        pdb_file = distance_interface_path + 'dist_regions_' + structure[0:4] + '.pdb' 

        print pdb_file
        print P1, P2
        
        # Parse the interface dictionary if it is a HET
        if complex_type == 'HET':
            # Retrieve regions from the PDB file
            pdb_dict, start_dict, end_dict, pdb_sequence_dict = region_parser(pdb_file, aa_dict)
            
            # Run first with the first paralog
            paralog_ID = P1
            paralog_ID_2 = P2
            
            # Get the phylo IDs to look at the PhylomeDB files
            phylo_ID = systematic2phylomeDB[paralog_ID]
            phylo_ID_2 = systematic2phylomeDB[paralog_ID_2]
            
            phylome_dict = fasta2dict(os.path.join(phylomedb_path, phylo_ID + '.raw.fasta'))
            
            full_region_dict, match_list = first_round_alignments2(phylome_dict, phylo_ID, paralog_ID, alignment_dict, chain_dict)
            
            record_list = fasta2list(os.path.join(phylomedb_path, phylo_ID + '.raw.fasta'))
                                      
            # Need to check if the paralog is in the phylome
            record_2 = phylome_dict.get(phylo_ID_2, -1)
            if type(record_2) != SeqRecord:
                record_2 = proteome_dict[paralog_ID_2]
            
                record_list.append(record_2)
            
                # Write the records and align
                outpath = output_009_phylo_PDB_aln + '/' + P1 + '_' + P2 + '.fasta'
                SeqIO.write(record_list, outpath, 'fasta')
                pair_align_dict = align_seqs(outpath)
                
                final_pair_dict = OrderedDict()
                
                for key, value in pair_align_dict.items():
                    if key in (phylo_ID, paralog_ID_2):
                        final_pair_dict[key] = value
                 
                out_files = save_final_sequences2(final_pair_dict, full_region_dict, paralog_ID_2, phylo_ID, aa2three_letter)
                
                same_phylome = 0
                
            else:
                final_pair_dict = OrderedDict()
            
                # Extract from that alignment the records for the paralogs 
                for key, value in phylome_dict.items():
                    if key in (phylo_ID, phylo_ID_2):
                        final_pair_dict[key] = value
            
                out_files = save_final_sequences2(final_pair_dict, full_region_dict, phylo_ID_2, phylo_ID, aa2three_letter)
            
                same_phylome = 1
            
            # Finally, read sequence similarity
            full_seq_ident = percentage_identity2(final_pair_dict)
            PDB_seq_ident = percentage_identity(in_file = out_files[0])
            interface_seq_ident = percentage_identity(in_file = out_files[1])
            non_interface_seq_ident = percentage_identity(in_file = out_files[2])

            # Check percentage identity
            print 'Full sequence identity:', full_seq_ident
            print "Whole PDB structure:", PDB_seq_ident
            print "Only interfaces:", interface_seq_ident
            print "Only non-interfaces:", non_interface_seq_ident
            print '---------'
            
            new_table_writer.writerow(line + ["HET_P1", full_seq_ident, PDB_seq_ident, interface_seq_ident, non_interface_seq_ident, same_phylome])
            
            # Run with the second one
            paralog_ID = P2
            paralog_ID_2 = P1
            
            # I should get the phylo IDs to look at the PhylomeDB files
            phylo_ID = systematic2phylomeDB[paralog_ID]
            phylo_ID_2 = systematic2phylomeDB[paralog_ID_2]
            
            phylome_dict = fasta2dict(os.path.join(phylomedb_path, phylo_ID + '.raw.fasta'))
            
            full_region_dict, match_list = first_round_alignments2(phylome_dict, phylo_ID, paralog_ID, alignment_dict, chain_dict)
            
            record_list = fasta2list(os.path.join(phylomedb_path, phylo_ID + '.raw.fasta'))
                                      
            # Need to check if the paralog is in the phylome
            record_2 = phylome_dict.get(phylo_ID_2, -1)
            if type(record_2) != SeqRecord:
                record_2 = proteome_dict[paralog_ID_2]
            
                record_list.append(record_2)
            
                # Write the records and align
                outpath = output_009_phylo_PDB_aln + '/' + P1 + '_' + P2 + '.fasta'
                SeqIO.write(record_list, outpath, 'fasta')
                pair_align_dict = align_seqs(outpath)
                
                final_pair_dict = OrderedDict()
                
                for key, value in pair_align_dict.items():
                    if key in (phylo_ID, paralog_ID_2):
                        final_pair_dict[key] = value
            
                out_files = save_final_sequences2(final_pair_dict, full_region_dict, paralog_ID_2, phylo_ID, aa2three_letter)
                
                same_phylome = 0
                
            else:
                final_pair_dict = OrderedDict()
            
                # Extract from that alignment the records for the two paralogs
                for key, value in phylome_dict.items():
                    if key in (phylo_ID, phylo_ID_2):
                        final_pair_dict[key] = value
            
                out_files = save_final_sequences2(final_pair_dict, full_region_dict, phylo_ID_2, phylo_ID, aa2three_letter)
            
                same_phylome = 1
            
            # Finally, read sequence similarity
            full_seq_ident = percentage_identity2(final_pair_dict)
            PDB_seq_ident = percentage_identity(in_file = out_files[0])
            interface_seq_ident = percentage_identity(in_file = out_files[1])
            non_interface_seq_ident = percentage_identity(in_file = out_files[2])

            # Check percentage identity
            print 'Full sequence identity:', full_seq_ident
            print "Whole PDB structure:", PDB_seq_ident
            print "Only interfaces:", interface_seq_ident
            print "Only non-interfaces:", non_interface_seq_ident
            print '---------'
            new_table_writer.writerow(line + ["HET_P2", full_seq_ident, PDB_seq_ident, interface_seq_ident, non_interface_seq_ident, same_phylome])
            
        else:
            # Retrieve regions from the PDB file
            pdb_dict, start_dict, end_dict, pdb_sequence_dict = region_parser(pdb_file, aa_dict)
            
            # As these are HM, get the sequence from any of them
            pdb_sequence = pdb_sequence_dict.values()[0]
            
            if complex_type == 'P1_HM':  
                paralog_ID = P1
                paralog_ID_2 = P2
            elif complex_type == 'P2_HM':
                paralog_ID = P2
                paralog_ID_2 = P1

            # Get the phylo IDs to look at the PhylomeDB files
            phylo_ID = systematic2phylomeDB[paralog_ID]
            phylo_ID_2 = systematic2phylomeDB[paralog_ID_2]
            
            phylome_dict = fasta2dict(os.path.join(phylomedb_path, phylo_ID + '.raw.fasta'))
            
            full_region_dict, match_list = first_round_alignments2(phylome_dict, phylo_ID, paralog_ID, alignment_dict, chain_dict)
            
            record_list = fasta2list(os.path.join(phylomedb_path, phylo_ID + '.raw.fasta'))
                                      
            # Need to check if the paralog is in the phylome
            record_2 = phylome_dict.get(phylo_ID_2, -1)
            if type(record_2) != SeqRecord:
                record_2 = proteome_dict[paralog_ID_2]
            
                record_list.append(record_2)
            
                # Write the records and align
                outpath = output_009_phylo_PDB_aln + '/' + P1 + '_' + P2 + '.fasta'
                SeqIO.write(record_list, outpath, 'fasta')
                pair_align_dict = align_seqs(outpath)
                
                final_pair_dict = OrderedDict()
                
                for key, value in pair_align_dict.items():
                    if key in (phylo_ID, paralog_ID_2):
                        final_pair_dict[key] = value
                        print key
                
                out_files = save_final_sequences2(final_pair_dict, full_region_dict, paralog_ID_2, phylo_ID, aa2three_letter)
                
                same_phylome = 0
                
            else:
                final_pair_dict = OrderedDict()
            
                # Extract from that alignment the records for the two paralogs
                for key, value in phylome_dict.items():
                    if key in (phylo_ID, phylo_ID_2):
                        final_pair_dict[key] = value
                        print key
            
                out_files = save_final_sequences2(final_pair_dict, full_region_dict, phylo_ID_2, phylo_ID, aa2three_letter)
                
                same_phylome = 1
                
            # Finally, read sequence similarity
            full_seq_ident = percentage_identity2(final_pair_dict)
            PDB_seq_ident = percentage_identity(in_file = out_files[0])
            interface_seq_ident = percentage_identity(in_file = out_files[1])
            non_interface_seq_ident = percentage_identity(in_file = out_files[2])
            
            # Check percentage identity
            print 'Full sequence identity:', full_seq_ident
            print "Whole PDB structure:", PDB_seq_ident
            print "Only interfaces:", interface_seq_ident
            print "Only non-interfaces:", non_interface_seq_ident
            print '---------'
            new_table_writer.writerow(line + [complex_type, full_seq_ident, PDB_seq_ident, interface_seq_ident, non_interface_seq_ident, same_phylome])
            
handle.close()
new_table.close()


# Define a function that gets the length of the matching chain in the PDB structure
def get_PDB_length(path_to_pdb_file, chain_query):
    '''
    A function designed to count the number of residues in a particular chain in a PDB structure.
    '''
    structure = parser.get_structure(path_to_pdb_file, path_to_pdb_file)
    
    res_no = 0
    
    for model in structure:
        for chain in model:
            if chain.get_id() == chain_query:
                for residue in chain.get_residues():
                    if residue.id[0] == ' ':
                        res_no +=1
                
    return(res_no)


parser = PDBParser()

handle = open(os.path.join(output_path, "Paralogs_sequence_identities_dist_2019-06-19.txt"), 'r')
reader = csv.reader(handle, delimiter = '\t')

out_handle = open("Data/Paralogs_sequence_identities_dist_completeness_2019-06-19.txt"), 'w')
writer = csv.writer(out_handle, delimiter = '\t')

for entry in reader:
    P1 = entry[0]
    P2 = entry[1]
    
    # Skip header
    if P1 == 'P1_ID':
        header = entry
        header.append('PDB_length')
        header.append('Full_protein_length')
        writer.writerow(header)
        continue
       
    compared_struc = entry[10]
    
    # Now that I have the dictionary, I can use it to select the best isoforms
    if compared_struc == 'P1_HM' or compared_struc == 'HET_P1':
    
        # Add the length of the protein chain from the PDB structure
        if compared_struc == 'P1_HM':
            pdb_structure = entry[8]
        elif 'HET_P1':
            pdb_structure = entry[5]
            
        matching_chain = alignment_dict[(P1, P2)][P1][pdb_structure[0:4]][0]
        
        protein_path = output_004_bio_assemblies + pdb_structure + '.pdb'
        
        pdb_length = get_PDB_length(protein_path, matching_chain)
        
        # Add the length of the full protein from the reference proteome
        full_protein_length = len(proteome_dict[P1].seq)
    
    elif compared_struc == 'P2_HM' or compared_struc == 'HET_P2':         
        # Add the length of the protein chain from the PDB structure
        if compared_struc == 'P2_HM':
            pdb_structure = entry[9]
        elif 'HET_P2':
            pdb_structure = entry[5]
            
        matching_chain = alignment_dict[(P1, P2)][P2][pdb_structure[0:4]][0]
        
        protein_path = output_004_bio_assemblies + pdb_structure + '.pdb'
        
        pdb_length = get_PDB_length(protein_path, matching_chain)
        
        # Add the length of the full protein from the reference proteome
        full_protein_length = len(proteome_dict[P2].seq)
    
    new_line = entry
    new_line.append(pdb_length)
    new_line.append(full_protein_length)
    writer.writerow(new_line)
                
            
handle.close()
out_handle.close()

