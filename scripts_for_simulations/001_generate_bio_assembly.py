import numpy as np
import re
from collections import OrderedDict
import os
import argparse
import csv

# This script will receive a PDB file and make the necessary calculations to obtain its biological assembly.

parser = argparse.ArgumentParser(description = 'This script will receive a PDB file and make the necessary calculations to obtain its biological assembly.')

parser.add_argument('-i', type = str, help = 'The path to the input PDB file', dest = 'infile')
parser.add_argument('-o', type = str, help = 'The path to the output PDB file', dest = 'outfile')
parser.add_argument('-n', type = str, help = 'The number of the biological assembly to retrieve', dest = 'bio_assembly')
parser.add_argument('-t', type = str, help = 'The output file where the info about the chain names will be saved', dest='out_table')

args = parser.parse_args()
pdb_file = args.infile
outfile = args.outfile
bio_assembly = args.bio_assembly
out_table = args.out_table

def get_coords_trans_matrix(pdb_file):
    '''This function will receive a pdb_file and return a dictionary with whose keys are the chains' IDs
    and whose values are numpy arrays with the coordinates for each atom and a dictionary with the transition
    matrices.
    '''
    # Get a handle for the file and loop through its lines, looking for the atoms and the transformation matrix
    handle = open(pdb_file, 'r')
    
    # As discussed in http://thread.gmane.org/gmane.comp.python.bio.general/8308 the transformation matrices
    # include rotation terms and translation terms
    
    coord_dict = OrderedDict()
    trans_matrix_dict = OrderedDict()

    modres = []
    all_chains = []
    seqres = OrderedDict()
    bio_bool = False

    for pdb_line in handle:
    
        # Check if there are modified residues
        if pdb_line.startswith('MODRES'):
            line_mod_res = pdb_line[12:15]
            if not line_mod_res in modres:
                modres.append(line_mod_res)
                
        # Save the sequence of the full protein chain.
        if pdb_line.startswith('SEQRES'):
            # Save this information in the seqres dictionary, with the chain ID as the key and a list
            # of all the SEQRES lines for that chain as the value 
            chain_ID = pdb_line[11]
            if seqres.get(chain_ID, -1) == -1:
                seqres[chain_ID] = []
            seqres[chain_ID].append(pdb_line)
            
        if pdb_line.startswith('ATOM') or (pdb_line.startswith('HETATM') and pdb_line[17:20] in modres):
            # Parse this as an atom line
            chain = pdb_line[21]
            # Add this chain to the all_chains list if it is not already there
            if not chain in all_chains:
                all_chains.append(chain)

            x = float(pdb_line[30:38].strip(' '))
            y = float(pdb_line[38:46].strip(' '))
            z = float(pdb_line[46:54].strip(' '))
    
            if coord_dict.get(chain, -1) == -1:
                coord_dict[chain] = [[x, y, z]]
            else:
                coord_dict[chain].append([x, y, z])
                    
        elif pdb_line.startswith('REMARK 350'):
            split_line = re.split('\s+', pdb_line)
            # Check if we have found the biological assembly we are looking for.
            if re.search("BIOMOLECULE:\s+" + bio_assembly, pdb_line):
                bio_bool = True
                
                chains = []
                
            # "BIOMOLECULE" indicates another biological assembly. Stop taking data.
            elif "BIOMOLECULE" in pdb_line:
                bio_bool = False 
                
            if bio_bool:
                if split_line[2] == 'APPLY' or "AND CHAINS:" in pdb_line:
                    
                    # Save the list of chains that will be transformed
                    chains = chains + re.search(': (.*)', pdb_line).group(1).strip().strip(',').split(', ')
                    
                elif split_line[2][0:-1] == 'BIOMT':
                    
                    for chain in chains:
                        if coord_dict.get(chain, -1) == -1:
                            coord_dict[chain] = []
                        # Save transformation matrices in a dictionary
                        if trans_matrix_dict.get(chain, -1) == -1:
                            trans_matrix_dict[chain] = OrderedDict()
                    
                    # Save the rows of the transformation matrix
                    for chain in chains:
                        if trans_matrix_dict[chain].get(split_line[3], -1) == -1:
                            trans_matrix_dict[chain][split_line[3]] = [[float(value) for value in split_line[4:8]]]
                        else:
                            trans_matrix_dict[chain][split_line[3]].append([float(value) for value in split_line[4:8]])

    # Convert the matrices into numpy objects
    for key, value in coord_dict.items():
        coord_dict[key] = np.array(value)

    # If there are no transformation matrices, trans_matrix will still be an empty OrderedDict.
    # In such a case, define it as an identity matrix for rotation and a vector of zeros for translation
    if trans_matrix_dict == OrderedDict():
        for chain in all_chains:
            trans_matrix_dict[chain] = OrderedDict()
            trans_matrix_dict[chain][1] = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]

    for key, dict_value in trans_matrix_dict.items():
        for key2, matrix in dict_value.items():
            trans_matrix_dict[key][key2] = np.array(matrix)
    
    return coord_dict, trans_matrix_dict, modres, seqres

##########################

def check_trans_matrices(trans_matrix_dict):
    '''This function will check if there are transition matrices that need to be kept or if the PDB file
    only contained identity matrices. It will return only the matrices that have a rotation that is not equal
    to 1 or whose translation vector contains something different from zeros.
    '''
    outdict = OrderedDict()
    # Assign new chain names to have consistent naming
    new_key = 'A'
    for key1, subdict in trans_matrix_dict.items():
        outdict[key1] = OrderedDict()
        for key2, matrix in subdict.items():
            outdict[key1][new_key] = matrix
            new_key = chr(ord(new_key) + 1)
    return outdict

##########################

def generate_new_chains(template_dict, trans_mat_dict, seqres):
    '''This function will take a dictionary with chains and a dictionary with the transformation matrices in
    order to get the new chains. It will also update the SEQRES dictionary.
    '''
    outdict = OrderedDict()
    new_seqres = OrderedDict()
    for chain, coords in template_dict.items():
        # Know which new chains will be daughters of the old chains, i.e.,
        # which chains are generated from which chains
        outdict[chain] = OrderedDict()

        # Check if the current chain is required in the biological assembly
        if chain in trans_mat_dict.keys():
    
            for new_chain, matrix in trans_mat_dict[chain].items():
                outdict[chain][new_chain] = np.dot(coords, matrix[0:3, 0:3].transpose()) + matrix[0:3, 3]
                outdict[chain][new_chain] = np.round(outdict[chain][new_chain], 3)
                # Update the seqres dictionary.
                new_seqres[new_chain] = []
                for i in range(0, len(seqres[chain])):
                    new_seqres[new_chain].append(seqres[chain][i][0:11] + new_chain + seqres[chain][i][12:])

    return outdict, new_seqres

##########################            

def numpy_list_generator(np_array):
    '''This function will receive a numpy array and return a generator of its rows as lists of strings.
    '''
    for row in np_array:
        row = row.tolist()
        for i in range(0, len(row)):
            row[i] = format(row[i], '.3f')
            row[i] = str(row[i])
        yield row

##########################

def write_biological_assembly(pdb_file, new_coords_dict, outfile, modres, seqres):
    '''This function receives a PDB file and writes the coordinates of the biological assembly. 
    new_coords is the dictionary that says which new chains are generated from which previous chains.
    '''
    handle = open(pdb_file, 'r')
    out_handle = open(outfile, 'w')

    # Write the SEQRES header for the chains in the structure
    for chain, seqres_lines in seqres.items():
        out_handle.writelines(seqres_lines)

    chain_dict = OrderedDict()
    generator_dict = OrderedDict()
    
    # Add the chain names we already know to the chain_dict
    for chain, new_chains in new_coords_dict.items():

        # Loop through the new chains and their coordinates
        for new_chain, new_coords in new_chains.items():
            chain_dict[new_chain] = []
            generator_dict[new_chain] = numpy_list_generator(new_coords)
    
    for line in handle:
        if line.startswith('ATOM') or (line.startswith('HETATM') and line[17:20] in modres):
            # Save the line's information to the dictionary
            chain = line[21]

            # Use the information for this line to write the information for its daughters
            for new_chain in new_coords_dict[chain].keys():
                # Get the new coordinates for that new_chain using its respective generator
                new_coords = generator_dict[new_chain].next()

                # Replace coordinates
                new_x = ' '*(len(line[30:38])-len(new_coords[0])) + new_coords[0]
                new_y = ' '*(len(line[38:46])-len(new_coords[1])) + new_coords[1]
                new_z = ' '*(len(line[46:54])-len(new_coords[2])) + new_coords[2]
                new_line = line[0:21] + new_chain + line[22:30] + new_x + new_y + new_z + line[54:]
                
                # Save the new coordinates to chain_dict
                chain_dict[new_chain].append(new_line)
    
    atom_counter = -1   
    # Write the new chains
    for chain, coord_list in chain_dict.items():
        # Format atom numbers so that they are continuous
        if atom_counter == -1:
            atom_counter = int(coord_list[-1][6:11].strip())
        else:
            for line_pos in range(0, len(coord_list)):
                atom_counter = atom_counter + 1
                atom_line = coord_list[line_pos]
                coord_list[line_pos] = atom_line.replace(atom_line[6:12],
                                                         ' '*(len(atom_line[6:11])-len(str(atom_counter))) + str(atom_counter) + ' ')
                
        
        out_handle.writelines(coord_list)
        out_handle.write('TER\n')
        
    out_handle.close()

##########################

# Get coordinates and transformation matrices
coords, trans_mat, modres, seqres = get_coords_trans_matrix(pdb_file)

# Organize the transformation matrices
trans_mat = check_trans_matrices(trans_mat)

# Use the transformation matrices to generate the new chains
new_coords_dict, seqres = generate_new_chains(coords, trans_mat, seqres)

# Save a dictionary to map chain IDs from the asymmetric unit to the ones in the
# biological assembly
out_handle = open(out_table, 'w')
writer = csv.writer(out_handle, delimiter = '\t')

for old_chain, new_chain_dict in new_coords_dict.items():
    new_chains = ','.join(new_chain_dict.keys())
    writer.writerow([old_chain, new_chains])

out_handle.close()

# Save the biological assembly
write_biological_assembly(pdb_file, new_coords_dict, outfile, modres, seqres)
