import argparse
import random
import sys
from collections import OrderedDict
import numpy

#### Some previously defined functions and variables ####

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

thorvaldsen_matrix_file = '/home/afcis2/FoldX_simulations/Subst_matrices/Thorvaldsen_genetic_code.txt'
zhu_matrix_file = '/home/afcis2/FoldX_simulations/Subst_matrices/clean_zhu_codon_aa.txt'

# thorvaldsen_matrix_file = '/Users/angelcisneros/Documents/Automne2017/Scripts/Python/Subst_matrices/Thorvaldsen_genetic_code.txt'
# zhu_matrix_file = '/Users/angelcisneros/Documents/Automne2018/Paralog_interference/Subs_matrices/clean_zhu_codon_aa.txt'

#########################################################

parser = argparse.ArgumentParser('This script will receive a PDB file for FoldX simulations and randomly choose amino acids from a specified region to be mutated to other random amino acids.')

parser.add_argument('-i', type = str, help = 'The path to the input PDB file', dest = 'infile')
parser.add_argument('-t', nargs='+', type = float, help = 'A series of float values separated by spaces that indicate the regions whose residues can be substituted (0.00 for the surface, 0.25 for the interior, 0.5 for the support, 0.75 for the rim, 1.00 for the core). Example: 0.75 1.00 will allow substitutions in the rim and the core', default = 0.00, dest = 'regions')
parser.add_argument('-c', type = str2bool, help = 'A boolean value that specifies whether the complex is a homodimer or not. If the complex is a homodimer, substitutions will be applied to both chains, default = True', dest = 'homodimer_check', default = True)
parser.add_argument('-p', type = int, help = 'An integer that specifies the program to be used (1: random substitutions, 2: exhaustive search of the one-substitution space)', dest = 'program', default = 1)
parser.add_argument('-m', type = int, help = 'An integer that specifies the substitution matrix to use (0 for a uniform probability of substitution, 1 for the Thorvaldsen matrix based on codons)', dest = 'matrix_type', default = 0)

args = parser.parse_args()

reference_file = args.infile
regions = args.regions
homodimer_check = args.homodimer_check
program = args.program
matrix_type = args.matrix_type

##### Tests for my matrix sampling #####

def sample_subs_matrix(subs_dict, curr_aa):

    if sum(list(subs_dict[curr_aa].values())) < 1:
        # If the probabilities don't add up to one (they are usually 0.999X), I add the option to resample
        candidate = numpy.random.choice(list(subs_dict[curr_aa].keys()) + ['X'], p = list(subs_dict[curr_aa].values()) + [1-sum(list(subs_dict[curr_aa].values()))])
    elif sum(list(subs_dict[curr_aa].values())) > 1:
        # If the probabilities add up to more than one (due to small rounding errors), I reduce them equally to make them become one
        for i in range(len(subs_dict[curr_aa].values())):
            subs_dict[curr_aa].values()[i] = subs_dict[curr_aa].values()[i] - (1-sum(list(subs_dict[curr_aa].values())))/20
        candidate = numpy.random.choice(list(subs_dict[curr_aa].keys()), p = list(subs_dict[curr_aa].values()))
    else:
        # If the probabilities add up to exactly one, I sample as is
        candidate = numpy.random.choice(list(subs_dict[curr_aa].keys()), p = list(subs_dict[curr_aa].values()))
    
    if candidate == 'X':
        return(sample_subs_matrix(subs_dict, curr_aa))
    else:
        return(candidate)

#########################

# Add a function that makes sure that the candidate is present on both chains

def residue_both_chains(selected_pos, all_res_pos):
    ''' This function will receive a selected residue for mutation and two dictionaries
        for the whole chain sequences. It will make sure that the selected position is
        present in both chains, keeping in mind that it may be a different amino acid
        in the case of the heteromeric simulation. It returns a boolean that will tell
        me if I can stop redrawing mutation candidates.
    '''
    stop_redrawing = True

    for chain, position_list in all_res_pos.items():
        # If it is not present on any of the two chains, then I must stop.
        if not selected_pos in position_list:
            stop_redrawing = False
            # print "Residue", selected_pos, "not present in chain", chain
            break

    return stop_redrawing

#########################

#########################

def parse_subs_matrix(subs_matrix):
    subs_dict = OrderedDict()

    with open(subs_matrix, 'r') as matrix:
        # Extract the first line
        header = matrix.readline().strip().split('\t')
        
        # I don't need the 'amino acid' tag
        # Add all the amino acids to the dictionary
        for aa in header[1:]:
            subs_dict[aa] = OrderedDict()
        
        # Read each line and load the substitution probabilities in the dictionary
        for line in matrix:
            line = line.strip().split('\t')
            curr_aa = line[0]
            
            for new_aa_pos in range(1,len(line)):
                new_aa = header[new_aa_pos]
                new_aa_prob = line[new_aa_pos]
                subs_dict[curr_aa][new_aa] = float(new_aa_prob)

    return subs_dict

#########################

def propose_mutation(resname, matrix_type, aa_dict):
    ''' This function proposes a mutation based on the selected transition matrix.
    '''
    if matrix_type == 0:
        # Uniform distribution
        aa_candidates = aa_dict.values()
        aa_candidates.remove(resname)
        new_aa = random.sample(aa_candidates, 1)[0]
    elif matrix_type == 1:
        # Thorvaldsen's substitution matrix
        subs_dict = parse_subs_matrix(thorvaldsen_matrix_file)
        new_aa = sample_subs_matrix(subs_dict, resname)
    elif matrix_type == 2:
        # The substitution matrix derived from Zhu's data and codon usage
        subs_dict = parse_subs_matrix(zhu_matrix_file)
        new_aa = sample_subs_matrix(subs_dict, resname)
    return new_aa

#########################

aa_list = [
    ('ALA','A'),
    ('ARG','R'),
    ('ASN', 'N'),
    ('ASP', 'D'),
    ('CYS','C'),
    ('GLU','E'),
    ('GLN','Q'),
    ('GLY','G'),
    ('HIS','H'),
    ('ILE','I'),
    ('LEU','L'),
    ('LYS','K'),
    ('MET','M'),
    ('PHE','F'),
    ('PRO','P'),
    ('SER','S'),
    ('THR','T'),
    ('TRP','W'),
    ('TYR','Y'),
    ('VAL','V')
]

aa_dict = OrderedDict(aa_list)

# Read the reference PDB file to check the amino acid positions I will be able to mutate
ref_handle = open(reference_file, 'r')

mutatable_res = {}
all_res_pos = OrderedDict()
chains = []
# Read through the file and return a list of all the positions that could be modified with the current regions
for line in ref_handle:
    # I am only interested in lines that clear the given regions and contain atoms
    if line.startswith('ATOM'):
        # I will save all the residue positions to the all_res_pos dictionary
        # I will only look at alpha carbons
        if line[12:16].strip() == 'CA':
            res_num = line[22:27].strip(' ')
            resname = aa_dict[line[17:20].strip(' ')]
            chain = line[21]
        
            # Register the chain if it is new
            if not chain in chains:
                chains.append(chain)
                all_res_pos[chain] = []
                mutatable_res[chain] = []
        
            # I only save the residue number here. This is all I need to make sure that a particular position is present on both chains
            all_res_pos[chain].append(res_num)
        
            if float(line[60:66].strip()) in regions:

                # Add the mutatable residue to the list for that chain
                mutatable_res[chain].append(resname + chain + res_num)

ref_handle.close() 

##### Check the chosen program (1 for random, 2 for exhaustive) #####

if program == 1:

    if homodimer_check:

        # Make sure that the selected position is in both chains
        stop_redrawing = False
        while not stop_redrawing:
            # Once I have the whole list of mutatable residues, I need to select one of them at random
            mutated_res = random.sample(mutatable_res[chains[0]] + mutatable_res[chains[1]], 1)[0]
            resname = mutated_res[0]
            respos = mutated_res[2:]
            # Check if I need to stop
            stop_redrawing = residue_both_chains(respos, all_res_pos)

        # Now that I am out of the loop, I can save the proposed mutations
        # Select a random value from the other amino acids for this substitution
        # This is where I can sample from either a uniform distribution or Thorvaldsen's substitution matrix
        new_aa = propose_mutation(resname, matrix_type, aa_dict)
	
        outfile = open('individual_list.txt', 'w')
        outfile.write(resname + chains[0] + respos +  new_aa + ',' + resname + chains[1] + respos + new_aa + ';' )
        outfile.close()

    else:
        
        # Make sure that the selected positions are in both chains
        stop_redrawing = False
        while not stop_redrawing:
            mutated_res1 = random.sample(mutatable_res[chains[0]], 1)[0]
            resname1 = mutated_res1[0]
            respos1 = mutated_res1[2:]
            stop_redrawing = residue_both_chains(respos1, all_res_pos)

        # Mutation one is now set. We must move on to mutation 2.
        stop_redrawing = False
        while not stop_redrawing:
            mutated_res2 = random.sample(mutatable_res[chains[1]], 1)[0]
            resname2 = mutated_res2[0]
            respos2 = mutated_res2[2:]
            stop_redrawing = residue_both_chains(respos2, all_res_pos)

        new_aa1 = propose_mutation(resname1, matrix_type, aa_dict)
        new_aa2 = propose_mutation(resname2, matrix_type, aa_dict)

        # Write the substitutions for both homodimers and the heterodimer
        outfile = open('homodimer_' + chains[0] + '/individual_list.txt', 'w')
        outfile.write(resname1 + chains[0] + respos1 + new_aa1 + ',' + resname1 + chains[1] + respos1 + new_aa1 + ';' )
        outfile.close()

        outfile = open('homodimer_' + chains[1] + '/individual_list.txt', 'w')
        outfile.write(resname2 + chains[0] + respos2 + new_aa2 + ',' + resname2 + chains[1] + respos2 + new_aa2 + ';' )
        outfile.close()

        outfile = open('heterodimer/individual_list.txt', 'w')
        outfile.write(mutated_res1 + new_aa1 + ',' + mutated_res2 + new_aa2 + ';')
        outfile.close()	

elif program == 2:
    # I can loop through the list of mutatable residues, select them in order and start substituting them
    # As some residues belong to the interface in both chains, for mutations that consider the complex a homodimer, I will keep track of the positions of the residues I have already substituted
    # This will avoid repeating substitutions in said homodimeric complexes.
    substituted_positions = {}
    for res in mutatable_res:
        # Extract its information
        resname = res[0]
        reschain = res[1]
        respos = res[2:]
        # Get the identities of the amino acids that could replace this residue
        # Uniform distribution
        if matrix_type == 0:
            aa_candidates = aa_dict.values()
            # All other amino acids
            aa_candidates.remove(resname)
        elif matrix_type == 1:
            # Restricted to the ones that are possible with Thorvaldsen's matrix
            aa_candidates = []
            for new_aa, prob in subs_dict[resname].items():
                if prob > 0:
                    aa_candidates.append(new_aa)

        if substituted_positions.get(respos, -1) == -1:
            # Write to the standard output the list of all the possible substitutions
            for subs in aa_candidates:
                # Options for homodimers and heterodimers
                if homodimer_check:
                    sys.stdout.write(resname + chains[0] + respos +  subs + ',' + resname + chains[1] + respos + subs + ';' )
                    substituted_positions[respos] = 1
                else:
                    sys.stdout.write(res + subs + ';')



