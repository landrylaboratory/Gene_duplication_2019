# This file contains all the helper functions and classes for the interface_identifier.py script

import re
import os
from collections import OrderedDict
import math
from Bio.PDB import *
from Bio import SeqIO
from Bio.Seq import *
from Bio.SeqRecord import *
import csv

##################################

def parse_pdb_line(pdb_line):
    '''This function will receive a line from a PDB file and parse it as a list. It will do so based on the
    PDB format explanation from this site:

    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

    I used to parse lines with re.split using whitespaces to separate the elements in each line. However, I ran
    into cases that don't have spaces, so I am probably better off separating them based on the characters' positions.
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

    # This list would have the following elements in these positions (with Python numbering):
    # 'ATOM' in position 0
    # atom number in position 1
    # atom type in position 2
    # residue name in position 3
    # chain in position 4
    # residue's number in the chain in position 5
    # coordinates in positions 6 to 8
    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z]

##################################

def parse_rsa_line(rsa_line):
    '''This function will do the same as the parse_pdb_line function but with FreeSasa's RSA
    files, as those also had problems with some residue numbers.
    '''
    res = rsa_line[0:3]
    resname = rsa_line[4:7]
    chain = rsa_line[8]
    res_num = rsa_line[9:13].strip(' ')
    abs_sasa = rsa_line[17:23].strip(' ')
    rel_sasa = rsa_line[24:29].strip(' ')

    # Position 1 is the residue's three-letter name
    # Position 2 is the chain ID
    # Position 3 is the residue's position in the chain
    # Position 5 is the relative SASA with respect to the amino acid of interest between two Ala residues
    return [res, resname, chain, res_num, abs_sasa, rel_sasa]

##################################

# Parse the FreeSasa RSA output
def parse_rsa(infile):
    '''This function receives a FreeSasa RSA output file and returns a dictionary that organizes the chains
    and their residues, such that one can recover the relative SASA for any particular aminoacid.'''
    handle = open(infile, 'r')

    out_dict = OrderedDict()

    # Get this file's chains' names
    for line in handle:
        if 'Chains' in line:
            dummy = re.search('\w+$', line.strip()).group(0)
            chains = []
            for chain in dummy:
                out_dict[chain] = OrderedDict()
            break

    # I don't need to check the top of the file again as the chains' information is always
    # before the residues' information
    for line in handle:
        # I will have found the residues' information once I find lines that start with "RES"
        if re.search('^RES', line):
            # Split the line with spaces.
            # Position 0 is always RES
            # Position 1 is the residue's three-letter name
            # Position 2 is the chain ID
            # Position 3 is the residue's position in the chain
            # Position 5 is the relative SASA with respect to the amino acid of interest between two Ala residues
            # parsed_line = re.split('\s+', line.strip())
            parsed_line = parse_rsa_line(line)

            chain = parsed_line[2]
            resid = parsed_line[3] + parsed_line[1]
            rel_sasa = float(parsed_line[5])

            out_dict[chain][resid] = rel_sasa

    return out_dict

##################################

def replace_b_factor(pdb_infile, in_dict, outfile):
    '''This function uses an input PDB file, selects the lines that correspond to the atoms, and replaces their
    b-factor values with the relative SASA from the dictionary. The output is written directly to a the outfile.'''
    in_handle = open(pdb_infile, 'r')
    out_handle = open(outfile, 'w')

    for line in in_handle:
        if line.startswith('ATOM'):

            # I will need to parse the line so that I can get the residue's chain, name, and position,
            # which I need to recover its relative SASA value from my dictionary
            # parsed_line = re.split('\s+', line.strip())
            parsed_line = parse_pdb_line(line)

            # This list would have the following elements in these positions (with Python numbering):
            # 'ATOM' in position 0
            # atom number in position 1
            # atom type in position 2
            # residue name in position 3
            # chain in position 4
            # residue's number in the chain in position 5
            chain = parsed_line[4]
            resid = parsed_line[5] + parsed_line[3]

            # Characters in positions 61-66 (60-66 with Python numbering) contain the b-value
            # These two lines allow me to recover the previously calculated SASA values and to format them
            # such that they respect the standard PDB format. The second line is particularly useful because
            # it adds any necessary spaces.
            dict_sasa = str(in_dict[chain][resid]) +'0'
            final_sasa = (6-len(dict_sasa))*' ' + dict_sasa

            final_line = line.replace(line[60:66], final_sasa)
            out_handle.write(final_line)

    # Once the loop has finished, I can close the outfile
    out_handle.close()

##################################

def get_sasa(infile, outfile):
    '''This function is just a wrapper for me to get the SASA of a structure with FreeSasa.
    It will just receive an input file and an output file.'''
    os.system('freesasa -n 100 --format=rsa ' + infile + ' > ' + outfile)

##################################

# Catherine thought I could use classes, so I'll try to do that
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
    '''This function will parse a pdb file using the PDB_coordinates class I defined.
    The "alpha_only" argument will help me decide if I only want alpha carbons ("CA")
    if it is True or all atoms if it is False.
    '''
    handle = open(infile, 'r')
    coord_dict = OrderedDict()

    for line in handle:
        if line.startswith('ATOM'):

            # atom_line = re.split('\s+', line.strip())
            atom_line = parse_pdb_line(line)

            if alpha_only:

                if atom_line[2] == 'CA':
                    # Then this line contains an alpha carbon and we are interested in preserving it
                    # If this is the first time we see that chain, initialize the list of alpha carbons
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
                # I will include the terminal oxygen right away. Its vdw radius does not seem to change
                # with regard to the aminoacid
                vdw_dict[aminoacid] = {'OXT' : 1.46}
            # Add the vdw radius for the current atomtype for this aminoacid if it has not been registered
            if vdw_dict[aminoacid].get(atomtype, -1) == -1:
                vdw_dict[aminoacid][atomtype] = vdw_radius
            # If the current vdw radius is different from the one I had previously registered for this
            # atomtype in this aminoacid, there should be something wrong, so I'll stop this and print
            elif vdw_dict[aminoacid][atomtype] != vdw_radius:
                print 'Conflicting vdw radius for ' + atomtype + ' in ' + aminoacid

    # Now that I finished, I can return the dictionary
    return vdw_dict

##################################
