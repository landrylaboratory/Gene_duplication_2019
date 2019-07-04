# This file contains all the helper functions and classes for the 004_call_interfaces.py script

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

    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html.
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

    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z]

##################################

def replace_b_factor(pdb_infile, in_dict, outfile):
    '''This function uses an input PDB file and replaces b-factors with values from
    a dictionary.'''
    in_handle = open(pdb_infile, 'r')
    out_handle = open(outfile, 'w')

    for line in in_handle:
        if line.startswith('ATOM'):

            # Parse the line
            parsed_line = parse_pdb_line(line)

            chain = parsed_line[4]
            resid = parsed_line[5] + parsed_line[3]

            # Replace the b-factor
            dict_sasa = str(in_dict[chain][resid]) +'0'
            final_sasa = (6-len(dict_sasa))*' ' + dict_sasa

            final_line = line.replace(line[60:66], final_sasa)
            out_handle.write(final_line)

    # Close the outfile
    out_handle.close()

##################################

# Use a class for PDB coordinates
class PDB_coordinates:
    def __init__(self, x_coord, y_coord, z_coord, residue, atomtype):
        '''The constructor for the PDB coordinates class'''
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
        self.residue = residue
        self.atomtype = atomtype

    def measure_distance(self, target):
        '''A function that measures the distance between two PDB_coordinates objects.'''
        x_dist = (self.x - target.x)**2
        y_dist = (self.y - target.y)**2
        z_dist = (self.z - target.z)**2
        return math.sqrt(x_dist + y_dist + z_dist)

##################################

def parse_coordinates(infile, alpha_only):
    '''This function will parse a pdb file using the PDB_coordinates class.
    If the "alpha_only" argument is set to True the function only reads alpha carbons ("CA").
    If it is set to False it will read all atoms.
    '''
    handle = open(infile, 'r')
    coord_dict = OrderedDict()

    for line in handle:
        if line.startswith('ATOM'):

            atom_line = parse_pdb_line(line)

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
    '''This function parses the FreeSasa file to use as a reference for the van der Waals radii.
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
                # Include the terminal oxygen
                vdw_dict[aminoacid] = {'OXT' : 1.46}
            # Add the vdw radius for the current atomtype for this aminoacid if it has not been registered
            if vdw_dict[aminoacid].get(atomtype, -1) == -1:
                vdw_dict[aminoacid][atomtype] = vdw_radius


    return vdw_dict

##################################
