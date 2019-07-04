import argparse
import sys
import re

# This code will receive an input PDB file and look for both insertion residues 
# (specified with a residue number and a letter, e.g. PHE184A) and alternative coordinates for residues 
# (specified by a residue name with a letter designation, e.g. ASER and BSER). 
# For insertion residues, it will renumber accordingly as to integrate them in 
# the numbering and make sure there are no repeated numbers. For alternative
# coordinates, it will keep the one with the highest occupancy because it is the one that would be found the most in the crystal.
# It will print messages on insertion residues and alternative coordinates to the standard output.

parser = argparse.ArgumentParser('This code will receive an input PDB file and look for both insertion residues (specified with a residue number and a letter, e.g. PHE184A) and alternative coordinates for residues (specified by a residue name with a letter designation, e.g. ASER and BSER). For insertion residues, it will renumber accordingly as to integrate them in the numbering and make sure there are no repeated numbers. For alternative coordinates, it will keep the one with the highest occupancy because it is the one that would be found the most in the crystal. It will print messages on insertion residues and alternative coordinates to the standard output.')
parser.add_argument('-i', help = 'Input PDB', dest = 'infile', type = str)
parser.add_argument('-o', help = 'Output PDB', dest = 'outfile', type = str)

args = parser.parse_args()
infile = args.infile
outfile = args.outfile

def check_alternative_coords(current_line, prev_line):
    '''This function will look for atoms that have more than one coordinate. It will keep the one
    with the highest occupancy since it would be the one that is better supported by the model. If the
    occupancies are the same, then it will keep the first orientation.
    '''
    # Save the occupancy values
    curr_occ = float(current_line[54:60].strip())
    prev_occ = float(prev_line[54:60].strip())
    
    atom_name = current_line[12:16].strip()
    resname = current_line[17:20].strip()
    resid = current_line[22:26].strip()
    
    sys.stdout.write('Alternative coordinates found for ' + atom_name + ' in ' + resname + resid + '.\n')

    # Compare the two occupancy values
    if curr_occ == prev_occ:
        sys.stdout.write('The occupancies are the same for ' + atom_name + ' in ' + resname + resid + '.\n')
        return prev_line[0:16] + ' ' + prev_line[17:]
    elif curr_occ > prev_occ:
        if curr_occ < 0.50:
            sys.stdout.write('The best occupancy for ' + atom_name + ' in ' + resname + resid + ' is smaller than 50.\n')
        return current_line[0:16] + ' ' + current_line[17:]
    else:
        if prev_occ < 0.50:
            sys.stdout.write('The best occupancy for ' + atom_name + ' in ' + resname + resid + ' is smaller than 50.\n')
        return prev_line[0:16] + ' ' + prev_line[17:]
    
#############

def shift_line(line, shift):
    '''This is a helper function for the renumber_insertions function. It will shift the residue number in a
    given line from a PDB file by <shift> positions.
    '''
    # Extract the residue's id
    resid = line[22:26].strip(' ')
    # Check how many digits the resid has
    resid_dig = len(resid)
    # Apply the shift 
    resid = str(int(resid) + shift)
    # Update the numbering
    new_line = line.replace(line[22:27], ' '*(4-resid_dig) + resid + ' ')
    return(new_line)

#############

# Open the files
in_handle = open(infile, 'r')
out_handle = open(outfile, 'w')

# This counter for insertion residues will determine the shift in numbering
chain = ''
shift = 0
insertion_residues = {}
prev_line = ''

for line in in_handle:
    # Write directly to the outfile all the lines that do not contain information for atoms
    if line.startswith('ATOM'):
        # Parse the line and check for insertion residues
        # Set the shift to 0 if this is a new chain
        if line[21] != chain:
            chain = line[21]
            shift = 0
            if insertion_residues.get(chain, -1) == -1:
                insertion_residues[chain] = []

        if line[16] != ' ':
            # This is an atom with more than one coordinates
            if prev_line == '':
                # This is the first of the multiple coordinates
                prev_line = line
            else:
                # Compare both lines to know which one has the better occupancy
                best_line = check_alternative_coords(line, prev_line)
                # Write the line with the best occupancy to the file
                out_handle.write(best_line)
                # Restart the prev_line variable
                prev_line = ''

        # Look for insertion residues with the regular expression
        elif re.search('\d+[A-Z]', line[22:27].strip(' ')):

            resid = line[22:27].strip(' ')
            # The shift must be added only once per residue, not once per atom
            if resid in insertion_residues[chain]:
                # A previous atom of this residue has been considered in the shift.
                # Apply the shift.
                new_line = shift_line(line, shift)
                out_handle.write(new_line)
            else:
                # Add one to the shift and add this residue to insertion_residues
                shift = shift + 1
                insertion_residues[chain].append(resid)
                new_line = shift_line(line, shift)
                out_handle.write(new_line)

        else:
            # Apply the shift to this residue
            new_line = shift_line(line, shift)
            out_handle.write(new_line)

    elif not line.startswith('ANISOU'):
        # Ignore ANISOU records
        out_handle.write(line)


in_handle.close()
out_handle.close()

# Print a list of insertion residues in this structure
if shift > 0:
    sys.stdout.write('Insertion residues: ')

    for key, val_list in insertion_residues.items():
	sys.stdout.write('chain' + key + '\n')
	for residue in val_list:
	    sys.stdout.write(residue + '\n')

