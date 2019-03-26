import argparse
import sys
import re

parser = argparse.ArgumentParser('This code will receive an input PDB file and look for both insertion residues (specified with a residue number and a letter, e.g. PHE184A) and alternative coordinates for residues (specified by a residue name with a letter designation, e.g. ASER and BSER). For insertion residues, it will renumber accordingly as to integrate them in the numbering and make sure there are no repeated numbers. For alternative coordinates, it will keep the one with the highest occupancy because it is the one that would be found the most in the crystal. I will print messages on insertion residues and alternative coordinates to the standard output.')
parser.add_argument('-i', help = 'Input PDB', dest = 'infile', type = str)
parser.add_argument('-o', help = 'Output PDB', dest = 'outfile', type = str)

args = parser.parse_args()
infile = args.infile
outfile = args.outfile

def check_alternative_coords(current_line, prev_line):
    '''This function will be on the lookout for atoms that have more than one coordinate. I will keep the one
    with the best quality occupancy, as it would be the one that is better supported by the model. If the
    occupancies are the same, then I will keep the first orientation.
    '''
    # The occupancy in the variable "line"
    curr_occ = float(current_line[54:60].strip())
    prev_occ = float(prev_line[54:60].strip())
    
    atom_name = current_line[12:16].strip()
    resname = current_line[17:20].strip()
    resid = current_line[22:26].strip()
    
    sys.stdout.write('Alternative coordinates found for ' + atom_name + ' in ' + resname + resid + '.\n')

    # Compare the two occupancies
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
    # I will extract the residue's id
    resid = line[22:26].strip(' ')
    # I will check how many digits the resid has
    resid_dig = len(resid)
    # resid is a string, so I need to change it to int to add the shift and back to str to write it 
    resid = str(int(resid) + shift)
    # Now that the resid has been shifted, I can modify the line and write it
    # I add spaces at the start to make sure I can maintain the PDB file format (columns 23-26, 
    # or 22-25 with Python numbering) are the residue id, which has at most four digits in the
    # current format
    # The concatenated space at the end would take the place of the insertion residue code
    new_line = line.replace(line[22:27], ' '*(4-resid_dig) + resid + ' ')
    return(new_line)

#############

# Open the files
in_handle = open(infile, 'r')
out_handle = open(outfile, 'w')

# This will be counter for insertion residues, which will tell me if I had to shift the numbering
# I will also keep track of insertion residues
chain = ''
shift = 0
insertion_residues = {}
prev_line = ''

for line in in_handle:
    # I will write directly to the outfile all the lines that do not contain information for atoms
    if line.startswith('ATOM'):
        # Then I have to parse the line and check for insertion residues
        # Before doing anything, I need to check the chain. If I find a new chain, I will set
        # the shift back to 0
        if line[21] != chain:
            chain = line[21]
            shift = 0
            if insertion_residues.get(chain, -1) == -1:
                insertion_residues[chain] = []

        if line[16] != ' ':
            # Then, I have found an atom that has more than one coordinate
            # I will compare this atom's occupancy with that of the previous line or
            # save this line as prev_line if there is nothing to compare
            if prev_line == '':
                # Then, there is nothing to compare
                prev_line = line
            else:
                # I need to compare both lines to know which one has the better occupancy and write it to the file
                best_line = check_alternative_coords(line, prev_line)
                # Write the line with the best occupancy to the file
                out_handle.write(best_line)
                # Restart the prev_line variable
                prev_line = ''

        # I will remove spaces and look for the following pattern: <one or more digits><a letter>
        # Only insertion residues will match that pattern
        elif re.search('\d+[A-Z]', line[22:27].strip(' ')):
            # Then I found an insertion residue
            resid = line[22:27].strip(' ')
            # Before adding it to the shift, I must check if I have already seen it. Remember that I will
            # find the insertion code in the lines for all the atoms that belong to the insertion residue
            if resid in insertion_residues[chain]:
                # Then I have already considered this residue. I only need to shift it
                new_line = shift_line(line, shift)
                out_handle.write(new_line)
            else:
                # Then I need to add one to the shift and add this residue to insertion_residues
                shift = shift + 1
                insertion_residues[chain].append(resid)
                new_line = shift_line(line, shift)
                out_handle.write(new_line)

        else:
            # Then it is not an insertion residue, but I must shift this line's residue id and write the line
            new_line = shift_line(line, shift)
            out_handle.write(new_line)

    elif not line.startswith('ANISOU'):
        # I will write everything to the outfile except for the ANISOU (anisotropic temperature) records
        # I do not believe I need said records
        out_handle.write(line)


# Once I'm finished, I can close both files
in_handle.close()
out_handle.close()

if shift > 0:
    sys.stdout.write('Insertion residues: ')

    for key, val_list in insertion_residues.items():
	sys.stdout.write('chain' + key + '\n')
	for residue in val_list:
	    sys.stdout.write(residue + '\n')

