#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# This script takes the FASTA formatted interface files and outputs per chain counts of total residues and interface residues.

import argparse
from collections import OrderedDict
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(description = 'This script will receive an input FASTA file and return counts of lowercase (interfaces) and total residues.')

parser.add_argument('-i', type = str, help = 'The input file to work with', dest = 'input')

args = parser.parse_args()
infile= args.input

def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

# Parse the file
records = SeqIO.parse(infile, 'fasta')

for record in records:

    sys.stdout.write(record.description+ '\n')
    sys.stdout.write(str(n_lower_chars(record.seq))+ '\n')
    sys.stdout.write(str(len(record.seq))+ '\n')
    sys.stdout.write('--------'+ '\n')

