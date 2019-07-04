#!/bin/bash

# This script downloads a file from the PDB. It is based on the thread at (https://nsaunders.wordpress.com/2008/01/14/rapid-command-line-access-to-the-pdb/). It receives the following argument:
# $1 = 4-letter ID for the file to download

wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${1}.ent.gz
gunzip pdb${1}.ent.gz
mv pdb${1}.ent ${1}.pdb

