#!/bin/bash

# This script provides the list of parameters for the admin script to use during
# the simulations.

# Current directory
curr_dir=$PWD

# the name of the queue in which to run the main loop
export queue=<queue>

# Choice of substitution matrix
# 0 for an unconstrained selection of substitutions
# 1 for the Thorvaldsen matrix (restricted by the genetic code)
export subs_matrix=1

# rotabase must contain the path to the rotabase.txt file that FoldX needs
export rotabase=</path/to/rotabase.txt>

# sub_gen must contain the path to the substitution generator script
export sub_gen=${curr_dir}/005.4_substitution_generator.py

# regions indicates which regions will be allowed for substitutions
# "0.75 1.00" for interfaces
export regions="0.75 1.00"

# mut_total should contain the total number of substitutions to be performed
export mut_total=200

# rep_total should contain the total number of replicates to be done
export rep_total=50

# folder_tag should contain a tag that will be used to identify the output folder
export folder_tag=$(date "+%Y-%m-%d-%H-%M-%S")

# ncores is the number of cores that will be used
# Parameter for SLURM
export ncores=10

# memory will contain the amount of memory that will be used for this run (in Mb)
# Parameter for SLURM
export memory=8000

# apply_selection will contain the path to the Rscript that will apply the selection criteria
export apply_selection=${curr_dir}/005.5_apply_selection_sims.R

# fitness_program will indicate how the thresholds for fitness will be used on homodimers
# 0 for selection based only on binding energy (Kachroo's low stability)
# 1 for selection based only on protein stability (Kachroo's non-bound)
# 2 for selection based on both protein stability and binding energy (Kachroo's wildtype)
export fitness_program=2

# scenario will indicate the selection scenario for the coevolution of paralogs:

# scenario = 1 will be selection on both HMs (the substitutions are passed to the next round if the two homodimers are accepted)
# scenario = 2 will be selection on the HET (the substitutions are passed to the next round if the heterodimer was accepted)
# scenario = 3 will be selection on homodimer AA (the substitutions are passed to the next round if the homodimer A is accepted)
# scenario = 4 will be selection on homodimer BB (the substitutions are passed to the next round if the homodimer B is accepted)
export scenario=4

# The population size parameter for the probability of acceptance function
export pop_size=1000

# The value of the fitness function to be used:
# 0 for the exponential function
# 1 for stabilizing selection with two exponential functions
export fit_func=0

# The beta value for the shape of the exponential function (if chosen)
# Neutral evolution is achieved by setting beta to 0
export beta=10


# The length of the plateau at the top of the stabilizing selection function, as to indicate a region of maximum 
export plateau_length=0.5



