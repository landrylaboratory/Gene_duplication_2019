#!/bin/bash

curr_dir=$PWD

# the name of the queue in which to run the main loop
# low-suspend was used before for manitou or ibismini for katak
# The current queue to which I send the simulations in Manitou is medium
export queue=medium

# subs_matrix can be 0 for an unconstrained selection of substitutions, 1 for the Thorvaldsen matrix, or 2 for the matrix with Zhu-s data and codon usage
export subs_matrix=1

# rotabase must contain the path to the rotabase.txt file
export rotabase=</path/to/rotabase.txt>

# sub_gen must contain the path to the substitution generator script
export sub_gen=${curr_dir}/005.4_substitution_generator.py

# Should be set to False always because the code will work with the heterodimer file as to generate all candidate substitutions and to explore the respective homodimers accordingly.
export homodimer_check=False

# regions indicates which regions will be allowed for substitutions
# "0.00 0.25 0.50" for non-interfaces or "0.75 1.00" for interfaces or "0.00 0.25 0.50 0.75 1.00" for all regions
export regions="0.75 1.00"

# mut_total should contain the total number of substitutions to be performed
export mut_total=200

# rep_total should contain the total number of replicates to be done
export rep_total=50

# folder_tag should contain a tag that will be used to identify the output folder
# I will save the tag based on the date because this file with all the parameters will be copied inside.
export folder_tag=$(date "+%Y-%m-%d-%H-%M-%S")

# ncores is the number of cores that will be used
export ncores=10

# memory will contain the amount of memory that will be used for this run (in Mb)
export memory=8000

# apply_selection will contain the path to the Rscript that will apply the selection criteria
export apply_selection=${curr_dir}/apply_selection_sims_2.R

# fitness_program will indicate how the thresholds for fitness will be used on homodimers
# 0 for selection based only on binding energy (Kachroo's low stability)
# 1 for selection based only on protein stability (Kachroo's non-bound)
# 2 for selection based on both protein stability and binding energy (Kachroo's wildtype)
export fitness_program=2

# scenario will indicate the selection scenario for the coevolution of paralogs:
# scenario = 0 will be redundance (the substitutions are passed to the next round if at least one of the three complexes was accepted)
# scenario = 1 will be selection on both HMs (the substitutions are passed to the next round if the two homodimers are accepted)
# scenario = 2 will be selection on the HET (the substitutions are passed to the next round if the heterodimer was accepted)
# scenario = 3 will be selection on homodimer AA (the substitutions are passed to the next round if the homodimer A is accepted)
# scenario = 4 will be selection on homodimer BB (the substitutions are passed to the next round if the homodimer A is accepted)
export scenario=4

# The population size parameter for the probability of acceptance function
export pop_size=1000

# The value of the fitness function to be used (0 for the exponential or 1 for the stabilizing selection)
# Simulations shown for the paper were run with the exponential fitness function
export fit_func=0

# The beta value for the shape of the exponential function (if chosen)
# Neutral evolution is achieved by setting beta to 0
export beta=10

# The shape (k) parameter for the gamma distribution that will be used to implement stabilizing selection (if chosen)
export shape=2

# The scale (theta) parameter for the gamma distribution that will be used to implement stabilizing selection (if chosen)
export scale=2

# The length of the plateau at the top of the gamma distribution, as to indicate a region of maximum 
export plateau_length=0.5



