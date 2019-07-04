##!/bin/bash

# This script will go to a folder with the output of the simulations and gather the fixed and proposed mutations in tables
# in a subfolder called "Final_results"

# $1 = path to the directory with results from simulations
# $2 = number of replicates in that simulation

cd $1
num_reps=$2

mkdir Final_results

# Gather the fixed mutations for all replicates
cat <(head -n 1 1/all_deltaGs_results_heterodimer_subs.tab) <(ls */all_deltaGs_results_heterodimer_subs.tab | sort -g | xargs -I {} tail -n +2 {}) > Final_results/all_deltaGs_heterodimer_all_reps.tab
cat <(head -n 1 1/all_deltaGs_results_homodimer_A_subs.tab) <(ls */all_deltaGs_results_homodimer_A_subs.tab | sort -g | xargs -I {} tail -n +2 {}) > Final_results/all_deltaGs_homodimer_A_all_reps.tab
cat <(head -n 1 1/all_deltaGs_results_homodimer_B_subs.tab) <(ls */all_deltaGs_results_homodimer_B_subs.tab | sort -g | xargs -I {} tail -n +2 {}) > Final_results/all_deltaGs_homodimer_B_all_reps.tab

# Gather the data for proposed mutations for each replicate
for replicate in $(seq 1 ${num_reps})
do
    cd $replicate
    cat <(cat ../original_deltaGs*) <(ls */heterodimer/deltaGs.tab | sort -g | xargs -I {} sed '2q;d' {}) > all_deltaGs_proposed_heterodimer.txt
    cat <(cat ../original_deltaGs*) <(ls */homodimer_A/deltaGs.tab | sort -g | xargs -I {} sed '2q;d' {}) > all_deltaGs_proposed_homodimer_A.txt
    cat <(cat ../original_deltaGs*) <(ls */homodimer_B/deltaGs.tab | sort -g | xargs -I {} sed '2q;d' {}) > all_deltaGs_proposed_homodimer_B.txt
    
    paste all_deltaGs_proposed_heterodimer.txt all_substitutions_heterodimer_header.tab > all_deltaGs_proposed_heterodimer_subs.tab
    paste all_deltaGs_proposed_homodimer_A.txt all_substitutions_homodimer_A_header.tab > all_deltaGs_proposed_homodimer_A_subs.tab
    paste all_deltaGs_proposed_homodimer_B.txt all_substitutions_homodimer_B_header.tab > all_deltaGs_proposed_homodimer_B_subs.tab
    
    cd ..
done



# Gather the proposed mutations for all replicates
cat <(head -n 1 1/all_deltaGs_proposed_heterodimer_subs.tab) <(ls */all_deltaGs_proposed_heterodimer_subs.tab | sort -g | xargs -I {} tail -n +2 {}) > Final_results/all_deltaGs_proposed_heterodimer_all_reps.tab
cat <(head -n 1 1/all_deltaGs_proposed_homodimer_A_subs.tab) <(ls */all_deltaGs_proposed_homodimer_A_subs.tab | sort -g | xargs -I {} tail -n +2 {}) > Final_results/all_deltaGs_proposed_homodimer_A_all_reps.tab
cat <(head -n 1 1/all_deltaGs_proposed_homodimer_B_subs.tab) <(ls */all_deltaGs_proposed_homodimer_B_subs.tab | sort -g | xargs -I {} tail -n +2 {}) > Final_results/all_deltaGs_proposed_homodimer_B_all_reps.tab


