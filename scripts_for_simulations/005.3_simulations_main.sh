#!/bin/bash

#### This is the main script for the simulations with paralog coevolution ####
# It applies mutations and obtains their energetic parameters, calls the script
# for selection, and manages the output files.


###############################
# Arguments:
# $1 = replicate number
# $2 = chain names
# $3 = main structure of the protein filename
# $4 = script used to source all the variables
###############################

rep=$1
chains=$2
prot=$3
config_file=`echo $4 | rev | cut -f 1 -d '/' | rev`

# Load the variables
source ./${config_file}

mkdir $rep
cd $rep
ln -s $rotabase rotabase.txt

# Add symlink to the starting file
ln -s ./../${prot}.pdb ${prot}.pdb

# Get chains
chainA=$(echo $chains | cut -c 1)
chainB=$(echo $chains | cut -c 3)

# Go through the different replicates
for mut in $(seq 1 ${mut_total})
do

	mkdir $mut
	cd $mut

	mkdir homodimer_${chainA} homodimer_${chainB} heterodimer

	# Generate the substitutions
	# First substitutions are generatedfrom the original structure
	# All the others are generated from the last accepted complex
	if [ $mut -eq 1 ]
	then
		ln -s ../../${prot}.pdb
		python $sub_gen -i ${prot}.pdb -t $regions -c False -p 1 -m ${subs_matrix}
		previous_deltaG=../../../original_deltaGs.tab
	else
		prev=$(($mut-1))

		# Take the last accepted substitution
		if [[ $accepted ]]
		then
			ln -s ../${prev}/heterodimer/Optimized_${prot}_${prev}.pdb ./
			python $sub_gen -i Optimized_${prot}_${prev}.pdb -t $regions -c False -p 1 -m ${subs_matrix}
		else
			ln ../${prev}/heterodimer/previous/previous.pdb ./
			python $sub_gen -i previous.pdb -t $regions -c False -p 1 -m ${subs_matrix}
		fi

	fi

	for subfolder in homodimer_${chainA} homodimer_${chainB} heterodimer
	do
		# Enter the subfolder
		cd $subfolder
	
		if [ $mut -eq 1 ]
		then
			# Save a copy of the previous file
			mkdir previous
			cd previous
			ln ../../../../${prot}.pdb ./previous.pdb
			cd ..

            # Bring the files to this folder
            ln -s ../../../${prot}.pdb ./
            ln -s ${rotabase} ./

			# Apply the substitutions
			FoldX --command=BuildModel --pdb=${prot}.pdb --mutant-file=individual_list.txt --ionStrength=0.05 --pH=7 --water=CRYSTAL --vdwDesign=2 --out-pdb=true --pdbHydrogens=false --numberOfRuns=1 > /dev/null

            # Remove unnecessary files to reduce disk usage
            rm ${prot}.pdb WT_${prot}_1.pdb

		else

			prev=$(($mut-1))

			# Bring the files from the previous substitution			
			if [[ $accepted ]]
			then

				ln -s ../../${prev}/${subfolder}/Optimized_${prot}_${prev}.pdb ./
				ln -s ${rotabase} ./
				mkdir previous
				cd previous
				ln ../../../${prev}/${subfolder}/Optimized_${prot}_${prev}.pdb ./previous.pdb
				cd ..
			else

				ln ../../${prev}/${subfolder}/previous/previous.pdb ./Optimized_${prot}_${prev}.pdb
				ln -s ${rotabase} ./
				mkdir previous
				cd previous
				ln ../../../${prev}/${subfolder}/previous/previous.pdb ./
				cd ..

			fi
			
			# Apply the substitution
			FoldX --command=BuildModel --pdb=Optimized_${prot}_${prev}.pdb --mutant-file=individual_list.txt --ionStrength=0.05 --pH=7 --water=CRYSTAL --vdwDesign=2 --out-pdb=true --pdbHydrogens=false --numberOfRuns=1 > /dev/null
			
			# Rename the files that BuildModel produces
			mv Optimized_${prot}_${prev}_1.pdb ${prot}_${mut}.pdb


            # Remove unnecessary files to reduce disk usage
            rm Optimized_${prot}_${prev}.pdb WT_Optimized_${prot}_${prev}_1.pdb

		fi

		# Optimize 
		FoldX --command=Optimize --pdb=${prot}_${mut}.pdb > /dev/null

		# Analyze complex
		FoldX --command=AnalyseComplex --pdb=Optimized_${prot}_${mut}.pdb --analyseComplexChains=$chains --complexWithDNA=False > /dev/null

		# Use the chains variable to know the names of the chains
		chainA=$(echo $chains | cut -c 1)
		chainB=$(echo $chains | cut -c 3)

		# Separate the two chans to calculate their stability with FoldX
		grep -E '.{21}'${chainA} Optimized_${prot}_${mut}.pdb > Optimized_${prot}_${mut}_${chainA}.pdb
		grep -E '.{21}'${chainB} Optimized_${prot}_${mut}.pdb > Optimized_${prot}_${mut}_${chainB}.pdb

		# Get the protein stability for both chains
		FoldX --command=Stability --pdb=Optimized_${prot}_${mut}_${chainA}.pdb > /dev/null
		FoldX --command=Stability --pdb=Optimized_${prot}_${mut}_${chainB}.pdb > /dev/null

		# Remove unnecessary files to reduce disk usage
		rm Optimized_${prot}_${mut}_${chainA}.pdb Optimized_${prot}_${mut}_${chainB}.pdb
 
		# Collect tbe data into a single file (separated by spaces)
		paste <(echo ChainA_stability) <(echo ChainB_stability) <(echo Binding_energy) > deltaGs.tab
		paste <(cut -f 2 Optimized_${prot}_${mut}_${chainA}_0_ST.fxout) <(cut -f 2 Optimized_${prot}_${mut}_${chainB}_0_ST.fxout) <(cat Interaction* | sed '10q;d' | cut -f 6) >> deltaGs.tab

		outfile=Selection_${mut}
		reference_file=../../../original_deltaGs.tab

        # If this is the first mutation, the previous_deltaG will be the original one
        if [ $mut -eq 1 ]
        then
            previous_deltaG=${reference_file}
        else
            previous_deltaG=../../${prev}/${subfolder}/deltaGs_after_selection_${prev}.tab
        fi

		# Apply selection
		if [ $fit_func -eq 0 ]
		then
			Rscript $apply_selection -c deltaGs.tab \
			-p $previous_deltaG -o $outfile -m $fitness_program -r $reference_file \
			-N $pop_size -f $fit_func -b $beta -s $threshold_factor
		elif [ $fit_func -eq 1 ]
		then
		Rscript $apply_selection -c deltaGs.tab \
                        -p $previous_deltaG -o $outfile -m $fitness_program -r $reference_file \
                        -N $pop_size -f $fit_func -b $beta -s $threshold_factor -l $plateau_length
		fi

        # Remove unnecessary files to reduce disk usage
        rm Summary* Average* Dif* Raw* Interface_Residues* OP*fxout PdbList* Indiv_energies*

		cd ..
	done

	# Decide if the mutations will be accepted or not
	
	# scenario = 1 will be selection on both HMs (the substitutions are passed to the next round if the two homodimers are accepted)
	# scenario = 2 will be selection on heterodimer AB (the substitutions are passed to the next round if the heterodimer was accepted)
	# scenario = 3 will be selection on homodimer AA (the substitutions are passed to the next round if the homodimer AA is accepted)
	# scenario = 4 will be selection on homodimer BB (the substitutions are passed to the next round if the homodimer BB is accepted)
	if [ ${scenario} -eq 0 ]
	then
		accepted=$(grep 'ACCEPTED' */${outfile}_verdict.txt)
	elif [ ${scenario} -eq 1 ]
	then
		check1=$(grep 'ACCEPTED' homodimer_${chainA}/${outfile}_verdict.txt)
		check2=$(grep 'ACCEPTED' homodimer_${chainB}/${outfile}_verdict.txt)

		if [[ $check1 ]] && [[ $check2 ]]
		then
			# Both were accepted
			accepted='ACCEPTED'
		else
			# Else, I save it as an empty variable
			accepted=
        fi
	elif [ ${scenario} -eq 2 ]
    then
		accepted=$(grep 'ACCEPTED' heterodimer/${outfile}_verdict.txt)
    elif [ ${scenario} -eq 3 ]
    then
        accepted=$(grep 'ACCEPTED' homodimer_${chainA}/${outfile}_verdict.txt)
    elif [ ${scenario} -eq 4 ]
    then
        accepted=$(grep 'ACCEPTED' homodimer_${chainB}/${outfile}_verdict.txt)
	fi

	# Keep track of the deltaGs after selection
	if [[ $accepted ]]
	then
		# Accepted substitution
		ln homodimer_${chainA}/deltaGs.tab homodimer_${chainA}/deltaGs_after_selection_${mut}.tab
		ln homodimer_${chainB}/deltaGs.tab homodimer_${chainB}/deltaGs_after_selection_${mut}.tab
		ln heterodimer/deltaGs.tab heterodimer/deltaGs_after_selection_${mut}.tab
	else
		# Rejected substitution, back to the last accepted one
		if [ $mut -eq 1 ]
		then
			ln ../../original_deltaGs.tab homodimer_${chainA}/deltaGs_after_selection_${mut}.tab
            ln ../../original_deltaGs.tab homodimer_${chainB}/deltaGs_after_selection_${mut}.tab
            ln ../../original_deltaGs.tab heterodimer/deltaGs_after_selection_${mut}.tab
		else
			# Then I need to save deltaGs_after_selection.tab as the old variables
			ln ../$prev/homodimer_${chainA}/deltaGs_after_selection_${prev}.tab homodimer_${chainA}/deltaGs_after_selection_${mut}.tab
			ln ../$prev/homodimer_${chainB}/deltaGs_after_selection_${prev}.tab homodimer_${chainB}/deltaGs_after_selection_${mut}.tab
			ln ../$prev/heterodimer/deltaGs_after_selection_${prev}.tab heterodimer/deltaGs_after_selection_${mut}.tab
		fi
	fi	

	cd ..

done

# Concatenate the deltaGs for each homodimer and heterodimer for every replicate.
for subfolder in heterodimer homodimer_${chainA} homodimer_${chainB}
do	
	# Concatenate the list of deltaGs for each type of complex
	cat <(cat ../original_deltaGs*) <(ls */${subfolder}/deltaGs_after_selection* | sort -g | xargs -I {} sed '2q;d' {}) > all_deltaGs_results_${subfolder}.txt 
 
	# Get the list of substitutions
	ls */${subfolder}/individual_list.txt | sort -g | xargs cat | sed -E 's/;/\n/g; s/,/\t/g' > all_substitutions_${subfolder}_extracted.tab

	# Format the list so that it will have a tabular format
	paste <(cut -f 1 all_substitutions_${subfolder}_extracted.tab | cut -c 1) <(cut -f 1 all_substitutions_${subfolder}_extracted.tab | rev | cut -c 1) <(cut -f 1 all_substitutions_${subfolder}_extracted.tab | grep -o -E '[0-9]+') \
	<(cut -f 2 all_substitutions_${subfolder}_extracted.tab | cut -c 1) <(cut -f 2 all_substitutions_${subfolder}_extracted.tab | rev | cut -c 1) \
	<(cut -f 2 all_substitutions_${subfolder}_extracted.tab | grep -o -E '[0-9]+') \
	<(ls */${subfolder}/*verdict.txt | sort -g | xargs -I {} tail -n 1 {}) > all_substitutions_${subfolder}_formatted.tab

	# Add a header and a filler row for the starting values
	cat <(echo "Chain_A_subs_initial,Chain_A_subs_final,Chain_A_subs_position,Chain_B_subs_initial,Chain_B_subs_final,Chain_B_subs_position,Verdict" | sed -E 's/,/\t/g') \
	<(echo "NA,NA,NA,NA,NA,NA,NA" | sed -E 's/,/\t/g') \
	all_substitutions_${subfolder}_formatted.tab > all_substitutions_${subfolder}_header.tab

	paste all_deltaGs_results_${subfolder}.txt all_substitutions_${subfolder}_header.tab > all_deltaGs_results_${subfolder}_subs.tab

done
