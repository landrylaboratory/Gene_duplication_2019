#!/bin/bash

#### This script will be used for the simulations on the coevolution of homodimers and heterodimers of paralogs. ####

################################################################
# Arguments:
# $1 = The path to the input pdb file's name without the '.pdb' extension
# $2 = The config file with all the other parameters that the simulations need
################################################################

infile=$1
config_file=$2

# Load the parameters from the config file
source ${config_file}

# Create a new folder for the results
prot=`echo $infile | rev | cut -d'/' -f 1 | rev`
mkdir -p ../007_coevolution/${prot}_${folder_tag}
cp $infile.pdb ../007_coevolution/${prot}_${folder_tag}/${prot}.pdb

# Save the config file with the list of parameters used
cp $config_file ../007_coevolution/${prot}_${folder_tag}/
cd ../007_coevolution/${prot}_${folder_tag}

#### Write a README.txt file with the command line I used for this simulation ####

cat > README.txt << EOF
# This code receives the following arguments:
# $1 = the path to the input pdb file's name without the '.pdb' extension
# $2 = the config file used for this run (copied to this directory)

The command line used for the simulations in this folder was:
$0 $1 $2

EOF

# Inside the new folder, create a symlink the rotabase.txt file
ln -s $rotabase rotabase.txt

# Get the binding energy of the starting complex
dum=`grep '^ATOM' ${prot}.pdb | cut -c 22 | sort | uniq`
chains=`echo $dum | sed 's/ /,/g'`
FoldX --command=AnalyseComplex --pdb=${prot}.pdb --analyseComplexChains=$chains --complexWithDNA=False

# Isolate the two chains
chainA=`echo $chains | cut -c 1`
chainB=`echo $chains | cut -c 3`

grep -E '.{21}'${chainA} ${prot}.pdb > ${prot}_${chainA}.pdb
grep -E '.{21}'${chainB} ${prot}.pdb > ${prot}_${chainB}.pdb

# Calculate the stability for the two chains
FoldX --command=Stability --pdb=${prot}_${chainA}.pdb > stability_chain${chainA}.log
FoldX --command=Stability --pdb=${prot}_${chainB}.pdb > stability_chain${chainB}.log

# Collect tbe data into a single file
paste <(echo ChainA_stability) <(echo ChainB_stability) <(echo Binding_energy) > original_deltaGs.tab
paste <(cut -f 2 ${prot}_${chainA}_0_ST.fxout) <(cut -f 2 ${prot}_${chainB}_0_ST.fxout) <(cat Interaction* | sed '10q;d' | cut -f 6) >> original_deltaGs.tab

# Save the main directory for this protein
export MAINDIR=$PWD

######### Write the SLURM script ##########

cat > ${prot}_sim.slurm << EOF
#!/bin/bash

#SBATCH -D $MAINDIR
#SBATCH -J ${prot}_sim
#SBATCH -o ${prot}_sim.out
#SBATCH -c $ncores
#SBATCH -p $queue
#SBATCH --time=4-00:00
#SBATCH --mem=$memory

# Define some parameters for the SLURM parallel run
# srun="srun -c1"
parallel="parallel -N 1 -j \$SLURM_CPUS_PER_TASK --joblog parallel_joblog --resume"

# Start the simulations with the repaired initial structure
# A loop for the replicates. 
reps=\`seq 1 ${rep_total}\`

date

# Parallelize with GNU parallel
\$parallel "./005.3_simulations_main.sh {1} ${chains} ${prot} ${config_file}" ::: \$reps

date

EOF

# Submit the job to SLURM
sbatch ${prot}_sim.slurm 
