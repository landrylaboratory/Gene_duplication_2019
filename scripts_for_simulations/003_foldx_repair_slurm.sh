#!/bin/bash

# This code receives the name of a file that should be looked inside the 002_ins_res_alt_crd folder and passes it on to the FoldX RepairPDB function.
# $1 = input PDB without the 'pdb' extension

cat > $1_repair.sbatch << EOF
#!/bin/bash

#SBATCH -D /home/afcis2/FoldX_simulations
#SBATCH -J $1_repair
#SBATCH -o $1_repair.out
#SBATCH -c 1
#SBATCH -p ibismini
#SBATCH --time=1-00:00
#SBATCH --mem=51200

cp 002_ins_res_alt_crd/$1.pdb 003_repair

cd 003_repair

ln -s `which rotabase.txt` rotabase.txt

FoldX --command=RepairPDB --pdb=$1.pdb --ionStrength=0.05 --pH=7 --water=CRYSTAL --vdwDesign=2 --out-pdb=true --pdbHydrogens=false > $1_Repair.log

mkdir $1_Repair

mv $1_Repair?* $1_Repair

rm $1.pdb
EOF

sbatch $1_repair.sbatch

