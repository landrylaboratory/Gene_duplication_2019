##########################################################################################
##########################################################################################
### STEP 1: FROM BAM TO FASTQ

module load samtools
my_bam=/home/chebe1/Transcriptomics_Paradoxus/Rawdata/Fastq_Transcriptomics_CE_2017/Chris3
my_dir=/home/chebe1/Transcriptomics_Paradoxus/Rawdata/Fastq_Transcriptomics_CE_2017
for sample in $(ls ${my_bam}/*.bam) ;
do
base=$(basename $sample)
name=$(echo $base | rev | cut -d”.” -f2 | rev)
echo $name
samtools bam2fq ${my_bam}/${name}.bam > ${my_dir}/${name}.fq
done

##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
### STEP 2: CUTADAPT TO REMOVE POOR QUALITY PARTS OF EACH SEQUENCE

# Load the program
module load cutadapt

# Defining working condition
sourcedir=/home/chebe1/Transcriptomics_Paradoxus/Rawdata/Fastq_Transcriptomics_CE_2017
wkdir=/home/chebe1/Transcriptomics_Paradoxus/Analysis/Analysis_Transcriptome_Evolution_48samples


#### CONDITIONS
for file in ${sourcedir}/*.fq ; 

do
base=$(basename $file)
name=$(echo $base | rev | cut -d”.” -f2 | rev)
echo $name 

# 1 # I remove the first 12 bp, as recommended by LEXOGEN
cutadapt -u 12 -o ${wkdir}/${name}_1.fq ${sourcedir}/${name}.fq

# 2 # Trimming the poly-A tail from the 3’ end of your reads, use the 3’ adapter type (-a) with an adapter sequence of many repeated A nucleotides; I use a sequence of 100, though even when shorterthe end 
cutadapt -a “A{100}” -o ${wkdir}/${name}_2.fq ${wkdir}/${name}_1.fq

# 3 # Quality trimming
cutadapt -q 15 -o ${wkdir}/${name}_3.fq ${wkdir}/${name}_2.fq

# 4 # Discard the reads that are smaller than 30 bp
cutadapt --minimum-length 30 -o ${wkdir}/${name}_4.fq ${wkdir}/${name}_3.fq ;

done

### REMOVING AND RENAMING
rm *_1.fq
rm *_2.fq
rm *_3.fq

for file in ${wkdir}/*_4.fq ; 
do
base=$(basename $file)
name=$(echo $base | rev | cut -d”.” -f2 | rev) # I cut the suffix off
echo $name ;
name_2=$(echo $name | rev | cut -d”_” -f2- | rev) # I cut the “_4” off
echo $name_2 ;
name_3=$(echo $name_2 | cut -d”_” -f3-) # I cut the the library info in the front
echo $name_3 ;
mv ${wkdir}/$name.fq ${wkdir}/${name_3}.fq ;
done
