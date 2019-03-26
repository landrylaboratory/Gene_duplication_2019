#!/bin/bash

####Alignment of RNAseq reads
#Indicate the address where the samples are to be aligned
adresse="/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/RNAseq"
#Indicate the address and the name of the reference genome file
genome="/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Genomes/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa"
annotation="/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Genomes/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"
#Indicate the different samples to be aligned without the "fq" ending
echant=("DMSO1cat" "DMSO2cat" "DMSO3cat" "MTX1cat" "MTX2cat" "MTX3cat")

#1- index the reference
bwa index $genome

#2 - alignment
#"BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, 
#while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features 
#such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended 
#for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than 
#BWA-backtrack for 70-100bp Illumina reads"
#--> We took mem option

#-c INT	Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]
#-M	Mark shorter split hits as secondary (for Picard compatibility).
for i in ${!echant[@]} ;
        do
bwa mem -M -t 2 $genome $adresse/${echant[$i]}.fq > $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C.sam
	done

#3 - Check the proportion of reads that have mapped to several places
#XA: alternative alignment proposed by bwa
for i in ${!echant[@]} ;
        do
echo ${echant[$i]}_S288C.sam
tot=$(grep -v "@SQ" $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C.sam | wc -l)
XA=$(grep -c "XA:" $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C.sam)
echo $tot
echo $XA
echo "(100*($XA/$tot))" | bc -l
	done


#4- BAM format sorting
for i in ${!echant[@]} ;
        do
#Deletion of the reads that have mapped to several places marked by '/XA:/d' in the sam files
sed '/XA:/d' $adresse/${echant[$i]}_S288C.sam > $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C_no_several_map.sam
samtools view -bS  $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C_no_several_map.sam > $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C_no_several_map.bam

#Sort alignment based on alignment position of reads
samtools sort $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C_no_several_map.bam $adresse/alignment_bwa_S288C_SGD/${echant[$i]}_S288C_no_several_map_sort
	done

#5- Run htseq-count:
for i in ${!echant[@]} ;
        do
htseq-count -f=bam -t=gene $adresse/$alignment/${echant[$i]}_no_several_map_sorted.bam $adresse/$gff > $adresse/$alignment/${echant[$i]}_count.table
        done