#!/bin/bash
# Laura Dean
# 2/9/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=extract_region_from_bam
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=20g
#SBATCH --time=6:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out


# activate software
source $HOME/.bash_profile
conda activate samtools1.22


# set variables
wkdir=/gpfs01/home/mbzlld/data/hagfish/OPR4
cd $wkdir
bam=all_reads_2_chrs_only_asm.bam


# extract region of interest from bam then sort and index it
# added 1 million BP to each end of the region so I can get a good look in IGV
samtools view -b $bam "cluster2:43976023-46132782" | samtools sort -o ${bam%.*}_OPR4_region.bam
samtools index ${bam%.*}_OPR4_region.bam


# deactivate software
conda deactivate

