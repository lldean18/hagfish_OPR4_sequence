#!/bin/bash
# Laura Dean
# 19/8/25
# script written for running on the UoN HPC Ada

# script to extract the region containing the OPR4 gene from our de novo hagfish assembly

#SBATCH --job-name=extract_gene_region
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10g
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/hagfish
genome=$wkdir/flye_1/assembly_ragtag/ragtag.scaffold_chrs_only.fasta

cd $wkdir/OPR4

# activate software
source $HOME/.bash_profile
#conda create --name samtools1.22 bioconda::samtools=1.22 bioconda::bcftools=1.22 -y
conda activate samtools1.22

# extract the gene region
samtools faidx --threads 15 $genome
samtools faidx --threads 15 $genome 'cluster2:44986194-45132782' > OPR4_region.fasta

# deactivate software
conda deactivate
