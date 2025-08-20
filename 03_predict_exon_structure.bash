#!/bin/bash
# Laura Dean
# 19/8/25
# script written for running on the UoN HPC Ada

# script to predict the exon strucuture of the OPR4 protein based on OPR4 sequences from
# closely related / other species and the DNA sequence of the region containing the gene in
# our de novo assembly

#SBATCH --job-name=predict_exon_structure
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/hagfish/OPR4
#gene_region=$wkdir/OPR4_region.fasta
gene_region=$wkdir/OPR4_AA_region.fasta

cd $wkdir

# activate software
source $HOME/.bash_profile
#conda create --name exonerate bioconda::exonerate -y
conda activate exonerate

# predict the exon structure
exonerate --model protein2genome \
          --showtargetgff yes \
          --showalignment no \
          --query OPR4_AA_protein.fasta \
          --target $gene_region > ${gene_region%.*}_exonerate.out

# deactivate software
conda deactivate

