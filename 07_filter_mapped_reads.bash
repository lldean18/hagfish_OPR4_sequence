#!/bin/bash
# Laura Dean
# 2/9/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=filter_mapped_reads
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=200g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out


# activate software
source $HOME/.bash_profile
conda activate samtools1.22


# set variables
wkdir=/gpfs01/home/mbzlld/data/hagfish/OPR4
cd $wkdir
bam=all_reads_2_chrs_only_asm.bam

# Keep only mapped, primary alignments with MAPQ 20 or more
samtools view -b \
  -q 20 \                # min MAPQ (adjust as needed)
  -F 0x900 \             # filter out secondary (0x100) and supplementary (0x800)
  $bam > ${bam%.*}_clean.bam

# Sort and index for downstream use
samtools sort -o filtered.sorted.bam filtered.bam
samtools index filtered.sorted.bam



# deactivate software
conda deactivate

