#!/bin/bash
# Laura Dean
# 20/8/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=map_reads_to_ref
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --mem=200g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out


# activate software
source $HOME/.bash_profile
conda activate minimap2


# set variables
wkdir=/gpfs01/home/mbzlld/data/hagfish/OPR4
cd $wkdir
reads=/gpfs01/home/mbzlld/data/hagfish/basecalls/native_and_pcr_calls.fastq.gz
assembly=/gpfs01/home/mbzlld/data/hagfish/flye_1/assembly_ragtag/ragtag.scaffold_chrs_only.fasta


# map the raw reads back to our assembly
minimap2 \
	-a \
	-x map-ont \
	--split-prefix temp_prefix \
	-t 32 \
	-o $wkdir/all_reads_2_chrs_only_asm.sam \
	$assembly $reads


# sort and index the sam file and convert to bam format
samtools sort \
	--threads 32 \
	--write-index \
	--output-fmt BAM \
	-o $wkdir/all_reads_2_chrs_only_asm.bam $wkdir/all_reads_2_chrs_only_asm.sam


# remove the intermediate sam file
rm $wkdir/all_reads_2_chrs_only_asm.sam



# deactivate software
conda deactivate

