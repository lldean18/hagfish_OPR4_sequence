#!/bin/bash
# Laura Dean
# 19/8/25
# script written for running on the UoN HPC Ada

# script to concatenate DNA exon sequences and convert them to protein sequences

#SBATCH --job-name=convert_annotation_to_protein
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/hagfish/OPR4
annotation=$wkdir/OPR4_AA_region_exonerate.out


cd $wkdir


# activate software
source $HOME/.bash_profile
conda activate bedtools 

# === extract exons ===
grep -P "\texon\t" ${annotation} > ${annotation%.*}_OPR4_exons.gff
# === get the strand ===
STRAND=$(awk '{print $7}' ${annotation%.*}_OPR4_exons.gff | uniq)
# === convert exons to fasta ===
bedtools getfasta -fi ${annotation%_*}.fasta -bed ${annotation%.*}_OPR4_exons.gff -s -name > ${annotation%.*}_OPR4_AA_exons.fasta
# === merge exons into single CDS FASTA ===
# Remove headers, flatten, add new header
seqkit seq ${annotation%.*}_OPR4_AA_exons.fasta | grep -v ">" | tr -d '\n' > cds.tmp
ID=$(basename $annotation)  # create variable with the individual ID
ID=${ID%%_*}_AA  # create variable with the individual ID
echo ">${ID}_cds" > ${annotation%_*}_cds.fasta
cat cds.tmp >> ${annotation%_*}_cds.fasta
rm cds.tmp
# === reverse-complement if minus strand ===
#if [ "$STRAND" = "-" ]; then
#  echo "Reverse-complementing because gene is on minus strand..."
#  seqkit seq -r -p -t DNA ${annotation%_*}_cds.fasta > ${annotation%_*}_cds_rc.fasta
#  mv ${annotation%_*}_cds_rc.fasta ${annotation%_*}_cds.fasta
#fi
# === translate cds to protein ===
seqkit translate ${annotation%_*}_cds.fasta > ${annotation%_*}_protein.fasta


# deactivate software
conda deactivate
