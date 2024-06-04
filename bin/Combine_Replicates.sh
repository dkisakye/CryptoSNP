#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH -t 60
#SBATCH -p agsmall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Combine_Replicates_%j.err
#SBATCH -o Combine_Replicates_%j.out

# Combine the replicates for the 500 strain GWAS project

# Path to the file with the replicated sample names
REPLICATE_SAMPLE_NAMES="/home/Cryptococcus_project/GitHub_Repo/Results/500_Strain_GWAS/Replicate_Strain_Names.txt"
# Path to where the raw reads are stored
READS_DIR="/home/Cryptococcus_project/data"
# Set the output directory
OUT_DIR="/home/Cn_GWAS/Raw_Reads"

mkdir -p "${OUT_DIR}"
# For each sample, find the reads and concatenate them
for samplename in $(cat "${REPLICATE_SAMPLE_NAMES}")
do
    cat $(find "${READS_DIR}" -type f -name "*${samplename}*R1*" | sort -V) > "${OUT_DIR}/${samplename}_R1.fastq.gz"
    cat $(find "${READS_DIR}" -type f -name "*${samplename}*R2*" | sort -V) > "${OUT_DIR}/${samplename}_R2.fastq.gz"
done
