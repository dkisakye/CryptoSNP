#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 12:00:00
#SBATCH -p agsmall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o Reads_QC_%j.out
#SBATCH -e Reads_QC_%j.err

# Screen the reads for duplication, adapter sequences, and low-quality bases.

module load fastqc/0.11.9

READS_DIR="/home/Cryptococcus_project/data"
OUT_DIR="/home/Cryptococcus_project/results/FastQC_Reports"

mkdir -p "${OUT_DIR}"

fastqc \
    -f fastq \
    -o "${OUT_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    $(find "${READS_DIR}" -type f -name '*.fastq.gz')
