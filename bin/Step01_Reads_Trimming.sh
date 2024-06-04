#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 16:00:00
#SBATCH -p agsmall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o Read_Trimming_%j.out
#SBATCH -e Read_Trimming_%j.err

# Clean the reads of adapter contamination, low-quality bases, and reads that
# are too short post-cleaning.
module load parallel/20210822
module load trimmomatic/0.39

# Do a little arithmetic. We will run trimmomatic jobs in parallel, but keep in
# mind that trimmomatic is multithreaded
TRIM_THREADS="4"
N_PAR=$(echo "scale=0; ${SLURM_CPUS_PER_TASK} / ${TRIM_THREADS} / 1" | bc -l)

# Path to input and output directories
READS_DIR="/home/Cryptococcus_project/data"
LOG_DIR="/home/Cryptococcus_project/results/Combined_Replicates/Trimmomatic_Summaries"
OUT_DIR="/scratch.global/konox006/Cn_GWAS/Trimmed_Reads"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# Path to adapters file
ADAPTERS="/home/Cryptococcus_project/GitHub_Repo/Resources/all_illumina_adapters.fa"

# Path to replicated sample names
REPLICATE_SAMPLE_NAMES="/home/Cryptococcus_project/500_Strain_GWAS/Replicate_Strain_Names.txt"

# For each R1, skipping the replicates:
for fwd_read in $(find "${READS_DIR}" -type f -name '*R1_001.fastq.gz' | grep -vf "${REPLICATE_SAMPLE_NAMES}")
do
    # This assumes that the reads are named in a standard Illumina format.
    # Adjust this to match your file naming scheme.
    sname=$(basename "${fwd_read}" | sed -r 's/_R1_001\.fastq\.gz//g')
    rev_read="$(dirname ${fwd_read})/${sname}_R2_001.fastq.gz"
    # Make a new variable that has the strain name without the _SXXX suffix
    real_sname=$(basename "${fwd_read}" | sed -r 's/_S[0-9]+_R1_001\.fastq\.gz//g')
    # Write the trimming command
    echo "java -jar ${TRIMMOMATIC}/trimmomatic-0.39.jar PE \
        -threads ${TRIM_THREADS} \
        -summary ${LOG_DIR}/${real_sname}_summary.txt \
        ${fwd_read} \
        ${rev_read} \
        ${OUT_DIR}/${real_sname}_1P.fastq.gz \
        ${OUT_DIR}/${real_sname}_1U.fastq.gz \
        ${OUT_DIR}/${real_sname}_2P.fastq.gz \
        ${OUT_DIR}/${real_sname}_2U.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18"
done | parallel -j "${N_PAR}" --joblog "${LOG_DIR}/UnRep_Read_Trimming.joblog"

# Trim the replicated samples, too:
REP_READS_DIR="/home/Cn_GWAS/Raw_Reads"
for fwd_read in $(find "${REP_READS_DIR}" -type f -name '*R1.fastq.gz')
do
    # The combined replicates do not have the _SXXX suffix on their names.
    sname=$(basename "${fwd_read}" | sed -e 's/_R1\.fastq\.gz//g')
    rev_read="$(dirname ${fwd_read})/${sname}_R2.fastq.gz"
    # Write the trimming command
    echo "java -jar ${TRIMMOMATIC}/trimmomatic-0.39.jar PE \
        -threads ${TRIM_THREADS} \
        -summary ${LOG_DIR}/${sname}_summary.txt \
        ${fwd_read} \
        ${rev_read} \
        ${OUT_DIR}/${sname}_1P.fastq.gz \
        ${OUT_DIR}/${sname}_1U.fastq.gz \
        ${OUT_DIR}/${sname}_2P.fastq.gz \
        ${OUT_DIR}/${sname}_2U.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18"
done | parallel -j "${N_PAR}" --joblog "${LOG_DIR}/Rep_Read_Trimming.joblog"
