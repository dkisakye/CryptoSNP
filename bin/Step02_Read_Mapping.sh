#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --tmp=400gb
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 12:00:00
#SBATCH -p agsmall
#SBATCH -A nielsenk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o Read_Mapping_%j.out
#SBATCH -e Read_Mapping_%j.err

# Map the reads to the C. neoformans genome with bwa-mem and process the
# BAM files with samtools.
module load parallel/20210822
module load samtools/1.14

# Path to the mapping software. We will use bwa-mem for this step
module load bwa/0.7.17

# Do a little arithmetic. We will run mapping tasks in parallel, but keep in
# mind that mapping software is multithreaded
export BWA_THREADS="8"
N_PAR=$(echo "scale=0; ${SLURM_CPUS_PER_TASK} / ${BWA_THREADS} / 1" | bc -l)

# Path to input and output directories
export READS_DIR="/home/Cn_GWAS/Trimmed_Reads"
export BAM_DIR="/home/Cn_GWAS/BAM.Raw"
export BAM_MAPQ="/home/Cn_GWAS/BAM.MAPQ20"
export UNMAPPED="/home/Cn_GWAS/Unmapped_Reads"
export STATS_DIR="/home/Cn_GWAS/BAM_Stats"
mkdir -p "${BAM_DIR}" "${BAM_MAPQ}" "${STATS_DIR}" "${UNMAPPED}"


# Path to the bwa-mem index
export GENOME_IDX="/home/Cryptococcus_project/Reference_Genome/bwa_idx"

# Define a function that does the mapping and algnment processing. Takes the
# forward read (R1) filename as an argument
function BWA_map () {
    fwd_read="${1}"
    sname=$(basename "${fwd_read}" | sed -e 's/_1P\.fastq\.gz//g')
    rev_read="${READS_DIR}/${sname}_2P.fastq.gz"
    # Run minimap and pipe it to samtools for processing
    # "${BWA}" \
    #     -ax sr \
    #     -t "${BWA_THREADS}" \
    #     "${GENOME_IDX}" \
    #     "${fwd_read}" "${rev_read}" \
    #     | samtools view -hbu \
    #     | samtools sort -o "${BAM_DIR}/${sname}.bam" -
    # Run bwa-mem and pipe to samtools for processing
    bwa mem \
        -t "${BWA_THREADS}" \
        -k 15 \
        -r 1.0 \
        -M \
        "${GENOME_IDX}" \
        "${fwd_read}" "${rev_read}" \
        | samtools view -hbu \
        | samtools sort -o "${BAM_DIR}/${sname}.bam" -
    # We also want to isolate the unmapped reads for follow-up studies. In
    # order for 'samtools fastq' to work, we will need to sort the file by
    # read name.
    samtools view -f 4 -hbu "${BAM_DIR}/${sname}.bam" \
        | samtools sort -@ "${BWA_THREADS}" -n -o "/home/Cn_GWAS/${sname}.unmapped.qsort.bam" -
    samtools fastq \
        -1 "${UNMAPPED}/${sname}_Unmapped_R1.fastq.gz" \
        -2 "${UNMAPPED}/${sname}_Unmapped_R2.fastq.gz" \
        -s "${UNMAPPED}/${sname}_Unmapped_U.fastq.gz" \
        -0 /dev/null \
        -@ "${BWA_THREADS}" \
        "/scratch.local/${sname}.unmapped.qsort.bam"
    # Finally, we will remove mapped reads that have a mapping quality of less
    # than 20.
    samtools view \
        -hb -q 20 \
        "${BAM_DIR}/${sname}.bam" \
        > "${BAM_MAPQ}/${sname}.mapq20.bam"
}
# Export the function so parallel can call it
export -f BWA_map
# Run the function with parallel
parallel -j "${N_PAR}" --joblog "${BAM_DIR}/Read_Mapping.joblog" BWA_map ::: $(find "${READS_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*1P.fastq.gz')

# Then, index the filtered BAMs
parallel samtools index ::: $(find "${BAM_MAPQ}" -mindepth 1 -maxdepth -type f -name '*.bam')

# Next, we will collect stats on the BAM files, both raw and MAPQ-filtered
# for bam in $(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.bam')
# do
#     sname=$(basename "${bam}" | sed -e 's/\.bam//g')
#     echo "samtools stats ${bam} > ${STATS_DIR}/${sname}_Stats.txt"
# done | parallel
