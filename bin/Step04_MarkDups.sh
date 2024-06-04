#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=15gb
#SBATCH --tmp=250gb
#SBATCH -t 8:00:00
#SBATCH -p ag2tb
#SBATCH -A nielsenk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Mark_Duplicates_%j.err
#SBATCH -o Mark_Duplicates_%j.out

# Remove duplicates from BAM files.
module load java/openjdk-17.0.2
module load parallel/20210822

PICARD="/home/nielsenk/shared/Cryptococcus_project/Software/picard-tools-3.0.0/picard.jar"
BAM_DIR="/scratch.global/konox006/Cn_GWAS/BAM.RGs"
BAM_OUT_DIR="/scratch.global/konox006/Cn_GWAS/BAM.Final"
mkdir -p "${BAM_OUT_DIR}"

# for each BAM, use Picard to remove duplicates
for bam in $(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.mapq20.RGs.bam')
do
    sname=$(basename "${bam}" | sed -e 's/\.mapq20\.RGs\.bam//g')
    echo "java -jar ${PICARD} \
        MarkDuplicates \
        -I ${bam} \
        -O ${BAM_OUT_DIR}/${sname}.mapq20.RGs.MarkDups.bam \
        -M ${BAM_OUT_DIR}/${sname}_MarkDups.metrics.txt \
        --REMOVE_DUPLICATES true \
        --CREATE_INDEX true"
done | parallel --joblog "${BAM_OUT_DIR}/MarkDups.joblog"
