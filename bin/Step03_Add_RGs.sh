#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 4:00:00
#SBATCH -p agsmall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Add_RGs_%j.err
#SBATCH -o Add_RGs_%j.out

# Add read groups to the files with picard and re-index for variant calling
module load java/openjdk-17.0.2
module load parallel/20210822

PICARD="/home/Cryptococcus_project/Software/picard-tools-3.0.0/picard.jar"
BAM_DIR="/home/Cn_GWAS/BAM.MAPQ20"
BAM_OUT_DIR="/home/Cn_GWAS/BAM.RGs"

mkdir -p "${BAM_OUT_DIR}"

# for each BAM, write the read group into the file. We will just use the
# samplename for the RG.
for bam in $(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*mapq20.bam')
do
    sname=$(basename "${bam}" | sed -e 's/\.mapq20\.bam//g')
    echo "java -jar ${PICARD} \
        AddOrReplaceReadGroups \
        -I ${bam} \
        -O ${BAM_OUT_DIR}/${sname}.mapq20.RGs.bam \
        --RGLB ${sname} \
        --RGPL ILLUMINA \
        --RGPU ${sname} \
        --RGSM ${sname} \
        --RGID ${sname} \
        --CREATE_INDEX true"
done | parallel --joblog "${BAM_OUT_DIR}/Add_RGs.joblog"
