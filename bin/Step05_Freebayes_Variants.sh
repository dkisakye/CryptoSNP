#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 80:00:00
#SBATCH -p agsmall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Freebayes_Variants_%j.err
#SBATCH -o Freebayes_Variants_%j.out

# Run Freebayes to call variants across the genome.

module load parallel/20210822
# Path to freebayes executable
FREEBAYES="/home/Cryptococcus_project/Software/freebayes-1.3.6/freebayes-1.3.6-linux-amd64-static"
# Directory where the final BAMs are stored
BAM_DIR="/home/Cn_GWAS/BAM.Final"
# File that contains a list of pathnames to each BAM file, for joint genotype
# calling with Freebayes. One path per line.
BAM_FOF="/home/Cryptococcus_project/GitHub_Repo/Resources/500_Strain_BAMs.fof"
# Path to the list of chromosomes to parallelize across
CTG_LIST="/home/Cryptococcus_project/GitHub_Repo/Resources/H99_Chromosome_List.txt"
# Path to the reference fasta
REF="/home/Cryptococcus_project/Reference_Genome/FungiDB-61_CneoformansH99_Genome.fasta"
# Path to output directory
OUT_DIR="/home/Cn_GWAS/Variant_Calls"

mkdir -p "${OUT_DIR}"

# Exclude these samples from variant calling. This is defined as an extended
# regex that we will use with grep
# EXCLUDE="(10WB)|(LS2)|(MO1)"

# Make the file-of-files for Freebayes input
# find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.bam' | grep -vE "${EXCLUDE}" | sort -V > "${BAM_FOF}"
find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.bam' > "${BAM_FOF}"

# Define some parameters of the Freebayes run:
# Adjust these for your data.
N_BEST_SNP_ALLELES="4"
MAX_COMPLEX_GAP="-1"
MIN_BASE_QUALITY="20"
MIN_ALLELE_QSUM="40"
MIN_COVERAGE="50"
LIMIT_COVERAGE="30000"

# Run the Freebayes, parallelizing across contigs:
for contig in $(cat "${CTG_LIST}")
do
    out_fname="${OUT_DIR}/Cryptococcus_neoformans_${contig}_Raw_Variants.vcf"
    echo "${FREEBAYES} \
        -f ${REF} \
        -L ${BAM_FOF} \
        -r ${contig} \
        -v ${out_fname} \
        --strict-vcf \
        --use-best-n-alleles ${N_BEST_SNP_ALLELES} \
        -E ${MAX_COMPLEX_GAP} \
        --haplotype-length ${MAX_COMPLEX_GAP} \
        -q ${MIN_BASE_QUALITY} \
        -R ${MIN_ALLELE_QSUM} \
        --min-coverage ${MIN_COVERAGE} \
        --limit-coverage ${LIMIT_COVERAGE} \
        --no-population-priors"
done | parallel -j "${SLURM_CPUS_PER_TASK}" --joblog "${OUT_DIR}/Freebayes_Variants.joblog"

