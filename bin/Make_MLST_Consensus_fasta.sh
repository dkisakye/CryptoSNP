#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH -t 6:00:00
#SBATCH -A nielsenk
#SBATCH -p amdsmall
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -e Make_MLST_Consensus_fasta_%j.err
#SBATCH -o Make_MLST_Consensus_fasta_%j.out

# Make a consensus fasta from the MLST loci for the ~500 crypto isolates. This
# depends on identification of the MLST regions for the H99 reference genome
# and a modern (>1.15) version of Samtools.

# Define path to BAM input directory
BAM_DIR="/home/nielsenk/shared/Cryptococcus_project/results/Combined_Replicates/BAM.Final"
# Define MLST regions. Note that samtools regions are 1-based!! The start
# coordinate is shifted "forward" by 1 relative to the BED (which is 0-based)
#   Note, based on dicussions with Carolina Firacative, we adjusted the SOD1
#   start to one base upstream from what BLAST searches would suggest.
CAP59="CP003820.1:1879601-1880160"
GPD1="CP003826.1:523911-524454"
IGS1="CP003821.1:280331-281053"
LAC1="CP003827.1:1022213-1022682"
PLB1="CP003831.1:260896-261428"
SOD1="CP003824.1:1416715-1417250"
URA5="CP003827.1:315161-315797"
# Define path to samtools. Needs to be version >1.15 for the 'consensus' command
SAMTOOLS="/home/nielsenk/shared/Cryptococcus_project/Software/samtools-1.17/samtools"
# Define output filename and path
OUT_FILE="/home/nielsenk/shared/Cryptococcus_project/results/MLST/500_Strain_Consensus.fasta"

# Print a header into the file
#echo "Strain,CAP59,GPD1,IGS1,LAC1,PLB1,SOD1,URA5" > "${OUT_FILE}"
# Run through the BAM files:
for bam in $(find "${BAM_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.mapq20.RGs.MarkDups.bam' | sort -V)
do
    # Get the strain name from the filename
    str_name=$(basename "${bam}" | cut -f 1 -d '.')
    # Get the consensus sequences for each locus
    # -a: output all bases in region (do not skip those with no coverage)
    # -m "simple": use a simple consensus algorithm rather than a Bayesian one
    # -c 0.5: require at 50% of reads to agree (majority rules consensus)
    # -H 0.5: require at least 50% of reads to show alternate allele to
    #         consider a heterozygous position
    # These will still have N characters for low coverage, but those can be
    # excluded from the MLST typing, perhaps.
    str_cap59=$("${SAMTOOLS}" consensus -r "${CAP59}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    str_gpd1=$("${SAMTOOLS}" consensus -r "${GPD1}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    str_igs1=$("${SAMTOOLS}" consensus -r "${IGS1}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    str_lac1=$("${SAMTOOLS}" consensus -r "${LAC1}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    str_plb1=$("${SAMTOOLS}" consensus -r "${PLB1}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    str_sod1=$("${SAMTOOLS}" consensus -r "${SOD1}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    str_ura5=$("${SAMTOOLS}" consensus -r "${URA5}" -m "simple" -c 0.5 -H 0.5 -a "${bam}" | tail -n +2 | tr -d '\n')
    
# Then print it out into a fasta

   {
    echo ">${str_name}" >>  "${OUT_FILE}" 
    echo "${str_cap59}${str_gpd1}${str_igs1}${str_lac1}${str_plb1}${str_sod1}${str_ura5}"  >> "${OUT_FILE}"
   }

done
