#!/bin/bash
# Silly little script to submit the data processing pipeline jobs

STEP_01=$(sbatch --parsable Step01_Reads_Trimming.sh)
STEP_02=$(sbatch --dependency=afterok:${STEP_01} --parsable Step02_Read_Mapping.sh)
STEP_03=$(sbatch --dependency=afterok:${STEP_02} --parsable Step03_Add_RGs.sh)
STEP_04=$(sbatch --dependency=afterok:${STEP_03} --parsable Step04_MarkDups.sh)
STEP_05=$(sbatch --dependency=afterok:${STEP_04} --parsable Step05_Freebayes_Variants.sh)

echo "Step 01: Read trimming is job ID ${STEP_01}"
echo "Step 02: Read mapping is job ID ${STEP_02}"
echo "Step 03: Add RGs is job ID ${STEP_03}"
echo "Step 04: Mark dups is job ID ${STEP_04}"
echo "Step 05: Freebayes is job ID ${STEP_05}"
