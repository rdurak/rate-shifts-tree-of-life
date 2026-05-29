#!/bin/bash
#SBATCH --job-name=mafft_array
#SBATCH --array=1-NUMFILES
#SBATCH --output=logs/mafft_%A_%a.out
#SBATCH --error=logs/mafft_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=4G

mkdir -p mafft_alignments
mkdir -p logs

FILES=(concatenated_faas/*.fa)

FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}

BASENAME=$(basename "$FILE" .fa)

mafft --auto "$FILE" > mafft_alignments/${BASENAME}.mafft.fa
