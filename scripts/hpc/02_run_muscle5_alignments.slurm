#!/bin/bash
#SBATCH --job-name=muscle5_array
#SBATCH --array=1-NUMFILES
#SBATCH --output=logs/muscle5_%A_%a.out
#SBATCH --error=logs/muscle5_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=4G

mkdir -p muscle5_alignments
mkdir -p logs

FILES=(concatenated_faas/*.fa)

FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}

BASENAME=$(basename "$FILE" .fa)

muscle -super5 "$FILE" \
  -output muscle5_alignments/${BASENAME}.muscle5.fa
