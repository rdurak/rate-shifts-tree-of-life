#!/bin/bash
#SBATCH --job-name=iqtree_array
#SBATCH --array=1-NUMFILES
#SBATCH --output=logs/iqtree_%A_%a.out
#SBATCH --error=logs/iqtree_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=5G

# Collect all alignment files
FILES=(results/*/*.filter2.fasta)

# Select the file corresponding to this array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID-1]}

# Run IQ-TREE
./iqtree2 -s "$FILE" -m LG4M
