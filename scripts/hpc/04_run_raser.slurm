#!/bin/bash

#SBATCH --job-name=raser_array
#SBATCH --array=1-NUMTASKS
#SBATCH --output=logs/raser_%A_%a.out
#SBATCH --error=logs/raser_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem=32G

# Directory containing per-gene RASER runs
RASER_RUN_DIR="raser_runs"

# Path to RASER executable
RASER_EXEC="./raser"

# Number of gene-family directories processed per array task
JOBS_PER_TASK=10

# Create logs directory
mkdir -p logs

# Collect subdirectories
mapfile -t SUBDIRS < <(
    find "$RASER_RUN_DIR" -mindepth 1 -maxdepth 1 -type d | sort
)

# Calculate directory range for this SLURM task
START_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) * JOBS_PER_TASK ))
END_INDEX=$(( START_INDEX + JOBS_PER_TASK - 1 ))

# Process assigned directories
for ((i=START_INDEX; i<=END_INDEX && i<${#SUBDIRS[@]}; i++)); do

    DIR="${SUBDIRS[$i]}"

    echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} processing: ${DIR}"

    cd "$DIR" || {
        echo "Failed to enter ${DIR}"
        continue
    }

    # Run RASER analyses
    "$RASER_EXEC" 01_raser.params > raser_01.log 2>&1
    "$RASER_EXEC" 02_null.params > raser_02.log 2>&1

    cd - > /dev/null
done

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} completed."
