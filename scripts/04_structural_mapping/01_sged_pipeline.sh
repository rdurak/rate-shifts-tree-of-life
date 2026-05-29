#!/usr/bin/env bash

# ============================================================
# SGED structural mapping pipeline
#
# For a single orthogroup:
#   1. convert RASER output to SGED format
#   2. translate filtered alignment coordinates
#   3. select experimentally resolved structures
#   4. create sequence-structure indices
#   5. map alignment coordinates to PDB coordinates
#   6. compute structural descriptors
#
# Usage:
#   ./sged_pipeline.sh ORTHOGROUP_ID
#
# Example:
#   ./sged_pipeline.sh DHVAKCW
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# input orthogroup
# ------------------------------------------------------------

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 ORTHOGROUP_ID" >&2
  exit 1
fi

OG_ID="$1"

# ------------------------------------------------------------
# directories
# ------------------------------------------------------------

RASER_DIR="raser_results"

FILTERED_ALIGNMENTS="filtered_alignments"

UNFILTERED_ALIGNMENTS="muscle5_alignments"

COORD_INDEX_DIR="coordinate_indices"

STRUCTURE_DIR="structures"

OUTPUT_DIR="sged_outputs/${OG_ID}"

mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------
# input files
# ------------------------------------------------------------

RASER_FILE="${RASER_DIR}/${OG_ID}_outresfile"

FILTERED_ALIGNMENT="${FILTERED_ALIGNMENTS}/${OG_ID}.filter2.fasta"

UNFILTERED_ALIGNMENT="${UNFILTERED_ALIGNMENTS}/${OG_ID}.fa"

COORD_INDEX="${COORD_INDEX_DIR}/${OG_ID}.cs_index.txt"

STRUCTURE_PATH="${STRUCTURE_DIR}/${OG_ID}"

# ------------------------------------------------------------
# checks
# ------------------------------------------------------------

required_files=(
  "$RASER_FILE"
  "$FILTERED_ALIGNMENT"
  "$UNFILTERED_ALIGNMENT"
  "$COORD_INDEX"
)

for file in "${required_files[@]}"; do

  if [ ! -f "$file" ]; then
    echo "ERROR: missing file: $file" >&2
    exit 1
  fi
done

if [ ! -d "$STRUCTURE_PATH" ]; then
  echo "ERROR: missing structure directory: $STRUCTURE_PATH" >&2
  exit 1
fi

# ------------------------------------------------------------
# SGED executables
# ------------------------------------------------------------

SGED_RASER2SGED="${SGED_RASER2SGED:-sged-raser2sged.py}"

SGED_TRANSLATE="${SGED_TRANSLATE:-sged-translate-coords.py}"

SGED_SELECT="${SGED_SELECT:-sged-select-structure.py}"

SGED_INDEX="${SGED_INDEX:-sged-create-structure-index.py}"

SGED_INFO="${SGED_INFO:-sged-structure-infos.py}"

# ------------------------------------------------------------
# Step 1
# ------------------------------------------------------------

echo "[$(date)] Step 1: RASER -> SGED"

python3 "$SGED_RASER2SGED" \
  --raser "$RASER_FILE" \
  --alignment "$FILTERED_ALIGNMENT" \
  --alignment-format fasta \
  --output "${OUTPUT_DIR}/sged_step1"

# ------------------------------------------------------------
# Step 2
# ------------------------------------------------------------

echo "[$(date)] Step 2: Coordinate translation"

python3 "$SGED_TRANSLATE" \
  --sged "${OUTPUT_DIR}/sged_step1" \
  --index "$COORD_INDEX" \
  --reverse \
  --output "${OUTPUT_DIR}/sged_step2"

# ------------------------------------------------------------
# Step 3
# ------------------------------------------------------------

echo "[$(date)] Step 3: Structure selection"

python3 "$SGED_SELECT" \
  --structures "${STRUCTURE_PATH}"/* \
  --format mmCif \
  --resolution-threshold 5 \
  --csv-report "${OUTPUT_DIR}/pdb_summary.csv"

# ------------------------------------------------------------
# retrieve selected structure
# ------------------------------------------------------------

BEST_PDB_LIST="best_pdb_list.txt"

if [ ! -f "$BEST_PDB_LIST" ]; then
  echo "ERROR: best_pdb_list.txt not generated" >&2
  exit 1
fi

BEST_PDB="$(head -n 1 "$BEST_PDB_LIST")"

if [ -z "$BEST_PDB" ]; then
  echo "ERROR: no structure selected for ${OG_ID}" >&2
  exit 1
fi

echo "Selected structure: $BEST_PDB"

# ------------------------------------------------------------
# Step 4
# ------------------------------------------------------------

echo "[$(date)] Step 4: Structure indexing"

python3 "$SGED_INDEX" \
  --pdb-list "$BEST_PDB_LIST" \
  --pdb-format mmCif \
  --alignment "$UNFILTERED_ALIGNMENT" \
  --alignment-format fasta \
  --output "${OUTPUT_DIR}/pdb_index" \
  --gap-open=-2 \
  --gap-extend=-0.5

# ------------------------------------------------------------
# Step 5
# ------------------------------------------------------------

echo "[$(date)] Step 5: Alignment -> structure coordinates"

python3 "$SGED_TRANSLATE" \
  --sged "${OUTPUT_DIR}/sged_step2" \
  --index "${OUTPUT_DIR}/pdb_index" \
  --name PDB \
  --output "${OUTPUT_DIR}/sged_step3"

# ------------------------------------------------------------
# Step 6
# ------------------------------------------------------------

echo "[$(date)] Step 6: Structural descriptors"

python3 "$SGED_INFO" \
  --sged "${OUTPUT_DIR}/sged_step3" \
  --pdb "$BEST_PDB" \
  --pdb-format mmCif \
  --group PDB \
  --measure Chain \
  --measure DSSP \
  --measure ResidueDepth \
  --measure SecondaryStructureLabel \
  --output "${OUTPUT_DIR}/sged_final"

# ------------------------------------------------------------
# done
# ------------------------------------------------------------

echo
echo "[$(date)] Pipeline completed for ${OG_ID}"
