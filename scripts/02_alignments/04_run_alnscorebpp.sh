#!/usr/bin/env bash

# ============================================================
# Run bppalnscore for all per-orthogroup alignment directories
#
# Inputs:
#   combined_alignments/<ORTHOGROUP_ID>/AlnScores.bpp
#
# Usage:
#   ./04_run_alnscorebpp.sh
#
#   or
#
#   ./04_run_alnscorebpp.sh combined_alignments
#
# Optional:
#   BPPALNSCORE=/path/to/bppalnscore ./04_run_alnscorebpp.sh
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# input directory
# ------------------------------------------------------------

alignment_dir_base="${1:-combined_alignments}"

# ------------------------------------------------------------
# executable
# ------------------------------------------------------------

BPPALNSCORE="${BPPALNSCORE:-bppalnscore}"

# ------------------------------------------------------------
# checks
# ------------------------------------------------------------

if [ ! -d "$alignment_dir_base" ]; then
  echo "ERROR: alignment directory not found: $alignment_dir_base" >&2
  exit 1
fi

if ! command -v "$BPPALNSCORE" >/dev/null 2>&1; then
  echo "ERROR: bppalnscore executable not found: $BPPALNSCORE" >&2
  exit 1
fi

# ------------------------------------------------------------
# process directories
# ------------------------------------------------------------

processed=0
failed=0

for alignment_dir in "$alignment_dir_base"/*; do

  if [ ! -d "$alignment_dir" ]; then
    continue
  fi

  echo "Processing: $alignment_dir"

  (
    cd "$alignment_dir"

    if [ ! -f "AlnScores.bpp" ]; then
      echo "  Missing AlnScores.bpp" >&2
      exit 1
    fi

    if "$BPPALNSCORE" param=AlnScores.bpp; then

      echo "  Completed successfully"

    else

      echo "  Failed" >&2

      exit 1
    fi
  )

  status=$?

  if [ "$status" -eq 0 ]; then
    processed=$((processed + 1))
  else
    failed=$((failed + 1))
  fi
done

# ------------------------------------------------------------
# summary
# ------------------------------------------------------------

echo
echo "Finished."
echo "Successful directories: $processed"
echo "Failed directories: $failed"
