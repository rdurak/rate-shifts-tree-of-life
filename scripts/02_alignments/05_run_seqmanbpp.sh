#!/usr/bin/env bash

# ============================================================
# Run bppseqman for all per-orthogroup alignment directories
#
# Inputs:
#   combined_alignments/<ORTHOGROUP_ID>/SeqMan.bpp
#
# Usage:
#   ./05_run_seqmanbpp.sh
#
#   or
#
#   ./05_run_seqmanbpp.sh combined_alignments
#
# Optional:
#   BPPSEQMAN=/path/to/bppseqman ./05_run_seqmanbpp.sh
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# input directory
# ------------------------------------------------------------

alignment_dir_base="${1:-combined_alignments}"

# ------------------------------------------------------------
# executable
# ------------------------------------------------------------

BPPSEQMAN="${BPPSEQMAN:-bppseqman}"

# ------------------------------------------------------------
# checks
# ------------------------------------------------------------

if [ ! -d "$alignment_dir_base" ]; then
  echo "ERROR: alignment directory not found: $alignment_dir_base" >&2
  exit 1
fi

if ! command -v "$BPPSEQMAN" >/dev/null 2>&1; then
  echo "ERROR: bppseqman executable not found: $BPPSEQMAN" >&2
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

    if [ ! -f "SeqMan.bpp" ]; then
      echo "  Missing SeqMan.bpp" >&2
      exit 1
    fi

    if "$BPPSEQMAN" param=SeqMan.bpp; then

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
