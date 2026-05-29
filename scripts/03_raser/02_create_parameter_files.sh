#!/usr/bin/env bash

# ============================================================
# Create RASER parameter files for all orthogroups
#
# For each orthogroup directory:
#   - generates 01_raser.params
#   - generates 02_null.params
#   - updates input/output paths
#   - creates result directories
#
# Inputs:
#   raser_run/<ORTHOGROUP_ID>/
#   01_raser.params
#   02_null.params
#
# Usage:
#   ./02_create_parameter_files.sh
#
# Optional:
#   ./02_create_parameter_files.sh raser_run
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# directories
# ------------------------------------------------------------

RASER_RUN_DIR="${1:-raser_run}"

PARAM_TEMPLATE_1="01_raser.params"

PARAM_TEMPLATE_2="02_null.params"

# ------------------------------------------------------------
# checks
# ------------------------------------------------------------

if [ ! -d "$RASER_RUN_DIR" ]; then
  echo "ERROR: directory not found: $RASER_RUN_DIR" >&2
  exit 1
fi

if [ ! -f "$PARAM_TEMPLATE_1" ]; then
  echo "ERROR: missing template: $PARAM_TEMPLATE_1" >&2
  exit 1
fi

if [ ! -f "$PARAM_TEMPLATE_2" ]; then
  echo "ERROR: missing template: $PARAM_TEMPLATE_2" >&2
  exit 1
fi

# ------------------------------------------------------------
# absolute path
# ------------------------------------------------------------

ABS_BASE_PATH="$(cd "$RASER_RUN_DIR" && pwd)"

processed=0

# ------------------------------------------------------------
# process orthogroups
# ------------------------------------------------------------

for DIR in "$RASER_RUN_DIR"/*; do

  if [ ! -d "$DIR" ]; then
    continue
  fi

  PREFIX="$(basename "$DIR")"

  ABS_DIR="$ABS_BASE_PATH/$PREFIX"

  PARAM_OUT_1="$DIR/01_raser.params"

  PARAM_OUT_2="$DIR/02_null.params"

  mkdir -p "$DIR/raser_results"

  mkdir -p "$DIR/null_results"

  echo "Processing: $PREFIX"

  # ==========================================================
  # 01_raser.params
  # ==========================================================

  awk -v prefix="$PREFIX" -v path="$ABS_DIR" '

  BEGIN {
      outTreeWritten = 0;
      outTreeNamesWritten = 0;
  }

  /^_outTreeFile[ \t]/ {

      if (!outTreeWritten) {

          print "_outTreeFile\t" path "/raser_results/r_" prefix "_out_tree.ph";

          outTreeWritten = 1;
      }

      next;
  }

  /^_outTreeFileWithBranchesNames[ \t]/ {

      if (!outTreeNamesWritten) {

          print "_outTreeFileWithBranchesNames\t" path "/raser_results/r_" prefix "_out_tree.namesBS.ph";

          outTreeNamesWritten = 1;
      }

      next;
  }

  {
      gsub(/_inSeqFile.*/, "_inSeqFile\t" path "/" prefix ".filter2.fasta");

      gsub(/_inTreeFile.*/, "_inTreeFile\t" path "/" prefix ".filter2.fasta.treefile");

      gsub(/_outResFile.*/, "_outResFile\t" path "/raser_results/r_" prefix "_outresfile");

      gsub(/_outNodesResFile.*/, "_outNodesResFile\t" path "/raser_results/r_" prefix "_nodes_results_file");

      gsub(/_logFile.*/, "_logFile\t" path "/raser_results/r_" prefix "_log_file");

      print;
  }

  END {

      if (!outTreeWritten)

          print "_outTreeFile\t" path "/raser_results/r_" prefix "_out_tree.ph";

      if (!outTreeNamesWritten)

          print "_outTreeFileWithBranchesNames\t" path "/raser_results/r_" prefix "_out_tree.namesBS.ph";
  }

  ' "$PARAM_TEMPLATE_1" > "$PARAM_OUT_1"

  # ==========================================================
  # 02_null.params
  # ==========================================================

  awk -v prefix="$PREFIX" -v path="$ABS_DIR" '

  BEGIN {
      outTreeWritten = 0;
  }

  /^_outTreeFile[ \t]/ {

      if (!outTreeWritten) {

          print "_outTreeFile\t" path "/null_results/n_" prefix "_out_tree.ph";

          outTreeWritten = 1;
      }

      next;
  }

  {
      gsub(/_inSeqFile.*/, "_inSeqFile\t" path "/" prefix ".filter2.fasta");

      gsub(/_inTreeFile.*/, "_inTreeFile\t" path "/" prefix ".filter2.fasta.treefile");

      gsub(/_outResFile.*/, "_outResFile\t" path "/null_results/n_" prefix "_results_file");

      gsub(/_outNodesResFile.*/, "_outNodesResFile\t" path "/null_results/n_" prefix "_nodes_results_file");

      gsub(/_logFile.*/, "_logFile\t" path "/null_results/n_" prefix "_log_file");

      print;
  }

  END {

      if (!outTreeWritten)

          print "_outTreeFile\t" path "/null_results/n_" prefix "_out_tree.ph";
  }

  ' "$PARAM_TEMPLATE_2" > "$PARAM_OUT_2"

  processed=$((processed + 1))

done

# ------------------------------------------------------------
# summary
# ------------------------------------------------------------

echo
echo "Finished."
echo "Parameter files generated: $processed"
