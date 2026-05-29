#!/usr/bin/env python3

"""
Prepare orthogroup datasets for downstream phylogenetic analyses.

For each orthogroup listed in the input table:
    - create a family-specific output directory
    - copy filtered alignment files
    - copy SGED files

Outputs:
    <OUT_DIR>/<Orthogroup>/

Usage:
    python3 prepare_phylo_test_dataset.py
"""

import os
import shutil

import pandas as pd


# ============================================================
# INPUT FILE
# ============================================================

SELECTED_FAMILIES = "selected_orthogroups_for_phylo_test.tsv"


# ============================================================
# INPUT DIRECTORIES
# ============================================================

ALIGNMENT_DIR = "filtered_alignments"

SGED_DIR = "sged_outputs"


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = "phylo_test_dataset"


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(
    SELECTED_FAMILIES,
    sep="\t",
)

print(
    f"{len(df)} orthogroups loaded"
)


# ============================================================
# PROCESS FAMILIES
# ============================================================

for orthogroup in df["Orthogroup"]:

    print(
        f"\nProcessing {orthogroup}"
    )

    # --------------------------------------------------------
    # Create family-specific directory
    # --------------------------------------------------------

    family_dir = os.path.join(
        OUTPUT_DIR,
        orthogroup,
    )

    os.makedirs(
        family_dir,
        exist_ok=True,
    )

    # ========================================================
    # FILTERED ALIGNMENT
    # ========================================================

    alignment_src = os.path.join(
        ALIGNMENT_DIR,
        f"{orthogroup}.filter2.fasta",
    )

    alignment_dst = os.path.join(
        family_dir,
        f"{orthogroup}.filter2.fasta",
    )

    if os.path.exists(alignment_src):

        shutil.copy2(
            alignment_src,
            alignment_dst,
        )

        print("  Alignment copied")

    else:

        print("  Alignment missing")

    # ========================================================
    # SGED FILE
    # ========================================================

    sged_src = os.path.join(
        SGED_DIR,
        orthogroup,
        "sged_1",
    )

    sged_dst = os.path.join(
        family_dir,
        "sged_1",
    )

    if os.path.exists(sged_src):

        shutil.copy2(
            sged_src,
            sged_dst,
        )

        print("  SGED file copied")

    else:

        print("  SGED file missing")


# ============================================================
# DONE
# ============================================================

print(
    "\nDataset preparation finished."
)
