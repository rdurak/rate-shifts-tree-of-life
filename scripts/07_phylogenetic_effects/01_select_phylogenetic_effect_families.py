#!/usr/bin/env python3

"""
Select orthogroups for phylogenetic effect analyses.

Filters orthogroups based on:
    - ancestral category
    - sequence count
    - alignment length
    - prevalence of rate-shift sites

A reproducible random subset of orthogroups is then
selected for downstream phylogenetic experiments.

Outputs:
    - selected orthogroup table
    - filtering summary

Example:
    python 01_select_phylogenetic_effect_families.py \
        --input sged_summary_per_orthogroup.with_tree_aln.tsv \
        --output selected_orthogroups.tsv
"""

import pandas as pd


# ============================================================
# PARAMETERS
# ============================================================

INPUT_FILE = "sged_summary_per_orthogroup.with_tree_aln.tsv"

OUTPUT_FILE = "selected_orthogroups_for_phylo_test.tsv"

# ------------------------------------------------------------
# Filtering thresholds
# ------------------------------------------------------------

MIN_SEQS = 150
MAX_SEQS = 500

MIN_ALIGNMENT = 300

MIN_RS_PROP = 0.15
MAX_RS_PROP = 0.35

# ------------------------------------------------------------
# Deep ancestry categories
# ------------------------------------------------------------

DEEP_LCA = [
    "LUCA",
    "Bacteria",
    "Eukaryota",
]

# ------------------------------------------------------------
# Number of families to sample
# ------------------------------------------------------------

N_SAMPLE = 50

RANDOM_SEED = 42


# ============================================================
# LOAD DATA
# ============================================================

print("Reading dataset...")

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)

print(
    f"Total orthogroups in dataset: "
    f"{len(df)}"
)


# ============================================================
# APPLY FILTERS
# ============================================================

filtered = df[
    df["Last_Common_Ancestor"].isin(DEEP_LCA)
    &
    (df["num_sequences"] >= MIN_SEQS)
    &
    (df["num_sequences"] <= MAX_SEQS)
    &
    (df["num_sites"] >= MIN_ALIGNMENT)
    &
    (df["prop_rate_shift"] >= MIN_RS_PROP)
    &
    (df["prop_rate_shift"] <= MAX_RS_PROP)
].copy()

print(
    f"Orthogroups passing filters: "
    f"{len(filtered)}"
)

if len(filtered) < N_SAMPLE:

    raise ValueError(
        "Not enough candidate families "
        "after filtering."
    )


# ============================================================
# RANDOMLY SELECT FAMILIES
# ============================================================

selected = filtered.sample(
    n=N_SAMPLE,
    random_state=RANDOM_SEED,
)

selected = selected.sort_values(
    "Orthogroup"
)


# ============================================================
# SAVE OUTPUT
# ============================================================

selected.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False,
)

print("\nSelected orthogroups:\n")

print(
    selected[
        [
            "Orthogroup",
            "num_sites",
            "num_rate_shift_sites",
            "prop_rate_shift",
            "num_sequences",
            "tree_length",
            "Last_Common_Ancestor",
        ]
    ]
)

print(
    f"\nSaved selected candidates to: "
    f"{OUTPUT_FILE}"
)


# ============================================================
# SUMMARY STATISTICS
# ============================================================

print(
    "\nSummary statistics of selected families:"
)

print(
    selected.describe(
        include="all"
    )
)
