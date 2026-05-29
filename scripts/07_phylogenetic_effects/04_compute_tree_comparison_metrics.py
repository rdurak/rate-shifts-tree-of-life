#!/usr/bin/env python3

# ============================================================
# 04_compute_tree_comparison_metrics.py
#
# Compute:
#   1. Robinson-Foulds distances
#   2. Branch length comparisons
#   3. Branch support comparisons
#
# Input:
#   Directory containing:
#       *.filter2.fasta.contree
#       *.filter2.nors.fasta.contree
#       *.filter2.rand*.fasta.contree
#
# Outputs:
#   rf_results.tsv
#   delta_tree_length.tsv
#   delta_support.tsv
#
#   delta_tree_length_distribution.{png,pdf}
#   delta_support_distribution.{png,pdf}
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ete3 import Tree

# ============================================================
# USER SETTINGS
# ============================================================

TREE_DIR = "Gene_Trees_50"

RF_OUTFILE = "rf_results.tsv"

TREE_LENGTH_OUTFILE = "delta_tree_length.tsv"

SUPPORT_OUTFILE = "delta_support.tsv"

N_RANDOM = 100

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_total_tree_length(tree_file):

    try:
        t = Tree(tree_file)
    except Exception:
        return None

    total_length = sum(
        node.dist for node in t.traverse()
        if node.dist is not None
    )

    return total_length


def get_mean_support(tree_file):

    try:
        t = Tree(tree_file)
    except Exception:
        return None

    supports = [
        n.support for n in t.traverse()
        if not n.is_leaf() and n.support is not None
    ]

    if len(supports) == 0:
        return None

    return np.mean(supports)

# ============================================================
# OUTPUT RECORDS
# ============================================================

rf_records = []

tree_length_records = []

support_records = []

# ============================================================
# MAIN LOOP
# ============================================================

for file in os.listdir(TREE_DIR):

    # --------------------------------------------------------
    # ORIGINAL TREES ONLY
    # --------------------------------------------------------

    if not file.endswith(".filter2.fasta.contree"):
        continue

    if ".rand" in file or ".nors" in file:
        continue

    # --------------------------------------------------------
    # ORTHOGROUP
    # --------------------------------------------------------

    og = file.split(".")[0]

    base = file.replace(
        ".filter2.fasta.contree",
        ""
    )

    print(f"Processing: {og}")

    # --------------------------------------------------------
    # REFERENCE TREE
    # --------------------------------------------------------

    ref_path = os.path.join(TREE_DIR, file)

    try:
        ref_tree = Tree(ref_path)
    except Exception:
        print(f"  Failed to load reference tree")
        continue

    # --------------------------------------------------------
    # NORS TREE
    # --------------------------------------------------------

    nors_file = f"{base}.filter2.nors.fasta.contree"

    nors_path = os.path.join(TREE_DIR, nors_file)

    if not os.path.exists(nors_path):

        print(f"  Missing nors tree")

        continue

    try:
        nors_tree = Tree(nors_path)
    except Exception:
        print(f"  Failed to load nors tree")
        continue

    # ========================================================
    # RF DISTANCE
    # ========================================================

    rf, rf_max = ref_tree.robinson_foulds(
        nors_tree,
        unrooted_trees=True
    )[:2]

    nrf = rf / rf_max

    rf_records.append([
        og,
        "nors",
        "NA",
        rf,
        nrf
    ])

    # ========================================================
    # TREE LENGTH
    # ========================================================

    nors_length = get_total_tree_length(nors_path)

    # ========================================================
    # SUPPORT
    # ========================================================

    nors_support = get_mean_support(nors_path)

    # --------------------------------------------------------
    # RANDOM TREES
    # --------------------------------------------------------

    rand_lengths = []

    rand_supports = []

    for i in range(1, N_RANDOM + 1):

        rand_file = (
            f"{base}.filter2.rand{i}.fasta.contree"
        )

        rand_path = os.path.join(
            TREE_DIR,
            rand_file
        )

        if not os.path.exists(rand_path):
            continue

        try:
            rand_tree = Tree(rand_path)
        except Exception:
            continue

        # ----------------------------------------------------
        # RF
        # ----------------------------------------------------

        rf, rf_max = ref_tree.robinson_foulds(
            rand_tree,
            unrooted_trees=True
        )[:2]

        nrf = rf / rf_max

        rf_records.append([
            og,
            "rand",
            i,
            rf,
            nrf
        ])

        # ----------------------------------------------------
        # TREE LENGTH
        # ----------------------------------------------------

        length = get_total_tree_length(rand_path)

        if length is not None:
            rand_lengths.append(length)

        # ----------------------------------------------------
        # SUPPORT
        # ----------------------------------------------------

        support = get_mean_support(rand_path)

        if support is not None:
            rand_supports.append(support)

    # ========================================================
    # TREE LENGTH DELTA
    # ========================================================

    if (
        nors_length is not None and
        len(rand_lengths) > 0
    ):

        mean_rand_length = np.mean(rand_lengths)

        delta_length = (
            nors_length - mean_rand_length
        )

        tree_length_records.append([
            og,
            nors_length,
            mean_rand_length,
            delta_length
        ])

    # ========================================================
    # SUPPORT DELTA
    # ========================================================

    if (
        nors_support is not None and
        len(rand_supports) > 0
    ):

        mean_rand_support = np.mean(rand_supports)

        for s in rand_supports:

            support_records.append([
                og,
                s - mean_rand_support
            ])

# ============================================================
# SAVE RF RESULTS
# ============================================================

with open(RF_OUTFILE, "w") as out:

    out.write(
        "Orthogroup\ttype\treplicate\tRF\tnRF\n"
    )

    for r in rf_records:

        out.write(
            f"{r[0]}\t{r[1]}\t{r[2]}"
            f"\t{r[3]}\t{r[4]}\n"
        )

print(f"Wrote: {RF_OUTFILE}")

# ============================================================
# SAVE TREE LENGTH RESULTS
# ============================================================

tree_df = pd.DataFrame(
    tree_length_records,
    columns=[
        "Orthogroup",
        "nors_tree_length",
        "mean_random_tree_length",
        "delta_tree_length"
    ]
)

tree_df.to_csv(
    TREE_LENGTH_OUTFILE,
    sep="\t",
    index=False
)

print(f"Wrote: {TREE_LENGTH_OUTFILE}")

# ============================================================
# SAVE SUPPORT RESULTS
# ============================================================

support_df = pd.DataFrame(
    support_records,
    columns=[
        "Orthogroup",
        "delta_support"
    ]
)

support_df.to_csv(
    SUPPORT_OUTFILE,
    sep="\t",
    index=False
)

print(f"Wrote: {SUPPORT_OUTFILE}")

# ============================================================
# TREE LENGTH PLOT
# ============================================================

plt.figure(figsize=(10, 5))

plt.hist(
    tree_df["delta_tree_length"],
    bins=50,
    edgecolor="black"
)

plt.axvline(
    0,
    linestyle="--",
    color="black"
)

plt.xlabel(
    "Δ Tree Length (nors − random)"
)

plt.ylabel(
    "Number of Orthogroups"
)

plt.tight_layout()

plt.savefig(
    "delta_tree_length_distribution.png",
    dpi=300
)

plt.savefig(
    "delta_tree_length_distribution.pdf"
)

plt.close()

# ============================================================
# SUPPORT PLOT
# ============================================================

plt.figure(figsize=(12, 5))

groups = support_df.groupby("Orthogroup")

data = [
    g["delta_support"].values
    for _, g in groups
]

labels = [
    name for name, _ in groups
]

plt.boxplot(
    data,
    labels=labels,
    patch_artist=True,
    boxprops=dict(
        facecolor="lightgray",
        color="black"
    ),
    whiskerprops=dict(color="black"),
    capprops=dict(color="black"),
    medianprops=dict(color="black"),
)

# jitter points
for i, (name, g) in enumerate(groups, start=1):

    x = np.random.normal(
        i,
        0.04,
        size=len(g)
    )

    plt.scatter(
        x,
        g["delta_support"],
        s=10,
        alpha=0.6,
        color="black"
    )

plt.axhline(
    0,
    linestyle="--",
    color="black",
    linewidth=1
)

plt.xticks(rotation=60)

plt.xlabel("Orthogroup")

plt.ylabel(
    "ΔSupport (nors − random)"
)

plt.tight_layout()

plt.savefig(
    "delta_support_distribution.png",
    dpi=300
)

plt.savefig(
    "delta_support_distribution.pdf"
)

plt.close()

# ============================================================
# SUMMARY
# ============================================================

print("\nSummary statistics:")
print(
    tree_df["delta_tree_length"].describe()
)

print("\nFinished.")
