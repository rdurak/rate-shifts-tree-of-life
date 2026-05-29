#!/usr/bin/env python3

"""
Generate alignment variants for phylogenetic sensitivity analyses.

For each orthogroup directory:
    - identify rate-shifting sites from SGED annotations
    - generate an alignment with all rate-shifting sites removed
    - generate 100 random-removal alignments with the same number of removed sites

Outputs per family:
    <family>.filter2.nors.fasta
    <family>.filter2.rand1.fasta
    ...
    <family>.filter2.rand100.fasta

Usage:
    python3 03_generate_modified_alignments.py
"""

import os
import random

from Bio import SeqIO


# ============================================================
# PARAMETERS
# ============================================================

BASE_DIR = "phylo_test_dataset"

N_RANDOM = 100

RANDOM_SEED = 42

random.seed(RANDOM_SEED)


# ============================================================
# PROCESS FAMILIES
# ============================================================

for family in os.listdir(BASE_DIR):

    family_dir = os.path.join(
        BASE_DIR,
        family,
    )

    if not os.path.isdir(family_dir):
        continue

    print(f"\nProcessing {family}")

    # ========================================================
    # INPUT FILES
    # ========================================================

    alignment_file = os.path.join(
        family_dir,
        f"{family}.filter2.fasta",
    )

    sged_file = os.path.join(
        family_dir,
        "sged_1",
    )

    if not os.path.exists(alignment_file):

        print("  Alignment missing")

        continue

    if not os.path.exists(sged_file):

        print("  SGED file missing")

        continue

    # ========================================================
    # LOAD ALIGNMENT
    # ========================================================

    records = list(
        SeqIO.parse(
            alignment_file,
            "fasta",
        )
    )

    if len(records) == 0:

        print("  Empty alignment")

        continue

    alignment_length = len(records[0].seq)

    print(
        f"  Alignment length: "
        f"{alignment_length}"
    )

    print(
        f"  Number of sequences: "
        f"{len(records)}"
    )

    # ========================================================
    # PARSE RATE-SHIFTING SITES
    # ========================================================

    rs_sites = []

    with open(sged_file) as handle:

        next(handle)

        for line in handle:

            cols = line.strip().split()

            if len(cols) < 4:
                continue

            # Convert 1-based SGED positions
            # to 0-based Python indices
            site = int(
                cols[0].strip("[]")
            ) - 1

            if cols[3] == "*":

                rs_sites.append(site)

    rs_sites = sorted(set(rs_sites))

    n_rs_sites = len(rs_sites)

    print(
        f"  Rate-shifting sites detected: "
        f"{n_rs_sites}"
    )

    if n_rs_sites == 0:

        print(
            "  WARNING: no rate-shifting "
            "sites detected"
        )

        continue

    # ========================================================
    # REMOVE ALIGNMENT COLUMNS
    # ========================================================

    def remove_sites(
        records,
        remove_set,
    ):

        keep = [
            i for i in range(alignment_length)
            if i not in remove_set
        ]

        new_records = []

        for record in records:

            seq = record.seq

            new_seq = "".join(
                seq[i]
                for i in keep
            )

            new_record = record[:]

            new_record.seq = type(seq)(new_seq)

            new_records.append(new_record)

        return new_records

    # ========================================================
    # REMOVE RATE-SHIFTING SITES
    # ========================================================

    nors_records = remove_sites(
        records,
        set(rs_sites),
    )

    nors_file = os.path.join(
        family_dir,
        f"{family}.filter2.nors.fasta",
    )

    SeqIO.write(
        nors_records,
        nors_file,
        "fasta",
    )

    print(
        "  Rate-shifting-site-removed "
        "alignment written"
    )

    # ========================================================
    # RANDOM SITE REMOVALS
    # ========================================================

    for i in range(1, N_RANDOM + 1):

        rand_sites = set(
            random.sample(
                range(alignment_length),
                n_rs_sites,
            )
        )

        rand_records = remove_sites(
            records,
            rand_sites,
        )

        rand_file = os.path.join(
            family_dir,
            f"{family}.filter2.rand{i}.fasta",
        )

        SeqIO.write(
            rand_records,
            rand_file,
            "fasta",
        )

    print(
        f"  Random alignments written "
        f"({N_RANDOM})"
    )


# ============================================================
# DONE
# ============================================================

print(
    "\nAlignment variant generation finished."
)
