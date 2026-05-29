#!/usr/bin/env python3

"""
Create per-orthogroup RASER run directories.

For each orthogroup:
  - creates a dedicated run directory
  - copies the corresponding filtered alignment
  - copies the corresponding phylogenetic tree
  - creates result directories for RASER analyses

Inputs:
  - trees_with_pdb/
  - alignments_with_pdb/

Outputs:
  raser_run/<ORTHOGROUP_ID>/

Example:
  python 01_create_raser_inputs.py \
    --tree-dir trees_with_pdb \
    --align-dir alignments_with_pdb \
    --out-dir raser_run
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Create per-orthogroup RASER "
            "input directories."
        )
    )

    ap.add_argument(
        "--tree-dir",
        default="trees_with_pdb",
        help="Directory containing tree files"
    )

    ap.add_argument(
        "--align-dir",
        default="alignments_with_pdb",
        help="Directory containing filtered alignments"
    )

    ap.add_argument(
        "--out-dir",
        default="raser_run",
        help="Output directory"
    )

    args = ap.parse_args()

    tree_dir = Path(args.tree_dir)

    align_dir = Path(args.align_dir)

    output_dir = Path(args.out_dir)

    # --------------------------------------------------------
    # input checks
    # --------------------------------------------------------

    if not tree_dir.exists():
        raise FileNotFoundError(
            f"Tree directory not found: "
            f"{tree_dir}"
        )

    if not align_dir.exists():
        raise FileNotFoundError(
            f"Alignment directory not found: "
            f"{align_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # identify orthogroups
    # --------------------------------------------------------

    align_files = sorted(
        align_dir.glob("*.fasta")
    )

    if len(align_files) == 0:

        print(
            "No alignment files found."
        )

        return

    prefixes = sorted(
        {
            f.name.replace(
                ".filter2.fasta",
                ""
            )
            for f in align_files
        }
    )

    created = 0

    missing_trees = 0

    # --------------------------------------------------------
    # process orthogroups
    # --------------------------------------------------------

    for prefix in prefixes:

        subdir = output_dir / prefix

        subdir.mkdir(
            parents=True,
            exist_ok=True
        )

        # ----------------------------------------------------
        # create result directories
        # ----------------------------------------------------

        (
            subdir / "raser_results"
        ).mkdir(
            exist_ok=True
        )

        (
            subdir / "null_results"
        ).mkdir(
            exist_ok=True
        )

        # ----------------------------------------------------
        # expected files
        # ----------------------------------------------------

        align_file = (
            f"{prefix}.filter2.fasta"
        )

        tree_file = (
            f"{align_file}.treefile"
        )

        align_path = (
            align_dir / align_file
        )

        tree_path = (
            tree_dir / tree_file
        )

        # ----------------------------------------------------
        # copy alignment
        # ----------------------------------------------------

        if align_path.exists():

            shutil.copy(
                align_path,
                subdir / align_file
            )

        else:

            print(
                f"Warning: alignment missing "
                f"for {prefix}"
            )

            continue

        # ----------------------------------------------------
        # copy tree
        # ----------------------------------------------------

        if tree_path.exists():

            shutil.copy(
                tree_path,
                subdir / tree_file
            )

        else:

            print(
                f"Warning: tree missing "
                f"for {prefix}"
            )

            missing_trees += 1

            continue

        created += 1

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    print("\nFinished.")

    print(
        f"RASER directories created: "
        f"{created}"
    )

    print(
        f"Missing trees: "
        f"{missing_trees}"
    )


if __name__ == "__main__":
    main()
