#!/usr/bin/env python3

"""
Create per-orthogroup alignment directories containing
corresponding MAFFT and MUSCLE5 alignments.

Inputs:
  - MAFFT alignment directory
  - MUSCLE5 alignment directory

Output:
  combined_alignments/<ORTHOGROUP_ID>/

Example:
  python 01_create_alignment_directories.py \
    --mafft-dir mafft_alignments \
    --muscle5-dir muscle5_alignments \
    --out-dir combined_alignments
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Create per-orthogroup directories and copy "
            "corresponding MAFFT and MUSCLE5 alignments."
        )
    )

    ap.add_argument(
        "--mafft-dir",
        default="mafft_alignments",
        help="Directory containing MAFFT alignments"
    )

    ap.add_argument(
        "--muscle5-dir",
        default="muscle5_alignments",
        help="Directory containing MUSCLE5 alignments"
    )

    ap.add_argument(
        "--out-dir",
        default="combined_alignments",
        help="Output directory"
    )

    args = ap.parse_args()

    mafft_dir = Path(args.mafft_dir)

    muscle5_dir = Path(args.muscle5_dir)

    output_dir = Path(args.out_dir)

    # --------------------------------------------------------
    # input checks
    # --------------------------------------------------------

    if not mafft_dir.exists():
        raise FileNotFoundError(
            f"MAFFT directory not found: {mafft_dir}"
        )

    if not muscle5_dir.exists():
        raise FileNotFoundError(
            f"MUSCLE5 directory not found: {muscle5_dir}"
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # process MAFFT alignments
    # --------------------------------------------------------

    mafft_files = sorted(
        mafft_dir.glob("*.fa")
    )

    if len(mafft_files) == 0:
        print("No MAFFT alignment files found.")
        return

    created = 0

    missing = 0

    for mafft_file in mafft_files:

        seq_id = (
            mafft_file.name
            .replace("mafft_", "")
            .replace(".fa", "")
        )

        muscle5_file = (
            muscle5_dir / f"{seq_id}.fa"
        )

        if not muscle5_file.exists():

            print(
                f"Warning: missing MUSCLE5 alignment "
                f"for {seq_id}"
            )

            missing += 1

            continue

        # ----------------------------------------------------
        # create orthogroup directory
        # ----------------------------------------------------

        alignment_dir = output_dir / seq_id

        alignment_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        # ----------------------------------------------------
        # copy files
        # ----------------------------------------------------

        shutil.copy(
            mafft_file,
            alignment_dir
        )

        shutil.copy(
            muscle5_file,
            alignment_dir
        )

        created += 1

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    print(
        f"Created alignment directories: {created}"
    )

    print(
        f"Missing MUSCLE5 alignments: {missing}"
    )


if __name__ == "__main__":
    main()
