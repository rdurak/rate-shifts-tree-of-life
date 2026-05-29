#!/usr/bin/env python3

"""
Create SeqMan.bpp configuration files for each
per-orthogroup alignment directory.

The script reads a template Bio++ configuration file
and generates one SeqMan.bpp file per orthogroup,
setting the DATA parameter to the orthogroup directory name.

Inputs:
  - combined_alignments/<ORTHOGROUP_ID>/
  - template SeqMan.bpp file

Outputs:
  - combined_alignments/<ORTHOGROUP_ID>/SeqMan.bpp

Example:
  python 03_create_seqmanbpp_config.py \
    --alignment-dir-base combined_alignments \
    --template config/SeqMan.bpp
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Create SeqMan.bpp configuration files "
            "for per-orthogroup alignment directories."
        )
    )

    ap.add_argument(
        "--alignment-dir-base",
        default="combined_alignments",
        help="Directory containing per-orthogroup folders"
    )

    ap.add_argument(
        "--template",
        default="config/SeqMan.bpp",
        help="Template SeqMan.bpp file"
    )

    args = ap.parse_args()

    alignment_dir_base = Path(
        args.alignment_dir_base
    )

    template_path = Path(args.template)

    # --------------------------------------------------------
    # input checks
    # --------------------------------------------------------

    if not alignment_dir_base.exists():
        raise FileNotFoundError(
            f"Alignment directory not found: "
            f"{alignment_dir_base}"
        )

    if not template_path.exists():
        raise FileNotFoundError(
            f"Template file not found: "
            f"{template_path}"
        )

    # --------------------------------------------------------
    # load template
    # --------------------------------------------------------

    with template_path.open("r") as template_file:

        seqman_template = (
            template_file.readlines()
        )

    if len(seqman_template) == 0:

        raise ValueError(
            "Template file is empty."
        )

    # --------------------------------------------------------
    # process orthogroup directories
    # --------------------------------------------------------

    alignment_dirs = sorted(
        [
            d for d in alignment_dir_base.iterdir()
            if d.is_dir()
        ]
    )

    if len(alignment_dirs) == 0:

        print(
            "No orthogroup directories found."
        )

        return

    created = 0

    for alignment_dir in alignment_dirs:

        orthogroup_id = alignment_dir.name

        output_path = (
            alignment_dir / "SeqMan.bpp"
        )

        # ----------------------------------------------------
        # modify template
        # ----------------------------------------------------

        tmp = list(seqman_template)

        tmp[0] = (
            f"DATA={orthogroup_id}\n"
        )

        # ----------------------------------------------------
        # write config file
        # ----------------------------------------------------

        with output_path.open("w") as out_file:

            out_file.writelines(tmp)

        print(
            f"Created: {output_path}"
        )

        created += 1

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    print(
        f"\nSeqMan.bpp files created: "
        f"{created}"
    )


if __name__ == "__main__":
    main()
