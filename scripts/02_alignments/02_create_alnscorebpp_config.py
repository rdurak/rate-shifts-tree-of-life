#!/usr/bin/env python3

"""
Create AlnScores.bpp configuration files for each
per-orthogroup alignment directory.

The script reads a template Bio++ configuration file
and generates one AlnScores.bpp file per orthogroup,
setting the DATA parameter to the orthogroup directory name.

Inputs:
  - combined_alignments/<ORTHOGROUP_ID>/
  - template AlnScores.bpp file

Outputs:
  - combined_alignments/<ORTHOGROUP_ID>/AlnScores.bpp

Example:
  python 02_create_alnscorebpp_config.py \
    --alignment-dir-base combined_alignments \
    --template config/AlnScores.bpp
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:

    ap = argparse.ArgumentParser(
        description=(
            "Create AlnScores.bpp configuration files "
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
        default="config/AlnScores.bpp",
        help="Template AlnScores.bpp file"
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

        aln_scores_template = (
            template_file.read()
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
            alignment_dir / "AlnScores.bpp"
        )

        # ----------------------------------------------------
        # set DATA parameter
        # ----------------------------------------------------

        aln_scores_content = (
            aln_scores_template.replace(
                "DATA=",
                f"DATA={orthogroup_id}"
            )
        )

        # ----------------------------------------------------
        # write config file
        # ----------------------------------------------------

        with output_path.open("w") as out_file:

            out_file.write(
                aln_scores_content
            )

        print(
            f"Created: {output_path}"
        )

        created += 1

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    print(
        f"\nAlnScores.bpp files created: "
        f"{created}"
    )


if __name__ == "__main__":
    main()
