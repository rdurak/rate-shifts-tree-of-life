#!/usr/bin/env python3
"""
Create one multi-FASTA per OMA group (fingerprint), for a selected set of fingerprints.

Inputs:
  - fingerprints.txt : one fingerprint per line to process
  - oma-groups.txt   : OMA groups table (tab-separated), where:
                       parts[1] is fingerprint and parts[2:] are sequence IDs
  - oma-seqs.fa      : FASTA containing all sequences

Output:
  - For each fingerprint, writes: <out_dir>/<fingerprint>/<fingerprint>.fa

Example:
  python scripts/01_create_fasta_files.py \
    --fingerprints data/fingerprints.txt \
    --groups data/oma-groups.txt \
    --fasta data/oma-seqs.fa \
    --out-dir Groups_Min_100
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Set

from Bio import SeqIO


def load_fingerprints(path: Path) -> Set[str]:
    fps: Set[str] = set()
    with path.open("r") as f:
        for line in f:
            s = line.strip()
            if s:
                fps.add(s)
    return fps


def parse_groups(groups_file: Path, selected_fps: Set[str]) -> Dict[str, List[str]]:
    """
    Returns mapping: fingerprint -> list of sequence IDs
    """
    fp_to_ids: Dict[str, List[str]] = {}
    with groups_file.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue  # skip malformed lines
            fp = parts[1]
            if fp in selected_fps:
                fp_to_ids[fp] = parts[2:]
    return fp_to_ids


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Write multi-FASTA files per selected OMA group fingerprint."
    )
    ap.add_argument("--fingerprints", required=True, type=Path, help="fingerprints.txt")
    ap.add_argument("--groups", required=True, type=Path, help="oma-groups.txt")
    ap.add_argument("--fasta", required=True, type=Path, help="oma-seqs.fa")
    ap.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help="Base output directory (e.g., Groups_Min_100)",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing <fingerprint>.fa files if present",
    )
    args = ap.parse_args()

    # Basic checks (fail fast)
    for p in [args.fingerprints, args.groups, args.fasta]:
        if not p.exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)

    selected_fps = load_fingerprints(args.fingerprints)
    if not selected_fps:
        print("ERROR: No fingerprints loaded (file empty?).", file=sys.stderr)
        return 2

    fp_to_ids = parse_groups(args.groups, selected_fps)
    if not fp_to_ids:
        print(
            "WARNING: No matching fingerprints found in groups file. Nothing to do.",
            file=sys.stderr,
        )
        return 0

    # Index FASTA for random access by record.id (more memory-friendly than loading all sequences)
    # Note: SeqIO.index builds an on-disk index mapping IDs -> file offsets.
    fasta_index = SeqIO.index(str(args.fasta), "fasta")

    missing_total = 0
    written_groups = 0

    for fp, ids in fp_to_ids.items():
        out_group_dir = args.out_dir / fp
        out_group_dir.mkdir(parents=True, exist_ok=True)

        out_fa = out_group_dir / f"{fp}.fa"
        if out_fa.exists() and not args.overwrite:
            print(f"SKIP (exists): {out_fa}")
            continue

        missing_here = 0
        written_here = 0

        with out_fa.open("w") as out_handle:
            for seq_id in ids:
                if seq_id in fasta_index:
                    SeqIO.write(fasta_index[seq_id], out_handle, "fasta")
                    written_here += 1
                else:
                    missing_here += 1

        missing_total += missing_here
        written_groups += 1
        print(
            f"{fp}: wrote {written_here} seqs -> {out_fa}"
            + (f" (missing {missing_here})" if missing_here else "")
        )

    fasta_index.close()

    print(
        f"DONE. Groups written: {written_groups}. "
        f"Total missing sequence IDs: {missing_total}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

