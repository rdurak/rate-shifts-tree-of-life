#!/usr/bin/env python3
"""
04_compute_contact_subgraphs_random.py

For each family that has:
  - rand_groups/<fam>.<suffix>.tsv   (default suffix: rand_groups_pdb.tsv)

Find the structure file from:
  raser_outputs/<fam>/pdb_index
    line like: "# SGED input PDB = /path/to/6ppp.cif"

Then run sged-structure-infos.py to compute ContactSubgraphs (default 8 Å) for:
  - randomized groups file -> contact_rand/<fam>.contact_rand.tsv

This script writes outputs + a CSV log.

Example:
  python3 compute_contactsubgraphs_rand_from_index.py \
    --rand-dir /path/to/rand_groups_pdb \
    --raser-outputs /media/rdurak/DATA/OMA/06_Raser_pass/all_results/sged_run/raser_outputs \
    --structure-script /path/to/sged-structure-infos.py \
    --out-dir /path/to/contact_outputs \
    --threshold 8 \
    --pdb-format mmCif \
    --suffix .rand_groups_pdb.tsv \
    --skip-existing
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path
from typing import Optional, Tuple

PDB_LINE_RE = re.compile(r"^\s*#\s*SGED input PDB\s*=\s*(.+?)\s*$")


def extract_pdb_path(pdb_index_path: Path) -> Optional[str]:
    try:
        with pdb_index_path.open("r") as fh:
            for line in fh:
                m = PDB_LINE_RE.match(line)
                if m:
                    return m.group(1).strip()
    except Exception:
        return None
    return None


def run_structure_infos(
    structure_script: Path,
    sged_file: Path,
    pdb_file: Path,
    pdb_format: str,
    threshold: float,
    out_file: Path,
) -> Tuple[int, str]:
    cmd = [
        "python3",
        str(structure_script),
        "--sged",
        str(sged_file),
        "--pdb",
        str(pdb_file),
        "--pdb-format",
        pdb_format,
        "--measure",
        "ContactSubgraphs",
        "--threshold",
        str(threshold),
        "--output",
        str(out_file),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    tail = (proc.stdout + "\n" + proc.stderr).strip()
    if len(tail) > 2500:
        tail = tail[-2500:]
    return proc.returncode, tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rand-dir", required=True, help="Directory containing randomized groups files")
    ap.add_argument("--raser-outputs", required=True, help="Directory containing <fam>/pdb_index")
    ap.add_argument("--structure-script", required=True, help="Path to sged-structure-infos.py")
    ap.add_argument("--out-dir", required=True, help="Output directory for contact results")
    ap.add_argument("--threshold", type=float, default=8.0, help="Contact distance threshold (Å)")
    ap.add_argument("--pdb-format", default="mmCif", help="PDB format argument for structure script (mmCif or PDB)")
    ap.add_argument("--suffix", default=".rand_groups_pdb.tsv",
                    help="Filename suffix to identify families (default: .rand_groups_pdb.tsv)")
    ap.add_argument("--skip-existing", action="store_true", help="Skip families whose outputs already exist")
    args = ap.parse_args()

    rand_dir = Path(args.rand_dir).expanduser().resolve()
    raser = Path(args.raser_outputs).expanduser().resolve()
    structure_script = Path(args.structure_script).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()

    if not rand_dir.is_dir():
        raise SystemExit(f"ERROR: rand-dir not found: {rand_dir}")
    if not raser.is_dir():
        raise SystemExit(f"ERROR: raser-outputs not found: {raser}")
    if not structure_script.is_file():
        raise SystemExit(f"ERROR: structure-script not found: {structure_script}")

    contact_rand_dir = out_dir / "contact_rand"
    log_dir = out_dir / "logs"
    for d in (contact_rand_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    log_csv = log_dir / "contactsubgraphs_rand_log.csv"

    rand_files = sorted(rand_dir.glob(f"*{args.suffix}"))
    if not rand_files:
        raise SystemExit(f"ERROR: no random group files matched '*{args.suffix}' in {rand_dir}")

    with log_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "family",
                "status",
                "rand_groups",
                "pdb_index",
                "structure_path",
                "contact_rand_out",
                "rc_rand",
                "message",
            ],
        )
        writer.writeheader()

        for rand_file in rand_files:
            fam = rand_file.name[: -len(args.suffix)]

            pdb_index = raser / fam / "pdb_index"
            if not pdb_index.is_file():
                writer.writerow(
                    {
                        "family": fam,
                        "status": "SKIP_NO_PDB_INDEX",
                        "rand_groups": str(rand_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": "",
                        "contact_rand_out": "",
                        "rc_rand": "",
                        "message": "Missing pdb_index in raser_outputs/<fam>/",
                    }
                )
                continue

            pdb_path_str = extract_pdb_path(pdb_index)
            if not pdb_path_str:
                writer.writerow(
                    {
                        "family": fam,
                        "status": "SKIP_NO_PDB_PATH",
                        "rand_groups": str(rand_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": "",
                        "contact_rand_out": "",
                        "rc_rand": "",
                        "message": "Could not parse '# SGED input PDB =' line",
                    }
                )
                continue

            pdb_path = Path(pdb_path_str)
            if not pdb_path.is_file():
                writer.writerow(
                    {
                        "family": fam,
                        "status": "SKIP_PDB_MISSING",
                        "rand_groups": str(rand_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": str(pdb_path),
                        "contact_rand_out": "",
                        "rc_rand": "",
                        "message": "Structure file path from pdb_index does not exist",
                    }
                )
                continue

            out_rand = contact_rand_dir / f"{fam}.contact_rand.tsv"

            if args.skip_existing and out_rand.exists():
                writer.writerow(
                    {
                        "family": fam,
                        "status": "SKIP_EXISTS",
                        "rand_groups": str(rand_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": str(pdb_path),
                        "contact_rand_out": str(out_rand),
                        "rc_rand": "",
                        "message": "Output already exists",
                    }
                )
                continue

            rc, msg = run_structure_infos(
                structure_script=structure_script,
                sged_file=rand_file,
                pdb_file=pdb_path,
                pdb_format=args.pdb_format,
                threshold=args.threshold,
                out_file=out_rand,
            )

            writer.writerow(
                {
                    "family": fam,
                    "status": "OK" if rc == 0 else "ERROR",
                    "rand_groups": str(rand_file),
                    "pdb_index": str(pdb_index),
                    "structure_path": str(pdb_path),
                    "contact_rand_out": str(out_rand) if rc == 0 else "",
                    "rc_rand": rc,
                    "message": msg.strip(),
                }
            )


if __name__ == "__main__":
    main()
