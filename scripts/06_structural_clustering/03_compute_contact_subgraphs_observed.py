#!/usr/bin/env python3
"""
03_compute_contact_subgraphs_observed.py

For each family that has:
  - obs_groups/<fam>.obs_group.tsv

Find the structure file from:
  raser_outputs/<fam>/pdb_index
    line like: "# SGED input PDB = /path/to/6ppp.cif"

Then run sged-structure-infos.py to compute ContactSubgraphs (default 8 Å) for:
  - observed group file  -> contact_obs/<fam>.contact_obs.tsv

This script does NOT change inputs. It only writes outputs + a CSV log.

Example:
  python3 compute_contactsubgraphs_obs_from_index.py \
    --obs-dir /path/to/output_dir/obs_groups \
    --raser-outputs /media/rdurak/DATA/OMA/06_Raser_pass/all_results/sged_run/raser_outputs \
    --structure-script /path/to/sged-structure-infos.py \
    --out-dir /path/to/contact_outputs \
    --threshold 8 \
    --pdb-format mmCif \
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
    """Return structure path from '# SGED input PDB = <path>' line."""
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
    ap.add_argument("--obs-dir", required=True, help="Directory containing <fam>.obs_group.tsv")
    ap.add_argument("--raser-outputs", required=True, help="Directory containing <fam>/pdb_index")
    ap.add_argument("--structure-script", required=True, help="Path to sged-structure-infos.py")
    ap.add_argument("--out-dir", required=True, help="Output directory for contact results")
    ap.add_argument("--threshold", type=float, default=8.0, help="Contact distance threshold (Å)")
    ap.add_argument("--pdb-format", default="mmCif", help="PDB format argument for structure script (mmCif or PDB)")
    ap.add_argument("--pattern", default="*.obs_group.tsv", help="Glob pattern for observed group files")
    ap.add_argument("--skip-existing", action="store_true", help="Skip families whose outputs already exist")
    args = ap.parse_args()

    obs_dir = Path(args.obs_dir).expanduser().resolve()
    raser = Path(args.raser_outputs).expanduser().resolve()
    structure_script = Path(args.structure_script).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()

    if not obs_dir.is_dir():
        raise SystemExit(f"ERROR: obs-dir not found: {obs_dir}")
    if not raser.is_dir():
        raise SystemExit(f"ERROR: raser-outputs not found: {raser}")
    if not structure_script.is_file():
        raise SystemExit(f"ERROR: structure-script not found: {structure_script}")

    contact_obs_dir = out_dir / "contact_obs"
    log_dir = out_dir / "logs"
    for d in (contact_obs_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    log_csv = log_dir / "contactsubgraphs_obs_log.csv"

    obs_files = sorted(obs_dir.glob(args.pattern))
    if not obs_files:
        raise SystemExit(f"ERROR: no observed group files matched {args.pattern} in {obs_dir}")

    with log_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "family",
                "status",
                "obs_group",
                "pdb_index",
                "structure_path",
                "contact_obs_out",
                "rc_obs",
                "message",
            ],
        )
        writer.writeheader()

        for obs_file in obs_files:
            fam = obs_file.name.replace(".obs_group.tsv", "")

            pdb_index = raser / fam / "pdb_index"
            if not pdb_index.is_file():
                writer.writerow(
                    {
                        "family": fam,
                        "status": "SKIP_NO_PDB_INDEX",
                        "obs_group": str(obs_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": "",
                        "contact_obs_out": "",
                        "rc_obs": "",
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
                        "obs_group": str(obs_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": "",
                        "contact_obs_out": "",
                        "rc_obs": "",
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
                        "obs_group": str(obs_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": str(pdb_path),
                        "contact_obs_out": "",
                        "rc_obs": "",
                        "message": "Structure file path from pdb_index does not exist",
                    }
                )
                continue

            out_obs = contact_obs_dir / f"{fam}.contact_obs.tsv"

            if args.skip_existing and out_obs.exists():
                writer.writerow(
                    {
                        "family": fam,
                        "status": "SKIP_EXISTS",
                        "obs_group": str(obs_file),
                        "pdb_index": str(pdb_index),
                        "structure_path": str(pdb_path),
                        "contact_obs_out": str(out_obs),
                        "rc_obs": "",
                        "message": "Output already exists",
                    }
                )
                continue

            rc1, msg1 = run_structure_infos(
                structure_script=structure_script,
                sged_file=obs_file,
                pdb_file=pdb_path,
                pdb_format=args.pdb_format,
                threshold=args.threshold,
                out_file=out_obs,
            )

            writer.writerow(
                {
                    "family": fam,
                    "status": "OK" if rc1 == 0 else "ERROR",
                    "obs_group": str(obs_file),
                    "pdb_index": str(pdb_index),
                    "structure_path": str(pdb_path),
                    "contact_obs_out": str(out_obs) if rc1 == 0 else "",
                    "rc_obs": rc1,
                    "message": msg1.strip(),
                }
            )


if __name__ == "__main__":
    main()

