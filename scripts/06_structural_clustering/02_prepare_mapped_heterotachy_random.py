#!/usr/bin/env python3
"""
02_prepare_mapped_heterotachy_random.py

Takes:
  - mapped_sites/<fam>.mapped_sites.tsv   (contains PDB column like [A:SER3])
  - obs_groups/<fam>.obs_group.tsv        (already in PDB coord space)
Creates:
  - mapped_sites_pdb/<fam>.mapped_sites_pdb.tsv  (Group := PDB)
  - rand_groups_pdb/<fam>.rand_groups_pdb.tsv    (random groups in PDB coords)

Then your contactsubgraphs script will work on rand_groups_pdb.
"""

import argparse
import csv
import subprocess
from pathlib import Path

import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapped-sites-dir", required=True)
    ap.add_argument("--obs-dir", required=True)
    ap.add_argument("--randomizer", required=True, help="Path to sged-randomize-groups.py")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--replicates", type=int, default=1000)
    ap.add_argument("--skip-existing", action="store_true")
    args = ap.parse_args()

    mapped_dir = Path(args.mapped_sites_dir).resolve()
    obs_dir = Path(args.obs_dir).resolve()
    randomizer = Path(args.randomizer).resolve()
    out_dir = Path(args.out_dir).resolve()

    pool_dir = out_dir / "mapped_sites_pdb"
    rand_dir = out_dir / "rand_groups_pdb"
    log_dir = out_dir / "logs"
    for d in (pool_dir, rand_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    log_csv = log_dir / "rebuild_rand_groups_pdb_log.csv"

    obs_files = sorted(obs_dir.glob("*.obs_group.tsv"))
    if not obs_files:
        raise SystemExit(f"No obs_group files found in {obs_dir}")

    with open(log_csv, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=[
            "family","status","mapped_in","pool_out","obs_in","rand_out","rc","message"
        ])
        w.writeheader()

        for obs in obs_files:
            fam = obs.name.replace(".obs_group.tsv", "")
            mapped_in = mapped_dir / f"{fam}.mapped_sites.tsv"
            if not mapped_in.exists():
                w.writerow({"family": fam, "status": "SKIP_NO_MAPPED",
                            "mapped_in": str(mapped_in), "pool_out": "", "obs_in": str(obs),
                            "rand_out": "", "rc": "", "message": "missing mapped_sites"})
                continue

            pool_out = pool_dir / f"{fam}.mapped_sites_pdb.tsv"
            rand_out = rand_dir / f"{fam}.rand_groups_pdb.tsv"

            if args.skip_existing and rand_out.exists():
                w.writerow({"family": fam, "status": "SKIP_EXISTS",
                            "mapped_in": str(mapped_in), "pool_out": str(pool_out), "obs_in": str(obs),
                            "rand_out": str(rand_out), "rc": "", "message": "rand exists"})
                continue

            # Build pool: Group := PDB
            df = pd.read_csv(mapped_in, sep="\t", dtype=str, keep_default_na=False)
            if "PDB" not in df.columns:
                w.writerow({"family": fam, "status": "ERROR_NO_PDB_COL",
                            "mapped_in": str(mapped_in), "pool_out": "", "obs_in": str(obs),
                            "rand_out": "", "rc": "", "message": "no PDB column"})
                continue

            df["Group"] = df["PDB"]
            df.to_csv(pool_out, sep="\t", index=False)

            # Run randomizer
            cmd = [
                "python3", str(randomizer),
                "-s", str(obs),
                "-r", str(pool_out),
                "-o", str(rand_out),
                "-n", str(args.replicates),
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            msg = (proc.stdout + "\n" + proc.stderr).strip()
            if len(msg) > 2000:
                msg = msg[-2000:]

            w.writerow({"family": fam,
                        "status": "OK" if proc.returncode == 0 else "ERROR_RANDOMIZER",
                        "mapped_in": str(mapped_in),
                        "pool_out": str(pool_out),
                        "obs_in": str(obs),
                        "rand_out": str(rand_out) if proc.returncode == 0 else "",
                        "rc": proc.returncode,
                        "message": msg})

if __name__ == "__main__":
    main()

