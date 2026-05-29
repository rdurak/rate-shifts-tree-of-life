#!/usr/bin/env python3
"""
01_prepare_mapped_heterotachy_observed.py

For each *_sged_4 file in an input directory:
  1) Create a mapped-sites pool (rows where PDB != [NA])
     -> <out>/mapped_sites/<fam>.mapped_sites.tsv
  2) Create ONE observed group per family from (heterotachy ∩ mapped) using PDB residue IDs
     where heterotachy is defined STRICTLY as Proba>0.95 == '*'
     -> <out>/obs_groups/<fam>.obs_group.tsv

Families are SKIPPED if:
  - no mapped sites
  - fewer than 2 mapped heterotachy sites

This script DOES NOT modify the original sged_4 files.
"""

from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path
from typing import List
import pandas as pd


def is_mapped_pdb(val: str) -> bool:
    if val is None:
        return False
    v = str(val).strip()
    return v != "" and v != "[NA]"


def strip_brackets(val: str) -> str:
    v = str(val).strip()
    if len(v) >= 2 and v[0] == "[" and v[-1] == "]":
        return v[1:-1].strip()
    return v


def build_obs_group_pdb(df_mapped: pd.DataFrame) -> List[str]:
    """Return list of mapped heterotachous PDB residue IDs."""
    flag = df_mapped["Proba>0.95"].fillna("").astype(str)
    het_df = df_mapped.loc[flag == "*"].copy()

    pdb_ids = []
    for v in het_df["PDB"].tolist():
        if is_mapped_pdb(v):
            pdb_ids.append(strip_brackets(v))
    return pdb_ids


def family_id_from_filename(path: Path) -> str:
    name = re.sub(r"\.(tsv|txt|csv)$", "", path.name, flags=re.IGNORECASE)
    if name.endswith("_sged_4"):
        return name[: -len("_sged_4")]
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sged4-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--pattern", default="*_sged_4*")
    args = ap.parse_args()

    sged4_dir = Path(args.sged4_dir).resolve()
    out_dir = Path(args.out_dir).resolve()

    mapped_dir = out_dir / "mapped_sites"
    obs_dir = out_dir / "obs_groups"
    log_dir = out_dir / "logs"

    for d in (mapped_dir, obs_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    log_csv = log_dir / "prep_obs_log.csv"

    files = sorted(sged4_dir.glob(args.pattern))
    if not files:
        raise SystemExit("No *_sged_4 files found")

    with open(log_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "family",
            "input_file",
            "status",
            "n_sites_total",
            "n_sites_mapped",
            "n_sites_het_mapped",
            "obs_group_file",
            "mapped_sites_file",
            "message",
        ])
        writer.writeheader()

        for f in files:
            fam = family_id_from_filename(f)
            obs_group_file = obs_dir / f"{fam}.obs_group.tsv"
            mapped_sites_file = mapped_dir / f"{fam}.mapped_sites.tsv"

            try:
                df = pd.read_csv(f, sep=args.sep, dtype=str, comment="#", keep_default_na=False)
            except Exception as e:
                writer.writerow({"family": fam, "input_file": str(f), "status": "ERROR_READ",
                                 "n_sites_total": "", "n_sites_mapped": "", "n_sites_het_mapped": "",
                                 "obs_group_file": "", "mapped_sites_file": "", "message": str(e)})
                continue

            if "PDB" not in df.columns or "Proba>0.95" not in df.columns:
                writer.writerow({"family": fam, "input_file": str(f), "status": "ERROR_COLUMNS",
                                 "n_sites_total": len(df), "n_sites_mapped": "", "n_sites_het_mapped": "",
                                 "obs_group_file": "", "mapped_sites_file": "", "message": "Missing required columns"})
                continue

            n_total = len(df)
            df_mapped = df[df["PDB"].apply(is_mapped_pdb)].copy()
            n_mapped = len(df_mapped)

            if n_mapped == 0:
                writer.writerow({"family": fam, "input_file": str(f), "status": "SKIP_NO_MAPPED",
                                 "n_sites_total": n_total, "n_sites_mapped": 0, "n_sites_het_mapped": 0,
                                 "obs_group_file": "", "mapped_sites_file": "", "message": "No mapped sites"})
                continue

            df_mapped.to_csv(mapped_sites_file, sep="\t", index=False)

            pdb_ids = build_obs_group_pdb(df_mapped)
            n_het_mapped = len(pdb_ids)

            if n_het_mapped < 2:
                writer.writerow({"family": fam, "input_file": str(f), "status": "SKIP_LT2_HET_MAPPED",
                                 "n_sites_total": n_total, "n_sites_mapped": n_mapped,
                                 "n_sites_het_mapped": n_het_mapped,
                                 "obs_group_file": "", "mapped_sites_file": str(mapped_sites_file),
                                 "message": "Less than 2 mapped heterotachous sites"})
                continue

            with open(obs_group_file, "w") as h:
                h.write("Group\n")
                h.write("[" + ";".join(pdb_ids) + "]\n")

            writer.writerow({"family": fam, "input_file": str(f), "status": "OK",
                             "n_sites_total": n_total, "n_sites_mapped": n_mapped,
                             "n_sites_het_mapped": n_het_mapped,
                             "obs_group_file": str(obs_group_file),
                             "mapped_sites_file": str(mapped_sites_file),
                             "message": ""})


if __name__ == "__main__":
    main()

