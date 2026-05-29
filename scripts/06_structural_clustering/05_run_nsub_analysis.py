#!/usr/bin/env python3
"""
05_run_nsub_analysis.py

Reads:
  Observed contact subgraphs:
    /media/rdurak/DATA/OMA/11_contactmap_rs/contact/contactsubgraphs_observed/contact_obs/*.contact_obs.tsv
  Random (1000 replicates per family) contact subgraphs:
    /media/rdurak/DATA/OMA/11_contactmap_rs/contact/contactsubgraphs_random/contact_rand/*.contact_rand.tsv

Computes per family:
  - m = NbMapped (from observed; should match random)
  - k_obs = NbSubgraphs (observed)
  - Nsub_obs = (k_obs - 1) / (m - 1)
  - For each replicate i:
      Nsub_rand_i = (k_i - 1) / (m - 1)
  - One-sided permutation p-value for clustering:
      p = (N + 1) / (R + 1)
    where N = count(Nsub_rand_i <= Nsub_obs), R = number of valid replicates.

Outputs:
  - master TSV with one row per family, plus summary stats for random distribution.

Usage:
  python3 make_master_nsub_pvalues.py \
    --obs-dir /media/rdurak/DATA/OMA/11_contactmap_rs/contact/contactsubgraphs_observed/contact_obs \
    --rand-dir /media/rdurak/DATA/OMA/11_contactmap_rs/contact/contactsubgraphs_random/contact_rand \
    --out /media/rdurak/DATA/OMA/11_contactmap_rs/contact/master_nsub_pvalues.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def fam_from_obs_filename(p: Path) -> str:
    # <fam>.contact_obs.tsv
    return re.sub(r"\.contact_obs\.tsv$", "", p.name)


def fam_from_rand_filename(p: Path) -> str:
    # <fam>.contact_rand.tsv
    return re.sub(r"\.contact_rand\.tsv$", "", p.name)


def compute_nsub(k: float, m: float) -> float:
    # Nsub = (k-1)/(m-1), defined only for m>=2 and k>=1
    return (k - 1.0) / (m - 1.0)


def first_valid_pair(df: pd.DataFrame, col_k: str, col_m: str):
    k = pd.to_numeric(df[col_k], errors="coerce")
    m = pd.to_numeric(df[col_m], errors="coerce")
    for i in range(len(df)):
        if np.isfinite(k.iloc[i]) and np.isfinite(m.iloc[i]):
            return float(k.iloc[i]), float(m.iloc[i])
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--obs-dir",
        default="/media/rdurak/DATA/OMA/11_contactmap_rs/contact/contactsubgraphs_observed/contact_obs/",
    )
    ap.add_argument(
        "--rand-dir",
        default="/media/rdurak/DATA/OMA/11_contactmap_rs/contact/contactsubgraphs_random/contact_rand/",
    )
    ap.add_argument(
        "--out",
        default="master_nsub_pvalues.tsv",
        help="Output TSV path",
    )
    ap.add_argument(
        "--min-valid-reps",
        type=int,
        default=100,
        help="Minimum number of valid random replicates required to compute p-value",
    )
    args = ap.parse_args()

    obs_dir = Path(args.obs_dir).expanduser().resolve()
    rand_dir = Path(args.rand_dir).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()

    if not obs_dir.is_dir():
        raise SystemExit(f"ERROR: obs-dir not found: {obs_dir}")
    if not rand_dir.is_dir():
        raise SystemExit(f"ERROR: rand-dir not found: {rand_dir}")

    obs_files = sorted(obs_dir.glob("*.contact_obs.tsv"))
    if not obs_files:
        raise SystemExit(f"ERROR: no observed files found in {obs_dir}")

    # Build quick lookup for random files
    rand_files = {fam_from_rand_filename(p): p for p in rand_dir.glob("*.contact_rand.tsv")}

    rows = []

    for obs_file in obs_files:
        fam = fam_from_obs_filename(obs_file)
        rand_file = rand_files.get(fam)

        base_row = {
            "family": fam,
            "obs_file": str(obs_file),
            "rand_file": str(rand_file) if rand_file else "",
        }

        # Read obs
        try:
            obs_df = pd.read_csv(obs_file, sep="\t", dtype=str, keep_default_na=False)
        except Exception as e:
            rows.append({**base_row, "status": "ERROR_OBS_READ", "message": str(e)})
            continue

        need_cols_obs = {"ContactSubgraphs.NbSubgraphs", "ContactSubgraphs.NbMapped"}
        if not need_cols_obs.issubset(obs_df.columns):
            rows.append({**base_row, "status": "ERROR_OBS_COLS",
                         "message": f"Missing columns: {sorted(need_cols_obs - set(obs_df.columns))}"})
            continue

        k_obs, m_obs = first_valid_pair(obs_df, "ContactSubgraphs.NbSubgraphs", "ContactSubgraphs.NbMapped")
        if k_obs is None or m_obs is None:
            rows.append({**base_row, "status": "SKIP_OBS_NAN", "message": "Obs k/m not finite"})
            continue

        if m_obs < 2:
            rows.append({**base_row, "status": "SKIP_OBS_M_LT2", "NbMapped": m_obs,
                         "message": "NbMapped < 2"})
            continue
        if k_obs < 1:
            rows.append({**base_row, "status": "SKIP_OBS_K_LT1", "NbSubgraphs": k_obs,
                         "message": "NbSubgraphs < 1"})
            continue

        nsub_obs = compute_nsub(k_obs, m_obs)

        # Need random file for p-values
        if rand_file is None or not Path(rand_file).is_file():
            rows.append({
                **base_row,
                "status": "SKIP_NO_RAND_FILE",
                "NbMapped": m_obs,
                "NbSubgraphs_obs": k_obs,
                "Nsub_obs": nsub_obs,
                "message": "Missing random contact file",
            })
            continue

        # Read random
        try:
            rand_df = pd.read_csv(rand_file, sep="\t", dtype=str, keep_default_na=False)
        except Exception as e:
            rows.append({**base_row, "status": "ERROR_RAND_READ",
                         "NbMapped": m_obs, "NbSubgraphs_obs": k_obs, "Nsub_obs": nsub_obs,
                         "message": str(e)})
            continue

        need_cols_rand = {"Replicate", "ContactSubgraphs.NbSubgraphs", "ContactSubgraphs.NbMapped"}
        if not need_cols_rand.issubset(rand_df.columns):
            rows.append({**base_row, "status": "ERROR_RAND_COLS",
                         "NbMapped": m_obs, "NbSubgraphs_obs": k_obs, "Nsub_obs": nsub_obs,
                         "message": f"Missing columns: {sorted(need_cols_rand - set(rand_df.columns))}"})
            continue

        k_rand = pd.to_numeric(rand_df["ContactSubgraphs.NbSubgraphs"], errors="coerce")
        m_rand = pd.to_numeric(rand_df["ContactSubgraphs.NbMapped"], errors="coerce")

        # Keep rows where k and m are finite and m>=2, k>=1
        valid = np.isfinite(k_rand) & np.isfinite(m_rand) & (m_rand >= 2) & (k_rand >= 1)

        if valid.sum() == 0:
            rows.append({**base_row, "status": "SKIP_NO_VALID_REPS",
                         "NbMapped": m_obs, "NbSubgraphs_obs": k_obs, "Nsub_obs": nsub_obs,
                         "message": "No valid random replicates"})
            continue

        # Prefer using observed m as denominator (should match by construction).
        # If some reps have different m, that's a warning but still computable.
        m_used = m_obs

        nsub_rand = compute_nsub(k_rand[valid].to_numpy(dtype=float), m_used)

        R = int(len(nsub_rand))
        if R < args.min_valid_reps:
            rows.append({**base_row, "status": "SKIP_TOO_FEW_REPS",
                         "NbMapped": m_obs, "NbSubgraphs_obs": k_obs, "Nsub_obs": nsub_obs,
                         "R_valid": R,
                         "message": f"Only {R} valid reps (< {args.min_valid_reps})"})
            continue

        # One-sided p-value: clustering => smaller Nsub than random
        N = int(np.sum(nsub_rand <= nsub_obs))
        p = (N + 1.0) / (R + 1.0)

        row = {
            **base_row,
            "status": "OK",
            "NbMapped": float(m_obs),
            "NbSubgraphs_obs": float(k_obs),
            "Nsub_obs": float(nsub_obs),
            "R_valid": R,
            "N_leq": N,
            "p_value": float(p),
            "Nsub_rand_mean": float(np.mean(nsub_rand)),
            "Nsub_rand_median": float(np.median(nsub_rand)),
            "Nsub_rand_sd": float(np.std(nsub_rand, ddof=1)) if R > 1 else np.nan,
            "Nsub_rand_q05": float(np.quantile(nsub_rand, 0.05)),
            "Nsub_rand_q95": float(np.quantile(nsub_rand, 0.95)),
            "delta_obs_minus_rand_mean": float(nsub_obs - np.mean(nsub_rand)),
            "message": "",
        }
        rows.append(row)

    out_df = pd.DataFrame(rows)

    # Sort: OK first, then by p-value
    if "p_value" in out_df.columns:
        out_df["p_value_sort"] = pd.to_numeric(out_df["p_value"], errors="coerce")
        out_df = out_df.sort_values(by=["status", "p_value_sort"], ascending=[True, True]).drop(columns=["p_value_sort"], errors="ignore")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)

    # Print small summary
    n_total = len(out_df)
    n_ok = int((out_df["status"] == "OK").sum()) if "status" in out_df.columns else 0
    print(f"Wrote: {out_path}")
    print(f"Families processed: {n_total}")
    print(f"OK (p-values computed): {n_ok}")
    if n_ok > 0:
        pvals = pd.to_numeric(out_df.loc[out_df["status"] == "OK", "p_value"], errors="coerce").dropna().to_numpy()
        print(f"p-value summary: n={len(pvals)}, min={pvals.min():.4g}, median={np.median(pvals):.4g}, mean={pvals.mean():.4g}, max={pvals.max():.4g}")


if __name__ == "__main__":
    main()
