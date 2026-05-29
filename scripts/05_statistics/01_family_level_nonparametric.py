#!/usr/bin/env python3
"""
Non-parametric tests on family-level rate-shifts summary.

Tests:
  1) Kendall tau: prop_rate_shift vs tree_length
  2) Kendall tau: prop_rate_shift vs num_sequences
  3) Kruskal-Wallis: prop_rate_shift across Last_Common_Ancestor groups
  4) Posthoc Dunn test (if KW significant), with BH-FDR correction

Outputs (in outdir):
  - kendall_results.tsv
  - kruskal_wallis.tsv
  - dunn_posthoc_FDR.tsv  (only if KW p <= alpha)

Notes:
  - Uses scipy.stats for kendalltau and kruskal.
  - Dunn test implemented internally (normal approximation on rank sums).
"""

import argparse
import os
import numpy as np
import pandas as pd
from math import erf, sqrt

MISSING_STRINGS = {"", "NA", "NaN", "NAN", "nan", ".", "None", "NONE", "null", "NULL"}

REQ_COLS = [
    "Orthogroup",
    "prop_rate_shift",
    "Last_Common_Ancestor",
    "tree_length",
    "num_sequences",
    "num_sites",
]


def clean_str(s: pd.Series) -> pd.Series:
    out = s.astype("string").str.strip()
    out = out.mask(out.isin(MISSING_STRINGS), pd.NA)
    return out


def coerce_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def norm_sf(z: float) -> float:
    """Survival function for standard normal using erf (no scipy dependency here)."""
    # sf(z) = 0.5 * erfc(z/sqrt(2)) = 0.5 * (1 - erf(z/sqrt(2)))
    return 0.5 * (1.0 - erf(z / sqrt(2.0)))


def two_sided_p_from_z(z: float) -> float:
    return 2.0 * norm_sf(abs(z))


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def dunn_posthoc(x: np.ndarray, groups: np.ndarray) -> pd.DataFrame:
    """
    Dunn's test on ranks with tie correction (approx normal).
    x: numeric array
    groups: array-like group labels (same length as x)

    Returns a tidy DataFrame with pairwise comparisons:
      group1, group2, z, p_raw
    """
    # Remove missing
    mask = np.isfinite(x) & pd.notna(groups)
    x = x[mask]
    groups = groups[mask]

    # Rank all observations (average ranks for ties)
    ranks = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    N = len(ranks)

    # Tie correction factor for variance (Kruskal-Wallis style)
    # T = sum(t^3 - t) over tie groups
    _, tie_counts = np.unique(x, return_counts=True)
    T = np.sum(tie_counts**3 - tie_counts)
    tie_corr = 1.0 - T / (N**3 - N) if N > 1 else 1.0

    # Group summaries
    df = pd.DataFrame({"rank": ranks, "group": groups})
    gstats = df.groupby("group")["rank"].agg(["count", "mean"]).reset_index()
    labels = gstats["group"].to_list()
    n_i = gstats["count"].to_numpy(dtype=float)
    mean_i = gstats["mean"].to_numpy(dtype=float)

    # Variance of mean rank difference
    # Var = (N*(N+1)/12) * (1/n_i + 1/n_j) / tie_corr
    base = (N * (N + 1) / 12.0)
    if tie_corr <= 0:
        tie_corr = 1.0

    out_rows = []
    for a in range(len(labels)):
        for b in range(a + 1, len(labels)):
            diff = mean_i[a] - mean_i[b]
            se = sqrt((base * (1.0 / n_i[a] + 1.0 / n_i[b])) / tie_corr)
            z = diff / se if se > 0 else np.nan
            p = two_sided_p_from_z(z) if np.isfinite(z) else np.nan
            out_rows.append((labels[a], labels[b], z, p))

    return pd.DataFrame(out_rows, columns=["group1", "group2", "z", "p_raw"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Family-level TSV file")
    ap.add_argument("-o", "--outdir", default="nonparam_family_out", help="Output directory")
    ap.add_argument("--min-sites", type=int, default=0, help="Filter to orthogroups with num_sites >= this")
    ap.add_argument("--min-per-group", type=int, default=10,
                    help="Drop LCA groups with fewer than this many orthogroups before KW/Dunn")
    ap.add_argument("--alpha", type=float, default=0.05, help="Alpha for deciding whether to run posthoc")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False)

    missing = [c for c in REQ_COLS if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    # Clean & coerce
    df["Orthogroup"] = clean_str(df["Orthogroup"])
    df["Last_Common_Ancestor"] = clean_str(df["Last_Common_Ancestor"])

    for c in ["prop_rate_shift", "tree_length", "num_sequences", "num_sites"]:
        df[c] = coerce_numeric(df[c])

    df = df.dropna(subset=["Orthogroup", "Last_Common_Ancestor", "prop_rate_shift", "tree_length", "num_sequences", "num_sites"])

    if args.min_sites > 0:
        df = df.loc[df["num_sites"] >= args.min_sites].copy()

    # Clamp proportion safely
    df = df.loc[(df["prop_rate_shift"] >= -1e-9) & (df["prop_rate_shift"] <= 1 + 1e-9)].copy()
    df["prop_rate_shift"] = df["prop_rate_shift"].clip(0, 1)

    # -------------------------
    # 1) Kendall correlations
    # -------------------------
    from scipy.stats import kendalltau, kruskal  # scipy is expected in your environment

    kendall_rows = []
    for xcol in ["tree_length", "num_sequences"]:
        x = df[xcol].to_numpy(float)
        y = df["prop_rate_shift"].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)
        tau, p = kendalltau(x[mask], y[mask], nan_policy="omit")
        kendall_rows.append({
            "x": xcol,
            "y": "prop_rate_shift",
            "n": int(np.sum(mask)),
            "kendall_tau": float(tau),
            "p_value": float(p),
        })

    kendall_df = pd.DataFrame(kendall_rows)
    kendall_df.to_csv(os.path.join(args.outdir, "kendall_results.tsv"), sep="\t", index=False)

    # -------------------------
    # 2) Kruskal-Wallis across LCA
    # -------------------------
    # Drop small groups
    counts = df["Last_Common_Ancestor"].value_counts()
    keep = counts[counts >= args.min_per_group].index
    df_kw = df.loc[df["Last_Common_Ancestor"].isin(keep)].copy()

    groups = []
    labels = []
    for lca, g in df_kw.groupby("Last_Common_Ancestor"):
        vals = g["prop_rate_shift"].to_numpy(float)
        if len(vals) >= args.min_per_group:
            groups.append(vals)
            labels.append(str(lca))

    kw_H, kw_p = kruskal(*groups) if len(groups) >= 2 else (np.nan, np.nan)

    kw_out = pd.DataFrame([{
        "test": "kruskal_wallis",
        "response": "prop_rate_shift",
        "factor": "Last_Common_Ancestor",
        "n_total": int(len(df_kw)),
        "n_groups": int(len(groups)),
        "min_per_group": args.min_per_group,
        "H": float(kw_H) if np.isfinite(kw_H) else np.nan,
        "p_value": float(kw_p) if np.isfinite(kw_p) else np.nan,
        "alpha": args.alpha,
        "posthoc_ran": bool(np.isfinite(kw_p) and kw_p <= args.alpha),
    }])
    kw_out.to_csv(os.path.join(args.outdir, "kruskal_wallis.tsv"), sep="\t", index=False)

    # -------------------------
    # 3) Dunn posthoc (if KW significant)
    # -------------------------
    if np.isfinite(kw_p) and kw_p <= args.alpha:
        dunn = dunn_posthoc(
            df_kw["prop_rate_shift"].to_numpy(float),
            df_kw["Last_Common_Ancestor"].to_numpy(object),
        )
        dunn["p_fdr_bh"] = bh_fdr(dunn["p_raw"].to_numpy(float))
        dunn = dunn.sort_values("p_fdr_bh", ascending=True)
        dunn.to_csv(os.path.join(args.outdir, "dunn_posthoc_FDR.tsv"), sep="\t", index=False)

    print("Wrote outputs to:", args.outdir)


if __name__ == "__main__":
    main()
