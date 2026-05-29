#!/usr/bin/env python3
"""
06_analyze_nsub_results.py

Input:
  master TSV with per-family Nsub and permutation p-values.

Outputs (in --out-dir):
  - master_nsub_pvalues.with_q.tsv           (master with BH q-values added)
  - summary_stats.tsv                        (one-row machine-readable summary)
  - summary_stats.txt                        (human-readable summary)
  - PNG plots (NO TITLES):
      pvalue_hist.png
      qvalue_hist.png
      delta_hist.png
      obs_vs_rand_mean_scatter.png
      volcano_delta_neglogp.png
      nsub_obs_vs_m.png
      delta_vs_m.png
  - optionally PDF versions of the plots (if --pdf)

Usage:
  python3 analyze_master_nsub_no_titles.py \
    --master /media/rdurak/DATA/OMA/11_contactmap_rs/contact/master_nsub_pvalues.tsv \
    --out-dir /media/rdurak/DATA/OMA/11_contactmap_rs/contact/summary_plots \
    --min_m 2 --bins 40
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR (q-values), returned in original order."""
    n = p.size
    order = np.argsort(p)
    p_sorted = p[order]
    ranks = np.arange(1, n + 1, dtype=float)

    q_sorted = p_sorted * n / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]  # monotone
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return q


def safe_neglog10(p: np.ndarray) -> np.ndarray:
    eps = np.nextafter(0, 1)
    return -np.log10(np.maximum(p, eps))


def savefig_all(figpath_png: Path, pdf: bool):
    plt.tight_layout()
    plt.savefig(figpath_png, dpi=300)
    if pdf:
        plt.savefig(figpath_png.with_suffix(".pdf"))
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--master",
        default="/media/rdurak/DATA/OMA/11_contactmap_rs/contact/master_nsub_pvalues.tsv",
        help="Path to master TSV",
    )
    ap.add_argument(
        "--out-dir",
        default="/media/rdurak/DATA/OMA/11_contactmap_rs/contact/summary_plots",
        help="Directory to write plots and outputs",
    )
    ap.add_argument("--bins", type=int, default=40, help="Bins for histograms")
    ap.add_argument("--min_m", type=int, default=2, help="Minimum NbMapped to include (default 2)")
    ap.add_argument("--pdf", action="store_true", help="Also save PDF versions of plots")
    args = ap.parse_args()

    master = Path(args.master).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not master.is_file():
        raise SystemExit(f"ERROR: master file not found: {master}")

    df = pd.read_csv(master, sep="\t", dtype=str, keep_default_na=False)

    required = [
        "family",
        "NbMapped",
        "Nsub_obs",
        "Nsub_rand_mean",
        "p_value",
        "delta_obs_minus_rand_mean",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: master TSV missing columns: {missing}")

    # numeric conversion
    for col in ["NbMapped", "Nsub_obs", "Nsub_rand_mean", "p_value", "delta_obs_minus_rand_mean"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # filter
    keep = (
        np.isfinite(df["NbMapped"])
        & np.isfinite(df["Nsub_obs"])
        & np.isfinite(df["Nsub_rand_mean"])
        & np.isfinite(df["p_value"])
        & np.isfinite(df["delta_obs_minus_rand_mean"])
        & (df["NbMapped"] >= args.min_m)
        & (df["p_value"] >= 0.0)
        & (df["p_value"] <= 1.0)
    )
    d = df.loc[keep].copy()
    if d.empty:
        raise SystemExit("ERROR: no valid rows after filtering.")

    # BH-FDR
    pvals = d["p_value"].to_numpy(dtype=float)
    d["q_value_bh"] = bh_fdr(pvals)

    # write master with q
    out_master = out_dir / "master_nsub_pvalues.with_q.tsv"
    d.to_csv(out_master, sep="\t", index=False)

    # summaries
    n = int(len(d))
    qvals = d["q_value_bh"].to_numpy(dtype=float)
    delta = d["delta_obs_minus_rand_mean"].to_numpy(dtype=float)

    n_fdr05 = int((qvals <= 0.05).sum())
    n_fdr10 = int((qvals <= 0.10).sum())

    # Wilcoxon signed-rank (Delta < 0)
    w_stat, w_p = wilcoxon(delta, alternative="less", zero_method="wilcox")

    # Extra helpful summary stats
    prop_delta_lt0 = float(np.mean(delta < 0))
    prop_delta_gt0 = float(np.mean(delta > 0))
    prop_delta_eq0 = float(np.mean(delta == 0))

    summary_row = {
        "input_master": str(master),
        "min_m": args.min_m,
        "n_families": n,

        "nsub_obs_min": float(d["Nsub_obs"].min()),
        "nsub_obs_median": float(d["Nsub_obs"].median()),
        "nsub_obs_mean": float(d["Nsub_obs"].mean()),
        "nsub_obs_max": float(d["Nsub_obs"].max()),

        "delta_min": float(np.min(delta)),
        "delta_median": float(np.median(delta)),
        "delta_mean": float(np.mean(delta)),
        "delta_max": float(np.max(delta)),
        "delta_prop_lt0": prop_delta_lt0,
        "delta_prop_gt0": prop_delta_gt0,
        "delta_prop_eq0": prop_delta_eq0,

        "p_min": float(np.min(pvals)),
        "p_median": float(np.median(pvals)),
        "p_mean": float(np.mean(pvals)),
        "p_max": float(np.max(pvals)),

        "bh_fdr_sig_q05": n_fdr05,
        "bh_fdr_sig_q10": n_fdr10,
        "bh_fdr_sig_q05_frac": float(n_fdr05 / n),
        "bh_fdr_sig_q10_frac": float(n_fdr10 / n),

        "wilcoxon_stat": float(w_stat),
        "wilcoxon_p_less": float(w_p),

        "output_master_with_q": str(out_master),
    }

    # write summary tsv
    summary_tsv = out_dir / "summary_stats.tsv"
    pd.DataFrame([summary_row]).to_csv(summary_tsv, sep="\t", index=False)

    # write summary txt
    summary_txt = out_dir / "summary_stats.txt"
    with summary_txt.open("w") as fh:
        fh.write("=== Summary (filtered) ===\n")
        fh.write(f"Input master: {master}\n")
        fh.write(f"Output dir:   {out_dir}\n")
        fh.write(f"Rows used:    {n}  (NbMapped >= {args.min_m})\n\n")

        fh.write("=== Nsub / Delta summaries ===\n")
        fh.write(
            "Nsub_obs:  "
            f"min={summary_row['nsub_obs_min']:.6g}, "
            f"median={summary_row['nsub_obs_median']:.6g}, "
            f"mean={summary_row['nsub_obs_mean']:.6g}, "
            f"max={summary_row['nsub_obs_max']:.6g}\n"
        )
        fh.write(
            "Delta (obs - rand_mean): "
            f"min={summary_row['delta_min']:.6g}, "
            f"median={summary_row['delta_median']:.6g}, "
            f"mean={summary_row['delta_mean']:.6g}, "
            f"max={summary_row['delta_max']:.6g}\n"
        )
        fh.write(
            "Delta sign: "
            f"prop(Δ<0)={summary_row['delta_prop_lt0']:.6g}, "
            f"prop(Δ>0)={summary_row['delta_prop_gt0']:.6g}, "
            f"prop(Δ=0)={summary_row['delta_prop_eq0']:.6g}\n\n"
        )

        fh.write("=== Permutation p-values ===\n")
        fh.write(
            "p: "
            f"min={summary_row['p_min']:.6g}, "
            f"median={summary_row['p_median']:.6g}, "
            f"mean={summary_row['p_mean']:.6g}, "
            f"max={summary_row['p_max']:.6g}\n"
        )
        fh.write(
            f"BH-FDR significant: q<=0.05: {n_fdr05}/{n} ({summary_row['bh_fdr_sig_q05_frac']:.2%})\n"
        )
        fh.write(
            f"BH-FDR significant: q<=0.10: {n_fdr10}/{n} ({summary_row['bh_fdr_sig_q10_frac']:.2%})\n\n"
        )

        fh.write("=== Global test (paired across families) ===\n")
        fh.write("Wilcoxon signed-rank test on Delta (alternative: Delta < 0):\n")
        fh.write(f"  statistic={w_stat:.6g}, p={w_p:.6g}\n\n")
        fh.write(f"Wrote updated master with q-values: {out_master}\n")

    # also print a short note
    print(f"Wrote: {out_master}")
    print(f"Wrote: {summary_tsv}")
    print(f"Wrote: {summary_txt}")

    # ---------- Plots (NO TITLES) ----------

    # 1) p-value histogram
    plt.figure()
    plt.hist(pvals, bins=args.bins)
    plt.xlabel("Permutation p-value")
    plt.ylabel("Number of protein families")
    savefig_all(out_dir / "pvalue_hist.png", args.pdf)

    # 2) q-value histogram
    plt.figure()
    plt.hist(qvals, bins=args.bins)
    plt.xlabel("BH FDR q-value")
    plt.ylabel("Number of protein families")
    savefig_all(out_dir / "qvalue_hist.png", args.pdf)

    # 3) Delta histogram
    plt.figure()
    plt.hist(delta, bins=args.bins)
    plt.xlabel(r"$\Delta = N_{\mathrm{sub}}^{obs} - \mathrm{mean}(N_{\mathrm{sub}}^{rand})$")
    plt.ylabel("Number of protein families")
    savefig_all(out_dir / "delta_hist.png", args.pdf)

    # 4) Observed vs Random mean scatter + y=x line
    x = d["Nsub_rand_mean"].to_numpy(dtype=float)
    y = d["Nsub_obs"].to_numpy(dtype=float)
    plt.figure()
    plt.scatter(x, y, s=8)
    lo = float(min(x.min(), y.min()))
    hi = float(max(x.max(), y.max()))
    plt.plot([lo, hi], [lo, hi])  # y=x
    plt.xlabel(r"mean($N_{\mathrm{sub}}^{rand}$)")
    plt.ylabel(r"$N_{\mathrm{sub}}^{obs}$")
    savefig_all(out_dir / "obs_vs_rand_mean_scatter.png", args.pdf)

    # 5) Volcano: Delta vs -log10(p)
    neglogp = safe_neglog10(pvals)
    plt.figure()
    plt.scatter(delta, neglogp, s=8)
    plt.xlabel(r"$\Delta = N_{\mathrm{sub}}^{obs} - \mathrm{mean}(N_{\mathrm{sub}}^{rand})$")
    plt.ylabel(r"$-\log_{10}(p)$")
    savefig_all(out_dir / "volcano_delta_neglogp.png", args.pdf)

    # 6) Nsub_obs vs NbMapped
    m = d["NbMapped"].to_numpy(dtype=float)
    plt.figure()
    plt.scatter(m, y, s=8)
    plt.xlabel("NbMapped (group size, m)")
    plt.ylabel(r"$N_{\mathrm{sub}}^{obs}$")
    savefig_all(out_dir / "nsub_obs_vs_m.png", args.pdf)

    # 7) Delta vs NbMapped
    plt.figure()
    plt.scatter(m, delta, s=8)
    plt.xlabel("NbMapped (group size, m)")
    plt.ylabel(r"$\Delta = N_{\mathrm{sub}}^{obs} - \mathrm{mean}(N_{\mathrm{sub}}^{rand})$")
    savefig_all(out_dir / "delta_vs_m.png", args.pdf)

    print(f"Wrote plots in: {out_dir}")


if __name__ == "__main__":
    main()

