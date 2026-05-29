#!/usr/bin/env python3
"""
Fit separate logistic mixed models (one per LCA_2 level) for site-level rate shifts.

For each LCA_2 in {Archaea, Bacteria, Eukaryota, LUCA}, fits:

  is_rate_shift ~ Rsa_z + Depth_z + C(SS3) + log_tree_z + log_numseq_z
  Random intercept: Orthogroup

Notes on scaling:
- Default: compute z-scores globally on the full (filtered) dataset, then subset by LCA_2.
  This makes coefficients comparable across LCA-specific models.
- Optional: --zscore-within-lca to compute z-scores separately within each LCA_2 subset.

Outputs:
  For each LCA_2, writes under --outdir/LCA2_<level>/:
    - model_summary.txt
    - fixed_effects_OR_pvalues.tsv
    - run_info.txt

Usage:
  python3 glmm_python_LCA2_separate_models.py \
      -i all_sged_with_lca_final.with_tree_seq.with_LCA2.tsv \
      -o py_glmm_out_LCA2_separate --depth CalphaDepth
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

MISSING_STRINGS = {"", "NA", "NaN", "NAN", "nan", ".", "None", "NONE", "null", "NULL"}
LCA_ORDER = ["Archaea", "Bacteria", "Eukaryota", "LUCA"]


def ss3_from_dssp(ss: pd.Series) -> pd.Series:
    """Map DSSP one-letter code to SS3: Helix/Strand/Coil."""
    s = ss.astype("string").str.strip()
    s = s.mask(s.isin(MISSING_STRINGS), pd.NA)

    helix = {"H", "G", "I"}
    strand = {"E", "B"}

    out = pd.Series(pd.NA, index=s.index, dtype="string")
    out[s.isin(list(helix))] = "Helix"
    out[s.isin(list(strand))] = "Strand"
    out[out.isna() & s.notna()] = "Coil"
    return out


def clean_cat(df: pd.DataFrame, col: str) -> None:
    df[col] = df[col].astype("string").str.strip()
    df[col] = df[col].mask(df[col].isin(MISSING_STRINGS), pd.NA)


def zscore_from_array(arr: np.ndarray) -> np.ndarray:
    """Z-score an array, ignoring NaNs. Returns z-scored array with NaNs preserved."""
    mu = np.nanmean(arr)
    sd = np.nanstd(arr, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(arr, np.nan, dtype=float)
    return (arr - mu) / sd


def log1p_numeric(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    return np.log1p(x)


def fit_and_write(df: pd.DataFrame, formula: str, vc_formulas: dict, outdir: str, run_meta: dict):
    """Fit BinomialBayesMixedGLM (VB) and write outputs."""
    from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
    from scipy.stats import norm

    os.makedirs(outdir, exist_ok=True)

    print(f"\n--- Fitting model in: {outdir} ---")
    print("Rows:", len(df))
    print("Formula:", formula)
    print("Random effects (vc_formulas):", vc_formulas)

    model = BinomialBayesMixedGLM.from_formula(formula, vc_formulas, df)
    res = model.fit_vb()

    # Save summary
    summary_path = os.path.join(outdir, "model_summary.txt")
    with open(summary_path, "w") as f:
        f.write(str(res.summary()))
        f.write("\n")
    print(f"Wrote model summary to: {summary_path}")

    # Fixed effects table: OR, CI, z, p (approx)
    fe_names = res.model.exog_names
    fe_mean = np.asarray(res.fe_mean)
    fe_sd = np.asarray(res.fe_sd)

    z = fe_mean / fe_sd
    p = 2 * (1 - norm.cdf(np.abs(z)))

    lo = fe_mean - 1.96 * fe_sd
    hi = fe_mean + 1.96 * fe_sd

    out = pd.DataFrame({
        "term": fe_names,
        "estimate_logit": fe_mean,
        "se": fe_sd,
        "z": z,
        "p_value": p,
        "OR": np.exp(fe_mean),
        "CI_low": np.exp(lo),
        "CI_high": np.exp(hi),
    })

    out_path = os.path.join(outdir, "fixed_effects_OR_pvalues.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote fixed effects table to: {out_path}")

    # Metadata
    meta_path = os.path.join(outdir, "run_info.txt")
    with open(meta_path, "w") as f:
        for k, v in run_meta.items():
            f.write(f"{k}\t{v}\n")
        f.write(f"formula\t{formula}\n")
        f.write(f"vc_formulas\t{vc_formulas}\n")
        f.write(f"rows_used\t{len(df)}\n")
        f.write(f"orthogroups\t{df['Orthogroup'].nunique()}\n")
    print(f"Wrote run metadata to: {meta_path}")

    return res


def main():
    ap = argparse.ArgumentParser(description="Fit separate logistic mixed models for each LCA_2.")
    ap.add_argument("-i", "--input", required=True, help="Path to site-level TSV (must include LCA_2)")
    ap.add_argument("-o", "--outdir", default="py_glmm_out_LCA2_separate", help="Base output directory")
    ap.add_argument("--depth", choices=["CalphaDepth", "ResidueDepth"], default="CalphaDepth",
                    help="Depth variable to use (default: CalphaDepth)")
    ap.add_argument("--sample-n", type=int, default=0, help="Optional subsample size per LCA (0 = use all)")
    ap.add_argument("--seed", type=int, default=1, help="Random seed for subsampling")
    ap.add_argument("--no-pdb-filter", action="store_true",
                    help="If set, do NOT filter to PDB-mapped sites.")
    ap.add_argument("--zscore-within-lca", action="store_true",
                    help="If set, compute z-scores separately within each LCA subset (not globally).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading {args.input} ...")
    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False, na_filter=True)

    required_cols = ["Proba>0.95", "Orthogroup", "SecondaryStructure", "Rsa", args.depth,
                     "tree_length", "num_sequences", "LCA_2"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: Missing required columns: {missing}")

    # Response
    df["is_rate_shift"] = (
        df["Proba>0.95"].astype("string").fillna("").str.strip().eq("*").astype("int8")
    )

    # Clean categorical
    clean_cat(df, "Orthogroup")
    clean_cat(df, "LCA_2")

    # Filter to structurally mapped unless disabled
    if not args.no_pdb_filter:
        if "PDB" not in df.columns:
            raise SystemExit("ERROR: PDB column missing, but PDB filter is enabled.")
        pdb = df["PDB"].astype("string").str.strip()
        pdb = pdb.mask(pdb.isin(MISSING_STRINGS), pd.NA)
        df = df.loc[pdb.notna()].copy()
        print(f"Mapped rows (PDB present): {len(df):,}")
    else:
        df = df.copy()
        print(f"Using ALL rows (no PDB filter): {len(df):,}")

    # Keep only expected LCA_2 levels
    df = df[df["LCA_2"].isin(LCA_ORDER)].copy()
    print("Rows after restricting to LCA_2 in", LCA_ORDER, ":", f"{len(df):,}")

    # Build predictors that don't depend on z-scoring choice
    df["SS3"] = ss3_from_dssp(df["SecondaryStructure"])
    df["SS3"] = pd.Categorical(df["SS3"], categories=["Coil", "Helix", "Strand"])

    df["log_tree_length"] = log1p_numeric(df["tree_length"])
    df["log_num_sequences"] = log1p_numeric(df["num_sequences"])

    # Coerce continuous to numeric once
    df["_Rsa"] = pd.to_numeric(df["Rsa"], errors="coerce")
    df["_Depth"] = pd.to_numeric(df[args.depth], errors="coerce")
    df["_log_tree"] = pd.to_numeric(df["log_tree_length"], errors="coerce")
    df["_log_numseq"] = pd.to_numeric(df["log_num_sequences"], errors="coerce")

    # Drop rows missing essentials before any scaling
    base_keep = ["is_rate_shift", "Orthogroup", "LCA_2", "SS3", "_Rsa", "_Depth", "_log_tree", "_log_numseq"]
    df = df[base_keep].dropna().copy()
    print(f"Rows after dropping missing essentials: {len(df):,}")
    print(f"Orthogroups retained: {df['Orthogroup'].nunique():,}")

    rng = np.random.default_rng(args.seed)

    # Global z-scoring (default) so coefficients are comparable across LCA-specific fits
    if not args.zscore_within_lca:
        df["Rsa_z"] = zscore_from_array(df["_Rsa"].to_numpy(dtype=float))
        df["Depth_z"] = zscore_from_array(df["_Depth"].to_numpy(dtype=float))
        df["log_tree_z"] = zscore_from_array(df["_log_tree"].to_numpy(dtype=float))
        df["log_numseq_z"] = zscore_from_array(df["_log_numseq"].to_numpy(dtype=float))

    vc = {"Orthogroup_RE": "0 + C(Orthogroup)"}
    formula = "is_rate_shift ~ Rsa_z + Depth_z + C(SS3) + log_tree_z + log_numseq_z"

    for lca in LCA_ORDER:
        sub = df[df["LCA_2"] == lca].copy()
        if sub.empty:
            print(f"WARNING: No rows for LCA_2={lca}; skipping.")
            continue

        # Optional per-LCA subsample
        if args.sample_n and len(sub) > args.sample_n:
            take = rng.choice(sub.index.to_numpy(), size=args.sample_n, replace=False)
            sub = sub.loc[take].copy()

        # If requested, z-score within each LCA subset instead of globally
        if args.zscore_within_lca:
            sub["Rsa_z"] = zscore_from_array(sub["_Rsa"].to_numpy(dtype=float))
            sub["Depth_z"] = zscore_from_array(sub["_Depth"].to_numpy(dtype=float))
            sub["log_tree_z"] = zscore_from_array(sub["_log_tree"].to_numpy(dtype=float))
            sub["log_numseq_z"] = zscore_from_array(sub["_log_numseq"].to_numpy(dtype=float))

        # Final drop in case z-scoring produced NaNs (e.g., zero variance)
        sub = sub.dropna(subset=["Rsa_z", "Depth_z", "log_tree_z", "log_numseq_z"]).copy()

        outdir = os.path.join(args.outdir, f"LCA2_{lca}")
        meta = {
            "input": args.input,
            "depth_variable": args.depth,
            "lca_subset": lca,
            "pdb_filter": (not args.no_pdb_filter),
            "sample_n": args.sample_n,
            "seed": args.seed,
            "zscore_within_lca": bool(args.zscore_within_lca),
            "random_effects": "Orthogroup intercept",
            "notes": "Separate model fit on subset defined by LCA_2.",
        }

        print(f"\n=== LCA_2 = {lca} ===")
        print(f"Rows used: {len(sub):,} | Orthogroups: {sub['Orthogroup'].nunique():,}")

        fit_and_write(sub, formula, vc, outdir, meta)

    print("\nDone. Wrote per-LCA model directories under:", args.outdir)


if __name__ == "__main__":
    main()

