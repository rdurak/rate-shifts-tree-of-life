#!/usr/bin/env python3
"""
GLMM (logistic mixed model) in Python for site-level rate shifts with domain interactions.

Assumes input TSV contains:
  - LCA_2  (values: LUCA, Bacteria, Archaea, Eukaryota)

Response:
  is_rate_shift = 1 if Proba>0.95 == "*", else 0

Default filter:
  structurally mapped sites (PDB present), unless --no-pdb-filter

Base fixed effects (in ALL models):
  C(LCA_2) + Rsa_z + Depth_z + C(SS3) + log_tree_z + log_numseq_z

Random effect (ALL models):
  Orthogroup random intercept

Fits 5 models into separate subdirectories under --outdir:
  M0_base_no_interactions
  M1_int_LCA2_x_SS3
  M2_int_LCA2_x_Rsa
  M3_int_LCA2_x_Depth
  M4_int_all (SS3 + Rsa + Depth interactions with LCA_2)

Implementation:
  statsmodels.genmod.bayes_mixed_glm.BinomialBayesMixedGLM (variational Bayes)

Outputs per model:
  - model_summary.txt
  - fixed_effects_OR_pvalues.tsv
  - run_info.txt

Usage:
  python3 glmm_LCA2_interactions.py \
      -i all_sged_with_lca_final.with_tree_seq.with_LCA2.tsv \
      -o py_glmm_out_LCA2_interactions \
      --depth CalphaDepth
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd

MISSING_STRINGS = {"", "NA", "NaN", "NAN", "nan", ".", "None", "NONE", "null", "NULL"}


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


def zscore(series: pd.Series) -> pd.Series:
    """Z-score numeric series; returns NaNs if variance is zero or all missing."""
    x = pd.to_numeric(series, errors="coerce")
    arr = x.to_numpy(dtype=float)
    mu = np.nanmean(arr)
    sd = np.nanstd(arr, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)
    return (x - mu) / sd


def log1p_numeric(series: pd.Series) -> pd.Series:
    """log1p on a numeric series with coercion."""
    x = pd.to_numeric(series, errors="coerce")
    return np.log1p(x)


def clean_cat(df: pd.DataFrame, col: str) -> None:
    df[col] = df[col].astype("string").str.strip()
    df[col] = df[col].mask(df[col].isin(MISSING_STRINGS), pd.NA)


def fit_and_write(df: pd.DataFrame, formula: str, vc_formulas: dict, outdir: str, run_meta: dict):
    """Fit BinomialBayesMixedGLM (VB) and write outputs."""
    import statsmodels.api as sm  # noqa: F401
    from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
    from scipy.stats import norm

    os.makedirs(outdir, exist_ok=True)

    print(f"\n--- Fitting model in: {outdir} ---")
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
        f.write(f"LCA2_levels\t{df['LCA_2'].nunique()}\n")
    print(f"Wrote run metadata to: {meta_path}")

    return res


def main():
    ap = argparse.ArgumentParser(description="Fit GLMMs with LCA_2 interactions (one-by-one and combined).")
    ap.add_argument("-i", "--input", required=True, help="Site-level TSV (must include LCA_2)")
    ap.add_argument("-o", "--outdir", default="py_glmm_out_LCA2_interactions", help="Base output directory")
    ap.add_argument("--depth", choices=["CalphaDepth", "ResidueDepth"], default="CalphaDepth",
                    help="Depth variable to use (default: CalphaDepth)")
    ap.add_argument("--sample-n", type=int, default=0, help="Optional subsample size (0 = use all mapped sites)")
    ap.add_argument("--seed", type=int, default=1, help="Random seed for subsampling")
    ap.add_argument("--no-pdb-filter", action="store_true",
                    help="If set, do NOT filter to PDB-mapped sites.")
    args = ap.parse_args()

    # Imports needed for modeling
    try:
        import statsmodels.api as sm  # noqa: F401
        from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM  # noqa: F401
    except Exception:
        print("ERROR: statsmodels is required. Install with: pip install statsmodels", file=sys.stderr)
        raise

    try:
        from scipy.stats import norm  # noqa: F401
    except Exception:
        print("ERROR: scipy is required. Install with: pip install scipy", file=sys.stderr)
        raise

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading {args.input} ...")
    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False, na_filter=True)

    # Validate required cols
    required_cols = [
        "Proba>0.95", "Orthogroup", "SecondaryStructure", "Rsa", args.depth,
        "tree_length", "num_sequences", "LCA_2"
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: Missing required columns: {missing}")

    # -----------------------------
    # Response: is_rate_shift
    # -----------------------------
    df["is_rate_shift"] = (
        df["Proba>0.95"]
        .astype("string")
        .fillna("")
        .str.strip()
        .eq("*")
        .astype("int8")
    )

    # -----------------------------
    # Clean key categorical fields
    # -----------------------------
    clean_cat(df, "Orthogroup")
    clean_cat(df, "LCA_2")

    # -----------------------------
    # Filter to structurally mapped (unless disabled)
    # -----------------------------
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

    # Optional subsample
    if args.sample_n and len(df) > args.sample_n:
        rng = np.random.default_rng(args.seed)
        take = rng.choice(df.index.to_numpy(), size=args.sample_n, replace=False)
        df = df.loc[take].copy()
        print(f"Subsampled to: {len(df):,} rows")

    # -----------------------------
    # Predictors
    # -----------------------------
    df["SS3"] = ss3_from_dssp(df["SecondaryStructure"])
    df["SS3"] = pd.Categorical(df["SS3"], categories=["Coil", "Helix", "Strand"])

    df["Rsa_z"] = zscore(df["Rsa"])
    df["Depth_z"] = zscore(df[args.depth])

    df["log_tree_length"] = log1p_numeric(df["tree_length"])
    df["log_num_sequences"] = log1p_numeric(df["num_sequences"])
    df["log_tree_z"] = zscore(df["log_tree_length"])
    df["log_numseq_z"] = zscore(df["log_num_sequences"])

    keep = [
        "is_rate_shift", "Orthogroup", "LCA_2", "SS3",
        "Rsa_z", "Depth_z", "log_tree_z", "log_numseq_z"
    ]
    df = df[keep].copy().dropna()

    print(f"Rows after dropping missing predictors: {len(df):,}")
    print(f"Orthogroups retained: {df['Orthogroup'].nunique():,}")
    print(f"LCA_2 levels retained: {df['LCA_2'].nunique():,}")

    # Random intercept per Orthogroup
    vc = {"Orthogroup_RE": "0 + C(Orthogroup)"}

    base_terms = [
        "C(LCA_2)",
        "Rsa_z",
        "Depth_z",
        "C(SS3)",
        "log_tree_z",
        "log_numseq_z",
    ]

    # Define model formulas (interactions one-by-one + combined)
    models = [
        ("M0_base_no_interactions",
         "is_rate_shift ~ " + " + ".join(base_terms),
         "Base model with LCA_2 fixed; no interactions."),

        ("M1_int_LCA2_x_SS3",
         "is_rate_shift ~ " + " + ".join(base_terms + ["C(LCA_2):C(SS3)"]),
         "Adds interaction: LCA_2 × SS3."),

        ("M2_int_LCA2_x_Rsa",
         "is_rate_shift ~ " + " + ".join(base_terms + ["C(LCA_2):Rsa_z"]),
         "Adds interaction: LCA_2 × Rsa_z."),

        ("M3_int_LCA2_x_Depth",
         "is_rate_shift ~ " + " + ".join(base_terms + ["C(LCA_2):Depth_z"]),
         "Adds interaction: LCA_2 × Depth_z."),

        ("M4_int_all",
         "is_rate_shift ~ " + " + ".join(base_terms + [
             "C(LCA_2):C(SS3)",
             "C(LCA_2):Rsa_z",
             "C(LCA_2):Depth_z",
         ]),
         "Adds interactions: LCA_2 × (SS3 + Rsa_z + Depth_z)."),
    ]

    meta_common = {
        "input": args.input,
        "depth_variable": args.depth,
        "pdb_filter": (not args.no_pdb_filter),
        "sample_n": args.sample_n,
        "seed": args.seed,
        "random_effects": "Orthogroup intercept",
    }

    for name, formula, notes in models:
        outdir = os.path.join(args.outdir, name)
        meta = meta_common.copy()
        meta.update({"model": name, "notes": notes})
        fit_and_write(df, formula, vc, outdir, meta)

    print("\nDone. Wrote model directories under:", args.outdir)


if __name__ == "__main__":
    main()

