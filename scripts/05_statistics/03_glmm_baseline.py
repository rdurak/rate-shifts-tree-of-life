#!/usr/bin/env python3

"""
Site-level logistic mixed model for evolutionary rate shifts.

Model:
    is_rate_shift ~
        C(Last_Common_Ancestor) +
        Rsa_z +
        Depth_z +
        C(SS3) +
        log_tree_z +
        log_numseq_z

Random effect:
    Orthogroup intercept

Implementation:
    statsmodels.genmod.bayes_mixed_glm.BinomialBayesMixedGLM
    (variational Bayes approximation)

Outputs:
    - model_summary.txt
    - fixed_effects_OR_pvalues.tsv
    - run_info.txt

Usage:
    python3 glmm_rate_shift_model.py \
        -i all_sged_with_lca_final.with_tree_seq.tsv \
        -o py_glmm_out \
        --depth CalphaDepth
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

MISSING_STRINGS = {
    "",
    "NA",
    "NaN",
    "NAN",
    "nan",
    ".",
    "None",
    "NONE",
    "null",
    "NULL",
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def ss3_from_dssp(ss: pd.Series) -> pd.Series:
    """
    Map DSSP one-letter secondary structure codes
    to three-state secondary structure classes.
    """

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
    """
    Z-score numeric series.
    """

    x = pd.to_numeric(series, errors="coerce")

    arr = x.to_numpy(dtype=float)

    mu = np.nanmean(arr)

    sd = np.nanstd(arr, ddof=0)

    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)

    return (x - mu) / sd


def log1p_numeric(series: pd.Series) -> pd.Series:
    """
    Safe log1p transform.
    """

    x = pd.to_numeric(series, errors="coerce")

    return np.log1p(x)


# ============================================================
# MODEL FITTING
# ============================================================

def fit_and_write(
    df: pd.DataFrame,
    formula: str,
    vc_formulas: dict,
    outdir: str,
    run_meta: dict,
):

    import statsmodels.api as sm  # noqa: F401

    from scipy.stats import norm

    from statsmodels.genmod.bayes_mixed_glm import (
        BinomialBayesMixedGLM,
    )

    os.makedirs(outdir, exist_ok=True)

    print("\n--- Fitting model ---")

    print("Formula:", formula)

    print("Random effects:", vc_formulas)

    model = BinomialBayesMixedGLM.from_formula(
        formula,
        vc_formulas,
        df,
    )

    res = model.fit_vb()

    # --------------------------------------------------------
    # Save model summary
    # --------------------------------------------------------

    summary_path = os.path.join(
        outdir,
        "model_summary.txt",
    )

    with open(summary_path, "w") as f:

        f.write(str(res.summary()))

        f.write("\n")

    print(f"Wrote model summary to: {summary_path}")

    # --------------------------------------------------------
    # Fixed effects table
    # --------------------------------------------------------

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

    out_path = os.path.join(
        outdir,
        "fixed_effects_OR_pvalues.tsv",
    )

    out.to_csv(
        out_path,
        sep="\t",
        index=False,
    )

    print(f"Wrote fixed effects table to: {out_path}")

    # --------------------------------------------------------
    # Run metadata
    # --------------------------------------------------------

    meta_path = os.path.join(
        outdir,
        "run_info.txt",
    )

    with open(meta_path, "w") as f:

        for k, v in run_meta.items():

            f.write(f"{k}\t{v}\n")

        f.write(f"formula\t{formula}\n")

        f.write(f"vc_formulas\t{vc_formulas}\n")

        f.write(f"rows_used\t{len(df)}\n")

        f.write(
            f"orthogroups\t{df['Orthogroup'].nunique()}\n"
        )

    print(f"Wrote run metadata to: {meta_path}")

    return res


# ============================================================
# MAIN
# ============================================================

def main():

    ap = argparse.ArgumentParser(
        description="Site-level logistic mixed model for rate shifts."
    )

    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input TSV file",
    )

    ap.add_argument(
        "-o",
        "--outdir",
        default="py_glmm_out",
        help="Output directory",
    )

    ap.add_argument(
        "--depth",
        choices=["CalphaDepth", "ResidueDepth"],
        default="CalphaDepth",
        help="Depth variable to use",
    )

    ap.add_argument(
        "--sample-n",
        type=int,
        default=0,
        help="Optional subsample size",
    )

    ap.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed",
    )

    ap.add_argument(
        "--no-pdb-filter",
        action="store_true",
        help="Use all rows instead of only PDB-mapped sites",
    )

    args = ap.parse_args()

    # --------------------------------------------------------
    # Dependency checks
    # --------------------------------------------------------

    try:

        import statsmodels.api as sm  # noqa: F401

        from statsmodels.genmod.bayes_mixed_glm import (
            BinomialBayesMixedGLM,
        )

    except Exception:

        print(
            "ERROR: statsmodels is required",
            file=sys.stderr,
        )

        raise

    try:

        from scipy.stats import norm  # noqa: F401

    except Exception:

        print(
            "ERROR: scipy is required",
            file=sys.stderr,
        )

        raise

    os.makedirs(args.outdir, exist_ok=True)

    # ========================================================
    # READ DATA
    # ========================================================

    print(f"Reading {args.input} ...")

    df = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
        low_memory=False,
        na_filter=True,
    )

    # ========================================================
    # RESPONSE VARIABLE
    # ========================================================

    df["is_rate_shift"] = (
        df["Proba>0.95"]
        .astype("string")
        .fillna("")
        .str.strip()
        .eq("*")
        .astype("int8")
    )

    # ========================================================
    # CLEAN CATEGORICAL VARIABLES
    # ========================================================

    df["Orthogroup"] = (
        df["Orthogroup"]
        .astype("string")
        .str.strip()
    )

    df["Orthogroup"] = df["Orthogroup"].mask(
        df["Orthogroup"].isin(MISSING_STRINGS),
        pd.NA,
    )

    df["Last_Common_Ancestor"] = (
        df["Last_Common_Ancestor"]
        .astype("string")
        .str.strip()
    )

    df["Last_Common_Ancestor"] = (
        df["Last_Common_Ancestor"].mask(
            df["Last_Common_Ancestor"].isin(
                MISSING_STRINGS
            ),
            pd.NA,
        )
    )

    # ========================================================
    # OPTIONAL PDB FILTER
    # ========================================================

    if not args.no_pdb_filter:

        pdb = (
            df["PDB"]
            .astype("string")
            .str.strip()
        )

        pdb = pdb.mask(
            pdb.isin(MISSING_STRINGS),
            pd.NA,
        )

        mapped = pdb.notna()

        df = df.loc[mapped].copy()

        print(f"Mapped rows: {len(df):,}")

    else:

        df = df.copy()

        print(f"Using all rows: {len(df):,}")

    # ========================================================
    # OPTIONAL SUBSAMPLING
    # ========================================================

    if args.sample_n and len(df) > args.sample_n:

        rng = np.random.default_rng(args.seed)

        take = rng.choice(
            df.index.to_numpy(),
            size=args.sample_n,
            replace=False,
        )

        df = df.loc[take].copy()

        print(f"Subsampled rows: {len(df):,}")

    # ========================================================
    # PREDICTORS
    # ========================================================

    df["SS3"] = ss3_from_dssp(
        df["SecondaryStructure"]
    )

    df["SS3"] = pd.Categorical(
        df["SS3"],
        categories=["Coil", "Helix", "Strand"],
    )

    df["Rsa_z"] = zscore(df["Rsa"])

    df["Depth_z"] = zscore(df[args.depth])

    df["log_tree_length"] = log1p_numeric(
        df["tree_length"]
    )

    df["log_num_sequences"] = log1p_numeric(
        df["num_sequences"]
    )

    df["log_tree_z"] = zscore(
        df["log_tree_length"]
    )

    df["log_numseq_z"] = zscore(
        df["log_num_sequences"]
    )

    # ========================================================
    # KEEP REQUIRED COLUMNS
    # ========================================================

    keep = [
        "is_rate_shift",
        "Orthogroup",
        "Last_Common_Ancestor",
        "SS3",
        "Rsa_z",
        "Depth_z",
        "log_tree_z",
        "log_numseq_z",
    ]

    df = df[keep].copy()

    df = df.dropna()

    print(
        f"Rows after filtering: {len(df):,}"
    )

    print(
        f"Orthogroups retained: "
        f"{df['Orthogroup'].nunique():,}"
    )

    print(
        f"LCA levels retained: "
        f"{df['Last_Common_Ancestor'].nunique():,}"
    )

    # ========================================================
    # MODEL FORMULA
    # ========================================================

    fixed_terms = [
        "C(Last_Common_Ancestor)",
        "Rsa_z",
        "Depth_z",
        "C(SS3)",
        "log_tree_z",
        "log_numseq_z",
    ]

    formula = (
        "is_rate_shift ~ " +
        " + ".join(fixed_terms)
    )

    vc_formulas = {
        "Orthogroup_RE":
            "0 + C(Orthogroup)"
    }

    run_meta = {
        "input": args.input,
        "depth_variable": args.depth,
        "pdb_filter":
            (not args.no_pdb_filter),
        "sample_n": args.sample_n,
        "seed": args.seed,
        "model":
            "LCAfixed_OGrandom",
        "random_effects":
            "Orthogroup intercept",
    }

    # ========================================================
    # FIT MODEL
    # ========================================================

    fit_and_write(
        df=df,
        formula=formula,
        vc_formulas=vc_formulas,
        outdir=args.outdir,
        run_meta=run_meta,
    )

    print(
        "\nDone. Results written to:",
        args.outdir,
    )


if __name__ == "__main__":
    main()
