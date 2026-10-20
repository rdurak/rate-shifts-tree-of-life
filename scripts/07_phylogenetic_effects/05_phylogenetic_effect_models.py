#!/usr/bin/env python3

# ============================================================
# lm_lrm_bl_bs.py
#
# Linear models for:
#   1. Lin-Rajan-Moret matching-distance differences
#   2. Tree-length differences
#   3. Bootstrap-support differences
#
# Inputs:
#   lrm_permutation_summary.tsv
#   delta_tree_length.tsv
#   delta_support.tsv
#   lm_input_50_families.tsv
#
# Outputs:
#   merged datasets
#   regression tables
#   diagnostic plots
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

# ============================================================
# INPUT FILES
# ============================================================

LRM_FILE = "lrm_permutation_summary.tsv"

TREE_LENGTH_FILE = "delta_tree_length.tsv"

SUPPORT_FILE = "delta_support.tsv"

META_FILE = "lm_input_50_families.tsv"

# ============================================================
# OUTPUT PREFIXES
# ============================================================

LRM_PREFIX = "lm_lrm"

TREE_PREFIX = "lm_tree_length"

SUPPORT_PREFIX = "lm_support"

# ============================================================
# LOAD METADATA
# ============================================================

meta = pd.read_csv(META_FILE, sep="\t")

meta = meta[[
    "Orthogroup",
    "num_sites",
    "num_sequences",
    "prop_rate_shift"
]]

# ============================================================
# LOG TRANSFORMS
# ============================================================

meta["log_num_sequences"] = np.log10(
    meta["num_sequences"]
)

meta["log_num_sites"] = np.log10(
    meta["num_sites"]
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def save_model_results(model, outfile):

    results = pd.DataFrame({
        "predictor": model.params.index,
        "estimate": model.params.values,
        "std_error": model.bse.values,
        "p_value": model.pvalues.values,
        "CI_lower": model.conf_int()[0].values,
        "CI_upper": model.conf_int()[1].values
    })

    results.to_csv(
        outfile,
        sep="\t",
        index=False
    )


def make_diagnostics(model, outfile):

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10, 4)
    )

    # --------------------------------------------------------
    # residuals vs fitted
    # --------------------------------------------------------

    axes[0].scatter(
        model.fittedvalues,
        model.resid,
        color="0.4",
        s=20
    )

    axes[0].axhline(
        0,
        linestyle="--",
        color="black"
    )

    axes[0].set_xlabel("Fitted values")

    axes[0].set_ylabel("Residuals")

    axes[0].set_title(
        "Residuals vs Fitted"
    )

    # --------------------------------------------------------
    # residual histogram
    # --------------------------------------------------------

    axes[1].hist(
        model.resid,
        bins=10,
        color="0.8",
        edgecolor="black"
    )

    axes[1].set_xlabel("Residuals")

    axes[1].set_title(
        "Residual distribution"
    )

    plt.tight_layout()

    plt.savefig(outfile, dpi=300)

    plt.close()

# ============================================================
# PART 1
# LRM DISTANCES: one response per family
# ============================================================

print("\nLRM DISTANCE MODEL")
lrm = pd.read_csv(LRM_FILE, sep="\t")
required = {"Orthogroup", "Family", "LRM_nors", "mean_random_LRM", "delta_LRM"}
missing = required.difference(lrm.columns)
if missing:
    raise ValueError(f"Missing columns: {sorted(missing)}")
if lrm.empty:
    raise ValueError("LRM summary is empty.")
if lrm["Orthogroup"].isna().any() or meta["Orthogroup"].isna().any():
    raise ValueError("Missing orthogroup identifiers.")
if lrm["Orthogroup"].duplicated().any() or meta["Orthogroup"].duplicated().any():
    raise ValueError("Expected exactly one row per orthogroup in summary and metadata.")
if not np.allclose(lrm["delta_LRM"], lrm["LRM_nors"] - lrm["mean_random_LRM"]):
    raise ValueError("delta_LRM does not equal LRM_nors minus mean_random_LRM.")

lrm_df = lrm.merge(meta, on="Orthogroup", how="left", validate="one_to_one")
model_columns = ["delta_LRM", "log_num_sequences", "log_num_sites", "prop_rate_shift"]
if not np.isfinite(lrm_df[model_columns].to_numpy(dtype=float)).all():
    raise ValueError("Missing or nonfinite response/predictors; check metadata matches and log inputs.")

lrm_df.to_csv(f"{LRM_PREFIX}_dataset.tsv", sep="\t", index=False)
lrm_model = smf.ols(
    "delta_LRM ~ log_num_sequences + log_num_sites + prop_rate_shift",
    data=lrm_df,
    missing="raise",
).fit()
print(lrm_model.summary())
save_model_results(lrm_model, f"{LRM_PREFIX}_results.tsv")
make_diagnostics(lrm_model, f"{LRM_PREFIX}_diagnostics.png")

# ============================================================
# PART 2
# TREE LENGTH
# ============================================================

print("\n================================================")
print("TREE LENGTH MODELS")
print("================================================")

tree = pd.read_csv(
    TREE_LENGTH_FILE,
    sep="\t"
)

tree = tree[[
    "Orthogroup",
    "delta_tree_length"
]]

tree_df = pd.merge(
    tree,
    meta,
    on="Orthogroup"
)

# ------------------------------------------------------------
# save merged dataset
# ------------------------------------------------------------

tree_dataset_out = (
    f"{TREE_PREFIX}_dataset.tsv"
)

tree_df.to_csv(
    tree_dataset_out,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# linear model
# ------------------------------------------------------------

tree_model = smf.ols(
    "delta_tree_length ~ "
    "log_num_sequences + "
    "log_num_sites + "
    "prop_rate_shift",
    data=tree_df
).fit()

print(tree_model.summary())

# ------------------------------------------------------------
# outputs
# ------------------------------------------------------------

save_model_results(
    tree_model,
    f"{TREE_PREFIX}_results.tsv"
)

make_diagnostics(
    tree_model,
    f"{TREE_PREFIX}_diagnostics.png"
)

# ============================================================
# PART 3
# BOOTSTRAP SUPPORT
# ============================================================

print("\n================================================")
print("BOOTSTRAP SUPPORT MODELS")
print("================================================")

support = pd.read_csv(
    SUPPORT_FILE,
    sep="\t"
)

support = support[[
    "Orthogroup",
    "delta_support"
]]

support_df = pd.merge(
    support,
    meta,
    on="Orthogroup"
)

# ------------------------------------------------------------
# save merged dataset
# ------------------------------------------------------------

support_dataset_out = (
    f"{SUPPORT_PREFIX}_dataset.tsv"
)

support_df.to_csv(
    support_dataset_out,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# linear model
# ------------------------------------------------------------

support_model = smf.ols(
    "delta_support ~ "
    "log_num_sequences + "
    "log_num_sites + "
    "prop_rate_shift",
    data=support_df
).fit()

print(support_model.summary())

# ------------------------------------------------------------
# outputs
# ------------------------------------------------------------

save_model_results(
    support_model,
    f"{SUPPORT_PREFIX}_results.tsv"
)

make_diagnostics(
    support_model,
    f"{SUPPORT_PREFIX}_diagnostics.png"
)

# ============================================================
# DONE
# ============================================================

print("\nFinished.")

