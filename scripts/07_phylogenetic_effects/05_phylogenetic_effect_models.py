#!/usr/bin/env python3

# ============================================================
# 05_phylogenetic_effect_models.py
#
# Linear models for:
#   1. Robinson-Foulds distances
#   2. Tree-length differences
#   3. Bootstrap-support differences
#
# Inputs:
#   rf_results.tsv
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

RF_FILE = "rf_results.tsv"

TREE_LENGTH_FILE = "delta_tree_length.tsv"

SUPPORT_FILE = "delta_support.tsv"

META_FILE = "lm_input_50_families.tsv"

# ============================================================
# OUTPUT PREFIXES
# ============================================================

RF_PREFIX = "lm_rf"

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
# RF DISTANCES
# ============================================================

print("\n================================================")
print("RF DISTANCE MODELS")
print("================================================")

rf = pd.read_csv(RF_FILE, sep="\t")

rand = rf[rf["type"] == "rand"]

nors = rf[rf["type"] == "nors"]

# ------------------------------------------------------------
# mean random RF
# ------------------------------------------------------------

rand_mean = (
    rand.groupby("Orthogroup")["nRF"]
    .mean()
    .reset_index()
)

rand_mean = rand_mean.rename(
    columns={"nRF": "mean_rand"}
)

# ------------------------------------------------------------
# nors RF
# ------------------------------------------------------------

nors = nors[[
    "Orthogroup",
    "nRF"
]].rename(
    columns={"nRF": "nRF_nors"}
)

# ------------------------------------------------------------
# merge
# ------------------------------------------------------------

rf_df = pd.merge(
    nors,
    rand_mean,
    on="Orthogroup"
)

rf_df["delta_RF"] = (
    rf_df["nRF_nors"] -
    rf_df["mean_rand"]
)

rf_df = pd.merge(
    rf_df,
    meta,
    on="Orthogroup"
)

# ------------------------------------------------------------
# save merged dataset
# ------------------------------------------------------------

rf_dataset_out = (
    f"{RF_PREFIX}_dataset.tsv"
)

rf_df.to_csv(
    rf_dataset_out,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# linear model
# ------------------------------------------------------------

rf_model = smf.ols(
    "delta_RF ~ "
    "log_num_sequences + "
    "log_num_sites + "
    "prop_rate_shift",
    data=rf_df
).fit()

print(rf_model.summary())

# ------------------------------------------------------------
# outputs
# ------------------------------------------------------------

save_model_results(
    rf_model,
    f"{RF_PREFIX}_results.tsv"
)

make_diagnostics(
    rf_model,
    f"{RF_PREFIX}_diagnostics.png"
)

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
