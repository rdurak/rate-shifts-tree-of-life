#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Residual-based partial correlation analyses
#
# Input:
#   sged_summary_per_orthogroup.with_tree_aln.tsv
#
# Analyses:
#   - Original models
#   - Log-transformed models
#   - Partial Kendall correlations
#   - Diagnostic plots
# ------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(readr)
}))

# ============================================================
# USER SETTINGS
# ============================================================

infile <- "sged_summary_per_orthogroup.with_tree_aln.tsv"

out_prefix <- "partial_correlation"

summary_out <- paste0(out_prefix, "_results.txt")

diagnostic_pdf_original <- paste0(
  out_prefix,
  "_diagnostics_original.pdf"
)

diagnostic_pdf_log <- paste0(
  out_prefix,
  "_diagnostics_log_transformed.pdf"
)

# ============================================================
# LOAD DATA
# ============================================================

df <- read_tsv(infile, show_col_types = FALSE)

# ============================================================
# REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  "prop_rate_shift",
  "num_sequences",
  "tree_length",
  "Last_Common_Ancestor"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ============================================================
# CLEAN DATA
# ============================================================

df$prop_rate_shift <- as.numeric(df$prop_rate_shift)
df$num_sequences   <- as.numeric(df$num_sequences)
df$tree_length     <- as.numeric(df$tree_length)

df <- df[
  is.finite(df$prop_rate_shift) &
  is.finite(df$num_sequences) &
  is.finite(df$tree_length),
]

df$Last_Common_Ancestor <- factor(df$Last_Common_Ancestor)

# ============================================================
# LOG TRANSFORM VARIABLES
# ============================================================

df$log_num_sequences <- log10(df$num_sequences + 1)
df$log_tree_length   <- log10(df$tree_length + 1)

# ============================================================
# ORIGINAL MODELS
# ============================================================

# -------------------------
# Model 1
# prop_rate_shift ~ LCA
# -------------------------

model1_orig <- lm(
  prop_rate_shift ~ Last_Common_Ancestor,
  data = df
)

# -------------------------
# Model 2
# num_sequences ~ LCA
# -------------------------

model2_orig <- lm(
  num_sequences ~ Last_Common_Ancestor,
  data = df
)

# -------------------------
# Model 3
# prop_rate_shift ~ LCA + num_sequences
# -------------------------

model3_orig <- lm(
  prop_rate_shift ~
    Last_Common_Ancestor +
    num_sequences,
  data = df
)

# -------------------------
# Model 4
# tree_length ~ LCA + num_sequences
# -------------------------

model4_orig <- lm(
  tree_length ~
    Last_Common_Ancestor +
    num_sequences,
  data = df
)

# ============================================================
# LOG-TRANSFORMED MODELS
# ============================================================

# -------------------------
# Model 1
# unchanged
# -------------------------

model1_log <- lm(
  prop_rate_shift ~ Last_Common_Ancestor,
  data = df
)

# -------------------------
# Model 2
# log_num_sequences ~ LCA
# -------------------------

model2_log <- lm(
  log_num_sequences ~ Last_Common_Ancestor,
  data = df
)

# -------------------------
# Model 3
# prop_rate_shift ~ LCA + log_num_sequences
# -------------------------

model3_log <- lm(
  prop_rate_shift ~
    Last_Common_Ancestor +
    log_num_sequences,
  data = df
)

# -------------------------
# Model 4
# log_tree_length ~ LCA + log_num_sequences
# -------------------------

model4_log <- lm(
  log_tree_length ~
    Last_Common_Ancestor +
    log_num_sequences,
  data = df
)

# ============================================================
# PARTIAL CORRELATION TESTS
# ============================================================

# ------------------------------------------------------------
# TEST 1
# Does number of sequences affect rate-shift prevalence
# independently of ancestry?
# ------------------------------------------------------------

resid_rs_1 <- residuals(model1_log)

resid_seq_1 <- residuals(model2_log)

test1 <- cor.test(
  resid_rs_1,
  resid_seq_1,
  method = "kendall"
)

# ------------------------------------------------------------
# TEST 2
# Does tree length affect rate-shift prevalence
# independently of ancestry and sampling?
# ------------------------------------------------------------

resid_rs_2 <- residuals(model3_log)

resid_tree_2 <- residuals(model4_log)

test2 <- cor.test(
  resid_rs_2,
  resid_tree_2,
  method = "kendall"
)

# ============================================================
# WRITE RESULTS
# ============================================================

sink(summary_out)

cat("==================================================\n")
cat("PARTIAL CORRELATION ANALYSES\n")
cat("==================================================\n\n")

cat("Input file:\n")
cat(infile, "\n\n")

cat("Orthogroups analyzed:\n")
cat(nrow(df), "\n\n")

# ------------------------------------------------------------
# TEST 1
# ------------------------------------------------------------

cat("==================================================\n")
cat("TEST 1\n")
cat("Rate-shift prevalence vs number of sequences\n")
cat("Correcting for Last_Common_Ancestor\n")
cat("==================================================\n\n")

print(test1)

cat("\n\n")

# ------------------------------------------------------------
# TEST 2
# ------------------------------------------------------------

cat("==================================================\n")
cat("TEST 2\n")
cat("Rate-shift prevalence vs tree length\n")
cat("Correcting for Last_Common_Ancestor and sampling\n")
cat("==================================================\n\n")

print(test2)

cat("\n\n")

# ------------------------------------------------------------
# MODEL SUMMARIES
# ------------------------------------------------------------

cat("==================================================\n")
cat("ORIGINAL MODELS\n")
cat("==================================================\n\n")

print(summary(model1_orig))
cat("\n\n")

print(summary(model2_orig))
cat("\n\n")

print(summary(model3_orig))
cat("\n\n")

print(summary(model4_orig))
cat("\n\n")

cat("==================================================\n")
cat("LOG-TRANSFORMED MODELS\n")
cat("==================================================\n\n")

print(summary(model1_log))
cat("\n\n")

print(summary(model2_log))
cat("\n\n")

print(summary(model3_log))
cat("\n\n")

print(summary(model4_log))
cat("\n\n")

sink()

cat("Wrote results to:", summary_out, "\n")

# ============================================================
# DIAGNOSTIC PLOTS
# ============================================================

# ------------------------------------------------------------
# ORIGINAL MODELS
# ------------------------------------------------------------

pdf(diagnostic_pdf_original, width = 9, height = 7)

par(mfrow = c(2, 2))

plot(model1_orig)
plot(model2_orig)
plot(model3_orig)
plot(model4_orig)

dev.off()

# ------------------------------------------------------------
# LOG-TRANSFORMED MODELS
# ------------------------------------------------------------

pdf(diagnostic_pdf_log, width = 9, height = 7)

par(mfrow = c(2, 2))

plot(model1_log)
plot(model2_log)
plot(model3_log)
plot(model4_log)

dev.off()

par(mfrow = c(1, 1))

cat("Wrote diagnostic plots\n")

# ============================================================
# DONE
# ============================================================

cat("\nDone.\n")
