# =============================================================================
# src/functions/central_functions.R
#
# ANALYSIS COMPONENT: DATA CLEANING AND TRANSFORMATION (function library; shared recodes, also used by the regression)
#
# Small shared helpers for the regression data preparation (step 2):
#   REGION_MAP          named vector: country -> WHO Africa sub-region
#   recode_hrp2_risk()  collapse HRP2 deletion-risk categories to the model grouping
#   binarise_peak()     binarise site_months_to_peak and report class balance
# =============================================================================

library(data.table)


# =============================================================================
# REGION_MAP
#
# Maps every study country to its WHO Africa sub-region. Used to assign
# `africa_region` in model and plot scripts. The full 32-country version is
# the canonical one; plot scripts that previously used a truncated 17-country
# version will silently return NA for countries not in their data, which is
# safe.
#
# NOTE: Two distinct hrp2_risk groupings are used in this codebase (see
# recode_hrp2_risk() below). Scripts 06 and 07 use the "model" scheme
# (Marginal + Moderate -> "Marginal/Moderate"); diagnostic plot scripts
# previously used the "plot" scheme (Slight + Marginal -> "Slight/Marginal").
# The difference changes the reference category. Confirm the intended grouping
# before running models.
# =============================================================================

REGION_MAP <- c(
  # Eastern and Southern Africa
  Kenya          = "Eastern and Southern Africa",
  Uganda         = "Eastern and Southern Africa",
  Rwanda         = "Eastern and Southern Africa",
  Burundi        = "Eastern and Southern Africa",
  Madagascar     = "Eastern and Southern Africa",
  Malawi         = "Eastern and Southern Africa",
  Angola         = "Eastern and Southern Africa",
  Mozambique     = "Eastern and Southern Africa",
  Zambia         = "Eastern and Southern Africa",
  Tanzania       = "Eastern and Southern Africa",
  Zimbabwe       = "Eastern and Southern Africa",
  Botswana       = "Eastern and Southern Africa",
  Ethiopia       = "Eastern and Southern Africa",
  Lesotho        = "Eastern and Southern Africa",
  "South Africa" = "Eastern and Southern Africa",
  # Western and Central Africa
  Nigeria                     = "Western and Central Africa",
  Ghana                       = "Western and Central Africa",
  Guinea                      = "Western and Central Africa",
  Liberia                     = "Western and Central Africa",
  "Sierra Leone"              = "Western and Central Africa",
  Togo                        = "Western and Central Africa",
  Mali                        = "Western and Central Africa",
  "Burkina Faso"              = "Western and Central Africa",
  Senegal                     = "Western and Central Africa",
  Niger                       = "Western and Central Africa",
  Cameroon                    = "Western and Central Africa",
  Benin                       = "Western and Central Africa",
  Chad                        = "Western and Central Africa",
  "Congo Democratic Republic" = "Western and Central Africa",
  "Cote d'Ivoire"             = "Western and Central Africa",
  Gabon                       = "Western and Central Africa",
  Gambia                      = "Western and Central Africa",
  Mauritania                  = "Western and Central Africa"
)


# =============================================================================
# load_merged()
#
# Loads dhs_pfhrp2_merged.csv and renames the RDT count columns to the generic
# names expected by the filtering and modelling pipeline.
#
# Returns a data.table with columns n_positive, n_tested, n_negative (in place
# of n_positive_rdt, n_tested_rdt, n_negative_rdt).
# =============================================================================


# =============================================================================
# filter_rdt_data()
#
# Standard two-stage filtering pipeline applied identically in scripts 06,
# plot_colinearity, plot_gap_diagnostics, and plot_loo_model1.
#
#   Stage 1 (rdt_base):  sub-regional rows with a valid RDT denominator
#   Stage 2 (rdt_micro): additionally requires a microscopy estimate
#   dat_full:            additionally requires fever prevalence (for Model 1)
#   dat_temporal:        rdt_micro without the fever requirement (for Model 2)
#
# year_c (survey_year centred on the mean of rdt_micro) is added in place.
#
# Returns a named list: list(full = dat_full, temporal = dat_temporal).
# For scripts that only need dat_full, use filter_rdt_data(merged)$full.
# =============================================================================


# =============================================================================
# recode_hrp2_risk()
#
# Collapses the four-level HRP2 risk variable into a three-level factor.
# Two grouping schemes are used in this codebase:
#
#   "model" (scripts 06, 07, plot_loo_model1):
#     Marginal + Moderate -> "Marginal/Moderate"
#     Levels: Marginal/Moderate (reference), Slight, High
#     Rationale: treats Slight as a distinct low-risk reference category.
#
#   "plot" (previously plot_colinearity, plot_gap_diagnostics):
#     Slight + Marginal -> "Slight/Marginal"
#     Levels: Slight/Marginal (reference), Moderate, High
#     Rationale: groups the two lowest-risk categories for visualisation.
#
# [UNCERTAINTY] The two schemes produce different reference categories and
# would give numerically different estimates if applied to model data. Check
# which scheme is intended before using in any modelling script.
# =============================================================================

recode_hrp2_risk <- function(x, scheme = c("model", "plot")) {
  scheme <- match.arg(scheme)
  if (scheme == "model") {
    factor(
      ifelse(x %in% c("Marginal", "Moderate"), "Marginal/Moderate", x),
      levels = c("Slight", "Marginal/Moderate", "High")
    )
  } else {
    factor(
      ifelse(x %in% c("Slight", "Marginal"), "Slight/Marginal", x),
      levels = c("Slight/Marginal", "Moderate", "High")
    )
  }
}


# =============================================================================
# extract_fe()
#
# Extracts fixed-effect odds ratios and Wald (or profile) CIs from a glmmTMB
# or lme4 model into a tidy data.table with standardised column names.
#
# Arguments:
#   fit          — a fitted glmmTMB or lme4 model object
#   conf.method  — passed to broom.mixed::tidy(); "Wald" (default, fast) or
#                  "profile" (slower, more accurate for small samples)
#   add_log_or   — if TRUE, appends a log_or column (used in script 06)
#
# Returns a data.table with columns: term, or, or_lo, or_hi, p_value.
# =============================================================================


# =============================================================================
# binarise_peak()
#
# Converts a continuous site_months_to_peak variable to a binary 0/1
# indicator (1 = survey during/after peak rainfall) and prints a class-balance
# check. Modifies dt in place by reference.
#
# Arguments:
#   dt          — a data.table (modified in place)
#   months_col  — name of the continuous source column (character)
#   new_col     — name of the new binary column to create (character)
#   threshold   — values <= threshold map to 1 (default 0)
#   verbose     — if TRUE, prints distribution table and balance warning
# =============================================================================

binarise_peak <- function(dt, months_col, new_col,
                          threshold = 0, verbose = TRUE) {
  dt[, (new_col) := as.integer(get(months_col) <= threshold)]
  if (verbose) {
    cat(sprintf("\n--- Distribution: %s (before binarising) ---\n", months_col))
    print(table(dt[[months_col]], useNA = "ifany"))
    cat("Quantiles:\n")
    print(quantile(dt[[months_col]], na.rm = TRUE))
    cat(sprintf("\n--- Seasonality: %s class balance ---\n", new_col))
    tab <- table(dt[[new_col]])
    print(tab)
    if (min(tab) < 0.05 * sum(tab))
      cat(sprintf(
        "WARNING: %s is highly imbalanced (< 5%% in one class) — the estimated effect may be unstable.\n",
        new_col
      ))
  }
  invisible(dt)
}
