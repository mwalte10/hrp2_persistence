# =============================================================================
# src/1_data_prep/02_survey_design_effects.R      (was src/02b_survey_design_effects.R)
#
# ANALYSIS COMPONENT: DATA CLEANING AND TRANSFORMATION (survey design effects)
#
# Step 1b. Runs between 01_format_dhs_mis_data.R and 03_structure_mis_data.R:
# 03 merges these design effects onto the site-level prevalence.
#
# Per-survey design effects (DEFF) for DHS/MIS RDT and microscopy prevalence,
# used by step 3a (grid_simulation.R Section F) to down-weight each survey's
# contribution to the observation likelihood (LIKELIHOOD = "binomial_deff" /
# "multinomial_deff"). MIS surveys sample children in clusters, and malaria is
# focal, so a cluster sample of n children carries the information of roughly
# n / DEFF independent children -- the original likelihood treated every survey
# as a simple random sample (DEFF = 1).
#
# Inputs (01_format_dhs_mis_data.R, one row per country x survey year x DHS region):
#   dhs_mis_rdt.csv / dhs_mis_microscopy.csv -- raw unweighted n_tested /
#     n_positive, plus the survey-weighted prevalence SE from .make_design()
#     (clusters hv001 + strata hv022 + weights hv005).
#   dhs_mis_concordance.csv -- RDT x microscopy 2x2 for children with both
#     tests, incl. Kish's effective size of the weights (n_eff_both_tested).
#
# Per region and test:
#   deff_total   = se^2 / (p(1-p)/n)    -- DEFF of the WEIGHTED estimator
#   deff_kish    = n_both / n_eff_both  -- part due to unequal weights alone
#   deff_cluster = deff_total / deff_kish
# Kish's decomposition (DEFF ~ DEFF_weighting x DEFF_clustering): the
# likelihood uses RAW UNWEIGHTED counts, so the weighting component doesn't
# apply to it -- only the clustering component does. RDT-only surveys have no
# concordance row and so no Kish factor; they get the median Kish factor
# across surveys (logged).
#
# Per survey: median of its regions' deff_cluster (regional estimates rest on
# few clusters and are noisy; the median is robust to the odd extreme
# region), floored at 1. Regions with p = 0 or 1 (DEFF undefined) or
# n < MIN_N are dropped first; national rows are excluded (they'd double-count
# the regions). deff_joint (for the 2x2 multinomial) = mean(deff_rdt, deff_lm).
#
# Writes DIR_STRUCTURE:
#   survey_design_effects.RDS            one row per survey (read by step 1c)
# =============================================================================
suppressMessages(library(data.table))
source("./src/functions/paths.R")

IN_DIR  <- DIR_PULL
OUT_DIR <- DIR_STRUCTURE
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MIN_N   <- 20L

rdt  <- fread(file.path(IN_DIR, "dhs_mis_rdt.csv"))
lm   <- fread(file.path(IN_DIR, "dhs_mis_microscopy.csv"))
conc <- fread(file.path(IN_DIR, "dhs_mis_concordance.csv"))

is_national <- function(x) is.na(x) | x == "" | tolower(x) == "national"

region_deff <- function(dt, prev_col, se_col, test) {
  dt <- dt[!is_national(sub_region)]
  dt[, .(country, survey_year, sub_region, test = test,
         n = n_tested, p = get(prev_col), se = get(se_col))]
}
reg <- rbind(region_deff(rdt, "prev_rdt", "se_prev_rdt", "rdt"),
             region_deff(lm,  "prev_micro", "se_prev_micro", "lm"))

kish <- conc[!is_national(sub_region) & n_eff_both_tested > 0,
             .(country, survey_year, sub_region, deff_kish = n_both_tested / n_eff_both_tested)]
reg <- merge(reg, kish, by = c("country", "survey_year", "sub_region"), all.x = TRUE)
median_kish <- median(kish$deff_kish, na.rm = TRUE)
n_kish_fill <- reg[is.na(deff_kish), uniqueN(paste(country, survey_year))]
reg[, kish_imputed := is.na(deff_kish)]
reg[is.na(deff_kish), deff_kish := median_kish]

reg[, usable := p > 0 & p < 1 & n >= MIN_N & !is.na(se) & se > 0]
reg[usable == TRUE, `:=`(deff_total   = se^2 / (p * (1 - p) / n))]
reg[usable == TRUE, `:=`(deff_cluster = deff_total / deff_kish)]

survey <- reg[usable == TRUE, .(
  deff = max(1, median(deff_cluster)),
  deff_raw_median = median(deff_cluster),
  n_regions = .N,
  deff_kish = median(deff_kish),
  kish_imputed = all(kish_imputed)
), by = .(country, survey_year, test)]
survey_w <- dcast(survey, country + survey_year ~ test,
                  value.var = c("deff", "deff_raw_median", "n_regions", "deff_kish", "kish_imputed"))
survey_w[, deff_joint := rowMeans(cbind(deff_rdt, deff_lm), na.rm = FALSE)]
setorder(survey_w, country, survey_year)

saveRDS(survey_w, file.path(OUT_DIR, "survey_design_effects.RDS"))

# --- Distribution report ---
q <- function(x) sprintf("median %.2f (IQR %.2f-%.2f, range %.2f-%.2f), n = %d",
                         median(x), quantile(x, 0.25), quantile(x, 0.75), min(x), max(x), length(x))
cat(sprintf("Regions usable: %d of %d (dropped p in {0,1}, n < %d, missing SE); national rows excluded.\n",
            reg[usable == TRUE, .N], nrow(reg), MIN_N))
cat(sprintf("Kish weighting factor (n / n_eff): %s\n", q(kish$deff_kish)))
cat(sprintf("  %d survey(s) without a concordance row given the median Kish factor %.2f.\n", n_kish_fill, median_kish))
cat(sprintf("Regional cluster DEFF, RDT: %s\n", q(reg[usable == TRUE & test == "rdt", deff_cluster])))
cat(sprintf("Regional cluster DEFF, LM:  %s\n", q(reg[usable == TRUE & test == "lm",  deff_cluster])))
cat(sprintf("Per-survey DEFF (floored at 1), RDT: %s\n", q(survey_w[!is.na(deff_rdt), deff_rdt])))
cat(sprintf("Per-survey DEFF (floored at 1), LM:  %s\n", q(survey_w[!is.na(deff_lm),  deff_lm])))
cat(sprintf("Per-survey DEFF, joint (mean of RDT & LM): %s\n", q(survey_w[!is.na(deff_joint), deff_joint])))
cat(sprintf("Surveys floored at 1: RDT %d, LM %d\n",
            survey_w[deff_raw_median_rdt < 1, .N], survey_w[deff_raw_median_lm < 1, .N]))

cat(sprintf("Saved survey_design_effects.RDS to %s\n", OUT_DIR))
