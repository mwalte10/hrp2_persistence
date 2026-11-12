# =============================================================================
# src/functions/brms_int_cluster_fns.R
#
# ANALYSIS COMPONENT: REGRESSION ANALYSIS (function library; also used by
# src/plots.R and src/results_numbers.Rmd)
#
# Data preparation, formula and priors for the step 2 regression: the joint
# microscopy x RDT outcome (p_p / n_p / p_n / n_n) of children tested by both,
#   outcome ~ sr_rdt_cat * hrp2_risk + has_smc + during_or_after_peak + age_c
#             + care_or_act + fever + (1 | survid),  family = categorical.
# Regularising priors (normal(0, 2.5) on fixed effects, normal(0, 5) on
# intercepts) keep the sparse sr_rdt_cat x hrp2_risk cells from separating.
# Also holds the extended formula and data used for fixed-effect selection.
# The fever-only formula and the hipercow chain functions are in
# archive/src/functions/brms_int_cluster_fns_removed.R.
# =============================================================================

library(data.table)
library(brms)
source("./src/functions/paths.R")
source("./src/functions/central_functions.R")

# --- Use cmdstanr, not rstan, as the Stan backend. The cluster's own task
# environment sets CMDSTAN (pointing at an existing CmdStan install) and
# CMDSTANR_USE_RTOOLS -- both cmdstanr-specific -- indicating the cluster is
# provisioned for cmdstanr rather than rstan. Confirmed necessary: an rstan-
# backend chain task failed with "Error: invalid connection" immediately
# after a successful-looking g++/make compile, a known class of Windows
# rstan issue tied to R's pipe()/connection handling under redirected
# (non-interactive/batch) stdout -- exactly how hipercow runs tasks.
# cmdstanr avoids that path entirely by compiling/running via the standalone
# CmdStan binary instead of an in-process R DLL compile. Guarded so this
# file still loads (e.g. for local interactive use, where cmdstanr isn't
# installed) if CMDSTAN isn't set or cmdstanr isn't available.
if (nzchar(Sys.getenv("CMDSTAN")) && requireNamespace("cmdstanr", quietly = TRUE)) {
  cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
}

# --- Model formula, shared by every phase/chain. Identical to the
# regularised-prior interaction test that converged cleanly on the 3000-row
# subsample (see test_brms_switch_int_regprior.R). ---
BRMS_INT_FORMULA <- outcome ~ sr_rdt_cat * hrp2_risk + has_smc + during_or_after_peak +
  age_c + care_or_act + fever + (1 | survid)

# Fixed-effect selection (2_regression/fixed_effect_selection.R; appendix
# Table tab:elpd): the "full" model adds African region, its interaction with
# deletion risk, and survey year. The main model is this formula less region,
# region x deletion risk and year.
BRMS_INT_FORMULA_EXT <- outcome ~ sr_rdt_cat * hrp2_risk + hrp2_risk * africa_region +
  has_smc + during_or_after_peak + age_c + care_or_act + fever + year_c + (1 | survid)

# =============================================================================
# prepare_brms_full_data() -- PHASE 1 (prep task)
#
# Replicates 03_alt.R's / 03_model_rdt_individual.R's standard data-prep
# pipeline (same steps as the subsample smoke tests this session), on the
# FULL dataset -- no subsampling. Returns the individual-level data.table
# ready to pass straight to brm().
# =============================================================================
prepare_brms_full_data <- function(mis_dir = DIR_STRUCTURE) {
  ind <- fread(file.path(mis_dir, "dhs_mis_individual_hrp2.csv"))
  include <- readRDS(file.path(mis_dir, "countries_included.RDS"))
  ind <- ind[country %in% include$country, ]
  ind <- ind[!is.na(rdt_result) & !is.na(micro_result) & !is.na(fever)]
  ind[, survid := paste0(country, "_", survey_year)]
  age_mean <- mean(ind$age_months)
  ind[, age_c := age_months - age_mean]
  ind[, sr_rdt_cat := cut(prev_rdt, breaks = c(-Inf, 0.10, 0.20, Inf),
                           labels = c("<10%", "10-20%", ">=20%"), right = FALSE)]
  ind[, sr_rdt_cat := relevel(sr_rdt_cat, ref = "<10%")]
  SMC_COV_THRESHOLD <- 0.01

  # Uganda-specific override, NOT a general NA -> 0 rule -- see
  # 03_model_rdt_individual.R's identical override for the full rationale:
  # 02_structure_mis_data.R deliberately leaves site_smc_cov NA for years
  # before a country's SMC programme existed in its own source data (a
  # considered design choice, left alone for every other country here).
  # Uganda is a known exception -- its SMC source data only has real
  # coverage values in 2024 (0-96.8%), and its real-world programme is
  # independently known to be recent and Karamoja-only (2023/2024 on),
  # so earlier years (here: 2018) genuinely had 0% coverage. Without this,
  # Uganda_2018's entire has_smc column is NA, and since has_smc is a
  # fixed-effect term in BRMS_INT_FORMULA/BRMS_FEVER_FORMULA, brm()'s own
  # complete-case handling silently drops all 6,430+ of its rows (confirmed
  # directly this session, via two full refits that both lost this survey
  # despite prepare_brms_full_data() itself returning it intact).
  ind[country == "Uganda" & survey_year < 2023 & is.na(site_smc_cov), site_smc_cov := 0]

  ind[, has_smc := as.integer(site_smc_cov > SMC_COV_THRESHOLD)]
  binarise_peak(ind, "site_months_to_peak", "during_or_after_peak")
  ind[, care_or_act := ifelse(formal_care == 1 | act_treatment == 1, 1, 0)]
  ind[, hrp2_risk := recode_hrp2_risk(hrp2_risk, "model")]
  ind[, micro_f := factor(micro_result, levels = c(0, 1), labels = c("n", "p"))]
  ind[, rdt_f := factor(rdt_result, levels = c(0, 1), labels = c("n", "p"))]
  ind_all_tx <- ind[!is.na(act_treatment) & !is.na(formal_care)]
  ind_all_tx[, outcome := paste0(micro_f, "_", rdt_f)]
  ind_all_tx[, outcome := factor(outcome, levels = c("p_p", "n_p", "p_n", "n_n"))]
  ind_all_tx
}

# =============================================================================
# build_brms_regularising_prior() -- PHASE 1 (prep task)
#
# Same override used in the validated subsample test: normal(0, 2.5) on
# every fixed-effect coefficient (class "b", all three dpars -- munp/mupn/
# munn given refcat = "p_p"), normal(0, 5) on each dpar's intercept. Weakly-
# informative on the logit scale (2.5 SD already covers ~12x odds ratios);
# purely there to give the sparse sr_rdt_cat x hrp2_risk interaction cells a
# proper posterior mode instead of running to +-Inf under separation.
# =============================================================================
build_brms_regularising_prior <- function(formula, data) {
  p <- get_prior(formula, data = data, family = categorical(refcat = "p_p"))
  p$prior[p$class == "b"] <- "normal(0, 2.5)"
  p$prior[p$class == "Intercept"] <- "normal(0, 5)"
  p
}

# =============================================================================
# prepare_brms_full_data_extended() -- prepare_brms_full_data() plus the two
# columns the fixed-effect selection's full model needs: africa_region
# (REGION_MAP, central_functions.R) and survey year centred on its mean.
# =============================================================================
prepare_brms_full_data_extended <- function(mis_dir = DIR_STRUCTURE) {
  ind <- prepare_brms_full_data(mis_dir)
  ind[, africa_region := REGION_MAP[country]]
  ind[, year_c := survey_year - mean(survey_year)]
  if (anyNA(ind$africa_region)) {
    stop(sprintf("prepare_brms_full_data_extended(): no REGION_MAP entry for %s",
                 paste(unique(ind[is.na(africa_region), country]), collapse = ", ")))
  }
  ind
}
