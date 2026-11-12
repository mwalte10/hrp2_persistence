# =============================================================================
# src/2_regression/fixed_effect_selection.R
#   (was archive/src/03j_brms_int_ext_full_fit.R, 03g_brms_int_loo_covariate_refit.R
#    and 03h_brms_int_loo_compare.R)
#
# ANALYSIS COMPONENT: REGRESSION ANALYSIS (fixed-effect selection)
#
# Leave-one-out cross-validation (LOO) used to choose the fixed effects of the
# main regression model (appendix Section "Fixed effect selection", Table
# tab:elpd). The full model (BRMS_INT_FORMULA_EXT, brms_int_cluster_fns.R)
# adds African region, region x deletion risk and survey year to the main
# model. Each candidate model drops one term (region is dropped together with
# its interaction, keeping marginality); "main_model" drops region,
# region x deletion risk and year, giving the main model's fixed effects.
#
# Every model uses the same children, priors, family and sampler settings
# (2 chains x 300 iterations, 150 warmup), so ELPD differences are comparable.
# The comparison fits use fewer iterations than the main fit (Section A of
# fit_fdr_fpr_regression.R) and are used only for model choice.
#
# Sections
#   A  fit each model in MODELS and save its loo() (hours per model; run one
#      model per R process to use more cores, e.g.
#        Rscript run_sections.R src/2_regression/fixed_effect_selection.R \
#          MODELS='"age_c"' COMPARE=FALSE)
#   B  compare each reduced model against the full model
#
# Outputs
#   DIR_REGRESSION/loo_covariate_selection/loo_brms_int_ext_<model>.RDS
#   RESULTS_DIR/table_loo_fixed_effect_selection.csv   (appendix Table tab:elpd)
# =============================================================================
library(data.table)
library(brms)
library(loo)
source("./src/functions/paths.R")
source("./src/functions/pipeline_fns.R")
source("./src/functions/brms_int_cluster_fns.R")

LOO_DIR <- file.path(DIR_REGRESSION, "loo_covariate_selection")
dir.create(LOO_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Candidate models: update() terms applied to BRMS_INT_FORMULA_EXT
MODEL_DROPS <- list(
  full                 = . ~ .,
  interaction          = . ~ . - sr_rdt_cat:hrp2_risk,
  africa_region_full   = . ~ . - africa_region - hrp2_risk:africa_region,
  year_c               = . ~ . - year_c,
  care_or_act          = . ~ . - care_or_act,
  age_c                = . ~ . - age_c,
  during_or_after_peak = . ~ . - during_or_after_peak,
  has_smc              = . ~ . - has_smc,
  fever                = . ~ . - fever,
  main_model           = . ~ . - africa_region - hrp2_risk:africa_region - year_c
)
MODEL_LABELS <- c(
  full                 = "Full model",
  year_c               = "- Survey year",
  africa_region_full   = "- African region and pfhrp2 risk x region",
  interaction          = "- Sub-regional prevalence x pfhrp2 risk",
  during_or_after_peak = "- During or after peak season",
  has_smc              = "- SMC",
  care_or_act          = "- Formal care or ACT",
  fever                = "- Fever",
  age_c                = "- Age",
  main_model           = "Main model (- region, - year)"
)

# Settings ======================================================================
OVR <- take_overrides()
MODELS  <- setting(OVR, "MODELS", names(MODEL_DROPS))   # Section A: models to fit
FIT     <- setting(OVR, "FIT", FALSE)                   # Section A: refit (hours per model)
COMPARE <- setting(OVR, "COMPARE", TRUE)                # Section B: ELPD table from saved loo() results
N_CHAINS       <- 2L
CHAIN_N_ITER   <- 300L
CHAIN_N_WARMUP <- 150L
SEED           <- 1L

# Fit ==========================================================================
# SECTION A -- ONE FIT AND loo() PER MODEL
if (FIT) {
  data_full <- prepare_brms_full_data_extended(DIR_STRUCTURE)
  cat(sprintf("Full data: %d rows, %d survid groups\n", nrow(data_full), uniqueN(data_full$survid)))
  for (m in MODELS) {
    if (!m %in% names(MODEL_DROPS)) stop(sprintf("Unknown model '%s'", m))
    form <- update(BRMS_INT_FORMULA_EXT, MODEL_DROPS[[m]])
    cat(sprintf("\n=== %s ===\n", m)); print(form)
    t0 <- Sys.time()
    fit <- brm(form, data = data_full, family = categorical(refcat = "p_p"),
               prior = build_brms_regularising_prior(form, data_full),
               chains = N_CHAINS, iter = CHAIN_N_ITER, warmup = CHAIN_N_WARMUP, cores = N_CHAINS,
               seed = SEED, backend = "rstan", refresh = 100)
    cat(sprintf("Fit time: %.2f hours\n", as.numeric(Sys.time() - t0, units = "hours")))
    loo_m <- loo(fit)
    print(loo_m)
    saveRDS(loo_m, file.path(LOO_DIR, sprintf("loo_brms_int_ext_%s.RDS", m)))
    rm(fit); gc()
  }
}

# Compare ======================================================================
# SECTION B -- ELPD AND DELTA-ELPD AGAINST THE FULL MODEL
# loo_compare() reports differences relative to whichever model ranks first,
# so the sign of (reduced - full) is resolved explicitly.
# File names used before November 2026: loo_brms_int_ext_drop_<term>.RDS for
# the reduced models and loo_brms_int_full.RDS for the main model.
loo_path <- function(m) {
  new <- file.path(LOO_DIR, sprintf("loo_brms_int_ext_%s.RDS", m))
  old <- if (m == "main_model") file.path(LOO_DIR, "loo_brms_int_full.RDS")
         else file.path(LOO_DIR, sprintf("loo_brms_int_ext_drop_%s.RDS", m))
  if (file.exists(new)) new else old
}
if (COMPARE) {
  if (!file.exists(loo_path("full"))) stop("No loo() result for the full model -- run Section A (FIT = TRUE).")
  loo_full <- readRDS(loo_path("full"))
  tab <- rbindlist(lapply(names(MODEL_LABELS), function(m) {
    if (!file.exists(loo_path(m))) { cat(sprintf("Missing loo() result: %s\n", m)); return(NULL) }
    loo_m <- readRDS(loo_path(m))
    if (m == "full") { d <- 0; se_d <- 0 } else {
      cmp <- loo_compare(list(full = loo_full, reduced = loo_m))
      if (rownames(cmp)[1] == "full") { d <- cmp["reduced", "elpd_diff"]; se_d <- cmp["reduced", "se_diff"] }
      else                            { d <- -cmp["full", "elpd_diff"];   se_d <- cmp["full", "se_diff"] }
    }
    data.table(model = m, label = MODEL_LABELS[[m]],
               elpd = loo_m$estimates["elpd_loo", "Estimate"], elpd_se = loo_m$estimates["elpd_loo", "SE"],
               delta_elpd = d, delta_elpd_se = se_d,
               max_pareto_k = max(loo_m$diagnostics$pareto_k, na.rm = TRUE))
  }))
  fwrite(tab, file.path(RESULTS_DIR, "table_loo_fixed_effect_selection.csv"))
  print(tab[, .(label, elpd = round(elpd, 1), elpd_se = round(elpd_se, 1),
                delta_elpd = round(delta_elpd, 1), delta_elpd_se = round(delta_elpd_se, 1))])
  cat(sprintf("Saved %s\n", file.path(RESULTS_DIR, "table_loo_fixed_effect_selection.csv")))
}
