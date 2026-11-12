# =============================================================================
# src/2_regression/fit_fdr_fpr_regression.R
#   (was src/03b_brms_int_local_fit.R + Section D of src/03b_brms_int_cluster_fit.R)
#
# ANALYSIS COMPONENT: REGRESSION ANALYSIS (main model)
#
# Step 2. Individual-level categorical regression of the joint microscopy x
# RDT result (p_p / n_p / p_n / n_n) for children tested by both, and the
# false-positive rate (FPR, P(RDT+ | LM-)) and false-discovery rate
# (FDR, P(LM- | RDT+)) by covariate, from the posterior draws.
#
#   outcome ~ sr_rdt_cat * hrp2_risk + has_smc + during_or_after_peak + age_c
#             + care_or_act + fever + (1 | survid),   family = categorical(refcat = "p_p")
#
# Regularising priors (normal(0, 2.5) on fixed effects, normal(0, 5) on
# intercepts; brms_int_cluster_fns.R) because the sparse sr_rdt_cat x
# hrp2_risk cells separate under flat priors (Rhat up to 1.64 on a 3,000-row
# subsample; all Rhat <= 1.013 with these priors).
#
# Section A fits the model locally with rstan (4 chains x 500 iterations, 200
# warmup, max_treedepth 12). 300 iterations on the full 50-survey data left the
# munn intercept and survey random effects at Rhat 1.05 and bulk ESS ~130,
# with 33 treedepth hits. The fit is run locally because on the DIDE cluster the shared
# CmdStan install is read-only to the worker account and rstan fails with
# "invalid connection" under hipercow's redirected stdout. 2,000 iterations
# would take 20-30 h per chain on the full ~98,700 rows.
#
# Section B computes, from the posterior draws, the FPR and FDR at each level
# of each covariate (all other covariates at the reference cell: SR prevalence
# 10-20%, medium deletion risk, binary covariates off, mean age) and across
# the full factorial covariate grid. The ratio is taken within each draw, so
# the credible intervals carry the posterior correlation between terms.
#
# Outputs (DIR_REGRESSION):
#   fit_brms_int_full_local.RDS     fitted brmsfit (Section A)
#   fpfn_by_covariate.RDS           paired on-vs-off differences per covariate (Section B)
#   fpfn_ame.RDS                    average marginal effect of each binary covariate / age on FPR and FDR (pp)
#   fpfn_full_grid.RDS              FPR/FDR for all 144 covariate combinations
# fit_brms_int_full_local.RDS embeds the individual-level DHS/MIS data, so it
# stays in OUTPUT_ROOT (git-ignored). Shareable posterior draws and summaries
# are written by src/export_posteriors.R to results/posteriors/regression/.
# The AME forest plot (fpfn_ame_forest.png) is drawn in src/plots.R.
#
# Fixed-effect selection (LOO, appendix Table tab:elpd) is in
# fixed_effect_selection.R. The no-deletion-risk sensitivity fit (not in the
# manuscript), the fever-only model, the cluster (hipercow) fitting path and
# the country check are in archive/src/.
# =============================================================================
library(data.table)
library(brms)
source("./src/functions/paths.R")
source("./src/functions/pipeline_fns.R")
source("./src/functions/brms_int_cluster_fns.R")
source("./src/functions/check_fns.R")

OUT_DIR <- DIR_REGRESSION
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OVR <- take_overrides()
FIT_MODEL      <- setting(OVR, "FIT_MODEL", FALSE)       # Section A: refit the model (several hours)
COVARIATE_FPFN <- setting(OVR, "COVARIATE_FPFN", TRUE)   # Section B: FPR/FDR summaries from the fit

# Fit ==========================================================================
# SECTION A -- LOCAL FIT
N_CHAINS       <- 4L
CHAIN_N_ITER   <- setting(OVR, "CHAIN_N_ITER", 500L)
CHAIN_N_WARMUP <- setting(OVR, "CHAIN_N_WARMUP", 200L)
MAX_TREEDEPTH  <- setting(OVR, "MAX_TREEDEPTH", 12L)
SEED           <- 1L

if (FIT_MODEL) {
  data_full  <- prepare_brms_full_data(DIR_STRUCTURE)
  cat(sprintf("Full data: %d rows, %d survid groups\n", nrow(data_full), uniqueN(data_full$survid)))
  # Surveys with paired RDT/LM results in the included countries, against the
  # complete cases brm() will actually use (it drops NA covariates silently)
  ind_check <- fread(file.path(DIR_STRUCTURE, "dhs_mis_individual_hrp2.csv"))
  ind_check <- ind_check[country %in% readRDS(file.path(DIR_STRUCTURE, "countries_included.RDS"))$country &
                         !is.na(rdt_result) & !is.na(micro_result)]
  check_model_surveys(data_full, BRMS_INT_FORMULA,
                      expected  = unique(paste0(ind_check$country, "_", ind_check$survey_year)),
                      min_share = setting(OVR, "MIN_SURVEY_SHARE", 0.9))
  rm(ind_check)
  prior_full <- build_brms_regularising_prior(BRMS_INT_FORMULA, data_full)

  t0 <- Sys.time()
  fit_brms_int_full <- brm(
    BRMS_INT_FORMULA,
    data    = data_full,
    family  = categorical(refcat = "p_p"),
    prior   = prior_full,
    chains  = N_CHAINS,
    iter    = CHAIN_N_ITER,
    warmup  = CHAIN_N_WARMUP,
    cores   = N_CHAINS,
    seed    = SEED,
    backend = "rstan",
    control = list(max_treedepth = MAX_TREEDEPTH),
    refresh = 100
  )
  cat(sprintf("Total fit time: %.2f hours\n", as.numeric(Sys.time() - t0, units = "hours")))
  saveRDS(fit_brms_int_full, file.path(OUT_DIR, "fit_brms_int_full_local.RDS"))

  s <- summary(fit_brms_int_full)
  print(s$fixed[, c("Estimate", "Est.Error", "Rhat", "Bulk_ESS", "Tail_ESS")])
  cat(sprintf("%d/%d fixed-effect terms have Rhat > 1.05\n", sum(s$fixed$Rhat > 1.05, na.rm = TRUE), nrow(s$fixed)))
}

# FPR / FDR by covariate =======================================================
# SECTION B -- TWO CONDITIONAL VIEWS OF RDT/MICROSCOPY DISCORDANCE
#   FPR = P(RDT+ | LM-) = P(n_p) / [P(n_p) + P(n_n)]  conditions on the true
#         (microscopy) state: what drives discordance (1 - specificity).
#   FDR = P(LM- | RDT+) = P(n_p) / [P(n_p) + P(p_p)]  conditions on the test
#         result: how to read routine RDT-positive case data (1 - PPV).
# Population-average predictions (re_formula = NA). sr_rdt_cat and hrp2_risk
# interact, so each prevalence comparison is repeated at every deletion-risk
# level and each deletion-risk comparison at every prevalence level.
if (COVARIATE_FPFN) {
  fit_path <- file.path(OUT_DIR, "fit_brms_int_full_local.RDS")
  if (!file.exists(fit_path)) stop("No fit_brms_int_full_local.RDS -- run Section A (FIT_MODEL <- TRUE) first.")
  cat(sprintf("Loading fitted model from %s\n", fit_path))
  fit_int_full <- readRDS(fit_path)

  REF_SR   <- "10-20%"
  REF_HRP2 <- "Marginal/Moderate"

  # Each comparison holds every OTHER covariate/stratifier at the reference
  # cell (REF_SR/REF_HRP2/other 0-1 covariates off/age_c = 0) and varies ONE
  # thing from an "off" to an "on" level -- the first 4 are 0/1 covariate
  # toggles, held at REF_SR/REF_HRP2 as before (the model has NO interaction
  # involving them, so their effect doesn't depend on which SR/HRP2 level is
  # held fixed).
  #
  # sr_rdt_cat/hrp2_risk are different: the fitted formula is
  # sr_rdt_cat * hrp2_risk (confirmed directly from the model -- 4 real
  # interaction coefficients per outcome equation), so the SR effect
  # genuinely differs depending on which HRP2 level it's measured at, and
  # vice versa. Rather than collapsing to one reference cell (as the
  # earlier version of this section did), each SR comparison is repeated at
  # EVERY hrp2_risk level (via `hold_var`/`hold_value` overriding the
  # default reference for that one comparison), and each HRP2 comparison at
  # every sr_rdt_cat level -- this is what actually exposes the interaction
  # rather than assuming it away.
  sr_levels   <- levels(fit_int_full$data$sr_rdt_cat)     # "<10%", "10-20%", ">=20%"
  hrp2_levels <- levels(fit_int_full$data$hrp2_risk)      # "Slight", "Marginal/Moderate", "High"

  COMPARISONS <- c(
    list(
      list(var = "has_smc",              off = 0, on = 1, label = "SMC",                      hold_var = NULL, hold_value = NULL),
      list(var = "during_or_after_peak", off = 0, on = 1, label = "During/after peak season", hold_var = NULL, hold_value = NULL),
      list(var = "care_or_act",          off = 0, on = 1, label = "Formal care/ACT",           hold_var = NULL, hold_value = NULL),
      list(var = "fever",                off = 0, on = 1, label = "Fever",                     hold_var = NULL, hold_value = NULL),
      # age_c is centred MONTHS (age_months - mean(age_months)), so one
      # additional YEAR of age is off = 0 (mean age) -> on = 12.
      list(var = "age_c",                off = 0, on = 12, label = "Age (+1 year)",            hold_var = NULL, hold_value = NULL)
    ),
    # SR effect, at each HRP2 level in turn
    unlist(lapply(c("10-20%", ">=20%"), function(sr_on) {
      lapply(hrp2_levels, function(h) list(
        var = "sr_rdt_cat", off = "<10%", on = sr_on,
        label = sprintf("Prevalence: %s vs <10%%\n(Reference deletion risk: %s)", sr_on, h),
        hold_var = "hrp2_risk", hold_value = h
      ))
    }), recursive = FALSE),
    # HRP2 effect, at each SR level in turn
    unlist(lapply(c("Marginal/Moderate", "High"), function(h_on) {
      lapply(sr_levels, function(s) list(
        var = "hrp2_risk", off = "Slight", on = h_on,
        label = sprintf("PfHRP2: %s vs Slight\n(Reference prevalence: %s)", h_on, s),
        hold_var = "sr_rdt_cat", hold_value = s
      ))
    }), recursive = FALSE)
  )

  build_ref_row <- function() {
    data.table::data.table(
      sr_rdt_cat           = factor(REF_SR, levels = levels(fit_int_full$data$sr_rdt_cat)),
      hrp2_risk            = factor(REF_HRP2, levels = levels(fit_int_full$data$hrp2_risk)),
      has_smc              = 0,
      during_or_after_peak = 0,
      care_or_act          = 0,
      fever                = 0,
      age_c                = 0,
      survid               = fit_int_full$data$survid[1]   # ignored: re_formula = NA below
    )
  }

  # sr_rdt_cat/hrp2_risk need a properly-levelled factor value (not a bare
  # string) to match the model's own training data; the four 0/1 covariates
  # are set directly.
  set_comparison_value <- function(row, var, value) {
    if (var %in% c("sr_rdt_cat", "hrp2_risk")) {
      data.table::set(row, j = var, value = factor(value, levels = levels(fit_int_full$data[[var]])))
    } else {
      data.table::set(row, j = var, value = value)
    }
    row
  }

  OUTCOME_SPECS <- list(
    list(level = "n_p", complement = "n_n", label = "P(RDT+ | microscopy-)"),
    list(level = "n_p", complement = "p_p", label = "P(microscopy- | RDT+)")
  )

  fpfn_raw <- lapply(COMPARISONS, function(cmp) {
    row_off <- build_ref_row()
    row_on  <- build_ref_row()
    if (!is.null(cmp$hold_var)) {
      row_off <- set_comparison_value(row_off, cmp$hold_var, cmp$hold_value)
      row_on  <- set_comparison_value(row_on,  cmp$hold_var, cmp$hold_value)
    }
    row_off <- set_comparison_value(row_off, cmp$var, cmp$off)
    row_on  <- set_comparison_value(row_on,  cmp$var, cmp$on)
    newdata <- data.table::rbindlist(list(row_off, row_on))
    epred <- brms::posterior_epred(fit_int_full, newdata = newdata, re_formula = NA)
    list(covariate = cmp$label, epred = epred)
  })

  # P(level | level-or-complement), within each posterior draw -- see
  # Section D header for why this (not the raw marginal epred[,,level])
  # is the right quantity.
  conditional_prob <- function(epred, row_idx, spec) {
    num <- epred[, row_idx, spec$level]
    num / (num + epred[, row_idx, spec$complement])
  }

  fpfn_by_covariate <- data.table::rbindlist(lapply(fpfn_raw, function(x) {
    data.table::rbindlist(lapply(OUTCOME_SPECS, function(spec) {
      p0   <- conditional_prob(x$epred, 1, spec)
      p1   <- conditional_prob(x$epred, 2, spec)
      diff <- p1 - p0
      data.table::data.table(
        covariate = x$covariate,
        outcome   = spec$label,
        off_mean  = mean(p0) * 100, off_lo = stats::quantile(p0, 0.025) * 100, off_hi = stats::quantile(p0, 0.975) * 100,
        on_mean   = mean(p1) * 100, on_lo  = stats::quantile(p1, 0.025) * 100, on_hi  = stats::quantile(p1, 0.975) * 100,
        diff_mean = mean(diff) * 100, diff_lo = stats::quantile(diff, 0.025) * 100, diff_hi = stats::quantile(diff, 0.975) * 100,
        p_on_greater = mean(p1 > p0)   # paired (same-draw) posterior probability that On > Off
      )
    }))
  }))

  print(fpfn_by_covariate)
  
  ##Manual fixes:
  fpfn_by_covariate[,covariate := gsub(pattern = '>=20%', replacement = '>20%', covariate)]
  fpfn_by_covariate[,covariate := gsub(pattern = 'Marginal/Moderate', replacement = 'Medium', covariate)]
  
  
  saveRDS(fpfn_by_covariate, file.path(OUT_DIR, "fpfn_by_covariate.RDS"))

  summ <- function(x) list(mean = mean(x) * 100, lo = stats::quantile(x, 0.025) * 100,
                           hi = stats::quantile(x, 0.975) * 100)

  # --- Average marginal effects (AME): each covariate switched off and on for
  # every child in the model data (everything else kept as observed), and the
  # FPR/FDR difference averaged over children within each draw. One number per
  # covariate on the map's percentage scale, without choosing a reference cell.
  # Children are collapsed to unique covariate rows with counts (~5,000 rows
  # instead of ~200,000). Formal care/ACT is averaged over febrile children
  # only, the group asked about care seeking. Age: each child one year older.
  ame_dt <- data.table::as.data.table(fit_int_full$data)[
    , .(n = .N), by = .(sr_rdt_cat, hrp2_risk, has_smc, during_or_after_peak, care_or_act, fever, age_c)]
  ame_dt[, survid := fit_int_full$data$survid[1]]   # ignored: re_formula = NA
  # Prevalence and deletion risk: every child moved from the reference level
  # to each other level, keeping their own level of the other factor, so the
  # interaction is averaged over the observed mix of that factor.
  ame_cat <- c(
    lapply(sr_levels[-1], function(s) list(var = "sr_rdt_cat", off = sr_levels[1], on = s,
                                           label = sprintf("Prevalence: %s vs %s", s, sr_levels[1]))),
    lapply(hrp2_levels[-1], function(h) list(var = "hrp2_risk", off = hrp2_levels[1], on = h,
                                             label = sprintf("PfHRP2: %s vs %s", h, hrp2_levels[1])))
  )
  fpfn_ame <- data.table::rbindlist(lapply(c(COMPARISONS[1:5], ame_cat), function(cmp) {
    rows <- if (cmp$var == "care_or_act") ame_dt[fever == 1] else ame_dt
    off <- data.table::copy(rows); on <- data.table::copy(rows)
    if (cmp$var == "age_c") {
      data.table::set(on, j = "age_c", value = on$age_c + cmp$on)
    } else if (cmp$var %in% c("sr_rdt_cat", "hrp2_risk")) {
      lv <- levels(fit_int_full$data[[cmp$var]])
      data.table::set(off, j = cmp$var, value = factor(rep(cmp$off, nrow(off)), levels = lv))
      data.table::set(on,  j = cmp$var, value = factor(rep(cmp$on,  nrow(on)),  levels = lv))
    } else {
      data.table::set(off, j = cmp$var, value = cmp$off)
      data.table::set(on,  j = cmp$var, value = cmp$on)
    }
    ep <- brms::posterior_epred(fit_int_full, newdata = rbind(off, on), re_formula = NA)
    n_rows <- nrow(rows); w <- rows$n / sum(rows$n)
    data.table::rbindlist(lapply(OUTCOME_SPECS, function(spec) {
      diff <- (conditional_prob(ep, n_rows + seq_len(n_rows), spec) - conditional_prob(ep, seq_len(n_rows), spec)) %*% w
      c(list(covariate = cmp$label, outcome = spec$label, n_children = sum(rows$n)), summ(as.vector(diff)))
    }))
  }))
  print(fpfn_ame)
  fpfn_ame[, covariate := sub(">=20%", ">20%", covariate, fixed = TRUE)]
  fpfn_ame[, covariate := sub("Marginal/Moderate", "Medium", covariate, fixed = TRUE)]
  saveRDS(fpfn_ame, file.path(OUT_DIR, "fpfn_ame.RDS"))

  # --- Full factorial grid: EVERY combination of covariate values (3 SR
  # levels x 3 HRP2 levels x 2^4 binary covariates = 144 combinations,
  # age_c held at 0/mean) -- this is the
  # full joint response surface. One posterior_epred() call over the whole
  # grid (not 144 separate calls).
  grid_dt <- data.table::CJ(
    sr_rdt_cat           = levels(fit_int_full$data$sr_rdt_cat),
    hrp2_risk            = levels(fit_int_full$data$hrp2_risk),
    has_smc              = c(0L, 1L),
    during_or_after_peak = c(0L, 1L),
    care_or_act          = c(0L, 1L),
    fever                = c(0L, 1L),
    sorted               = FALSE
  )
  grid_dt[, sr_rdt_cat := factor(sr_rdt_cat, levels = levels(fit_int_full$data$sr_rdt_cat))]
  grid_dt[, hrp2_risk  := factor(hrp2_risk,  levels = levels(fit_int_full$data$hrp2_risk))]
  grid_dt[, age_c      := 0]
  grid_dt[, survid     := fit_int_full$data$survid[1]]   # ignored: re_formula = NA below

  grid_epred <- brms::posterior_epred(fit_int_full, newdata = grid_dt, re_formula = NA)

  grid_results <- data.table::rbindlist(lapply(seq_len(nrow(grid_dt)), function(i) {
    data.table::rbindlist(lapply(OUTCOME_SPECS, function(spec) {
      p <- conditional_prob(grid_epred, i, spec)
      data.table::data.table(
        sr_rdt_cat            = as.character(grid_dt$sr_rdt_cat[i]),
        hrp2_risk             = as.character(grid_dt$hrp2_risk[i]),
        has_smc               = grid_dt$has_smc[i],
        during_or_after_peak  = grid_dt$during_or_after_peak[i],
        care_or_act           = grid_dt$care_or_act[i],
        fever                 = grid_dt$fever[i],
        outcome               = spec$label,
        p_mean = mean(p) * 100, p_lo = stats::quantile(p, 0.025) * 100, p_hi = stats::quantile(p, 0.975) * 100
      )
    }))
  }))

  ## Manual fixes (same convention as fpfn_by_covariate above):
  grid_results[, sr_rdt_cat := gsub(">=20%", ">20%", sr_rdt_cat, fixed = TRUE)]
  grid_results[, hrp2_risk  := gsub("Marginal/Moderate", "Medium", hrp2_risk, fixed = TRUE)]

  saveRDS(grid_results, file.path(OUT_DIR, "fpfn_full_grid.RDS"))

  cat(sprintf("Full factorial grid: %d covariate combinations.\n", nrow(grid_dt)))
}
