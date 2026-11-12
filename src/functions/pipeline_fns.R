# =============================================================================
# src/functions/pipeline_fns.R
#
# ANALYSIS COMPONENT: PIPELINE CONTROL (settings overrides, cluster helpers, MCMC fit variants)
#
# Driver helpers shared by the step scripts (never sourced on the cluster):
#   - take_overrides()/setting(): let a sensitivity script or run_sections.R
#     run a step script with different settings without editing it.
#   - collect_bundle(): the wait / retry / collect pattern every cluster step uses.
#   - MCMC variants: the output folder and priors for each fit (main analysis,
#     RDT sensitivity analysis, sigma_d prior-SD analysis, mu_d prior-centre
#     analysis).
# =============================================================================

source("./src/functions/paths.R")

# Settings overrides ===========================================================
# A step script starts with OVR <- take_overrides(). If the caller defined
# STEP_OVERRIDES (a named list), the script uses those values and removes the
# list straight away, so a later interactive run of the same script falls back
# to the defaults written in the script.
take_overrides <- function() {
  if (!exists("STEP_OVERRIDES", envir = globalenv(), inherits = FALSE)) return(list())
  ovr <- get("STEP_OVERRIDES", envir = globalenv())
  rm("STEP_OVERRIDES", envir = globalenv())
  ovr
}
setting <- function(ovr, name, default) if (is.null(ovr[[name]])) default else ovr[[name]]

# Run a step script once per element of `overrides` (a list of named lists).
run_step <- function(script, overrides) {
  for (o in overrides) {
    cat(sprintf("\n######## %s | %s ########\n", script,
                paste(sprintf("%s = %s", names(o), vapply(o, function(x) paste(format(x), collapse = ","), "")), collapse = "; ")))
    assign("STEP_OVERRIDES", o, envir = globalenv())
    source(script, local = globalenv())
  }
  invisible(NULL)
}

# Cluster bundles ==============================================================
# Waits for a hipercow bundle, retries failed/cancelled tasks up to
# max_retries times (transient node failures are routine at this scale), and
# returns the successful results. allow_failures = FALSE stops instead of
# dropping tasks (used for MCMC chains: R-hat needs every chain).
collect_bundle <- function(bundle, max_retries = 3L, allow_failures = TRUE, label = "task") {
  hipercow::hipercow_bundle_wait(bundle, timeout = Inf)
  for (attempt in seq_len(max_retries)) {
    statuses <- hipercow::hipercow_bundle_status(bundle)
    n_failed <- sum(statuses != "success")
    if (n_failed == 0L) break
    cat(sprintf("Retry %d/%d: %d %s(s) not successful, retrying...\n", attempt, max_retries, n_failed, label))
    hipercow::hipercow_bundle_retry(bundle, if_status_in = c("failure", "cancelled"))
    hipercow::hipercow_bundle_wait(bundle, timeout = Inf)
  }
  statuses <- hipercow::hipercow_bundle_status(bundle)
  n_failed <- sum(statuses != "success")
  if (n_failed > 0L && !allow_failures) {
    stop(sprintf("%d/%d %s(s) still failed after %d retries -- fix and resubmit.", n_failed, length(statuses), label, max_retries))
  }
  if (n_failed > 0L) cat(sprintf("%d/%d %s(s) still failed after %d retries -- dropped.\n", n_failed, length(statuses), label, max_retries))
  # hipercow_bundle_result() returns the condition object in a failed task's slot.
  raw <- hipercow::hipercow_bundle_result(bundle)
  ok  <- !vapply(raw, inherits, logical(1), what = "condition")
  list(results = raw[ok], ok = ok, raw = raw)
}

# MCMC variants ================================================================
# Main analysis: RDT sensitivity 0.9, sigma_d ~ Normal(0.3, 0.01).
# RDT sensitivity analysis: 0.8 and 1.0 (sigma_d prior SD 0.01).
# sigma_d prior-SD analysis: 0.05 and 0.1 (RDT sensitivity 0.9 only).
SIGMA_D_PRIOR_SD_MAIN    <- 0.01
SIGMA_D_PRIOR_SD_OPTIONS <- c(0.05, 0.1)

is_main_sigma_sd <- function(sigma_sd) abs(sigma_sd - SIGMA_D_PRIOR_SD_MAIN) < 1e-9   # vectorised

# mu_d prior centre ------------------------------------------------------------
# Main fit: the 15-day median time to HRP2-negative (Dalrymple et al. 2018) is
# the median of the S-stage clearance alone, Exp(mean d), so d = 15 / log 2
# (21.6 days). In the malariasimulation fork the clearance clock starts on
# entry to S (create_hrp2_clearance_process() only clears HRP2-positive people
# in S), and a treated case reaches S after Tr ~ Exp(mean dt). If the 15 days
# instead run from treatment start, the time to HRP2-negative is
# Exp(dt) + Exp(d) (hypoexponential), and d is the value that puts that sum's
# median at 15 days (about 13.8 days, half-life about 9.5 days, for dt = 5).
# Prior-centre sensitivity analysis: mu_d centred on that d, with logit-normal
# priors, RDT sensitivity 0.8, 0.9, 1.0 (archived: archive/src/3c_calibration_mcmc/sensitivity/mu_d_prior_mcmc.R). Not part
# of mcmc_variants().
HRP2_TREATMENT_DT <- 5   # Tr -> S mean (days); $dt in all 249 cached grid base_pars (checked 2026-10-19)

# CDF of Exp(mean d) + Exp(mean dt) at t; Erlang(2) when d == dt.
hypoexp_cdf <- function(t, d, dt) {
  if (abs(d - dt) < 1e-6) return(1 - exp(-t / d) * (1 + t / d))
  1 - (d * exp(-t / d) - dt * exp(-t / dt)) / (d - dt)
}
# Mean S-stage clearance time d whose post-treatment median is median_days.
# No root (uniroot() error) when Exp(dt) alone has median >= median_days.
d_from_median_post_treatment <- function(median_days = 15, dt = HRP2_TREATMENT_DT) {
  uniroot(function(d) hypoexp_cdf(median_days, d, dt) - 0.5,
          lower = 1e-6, upper = median_days / log(2), tol = 1e-10)$root
}
stopifnot(abs(d_from_median_post_treatment(15, 5) - 13.8) < 0.1)

MU_D_PRIOR_MEAN_MAIN    <- log(15 / log(2))
MU_D_PRIOR_MEAN_OPTIONS <- c(post_treatment = log(d_from_median_post_treatment(15, HRP2_TREATMENT_DT)))

is_main_mu_d <- function(mu_d_mean) abs(mu_d_mean - MU_D_PRIOR_MEAN_MAIN) < 1e-9   # vectorised
# Folder/bundle label of a prior centre: its name in MU_D_PRIOR_MEAN_OPTIONS,
# else the value itself.
mu_d_label <- function(mu_d_mean) {
  hit <- which(abs(MU_D_PRIOR_MEAN_OPTIONS - mu_d_mean) < 1e-9)
  if (length(hit)) names(MU_D_PRIOR_MEAN_OPTIONS)[hit[1]] else sprintf("%.3f", mu_d_mean)
}

# Prior family (mu_d and site log EIR) -----------------------------------------
# The emulators only cover the grid: half-life 5-365 days (GRID_HALFLIFE_RANGE)
# and EIR GRID_EIR_RANGE, so both priors are cut to those ranges.
#   "normal":       mu_d ~ Normal(mu_d_mean, mu_d_sd) on log d, unbounded, with
#                   7% (main centre) or 19% (post-treatment centre) of its mass
#                   below the 5-day floor, where no site can follow.
#                   log EIR_i ~ Normal(eir_mu, eir_sd), hard-truncated to the
#                   grid (2.4% of its mass outside).
#   "logit_normal": each prior is Normal on logit((x - L) / (U - L)), [L, U] the
#                   grid range (identical to the emulators' log_d_bounds and
#                   log_eir_bounds), so nothing is truncated. Each is matched to
#                   its normal's median and one upper quantile:
#                     mu_d: the 97.5% point, so the upper tail is unchanged
#                       and the lower tail is squeezed above the 5-day floor.
#                       Half-life 15.0 (6.9-65.2) days at the main centre,
#                       9.5 (5.7-41.5) at the post-treatment centre.
#                     EIR: the 75% point (so the same interquartile range),
#                       centred on the grid's geometric midpoint; logit SD
#                       ~0.92. Keeps the normal's peak (matching the 97.5%
#                       point instead gave logit SD ~1.37, nearly flat across
#                       EIR 1-1,000, which weakens the EIR prior's role in
#                       separating decay from EIR at each site).
# PLANNED: "logit_normal" will become the main analysis. Until then its fits
# are sensitivity fits (archived: archive/src/3c_calibration_mcmc/sensitivity/mu_d_prior_mcmc.R). To switch, set
# PRIOR_FAMILY_MAIN <- "logit_normal" -- but first move the current main-folder
# outputs (DIR_MCMC, sensitivity_rdt_se_*/, sensitivity_sigma_d_sd_*/, all
# fitted with the normal priors) aside, because those folders will then mean
# logit-normal fits, and re-run them.
PRIOR_FAMILY_MAIN    <- "normal"
PRIOR_FAMILY_OPTIONS <- c("normal", "logit_normal")

is_main_prior_family <- function(prior_family) identical(prior_family, PRIOR_FAMILY_MAIN)

# Logit-scale mean and SD of the logit-normal on `bounds` matched to a
# Normal(mean, sd)'s median and its upper quantile q (q > 0.5).
logit_normal_pars <- function(mean, sd, bounds, q = 0.975) {
  p_med <- (mean - bounds[1]) / diff(bounds)
  p_hi  <- (mean + qnorm(q) * sd - bounds[1]) / diff(bounds)
  if (!(p_med > 0 && p_hi < 1)) stop(sprintf("logit_normal_pars(): median or %g quantile outside the bounds.", q))
  c(mean = qlogis(p_med), sd = (qlogis(p_hi) - qlogis(p_med)) / qnorm(q))
}

# Folder/bundle label of a non-main prior (mu_d centre and/or prior family);
# NULL for the main priors.
prior_variant_label <- function(mu_d_mean, prior_family) {
  parts <- c(if (!isTRUE(is_main_mu_d(mu_d_mean))) sprintf("mu_d_%s", mu_d_label(mu_d_mean)),
             if (!is_main_prior_family(prior_family)) sprintf("%s_priors", prior_family))
  if (length(parts)) paste(parts, collapse = "_") else NULL
}

# Folder of one MCMC fit (and its posterior predictive check). Main and RDT
# sensitivity fits use cluster_fns.R::rdt_sensitivity_dir(); the sigma_d
# prior-SD and mu_d/prior-family fits use their own folders, tagged with the
# RDT sensitivity.
mcmc_variant_dir <- function(rdt_sensitivity = RDT_SENSITIVITY_MAIN, sigma_d_prior_sd = SIGMA_D_PRIOR_SD_MAIN,
                             base_dir = DIR_MCMC, mu_d_mean = MU_D_PRIOR_MEAN_MAIN,
                             prior_family = PRIOR_FAMILY_MAIN) {
  prior_part <- prior_variant_label(mu_d_mean, prior_family)
  if (!is.null(prior_part)) {
    sd_part <- if (isTRUE(is_main_sigma_sd(sigma_d_prior_sd))) "" else sprintf("_sigma_d_sd_%s", format(sigma_d_prior_sd))
    return(file.path(base_dir, sprintf("sensitivity_%s%s_%s", prior_part, sd_part,
                                       rdt_sensitivity_tag(rdt_sensitivity))))
  }
  if (isTRUE(is_main_sigma_sd(sigma_d_prior_sd))) return(rdt_sensitivity_dir(base_dir, rdt_sensitivity))
  file.path(base_dir, sprintf("sensitivity_sigma_d_sd_%s_%s", format(sigma_d_prior_sd), rdt_sensitivity_tag(rdt_sensitivity)))
}
mcmc_variant_tag <- function(rdt_sensitivity = RDT_SENSITIVITY_MAIN, sigma_d_prior_sd = SIGMA_D_PRIOR_SD_MAIN,
                             mu_d_mean = MU_D_PRIOR_MEAN_MAIN, prior_family = PRIOR_FAMILY_MAIN) {
  tag <- sprintf("%s_sigma_sd_%s", rdt_sensitivity_tag(rdt_sensitivity), format(sigma_d_prior_sd))
  prior_part <- prior_variant_label(mu_d_mean, prior_family)
  if (!is.null(prior_part)) tag <- sprintf("%s_%s", tag, prior_part)
  tag
}

# Every fit, one row each. `analysis` lists the sensitivity analyses a fit
# belongs to (the main fit is the reference in both).
mcmc_variants <- function() {
  data.table::data.table(
    rdt_sensitivity  = c(RDT_SENSITIVITY_OPTIONS, rep(RDT_SENSITIVITY_MAIN, length(SIGMA_D_PRIOR_SD_OPTIONS))),
    sigma_d_prior_sd = c(rep(SIGMA_D_PRIOR_SD_MAIN, length(RDT_SENSITIVITY_OPTIONS)), SIGMA_D_PRIOR_SD_OPTIONS)
  )
}

# Priors for one fit. sigma_d_prior_sd, mu_d_mean and prior_family differ between variants.
#   mu_d: on log(d), median mu_d_mean, Normal(mu_d_mean, 0.75) or its bounded
#     logit-normal match (see "Prior family" above). The main mu_d_mean,
#     log(15 / log 2), centres it on a 15-day HRP2 half-life (Dalrymple et al.
#     2018); 95% prior range about 3.4-65 days under the normal. See "mu_d prior
#     centre" above for the alternative centre.
#   sigma_d: Normal(0.3, sigma_d_prior_sd), truncated at zero, on its natural
#     scale. 0.3 is a ~0.55x-1.8x multiplicative spread around the population
#     median. The main fit's SD = 0.01 is deliberately tight: with no country
#     layer, it stops outlying countries inflating the population spread.
#   EIR: on log EIR, the same for every site; Normal(eir_mu, eir_sd) centred on
#     the grid's geometric midpoint with the grid ends at about +-2.25 SD, or
#     its bounded logit-normal match.
# Grid bounds are read at call time (GRID_HALFLIFE_RANGE, GRID_EIR_RANGE are in cluster_fns.R).
build_mcmc_priors <- function(sigma_d_prior_sd = SIGMA_D_PRIOR_SD_MAIN, mu_d_mean = MU_D_PRIOR_MEAN_MAIN,
                              prior_family = PRIOR_FAMILY_MAIN) {
  prior_family <- match.arg(prior_family, PRIOR_FAMILY_OPTIONS)
  priors <- list(
    mu_d_mean          = mu_d_mean,
    mu_d_sd            = 0.75,
    sigma_d_prior_mean = 0.3,
    sigma_d_prior_sd   = sigma_d_prior_sd,
    eir_mu             = mean(log(GRID_EIR_RANGE)),
    eir_sd             = diff(log(GRID_EIR_RANGE)) / 4.5,
    decay_family       = "normal",
    decay_t_df         = 4,
    mu_d_family        = prior_family,
    eir_family         = prior_family
  )
  if (prior_family == "logit_normal") {
    priors$mu_d_bounds <- log(GRID_HALFLIFE_RANGE / log(2))
    lp <- logit_normal_pars(mu_d_mean, priors$mu_d_sd, priors$mu_d_bounds)
    priors$mu_d_logit_mean <- lp[["mean"]]
    priors$mu_d_logit_sd   <- lp[["sd"]]
    priors$eir_bounds <- log(GRID_EIR_RANGE)
    lp <- logit_normal_pars(priors$eir_mu, priors$eir_sd, priors$eir_bounds, q = 0.75)   # same IQR: keeps the peak
    priors$eir_logit_mean <- lp[["mean"]]
    priors$eir_logit_sd   <- lp[["sd"]]
  }
  priors
}


# Misc =========================================================================
sanitise_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
# hipercow bundle names allow letters, numbers, "-" and "_" only (no "." as in "sigma_sd_0.01").
safe_bundle_name <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

# Fixed per-country deletion prevalence (data/hrp2_deletion.csv central value).
deletion_by_iso <- function(isos) {
  deletion_dt <- data.table::fread(HRP2_DELETION_CSV)
  data.table::rbindlist(lapply(isos, function(i) {
    m <- get_hrp2_deletion_mean(i, deletion_dt)
    data.table::data.table(iso = i, deletion_pct = m$mean_pct, source = m$source)
  }))
}
