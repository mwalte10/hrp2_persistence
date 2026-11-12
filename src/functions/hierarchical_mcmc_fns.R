# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: hierarchical MCMC)
#
# src/functions/hierarchical_mcmc_fns.R
# Pure function library for step 3c's (3c_calibration_mcmc/hierarchical_mcmc.R) blocked
# Metropolis-Hastings sampler for Bayesian inference of HRP2 decay and
# site-specific EIR using GP likelihood emulators.
#
# MODEL: every site's decay prior shares ONE population-level sigma_d --
#
#   Population-level (hyper)parameters:
#     mu_d    ~ Normal(priors$mu_d_mean, priors$mu_d_sd), or its logit-normal
#               match on log_d_bounds (priors$mu_d_family; see log_prior_mu_d())
#     sigma_d ~ Normal(priors$sigma_d_prior_mean, priors$sigma_d_prior_sd),
#               truncated at 0
#   Site-specific, partial pooling around the SINGLE population mean AND
#   SINGLE population spread (no country-specific layer at all):
#     log(d_i) ~ Normal/StudentT(mu_d, sigma_d), truncated to log_d_bounds
#   Site-specific, no pooling:
#     log(eir_i) ~ Normal(priors$eir_mu, priors$eir_sd), truncated to
#       [log(EIR_GRID[1]), log(EIR_GRID[end])], or its logit-normal match on
#       that range (priors$eir_family; see log_prior_site_eir())
#
# This replaced an earlier 4-level hierarchy with a country-specific decay
# spread (sigma_d_c[c]); that let Cote d'Ivoire, Ghana and Nigeria's sites
# float almost unconstrained from the population mean. The retired sampler is
# in archive/src/functions/hierarchical_mcmc_fns_country_spread.R.
#
# Sampling blocks per sweep (this primary model):
#   Block 1: mu_d      population grand mean decay -- no emulator calls
#   Block 2: sigma_d    population decay spread -- no emulator calls
#   Block 3: log_d_i    per-site decay, individually (partial pooling)
#   Block 4: log_eir_i  per-site EIR, individually (no pooling)
#
# Usage: see run_mcmc() and initialise_theta() below.
#
# IMPORTANT: this file must contain ONLY function definitions (plus the
# library()/source() calls they need) -- no top-level submit/wait/collect
# side effects. It's sourced by hipercow_environment_create() (the
# "hierarchical_mcmc" environment), which means every worker task loading
# that environment re-executes this file's top level. The driver script
# (06_hierarchical_mcmc.R) sources this file and keeps all the
# if(TRUE)/if(FALSE)-gated submit/wait/collect sections for interactive use.

library(data.table)
library(hipercow)

# predict_emulator_ll() (used by compute_log_posterior()/run_mcmc() below) is
# defined in emulator_fns.R, not redefined here.
source("./src/functions/emulator_fns.R")


# --- Log-priors ---

# Log-density of a logit-normal on [bounds[1], bounds[2]]: Normal(m, s) on
# logit((x - L) / (U - L)), with the log-Jacobian of that transform; -Inf
# outside (L, U). Used by the "logit_normal" prior family (pipeline_fns.R
# "Prior family"), which needs no truncation.
log_dlogitnorm_bounded <- function(x, bounds, m, s) {
  width <- diff(bounds)
  p <- (x - bounds[1]) / width
  if (p <= 0 || p >= 1) return(-Inf)
  dnorm(qlogis(p), mean = m, sd = s, log = TRUE) - log(p) - log1p(-p) - log(width)
}

log_prior_mu_d <- function(mu_d, priors) {
  # Population-level GRAND mean decay (across all countries). Centred
  # (06's PRIORS$mu_d_mean, normally log(15 / log(2))) on a reported HRP2
  # half-life of 15 days: "half of RDTs that detect HRP2 are still positive
  # 15 (5-32) days post-treatment" -- median 15d, loosely informed by the
  # reported 5-32d range. log(15 / log(2)) converts that half-life into
  # log(d), d being the mean-duration "scale" parameter malariasimulation
  # actually uses for its exponential clearance process
  # (half_life = d * log(2); see cluster_fns.R::halflife_to_scale()).
  # This reflects findings of the Dalrymple paper
  # https://link.springer.com/article/10.1186/s12936-018-2371-9
  # See 06_hierarchical_mcmc.R's own PRIORS$mu_d_mean/mu_d_sd for the actual
  # numbers and their calibration history.
  #
  # priors$mu_d_family (pipeline_fns.R "Prior family"; absent = "normal"):
  #   "normal":       Normal(mu_d_mean, mu_d_sd).
  #   "logit_normal": logit-normal on priors$mu_d_bounds (log_dlogitnorm_bounded()).
  if (!identical(priors$mu_d_family, "logit_normal")) {
    return(dnorm(mu_d, mean = priors$mu_d_mean, sd = priors$mu_d_sd, log = TRUE))
  }
  log_dlogitnorm_bounded(mu_d, priors$mu_d_bounds, priors$mu_d_logit_mean, priors$mu_d_logit_sd)
}


log_prior_sigma_d <- function(sigma_d, mean, sd) {
  # Population-level WITHIN-country decay spread -- the TRUE population
  # distribution's own SD (see this file's own REDESIGN 2 note). ESTIMATED
  # (2026-09-17; previously a FIXED constant, priors$sigma_d_mean), sampled
  # DIRECTLY on its natural (non-negative) scale via a plain Normal,
  # truncated at zero -- deliberately NOT log-transformed and NOT
  # Half-Normal. This is a real, substantive choice, not an oversight:
  # sigma_d's prior MODE is at `mean` (06's PRIORS$sigma_d_prior_mean,
  # currently 0.3) -- a genuine a priori expectation about typical
  # within-country spread -- unlike tau_sigma_d (log_prior_tau_sigma_d()
  # below), whose mode is deliberately at zero ("no additional spread" is
  # the default there). A Half-Normal centred at zero could not express "0.3
  # is the expected value," so a plain Normal is used instead; the
  # truncation at zero is a technical necessity (sigma_d is a SD, must stay
  # positive) rather than a modelling statement, similar in spirit to how
  # log_prior_site_decay()'s own truncation doesn't rigorously renormalize
  # for the same reason (see that function's own note) -- since `mean` sits
  # several `sd`s away from zero (e.g. 0.3 vs sd = 0.1, three SDs), the
  # truncated probability mass is negligible in practice and this
  # approximation (plain dnorm(), no renormalizing constant) is a reasonable
  # trade-off rather than a rigorous truncated-Normal treatment.
  #
  # mean/sd: 06's PRIORS$sigma_d_prior_mean/sigma_d_prior_sd -- see those
  # constants' own comments in 06_hierarchical_mcmc.R for the chosen values.
  if (sigma_d < 0) return(-Inf)
  dnorm(sigma_d, mean = mean, sd = sd, log = TRUE)
}


log_prior_site_decay <- function(log_d, mu_d, sigma_d_c, bounds, family, df) {
  # Site-level decay prior, log(d_i) | mu_d, sigma_d_c -- CONDITIONAL on the
  # SINGLE population mean (mu_d, shared by every site in every country) and
  # that site's OWN country's spread (sigma_d_c, from Block 3 below) -- see
  # this file's own header for the full hierarchy: mu_d -> tau_sigma_d ->
  # sigma_d_c[c] -> log_d_i. TWO interchangeable families, toggled by
  # `family` (06_hierarchical_mcmc.R's own PRIORS$decay_family, the value
  # every caller must now pass explicitly -- no module-level fallback
  # default any more, see this file's header for why):
  #   "normal": Normal(mu_d, sigma_d_c) -- the original, standard choice.
  #   "t":      Student-t(df, location = mu_d, scale = sigma_d_c) --
  #             heavier tails, tried 2026-09-14 (under the earlier
  #             mean-varying design) after a real fit showed the Normal
  #             version let a handful of sites with poor grid/emulator
  #             coverage on the decay axis (13/14 Côte d'Ivoire regions,
  #             most visibly CIV_Abidjan at a ~1123-day posterior mean
  #             half-life) inflate the (then-single, global) sigma_d to
  #             ~4.3 -- about 14x the prior's own expectation -- weakening
  #             the shrinkage every OTHER, well-behaved site received from
  #             partial pooling.
  # Kept at "normal" as the default through the 2026-09-16 country-spread
  # redesign -- the heavy-tailed option is kept available (not deleted) as a
  # general safety margin against sites that turn out to be inflated but
  # aren't (yet) well-explained by the country-specific spread layer alone
  # -- flip PRIORS$decay_family back to "t" in 06 if that turns out to still
  # matter empirically. Either way, the population-level hyperpriors
  # (log_prior_mu_d()/log_prior_tau_sigma_d()/log_prior_log_sigma_d_c()
  # above) are Normal regardless -- only this conditional, site-level
  # layer's family is toggled.
  #
  # log-density: dnorm() directly for "normal"; for "t", dt() gives the
  # standard (location 0, scale 1) t-density, so the location-scale version
  # needs the usual -log(scale) Jacobian term.
  #
  # bounds hard-truncates outside the emulator's own trained log_d range
  # (further capped by 06's own MAX_DECAY_HALFLIFE_DAYS, see
  # build_and_tune_mcmc_inputs()), for the identical reason
  # log_prior_site_eir() truncates eir -- the GP was never trained out
  # there. This applies under EITHER family.
  #
  # NOTE on the truncation's own normalizing constant: unlike log_eir
  # (whose Normal(priors$eir_mu, priors$eir_sd) has FIXED parameters, so
  # truncating it just drops an unimportant additive constant from the
  # log-posterior), this prior's mu_d/sigma_d_c are themselves being
  # SAMPLED -- the truncated distribution's true normalizing constant
  # technically depends on their current values, which this function does
  # NOT account for (same simple hard -Inf-outside-bounds treatment as EIR
  # gets, under either family). That's an approximation, not exact: several
  # blocks' acceptance ratios in run_mcmc() are very slightly biased by
  # ignoring how the truncated normalizing constant shifts with each
  # proposal. Accepted as a reasonable trade-off rather than building a
  # fully rigorous truncated-hierarchical MH step -- revisit if mu_d/
  # sigma_d_c's own posteriors ever look suspiciously pinned near the
  # truncation bounds themselves.
  if (log_d < bounds[1] || log_d > bounds[2]) return(-Inf)
  if (family == "normal") {
    dnorm(log_d, mean = mu_d, sd = sigma_d_c, log = TRUE)
  } else if (family == "t") {
    stats::dt((log_d - mu_d) / sigma_d_c, df = df, log = TRUE) - log(sigma_d_c)
  } else {
    stop(sprintf("log_prior_site_decay(): unknown family '%s' -- must be 'normal' or 't'.", family))
  }
}

# Vectorised-over-log_d twin of log_prior_site_decay() above -- SAME density
# formula and SAME hard -Inf-outside-bounds truncation, just usable on a
# whole grid at once (log_prior_site_decay() itself uses a scalar `if`, so
# can't take a vector). Kept as a separate function rather than rewriting
# log_prior_site_decay() itself, since that one is the real sampler's own
# per-iteration MH step and changing its scalar `if` to a vectorised
# ifelse() would very slightly slow down the hot path for no benefit there.
log_prior_site_decay_vec <- function(log_d_grid, mu_d, sigma_d_c, bounds, family, df) {
  ll <- if (family == "normal") {
    dnorm(log_d_grid, mean = mu_d, sd = sigma_d_c, log = TRUE)
  } else if (family == "t") {
    stats::dt((log_d_grid - mu_d) / sigma_d_c, df = df, log = TRUE) - log(sigma_d_c)
  } else {
    stop(sprintf("log_prior_site_decay_vec(): unknown family '%s' -- must be 'normal' or 't'.", family))
  }
  ifelse(log_d_grid < bounds[1] | log_d_grid > bounds[2], -Inf, ll)
}

# Recompute each site's site-level decay posterior under a COUNTERFACTUAL
# where every site is pooled with the POPULATION-level sigma_d instead of
# its own country's sigma_d_c[c] -- i.e. as if the country-spread layer
# (sigma_d_c/tau_sigma_d, see this file's own REDESIGN/REDESIGN 2 notes)
# didn't exist and every site shared ONE population-wide within-country
# spread. Used by 06_hierarchical_mcmc.R's own forest-plot comparison
# ("what would each country's half-life look like if only the
# population-level spread had been used").
#
# This does NOT refit mu_d/sigma_d/eir_i themselves -- their marginal
# posteriors are reused exactly AS ALREADY FITTED (a reasonable
# approximation: removing one noisy country's outsized sigma_d_c mainly
# changes how much THAT country's own sites shrink toward mu_d, only a
# second-order effect on mu_d/sigma_d's own population-level estimates).
# Only log_d_i's own CONDITIONAL posterior is genuinely re-derived -- EXACTLY
# (not approximated), via grid quadrature over log_d_bounds using the SAME
# emulator each site's actual posterior was fit against, for a THINNED
# subset of the real joint-posterior draws (mu_d/sigma_d/eir_i held fixed at
# that draw's own actual values; consecutive retained MCMC draws are highly
# autocorrelated anyway, so thinning loses little for a summary quantity
# like this, at a large compute saving -- n_sites x thin_draws emulator
# calls instead of n_sites x nrow(pooled_samples)).
#
# emulators: named list of GP emulators (one per site_name) -- the SAME
#   object run_mcmc() itself was fit against.
# pooled_samples: mcmc_combined.RDS's own pooled_samples data.table -- must
#   have `mu_d`, `sigma_d` (both log/natural as appropriate -- mu_d is
#   log-scale, sigma_d is already natural-scale, see this file's header) and
#   one `eir_<i>` column (natural/absolute) per site.
# site_names: mcmc_inputs.RDS's own site_names, in the SAME order pooled_
#   samples' own eir_<i>/d_<i> column suffixes use.
# log_d_bounds: mcmc_inputs.RDS's own log_d_bounds (single shared numeric(2)
#   -- see 06's own absolute-grid note for why this is shared, not per-site).
# priors: 06's own PRIORS list (only decay_family/decay_t_df are used here --
#   mu_d/sigma_d themselves come from the ACTUAL fitted draws, not resimulated).
# n_grid: log_d quadrature grid resolution -- 200 is enough for a smooth 1D
#   marginal here (checked directly: doubling barely moves the mean/quantiles).
# thin_draws: how many of the real posterior draws to recompute this
#   counterfactual for.
# Returns list(draw_idx, d) -- d is a (length(draw_idx) x n_sites) matrix on
# the NATURAL (not log) scale, same shape/units as the d_<i> columns
# pooled_samples already carries, so a caller can plug it straight into the
# SAME country-aggregation code (rowMeans()-per-country, then quantiles
# across draws) 06_hierarchical_mcmc.R's own Section Q already uses.
log_prior_site_eir <- function(log_eir, bounds, priors) {
  # Normal(mu, sd) on the log (absolute) EIR scale -- a single FIXED prior,
  # identical for every site, rather than a hierarchical/pooled one: there's
  # no assumption that sites share information about each other's
  # transmission intensity (unlike decay's mu_d/sigma_d), so mu/sd are
  # constants reflecting the shared EIR_GRID's own span (06's own
  # PRIORS$eir_mu/eir_sd), not population-level hyperparameters fit from
  # data.
  #
  # This replaces an earlier version of this prior that was instead centred
  # PER SITE at that site's own external eir_in baseline (log_eir_rel = 0,
  # i.e. eir_rel = 1) -- not just a flat prior across the full grid range,
  # because a flat EIR prior did nothing to break a decay/EIR
  # non-identifiability that showed up strongly in the real fitted chains
  # under that relative-EIR model: 143/223 sites showed
  # |corr(log_d_i, log_eir_rel_i)| > 0.3 in their posterior samples (several
  # sites > 0.9), and the handful of sites with the most extreme, least-
  # constrained decay estimates (e.g. implausible 100-360 day half-lives)
  # were exactly the ones with the strongest such correlation. The per-site
  # eir_in anchor gave each site a DIFFERENT informative centre, which was
  # part of what broke that ridge -- this fixed, same-for-every-site prior
  # provides less site-distinguishing information than that did, so this is
  # a real risk to re-check post-fit (the |corr(log_d_i, log_eir_i)|
  # diagnostic), not something assumed to still be fixed.
  #
  # bounds STILL hard-truncates outside the emulator's own training range
  # (bounds: numeric(2) = c(log(EIR_GRID[1]), log(EIR_GRID[end]))) -- unlike
  # the Normal(mu, sd) term itself, this isn't a modelling choice, it's a
  # technical necessity: the GP emulator was never trained out there and
  # predict_emulator_ll() would be extrapolating blindly if allowed to
  # evaluate a point outside its own design.
  #
  # priors$eir_family (pipeline_fns.R "Prior family"; absent = "normal"):
  #   "normal":       Normal(eir_mu, eir_sd), hard-truncated to bounds as above.
  #   "logit_normal": logit-normal on priors$eir_bounds (= bounds; checked in
  #                   build_and_tune_mcmc_inputs()), no truncation needed.
  if (identical(priors$eir_family, "logit_normal")) {
    return(log_dlogitnorm_bounded(log_eir, priors$eir_bounds, priors$eir_logit_mean, priors$eir_logit_sd))
  }
  if (log_eir < bounds[1] || log_eir > bounds[2]) return(-Inf)
  dnorm(log_eir, mean = priors$eir_mu, sd = priors$eir_sd, log = TRUE)
}



compute_log_posterior <- function(theta, emulators, log_eir_bounds, log_d_bounds, priors) {
  n_sites    <- length(emulators)
  site_names <- names(emulators)

  mu_d    <- theta["mu_d"]
  sigma_d <- theta["sigma_d"]
  log_d   <- theta[paste0("log_d_", seq_len(n_sites))]
  log_eir <- theta[paste0("log_eir_", seq_len(n_sites))]

  lp <- log_prior_mu_d(mu_d, priors) +
    log_prior_sigma_d(sigma_d, priors$sigma_d_prior_mean, priors$sigma_d_prior_sd)
  if (!is.finite(lp)) return(-Inf)

  for (i in seq_len(n_sites)) {
    sname <- site_names[i]
    lp_decay <- log_prior_site_decay(log_d[i], mu_d, sigma_d, log_d_bounds,
                                      family = priors$decay_family, df = priors$decay_t_df)
    if (!is.finite(lp_decay)) return(-Inf)
    lp_eir <- log_prior_site_eir(log_eir[i], log_eir_bounds, priors)
    if (!is.finite(lp_eir)) return(-Inf)
    ll_site <- predict_emulator_ll(emulators[[sname]], log_d[i], log_eir[i])
    lp <- lp + lp_decay + lp_eir + ll_site
  }
  lp
}


run_mcmc <- function(
  emulators, log_eir_bounds, log_d_bounds,
  init_theta, proposal_sds, priors,
  n_iter = 10000, n_warmup = 5000, thin = 1, seed = 123
) {
  set.seed(seed)
  n_sites    <- length(emulators)
  site_names <- names(emulators)
  n_total    <- n_iter + n_warmup
  n_retain   <- ceiling(n_iter / thin)

  lp_curr <- compute_log_posterior(init_theta, emulators, log_eir_bounds, log_d_bounds, priors)
  if (!is.finite(lp_curr)) {
    stop("Initial theta has -Inf log-posterior. Check init_theta and log_eir_bounds/log_d_bounds.")
  }

  mu_d    <- unname(init_theta["mu_d"])
  sigma_d <- unname(init_theta["sigma_d"])
  log_d   <- unname(init_theta[paste0("log_d_", seq_len(n_sites))])
  log_eir <- unname(init_theta[paste0("log_eir_", seq_len(n_sites))])

  lp_mu_d    <- log_prior_mu_d(mu_d, priors)
  lp_sigma_d <- log_prior_sigma_d(sigma_d, priors$sigma_d_prior_mean, priors$sigma_d_prior_sd)
  lp_decay <- vapply(seq_len(n_sites), function(i) {
    log_prior_site_decay(log_d[i], mu_d, sigma_d, log_d_bounds, family = priors$decay_family, df = priors$decay_t_df)
  }, numeric(1))
  lp_eir <- vapply(seq_len(n_sites), function(i) {
    log_prior_site_eir(log_eir[i], log_eir_bounds, priors)
  }, numeric(1))
  ll_site <- vapply(seq_len(n_sites), function(i) {
    predict_emulator_ll(emulators[[site_names[i]]], log_d[i], log_eir[i])
  }, numeric(1))

  samples <- matrix(NA_real_, nrow = n_retain, ncol = length(init_theta) + 1L)
  colnames(samples) <- c(names(init_theta), "log_posterior")

  accept_mu_d    <- 0L
  accept_sigma_d <- 0L
  accept_decay   <- integer(n_sites)
  accept_eir     <- integer(n_sites)
  sample_idx     <- 0L

  mu_d_sd    <- unname(proposal_sds["mu_d"])
  sigma_d_sd <- unname(proposal_sds["sigma_d"])
  decay_sds  <- unname(proposal_sds[paste0("log_d_", seq_len(n_sites))])
  eir_sds    <- unname(proposal_sds[paste0("log_eir_", seq_len(n_sites))])

  for (iter in seq_len(n_total)) {

    # Block 1: population grand mean decay (mu_d) -- touches EVERY site's own
    # decay prior directly (Normal(mu_d, sigma_d)); zero emulator calls.
    mu_d_prop    <- mu_d + rnorm(1, 0, mu_d_sd)
    lp_mu_d_prop <- log_prior_mu_d(mu_d_prop, priors)
    if (is.finite(lp_mu_d_prop)) {
      lp_decay_prop <- vapply(seq_len(n_sites), function(i) {
        log_prior_site_decay(log_d[i], mu_d_prop, sigma_d, log_d_bounds, family = priors$decay_family, df = priors$decay_t_df)
      }, numeric(1))
      lp_prop <- lp_curr + (lp_mu_d_prop - lp_mu_d) + sum(lp_decay_prop - lp_decay)
      if (log(runif(1)) < lp_prop - lp_curr) {
        lp_curr  <- lp_prop
        mu_d     <- mu_d_prop
        lp_mu_d  <- lp_mu_d_prop
        lp_decay <- lp_decay_prop
        if (iter > n_warmup) accept_mu_d <- accept_mu_d + 1L
      }
    }

    # Block 2: population decay spread (sigma_d) -- the ONLY spread
    # parameter in this variant (no country layer); touches EVERY site's own
    # decay prior directly. Zero emulator calls.
    sigma_d_prop    <- sigma_d + rnorm(1, 0, sigma_d_sd)
    lp_sigma_d_prop <- log_prior_sigma_d(sigma_d_prop, priors$sigma_d_prior_mean, priors$sigma_d_prior_sd)
    if (is.finite(lp_sigma_d_prop)) {
      lp_decay_prop <- vapply(seq_len(n_sites), function(i) {
        log_prior_site_decay(log_d[i], mu_d, sigma_d_prop, log_d_bounds, family = priors$decay_family, df = priors$decay_t_df)
      }, numeric(1))
      lp_prop <- lp_curr + (lp_sigma_d_prop - lp_sigma_d) + sum(lp_decay_prop - lp_decay)
      if (log(runif(1)) < lp_prop - lp_curr) {
        lp_curr    <- lp_prop
        sigma_d    <- sigma_d_prop
        lp_sigma_d <- lp_sigma_d_prop
        lp_decay   <- lp_decay_prop
        if (iter > n_warmup) accept_sigma_d <- accept_sigma_d + 1L
      }
    }

    # Block 3: per-site decay (log_d_i) -- partial pooling around mu_d/
    # sigma_d directly (no country-specific spread). One emulator call/site.
    for (i in seq_len(n_sites)) {
      log_d_i_prop    <- log_d[i] + rnorm(1, 0, decay_sds[i])
      lp_decay_i_prop <- log_prior_site_decay(log_d_i_prop, mu_d, sigma_d, log_d_bounds, family = priors$decay_family, df = priors$decay_t_df)
      if (is.finite(lp_decay_i_prop)) {
        ll_site_i_prop <- predict_emulator_ll(emulators[[site_names[i]]], log_d_i_prop, log_eir[i])
        lp_prop <- lp_curr + (lp_decay_i_prop - lp_decay[i]) + (ll_site_i_prop - ll_site[i])
        if (log(runif(1)) < lp_prop - lp_curr) {
          log_d[i]    <- log_d_i_prop
          lp_decay[i] <- lp_decay_i_prop
          ll_site[i]  <- ll_site_i_prop
          lp_curr     <- lp_prop
          if (iter > n_warmup) accept_decay[i] <- accept_decay[i] + 1L
        }
      }
    }

    # Block 4: per-site EIR -- unchanged from the original model.
    for (i in seq_len(n_sites)) {
      log_eir_prop <- log_eir[i] + rnorm(1, 0, eir_sds[i])
      lp_eir_prop  <- log_prior_site_eir(log_eir_prop, log_eir_bounds, priors)
      if (is.finite(lp_eir_prop)) {
        ll_site_i_prop <- predict_emulator_ll(emulators[[site_names[i]]], log_d[i], log_eir_prop)
        lp_prop <- lp_curr + (lp_eir_prop - lp_eir[i]) + (ll_site_i_prop - ll_site[i])
        if (log(runif(1)) < lp_prop - lp_curr) {
          log_eir[i] <- log_eir_prop
          lp_eir[i]  <- lp_eir_prop
          ll_site[i] <- ll_site_i_prop
          lp_curr    <- lp_prop
          if (iter > n_warmup) accept_eir[i] <- accept_eir[i] + 1L
        }
      }
    }

    if (iter > n_warmup && (iter - n_warmup) %% thin == 0) {
      sample_idx <- sample_idx + 1L
      samples[sample_idx, ] <- c(mu_d, sigma_d, log_d, log_eir, lp_curr)
    }

    if (iter %% 1000 == 0) {
      message(sprintf("Iter %d/%d | lp = %.2f", iter, n_total, lp_curr))
    }
  }

  samples_dt <- as.data.table(samples[seq_len(sample_idx), , drop = FALSE])
  samples_dt[, mu_d_natural := exp(mu_d)]
  for (i in seq_len(n_sites)) {
    samples_dt[, paste0("d_", i)   := exp(get(paste0("log_d_", i)))]
    samples_dt[, paste0("eir_", i) := exp(get(paste0("log_eir_", i)))]
  }

  list(
    samples      = samples_dt,
    site_names   = site_names,
    accept_rates = list(
      mu_d    = accept_mu_d / n_iter,
      sigma_d = accept_sigma_d / n_iter,
      decay   = setNames(accept_decay / n_iter, site_names),
      eir     = setNames(accept_eir / n_iter, site_names)
    ),
    n_warmup = n_warmup,
    n_iter   = n_iter,
    thin     = thin
  )
}


initialise_theta <- function(best_by_name1, site_names, log_eir_bounds, log_d_bounds) {
  n_sites <- length(site_names)
  log_d_init <- log(best_by_name1$hrp2_scale)
  mu_d_init <- mean(log_d_init)
  sigma_d_init <- sd(log_d_init)
  if (!is.finite(sigma_d_init) || sigma_d_init <= 0) sigma_d_init <- 0.3

  log_eir_init <- setNames(log(best_by_name1$eir), paste0("log_eir_", seq_along(site_names)))

  JITTER_OFF_GRID <- 1e-3
  log_d_init   <- pmin(pmax(log_d_init + JITTER_OFF_GRID, log_d_bounds[1] + JITTER_OFF_GRID), log_d_bounds[2] - JITTER_OFF_GRID)
  log_eir_init <- pmin(pmax(log_eir_init + JITTER_OFF_GRID, log_eir_bounds[1] + JITTER_OFF_GRID), log_eir_bounds[2] - JITTER_OFF_GRID)
  names(log_d_init) <- paste0("log_d_", seq_len(n_sites))

  c(mu_d = mu_d_init, sigma_d = sigma_d_init, log_d_init, log_eir_init)
}


tune_proposal_sds <- function(emulators, log_eir_bounds, log_d_bounds, init_theta, priors,
                                              n_tune = 2000, target = 0.30, seed = 123) {
  set.seed(seed)
  site_names <- names(emulators)
  n_sites    <- length(site_names)
  sds <- setNames(rep(0.1, length(init_theta)), names(init_theta))

  message("Tuning proposal SDs (shared-sigma_d variant) over ", n_tune, " iterations (3 rounds)...")
  for (round in seq_len(3)) {
    res <- run_mcmc(
      emulators = emulators, log_eir_bounds = log_eir_bounds, log_d_bounds = log_d_bounds,
      init_theta = init_theta, proposal_sds = sds, priors = priors,
      n_iter = n_tune, n_warmup = 0L, seed = seed + round
    )
    ar <- res$accept_rates
    sds["mu_d"]    <- sds["mu_d"]    * ifelse(ar$mu_d    > target, 2, 0.5)
    sds["sigma_d"] <- sds["sigma_d"] * ifelse(ar$sigma_d > target, 2, 0.5)
    for (i in seq_len(n_sites)) {
      sds[paste0("log_d_", i)]   <- sds[paste0("log_d_", i)]   * ifelse(ar$decay[i] > target, 2, 0.5)
      sds[paste0("log_eir_", i)] <- sds[paste0("log_eir_", i)] * ifelse(ar$eir[i]   > target, 2, 0.5)
    }
    message(sprintf("  Round %d | mu_d AR = %.2f | sigma_d AR = %.2f | mean decay AR = %.2f | mean EIR AR = %.2f",
                     round, ar$mu_d, ar$sigma_d, mean(ar$decay), mean(ar$eir)))
  }
  sds
}


build_and_tune_mcmc_inputs <- function(grid_out_dir, emulator_out_dir,
                                                       n_tune = 2000L, target = 0.30, seed = 123,
                                                       priors, max_decay_halflife_days = NULL,
                                                       log_lik_file = "hrp2_grid_log_lik.RDS") {
  emulators  <- readRDS(file.path(emulator_out_dir, "emulators.rds"))
  site_names <- names(emulators)

  # log_lik_file: the grid log-likelihood of the same RDT sensitivity as the
  # emulators, so the starting values come from that variant's best grid points.
  lik_obj <- readRDS(file.path(grid_out_dir, log_lik_file))
  best_by_name1 <- data.table::data.table(lik_obj$best_by_name1)
  best_by_name1[, site_name := sprintf("%s_%s", iso, name_1)]
  best_by_name1 <- best_by_name1[match(site_names, site_name)]
  if (anyNA(best_by_name1$site_name) || !identical(best_by_name1$site_name, site_names)) {
    stop("best_by_name1 doesn't cover every emulator site (or ordering mismatch).")
  }

  log_eir_bounds <- range(unlist(lapply(emulators, function(em) range(em@X[, "log_eir"]))))
  log_d_bounds   <- range(unlist(lapply(emulators, function(em) range(em@X[, "log_d"]))))
  if (!is.null(max_decay_halflife_days)) {
    log_d_bounds[2] <- min(log_d_bounds[2], log(max_decay_halflife_days / log(2)))
  }

  # Logit-normal priors are defined on the grid's decay and EIR ranges, which
  # must be the ranges the emulators were trained on.
  check_bounds <- function(family, prior_bounds, em_bounds, label) {
    if (identical(family, "logit_normal") && !isTRUE(all.equal(prior_bounds, em_bounds, tolerance = 1e-6))) {
      stop(sprintf("priors$%s [%.4f, %.4f] != emulator bounds [%.4f, %.4f].",
                   label, prior_bounds[1], prior_bounds[2], em_bounds[1], em_bounds[2]))
    }
  }
  check_bounds(priors$mu_d_family, priors$mu_d_bounds, log_d_bounds, "mu_d_bounds")
  check_bounds(priors$eir_family, priors$eir_bounds, log_eir_bounds, "eir_bounds")

  init_theta <- initialise_theta(best_by_name1, site_names, log_eir_bounds, log_d_bounds)
  if (identical(priors$mu_d_family, "logit_normal")) {   # start strictly inside the prior's support
    init_theta["mu_d"] <- min(max(init_theta["mu_d"], priors$mu_d_bounds[1] + 1e-3), priors$mu_d_bounds[2] - 1e-3)
  }
  proposal_sds <- tune_proposal_sds(
    emulators, log_eir_bounds, log_d_bounds, init_theta, priors,
    n_tune = n_tune, target = target, seed = seed
  )

  list(
    site_names     = site_names,
    log_eir_bounds = log_eir_bounds,
    log_d_bounds   = log_d_bounds,
    priors         = priors,
    max_decay_halflife_days = max_decay_halflife_days,
    init_theta     = init_theta,
    proposal_sds   = proposal_sds
  )
}


run_one_chain <- function(emulator_out_dir, log_eir_bounds, log_d_bounds,
                                          priors, init_theta, proposal_sds, n_iter, n_warmup, thin, seed) {
  emulators <- readRDS(file.path(emulator_out_dir, "emulators.rds"))
  run_mcmc(
    emulators, log_eir_bounds, log_d_bounds,
    init_theta, proposal_sds, priors,
    n_iter = n_iter, n_warmup = n_warmup, thin = thin, seed = seed
  )
}


# --- Posterior summary ---

summarise_posterior <- function(mcmc_result, flagged_sites = character(0)) {
  samps      <- mcmc_result$samples
  site_names <- mcmc_result$site_names
  n_sites    <- length(site_names)

  summarise_col <- function(col) {
    vals <- samps[[col]]
    data.table(
      parameter = col,
      mean      = mean(vals),
      sd        = sd(vals),
      q025      = quantile(vals, 0.025),
      q975      = quantile(vals, 0.975)
    )
  }

  # mu_d_natural/sigma_d: population-level GRAND mean half-life scale and
  # the SINGLE population-level decay spread every site shares (no
  # country-specific layer under this model -- see this file's own header).
  # d_i: per-site decay estimate, shrunk towards mu_d by an amount governed
  # by sigma_d relative to that site's own likelihood precision. eir_i:
  # per-site absolute EIR estimate -- no pooling, so no equivalent
  # population-level shrinkage.
  shared_cols <- c("mu_d_natural", "sigma_d")
  site_cols   <- c(paste0("d_", seq_len(n_sites)), paste0("eir_", seq_len(n_sites)))

  summary_dt <- rbindlist(lapply(c(shared_cols, site_cols), summarise_col))

  # Attach site labels
  summary_dt[grepl("^(d|eir)_", parameter),
             site := site_names[as.integer(sub(".*_", "", parameter))]]

  # flagged: TRUE for site rows (d_i/eir_i) matching flagged_sites -- kept
  # for forward compatibility with any future exclusion mechanism; this
  # model doesn't currently implement one (see this file's header), so
  # flagged_sites is empty by default and every row comes back FALSE.
  summary_dt[, flagged := !is.na(site) & site %in% flagged_sites]

  summary_dt[]
}

compute_rhat <- function(chain_samples) {
  # chain_samples: list of numeric vectors, one per chain, all the SAME
  # length (post-warmup, already-thinned samples of a single parameter --
  # run_mcmc() only ever returns post-warmup draws, so no extra trimming
  # needed here). Standard (non-split) R-hat (Gelman & Rubin 1992): compares
  # between-chain and within-chain variance. Values notably above ~1.1
  # suggest the chains haven't mixed to the same distribution -- run longer
  # or re-tune proposal_sds if so.
  m <- length(chain_samples)
  n <- length(chain_samples[[1]])
  stopifnot(m >= 2, all(vapply(chain_samples, length, integer(1)) == n))

  chain_means <- vapply(chain_samples, mean, numeric(1))
  chain_vars  <- vapply(chain_samples, var,  numeric(1))

  grand_mean <- mean(chain_means)
  B <- n / (m - 1) * sum((chain_means - grand_mean)^2)
  W <- mean(chain_vars)

  var_hat <- ((n - 1) / n) * W + B / n
  sqrt(var_hat / W)
}
