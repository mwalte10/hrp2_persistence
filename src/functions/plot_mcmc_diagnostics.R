# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: prior simulation and posterior predictive helpers)
#
# plot_mcmc_diagnostics.R
# Per-site vetting page for 06_hierarchical_mcmc.R's chain output: for each
# of the model's 4 parameters (2 population-level, 2 site-level -- see
# hierarchical_mcmc_fns.R's own header) show (a) the prior this model
# actually placed on that parameter and (b) the pooled posterior next to
# it, PLUS (c) a trace plot (one line per chain) so a reviewer can see
# whether any chain is stuck/not mixing. Every plotted quantity is on its
# NATURAL (not log) scale -- d (days-scale HRP2 clearance), eir (absolute),
# mu_d_natural, sigma_d -- since that's what a domain reviewer actually
# reasons about; the model itself samples the log-scale versions where
# relevant (see hierarchical_mcmc_fns.R), so priors are simulated on the
# log scale and exponentiated (Monte Carlo naturally handles the
# change-of-variables correctly, no Jacobian needed) -- EXCEPT sigma_d,
# which the model already samples on its own natural scale directly (see
# log_prior_sigma_d()'s own header for why), so it needs no transform at
# all here.
#
# PROMOTED 2026-09-30 from this session's "shared_sigma_d" sensitivity
# variant: every site shares ONE population sigma_d directly (no
# country-specific spread layer at all) -- see hierarchical_mcmc_fns.R's
# own header for the full model and the rationale for retiring the earlier
# 4-level (mu_d -> sigma_d -> tau_sigma_d -> sigma_d_c[c] -> log_d_i)
# hierarchy. mu_d/sigma_d are population-level -- their posterior is
# IDENTICAL on every site's page (not a bug: it's the same pooled draws
# regardless of which site's page is being built). d_i and eir_i are
# site-specific. The retired country-spread version of this file's own
# logic (6-panel layout, country_idx-indexed sigma_d_c) is preserved,
# unchanged, wherever the retired hierarchical_mcmc_fns_country_spread.R
# model is still fit (09d) -- not duplicated here.
#
# The d_i PRIOR shown is the MARGINAL (hierarchical) prior -- this means
# simulating mu_d/sigma_d from THEIR OWN priors first, then log_d_i | mu_d,
# sigma_d ~ Normal OR StudentT(df) location-scale, truncated to
# log_d_bounds, i.e. integrating over the full hyperprior -- the correct
# prior to compare a single site's posterior against, since log_d_i's prior
# is genuinely conditional on mu_d/sigma_d, not a fixed distribution on its
# own. WHICH family + truncation must match hierarchical_mcmc_fns.R::
# log_prior_site_decay() exactly -- see that function's own header for the
# full history. simulate_hierarchical_priors() below takes `priors` (06's
# own PRIORS list) as a required argument for exactly this reason -- there's
# no safe default to fall back on that's guaranteed to match a given run
# (see hierarchical_mcmc_fns.R's own header for why prior hyperparameters
# live in 06, not as constants duplicated across files).
# eir_i's prior has no hierarchy/pooling at all (log_prior_site_eir() is
# identical, fixed Normal(priors$eir_mu, priors$eir_sd) truncated to
# log_eir_bounds, for every site) but is still simulated rather than solved
# analytically, for the same reason -- one consistent mechanism (Monte Carlo
# + exp()) for every panel rather than several different ones.

library(ggplot2)
library(data.table)
library(patchwork)

`%||%` <- function(a, b) if (is.null(a)) b else a

# Panel titles here are long, descriptive parameter labels (e.g. "mu_d
# (population grand mean decay, shared by every site)") squeezed above
# narrow columns (4-5 panels across a page), wrapping to a fixed character
# width, rather than letting ggplot lay the string out on one line, is what
# stops adjacent panels' titles from running into each other.
wrap_label <- function(x, width = 20) paste(strwrap(x, width = width), collapse = "\n")

MCMC_CHAIN_COLOURS <- c("#2a78d6", "#e67e22", "#27ae60", "#8e44ad")
MCMC_PRIOR_POSTERIOR_COLOURS <- c(Prior = "#b8b6ae", Posterior = "#2a78d6")
MCMC_POSTPRED_COLOUR <- "#2a78d6"
MCMC_OBS_COLOUR       <- "#c0392b"

mcmc_vetting_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
    plot.title        = element_text(colour = "#0b0b0b", size = 8, lineheight = 0.9),
    plot.subtitle      = element_text(colour = "#52514e", size = 7.5),
    axis.text          = element_text(colour = "#898781", size = 7),
    axis.title          = element_text(colour = "#52514e", size = 8, lineheight = 0.9),
    legend.position     = "bottom",
    legend.text         = element_text(size = 7)
  )

# Simulate the hierarchical model's OWN priors, on the natural scale, once
# (shared across every site's page -- see this file's own header for why
# d_i's marginal prior is the SAME distribution for every site before
# conditioning on that site's own data). n_sim: Monte Carlo draw count:
# large enough for a smooth density curve without being slow to
# simulate/plot. priors: 06_hierarchical_mcmc.R's own PRIORS list
# (mu_d_mean/sd, sigma_d_prior_mean/sigma_d_prior_sd, eir_mu/sd,
# decay_family/decay_t_df) -- required, no default, for the identical
# reason family/df used to be required directly: there's no safe fallback
# guaranteed to match whatever a given run was actually fit with.
# mu_d and EIR draws follow priors$mu_d_family / priors$eir_family, as
# hierarchical_mcmc_fns.R::log_prior_mu_d()/log_prior_site_eir()
# (absent = "normal"; "logit_normal" = Normal on logit((x - L) / (U - L))).
r_mu_d_prior <- function(n, priors) {
  if (!identical(priors$mu_d_family, "logit_normal")) return(rnorm(n, mean = priors$mu_d_mean, sd = priors$mu_d_sd))
  priors$mu_d_bounds[1] + diff(priors$mu_d_bounds) * plogis(rnorm(n, mean = priors$mu_d_logit_mean, sd = priors$mu_d_logit_sd))
}
simulate_hierarchical_priors <- function(log_eir_bounds, log_d_bounds, priors,
                                          n_sim = 20000L, seed = 123L) {
  set.seed(seed)
  mu_d_draws <- r_mu_d_prior(n_sim, priors)
  # sigma_d is sampled DIRECTLY on its natural scale by the model itself
  # (truncated Normal, mode at priors$sigma_d_prior_mean -- see
  # log_prior_sigma_d()'s own header for why NOT log-transformed/Half-Normal),
  # so simulating it here is just a truncated-at-zero Normal draw via
  # rejection sampling (not abs() -- sigma_d's mode is away from zero, so
  # that trick doesn't apply the way it would for a mode-at-zero prior).
  sigma_d_raw <- rnorm(n_sim * 3L, mean = priors$sigma_d_prior_mean, sd = priors$sigma_d_prior_sd)
  sigma_d_in_bounds <- sigma_d_raw[sigma_d_raw >= 0]
  if (length(sigma_d_in_bounds) < n_sim) {
    stop("simulate_hierarchical_priors(): sigma_d rejection sampling under-filled -- widen the oversample factor.")
  }
  sigma_d_draws <- sigma_d_in_bounds[seq_len(n_sim)]

  # log_d_i | mu_d, sigma_d matches hierarchical_mcmc_fns.R::
  # log_prior_site_decay() exactly -- hard-truncated to log_d_bounds.
  # Rejection sampling needs FRESH mu_d/sigma_d draws at oversampled length,
  # not the n_sim-length draws above -- those must stay exactly n_sim long
  # since they're returned as-is for their own marginal-prior panels.
  n_oversample     <- n_sim * 3L
  mu_d_over        <- r_mu_d_prior(n_oversample, priors)
  sigma_d_over_raw <- rnorm(n_oversample, mean = priors$sigma_d_prior_mean, sd = priors$sigma_d_prior_sd)
  sigma_d_over     <- sample(sigma_d_over_raw[sigma_d_over_raw >= 0], n_oversample, replace = TRUE)
  raw_log_d <- if (priors$decay_family == "t") {
    mu_d_over + sigma_d_over * rt(n_oversample, df = priors$decay_t_df)
  } else if (priors$decay_family == "normal") {
    rnorm(n_oversample, mean = mu_d_over, sd = sigma_d_over)
  } else {
    stop(sprintf("simulate_hierarchical_priors(): unknown family '%s' -- must be 'normal' or 't'.", priors$decay_family))
  }
  in_bounds_d  <- raw_log_d[raw_log_d >= log_d_bounds[1] & raw_log_d <= log_d_bounds[2]]
  if (length(in_bounds_d) < n_sim) {
    stop("simulate_hierarchical_priors(): decay rejection sampling under-filled -- widen the oversample factor.")
  }
  log_d_draws <- in_bounds_d[seq_len(n_sim)]

  # Rejection sampling for the truncated Normal(priors$eir_mu, priors$eir_sd)
  # -- oversample then trim to bounds and take the first n_sim. The
  # "logit_normal" family (priors$eir_family) needs no truncation.
  if (identical(priors$eir_family, "logit_normal")) {
    log_eir_draws <- priors$eir_bounds[1] + diff(priors$eir_bounds) *
      plogis(rnorm(n_sim, mean = priors$eir_logit_mean, sd = priors$eir_logit_sd))
  } else {
    raw_eir <- rnorm(n_sim * 3L, mean = priors$eir_mu, sd = priors$eir_sd)
    in_bounds <- raw_eir[raw_eir >= log_eir_bounds[1] & raw_eir <= log_eir_bounds[2]]
    if (length(in_bounds) < n_sim) {
      stop("simulate_hierarchical_priors(): rejection sampling under-filled -- widen the oversample factor.")
    }
    log_eir_draws <- in_bounds[seq_len(n_sim)]
  }

  list(
    mu_d_natural = exp(mu_d_draws),
    sigma_d      = sigma_d_draws,
    d            = exp(log_d_draws),
    eir          = exp(log_eir_draws)
  )
}

# One (trace, prior-vs-posterior) panel PAIR for a single parameter.
# posterior_dt: data.table with columns `chain`, `idx` (retained-sample
# order within its own chain -- not a real MCMC iteration count, since only
# post-warmup/thinned rows are ever stored, but monotonic order is all a
# stuck/non-mixing check needs), `value`.
# prior_draws: numeric vector, this parameter's own simulated prior (natural
# scale) from simulate_hierarchical_priors() above.
build_mcmc_param_panels <- function(posterior_dt, prior_draws, param_label, x_label = NULL) {
  # All four of this model's parameters are log-normal-shaped even on their
  # "natural" scale (mu_d_natural = exp(mu_d), sigma_d = exp(log_sigma_d), d
  # and eir likewise) -- a linear axis lets rare, large prior draws (e.g.
  # sigma_d's own prior has a long right tail) squash the bulk of both prior
  # and posterior into an unreadable spike near zero (confirmed directly:
  # d's prior draws alone range from ~1e-7 to ~1e8). log10 on both the
  # trace's y-axis and the density's x-axis fixes this consistently.
  # scale_x_log10()/scale_y_log10() transform the DATA before stat_density()/
  # the line geometry are computed (not just the rendered axis), so this is
  # the statistically correct way to get a log-scale density, not a
  # post-hoc-relabelled linear one.
  p_trace <- ggplot(posterior_dt, aes(idx, value, colour = factor(chain))) +
    geom_line(linewidth = 0.3, alpha = 0.8) +
    scale_colour_manual(values = MCMC_CHAIN_COLOURS, name = "Chain") +
    scale_y_log10() +
    labs(title = wrap_label(param_label), x = "Retained sample", y = x_label %||% param_label) +
    mcmc_vetting_theme

  prior_post_dt <- data.table::rbindlist(list(
    data.table::data.table(value = prior_draws, series = "Prior"),
    data.table::data.table(value = posterior_dt$value, series = "Posterior")
  ))
  prior_post_dt[, series := factor(series, levels = c("Prior", "Posterior"))]

  # Zoomed to the combined prior+posterior's own 0.1st-99.9th percentile
  # range -- some priors' extreme tails (e.g. d's, driven by rare large
  # sigma_d hyperprior draws) otherwise stretch the axis out so far that both
  # curves' actual bulk collapses into an unreadable sliver. coord_cartesian()
  # (not scale_x_log10(limits=...)) zooms the VIEW without dropping any of
  # the density estimate's own underlying data.
  zoom_range <- stats::quantile(prior_post_dt$value, c(0.001, 0.999))

  p_density <- ggplot(prior_post_dt, aes(value, fill = series)) +
    geom_density(alpha = 0.45, colour = NA) +
    scale_fill_manual(values = MCMC_PRIOR_POSTERIOR_COLOURS, name = NULL) +
    scale_x_log10() +
    coord_cartesian(xlim = zoom_range) +
    labs(title = NULL, x = wrap_label(x_label %||% param_label, width = 24), y = "Density") +
    mcmc_vetting_theme

  list(trace = p_trace, density = p_density)
}

# Binomial-resampling CI on one site's own OBSERVED survey count: draw
# n_resample binomial samples at that survey's own (p, n) -- i.e. "what
# would repeating this exact survey look like" -- then take the empirical
# 2.5th/97.5th percentile of the resampled proportions. Same "resample the
# binomial noise a survey of this size carries" logic as 06_hierarchical_
# mcmc.R's own binom_resampled_quantiles() (see that function's own
# header), just applied here to a single OBSERVED (p, n) pair rather than a
# vector of the MODEL's own posterior draws. Vectorised over p/n (one
# resample-and-quantile per row, via Map()) since data.table's `:=` passes
# the whole filtered column at once, not row-by-row.
# `deff` (scalar or per-row): the survey's design effect -- resample at the
# effective size n / deff (see binom_resampled_quantiles()'s header); 1 =
# original behaviour.
binom_resampled_ci <- function(p, n, n_resample = 1000L, deff = 1) {
  deff <- rep_len(deff, length(n))
  deff[is.na(deff)] <- 1
  qs <- Map(function(p_i, n_i, d_i) {
    n_eff <- max(1, round(n_i / d_i))
    resampled <- rbinom(n_resample, size = n_eff, prob = p_i) / n_eff
    stats::quantile(resampled, c(0.025, 0.975))
  }, p, n, deff)
  list(lo = vapply(qs, `[[`, numeric(1), 1), hi = vapply(qs, `[[`, numeric(1), 2))
}


# --- Posterior predictive intervals for the joint RDT x LM cells, all sites --
# joint_ppi_check(): observed TP/FP/FN/TN shares (of children tested by both
# RDT and LM, from the MIS joint counts) merged onto a posterior predictive
# summary's tp/fp/fn/tn intervals (06 Section I), one row per site-year x cell.
joint_ppi_check <- function(ppi_summary, mis) {
  obs <- unique(data.table::as.data.table(mis)[!is.na(joint_n) & joint_n > 0,
           .(iso = iso3c, name_1, year, joint_n, TP = joint_tp / joint_n, FP = joint_fp / joint_n,
             FN = joint_fn / joint_n, TN = joint_tn / joint_n)])
  obs <- data.table::melt(obs, id.vars = c("iso", "name_1", "year", "joint_n"), variable.name = "cell", value.name = "observed")
  ppi <- data.table::as.data.table(ppi_summary)
  pred <- data.table::rbindlist(lapply(c("TP", "FP", "FN", "TN"), function(cl) {
    k <- tolower(cl)
    ppi[, .(iso, name_1, year, cell = cl, mean = get(paste0(k, "_mean")), lo = get(paste0(k, "_lo")), hi = get(paste0(k, "_hi")))]
  }))
  out <- merge(obs, pred, by = c("iso", "name_1", "year", "cell"))
  out[, inside := observed >= lo & observed <= hi]
  out[, cell := factor(cell, levels = c("TP", "FP", "FN", "TN"))]
  out[]
}

# One figure across all sites: observed share (x) vs posterior predictive mean
# and 95% interval (y), one panel per cell (columns) and RDT sensitivity (rows
# when more than one). Each panel has equal x and y ranges, so the dashed
# line is y = x; points inside their interval are blue, outside gold, and each
# panel states its coverage.
plot_joint_ppi_all_sites <- function(chk) {
  cell_lab <- c(TP = "TP (RDT+/LM+)", FP = "FP (RDT+/LM-)", FN = "FN (RDT-/LM+)", TN = "TN (RDT-/LM-)")
  chk <- data.table::copy(chk)[, cell_label := factor(cell_lab[as.character(cell)], levels = cell_lab)]
  chk[, status := factor(data.table::fifelse(inside, "Inside 95% interval", "Outside 95% interval"),
                         levels = c("Inside 95% interval", "Outside 95% interval"))]
  multi <- data.table::uniqueN(chk$rdt_sensitivity) > 1
  facet_vars <- if (multi) c("rdt_sensitivity", "cell_label") else "cell_label"
  # Equal x/y range within each panel, and the same range for a cell in every
  # sensitivity row (invisible points at the cell's min/max across all rows),
  # so a column can be compared top to bottom.
  lims <- chk[, .(lo = min(c(observed, lo)), hi = max(c(observed, hi))), by = cell_label]
  if (multi) lims <- lims[data.table::CJ(cell_label = lims$cell_label, rdt_sensitivity = unique(chk$rdt_sensitivity)), on = "cell_label"]
  blank <- data.table::rbindlist(list(lims[, `:=`(x = lo, y = lo)][], data.table::copy(lims)[, `:=`(x = hi, y = hi)][]), fill = TRUE)
  cov <- chk[, .(label = sprintf(" Coverage %s%%\n %s site-years", formatC(100 * mean(inside), format = "f", digits = 0),
                                 formatC(.N, big.mark = ","))), by = facet_vars]
  ggplot2::ggplot(chk, ggplot2::aes(x = observed, y = mean)) +
    ggplot2::geom_blank(data = blank, ggplot2::aes(x = x, y = y), inherit.aes = FALSE) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi, colour = status), alpha = 0.35, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = status), size = 0.8, alpha = 0.7) +
    ggplot2::geom_text(data = cov, ggplot2::aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1.15,
                       size = 3.3, lineheight = 0.95, inherit.aes = FALSE) +
    {if (multi) ggplot2::facet_wrap(ggplot2::vars(rdt_sensitivity, cell_label), ncol = 4, scales = "free",
                                    labeller = ggplot2::labeller(.multi_line = FALSE))
     else ggplot2::facet_wrap(ggplot2::vars(cell_label), nrow = 1, scales = "free")} +
    ggplot2::scale_colour_manual(values = c("Inside 95% interval" = "#156082", "Outside 95% interval" = "#E1AF00"), name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(x = "Observed share of children tested by both (MIS)",
                  y = "Posterior predictive share\n(mean and 95% interval)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
                   panel.spacing.x = ggplot2::unit(1.2, "lines"),
                   legend.position = "bottom")   # TODO: confirm legend position with Maggie
}

