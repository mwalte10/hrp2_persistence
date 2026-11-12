# =============================================================================
# src/export_posteriors.R
#
# ANALYSIS COMPONENT: OUTPUTS (posterior estimates; regression + calibration)
#
# Writes every posterior estimate reported in the manuscript to plain CSV in
# POSTERIOR_DIR (results/posteriors/, tracked in git). Reads the fitted models in
# OUTPUT_ROOT only; nothing is refitted. The exported files hold model
# parameters and model-derived quantities only, never survey microdata or
# site-file inputs, so they can be shared without DHS or site-file access.
#
# Section A, regression analysis (step 2, fit_brms_int_full_local.RDS):
#   results/posteriors/regression/
#     parameter_draws.csv.gz           every posterior draw of every fixed effect,
#                                      random-effect SD and survey random effect
#                                      (one row per draw; chain, iteration)
#     parameter_summary.csv            mean, SD, median, 2.5%, 97.5%, R-hat,
#                                      bulk and tail ESS for each parameter
#     fpr_fdr_by_covariate.csv         derived FPR/FDR summaries from step 2 Section B
#     fpr_fdr_average_marginal_effects.csv
#     fpr_fdr_full_covariate_grid.csv
#     loo_fixed_effect_selection.csv   ELPD table (fixed_effect_selection.R)
#
# Section B, calibration analysis (steps 3c and 4), one folder per fit:
#   results/posteriors/calibration/<fit>/, <fit> one of
#     main_rdt_se_090                  main analysis (RDT sensitivity 0.9, sigma_P prior SD 0.01)
#     sensitivity_rdt_se_080           RDT sensitivity 0.8
#     sensitivity_rdt_se_100           RDT sensitivity 1.0
#     sensitivity_sigma_sd_0.05        sigma_P prior SD 0.05 (RDT sensitivity 0.9)
#     sensitivity_sigma_sd_0.1         sigma_P prior SD 0.1  (RDT sensitivity 0.9)
#   each holding
#     population_draws.csv.gz          every pooled draw (chain, iteration) of mu_P
#                                      (log d), sigma_P, population median half-life
#     population_summary.csv           summaries of the above (get_population_parameters())
#     unobserved_site_summary.csv      predicted half-life for a site with no data
#     site_draws_thinned.csv.gz        N_SITE_DRAWS evenly thinned joint draws of each
#                                      site's d, half-life and EIR (long format)
#     site_summary.csv                 half-life and EIR per site (all draws)
#     country_summary.csv              half-life and EIR per country (sites weighted
#                                      by population at risk; needs the site files)
#     all_sites_summary.csv            the same across all sites
#     convergence.csv                  split R-hat and bulk/tail ESS per parameter
#     priors.csv                       prior settings of the fit
#   Full (unthinned) site-level draws stay in OUTPUT_ROOT (mcmc_combined.RDS).
#
# Scales: d is the mean duration of HRP2 positivity after clearance (days);
# half-life = d x ln 2. mu_P and sigma_P are on log(d) (manuscript notation;
# mu_d and sigma_d in the code). EIR is annual entomological inoculation rate.
#
# Run from the project root after steps 2-4:  Rscript src/export_posteriors.R
# =============================================================================
library(data.table)
source("./src/functions/paths.R")
source("./src/functions/pipeline_fns.R")
source("./src/functions/cluster_fns.R")             # RDT sensitivity settings, rdt_sensitivity_dir()
source("./src/functions/posterior_summary_fns.R")

OVR <- take_overrides()
EXPORT_REGRESSION  <- setting(OVR, "EXPORT_REGRESSION", TRUE)
EXPORT_CALIBRATION <- setting(OVR, "EXPORT_CALIBRATION", TRUE)
N_SITE_DRAWS       <- setting(OVR, "N_SITE_DRAWS", 1000L)   # thinned joint site draws per fit
AGG_WEIGHTS        <- "par_pf"                              # country summaries, as in plots.R / results_numbers.Rmd

dir.create(POSTERIOR_DIR, showWarnings = FALSE, recursive = TRUE)
write_csv <- function(x, dir, file) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fwrite(x, file.path(dir, file))
  cat(sprintf("  wrote %s (%d rows)\n", file.path(dir, file), nrow(x)))
}
q <- function(x, p) stats::quantile(x, p, names = FALSE)

# Regression ===================================================================
# SECTION A -- REGRESSION ANALYSIS POSTERIOR
if (EXPORT_REGRESSION) {
  cat("\n== Regression analysis ==\n")
  reg_out  <- file.path(POSTERIOR_DIR, "regression")
  unlink(reg_out, recursive = TRUE)   # the folder holds exactly what this run writes
  fit_path <- file.path(DIR_REGRESSION, "fit_brms_int_full_local.RDS")
  if (!file.exists(fit_path)) stop("No ", fit_path, " -- run step 2 first.")
  fit <- readRDS(fit_path)

  # Draws straight from the stanfit, so brms need not match the fitting version.
  # Keep model parameters only: fixed effects (b_), random-effect SDs (sd_) and
  # survey random effects (r_); drop lp__, lprior and brms internals (z_, Intercept copies).
  if (inherits(fit$fit, "stanfit")) loadNamespace("rstan")    # as.array() method for stanfit
  arr  <- as.array(fit$fit)                                    # iterations x chains x parameters
  keep <- grep("^(b_|sd_|r_)", dimnames(arr)[[3]], value = TRUE)
  arr  <- arr[, , keep, drop = FALSE]
  draws <- rbindlist(lapply(seq_len(dim(arr)[2]), function(ch) {
    d <- as.data.table(arr[, ch, ])
    d[, `:=`(chain = ch, iteration = seq_len(.N))]
  }))
  setcolorder(draws, c("chain", "iteration"))
  dir.create(reg_out, showWarnings = FALSE, recursive = TRUE)
  fwrite(draws, file.path(reg_out, "parameter_draws.csv.gz"))
  cat(sprintf("  wrote %s (%d draws x %d parameters)\n", file.path(reg_out, "parameter_draws.csv.gz"), nrow(draws), length(keep)))

  reg_summary <- rbindlist(lapply(keep, function(p) {
    m <- arr[, , p]
    data.table(parameter = p, mean = mean(m), sd = sd(m), median = median(m),
               q025 = q(m, 0.025), q975 = q(m, 0.975),
               rhat = posterior::rhat(m), ess_bulk = posterior::ess_bulk(m), ess_tail = posterior::ess_tail(m))
  }))
  write_csv(reg_summary, reg_out, "parameter_summary.csv")

  # Derived FPR/FDR quantities from step 2 Section B (already posterior summaries)
  derived <- c(fpfn_by_covariate.RDS = "fpr_fdr_by_covariate.csv",
               fpfn_ame.RDS          = "fpr_fdr_average_marginal_effects.csv",
               fpfn_full_grid.RDS    = "fpr_fdr_full_covariate_grid.csv")
  for (f in names(derived)) {
    p <- file.path(DIR_REGRESSION, f)
    if (!file.exists(p)) { cat(sprintf("  missing %s -- run step 2 Section B\n", f)); next }
    x <- readRDS(p)
    if (!is.data.frame(x)) x <- rbindlist(Filter(is.data.frame, x), idcol = "component", fill = TRUE)
    x <- as.data.table(x)
    list_cols <- names(x)[vapply(x, is.list, logical(1))]
    if (length(list_cols)) x[, (list_cols) := NULL]   # per-draw vectors stay in the RDS
    write_csv(x, reg_out, derived[[f]])
  }
  loo_tab <- file.path(RESULTS_DIR, "table_loo_fixed_effect_selection.csv")
  if (file.exists(loo_tab)) write_csv(fread(loo_tab), reg_out, "loo_fixed_effect_selection.csv")
  rm(fit, arr, draws); gc()
}

# Calibration ==================================================================
# SECTION B -- CALIBRATION ANALYSIS POSTERIOR, MAIN FIT AND SENSITIVITY FITS
if (EXPORT_CALIBRATION) {
  fits <- unique(rbind(
    data.table(rdt_sensitivity = RDT_SENSITIVITY_MAIN, sigma_d_prior_sd = SIGMA_D_PRIOR_SD_MAIN),
    mcmc_variants()))
  fits[, label := mapply(function(se, sd) {
    if (!is_main_sigma_sd(sd)) return(sprintf("sensitivity_sigma_sd_%s", format(sd)))
    sprintf("%s_%s", if (abs(se - RDT_SENSITIVITY_MAIN) < 1e-9) "main" else "sensitivity", rdt_sensitivity_tag(se))
  }, rdt_sensitivity, sigma_d_prior_sd)]

  unlink(file.path(POSTERIOR_DIR, "calibration"), recursive = TRUE)
  for (i in seq_len(nrow(fits))) {
    f <- fits[i]
    mcmc_dir <- mcmc_variant_dir(f$rdt_sensitivity, f$sigma_d_prior_sd)
    out <- file.path(POSTERIOR_DIR, "calibration", f$label)
    cat(sprintf("\n== Calibration: %s (%s) ==\n", f$label, mcmc_dir))
    if (!file.exists(file.path(mcmc_dir, "mcmc_combined.RDS"))) { cat("  no mcmc_combined.RDS -- skipped\n"); next }

    mc   <- readRDS(file.path(mcmc_dir, "mcmc_combined.RDS"))
    post <- load_hierarchical_posterior(mcmc_dir, DIR_GRID, mcmc_combined = mc)
    ps   <- post$pooled_samples
    if (!"iteration" %in% names(ps)) ps[, iteration := seq_len(.N), by = chain]

    # Population level: every draw
    pop <- ps[, .(chain, iteration, mu_P_log_d = mu_d, sigma_P = sigma_d,
                  halflife_population_median = log_d_to_halflife(mu_d))]
    dir.create(out, showWarnings = FALSE, recursive = TRUE)
    fwrite(pop, file.path(out, "population_draws.csv.gz"))
    cat(sprintf("  wrote %s (%d draws)\n", file.path(out, "population_draws.csv.gz"), nrow(pop)))
    write_csv(get_population_parameters(post), out, "population_summary.csv")
    write_csv(predict_unobserved_site(post)$summary, out, "unobserved_site_summary.csv")

    # Site level: thinned joint draws (long) and summaries from all draws
    keep_rows <- unique(round(seq(1, nrow(ps), length.out = min(N_SITE_DRAWS, nrow(ps)))))
    d_mat   <- exp(site_param_draws(post, "log_d"))[keep_rows, , drop = FALSE]
    eir_mat <- exp(site_param_draws(post, "log_eir"))[keep_rows, , drop = FALSE]
    site_draws <- data.table(draw = rep(seq_along(keep_rows), times = ncol(d_mat)),
                             chain = rep(ps$chain[keep_rows], times = ncol(d_mat)),
                             iteration = rep(ps$iteration[keep_rows], times = ncol(d_mat)),
                             site = rep(post$site_names, each = length(keep_rows)),
                             d = as.vector(d_mat), eir = as.vector(eir_mat))
    site_draws[, `:=`(iso = sub("_.*$", "", site), halflife = d_to_halflife(d))]
    setcolorder(site_draws, c("draw", "chain", "iteration", "site", "iso", "d", "halflife", "eir"))
    fwrite(site_draws, file.path(out, "site_draws_thinned.csv.gz"))
    cat(sprintf("  wrote %s (%d draws x %d sites)\n", file.path(out, "site_draws_thinned.csv.gz"),
                length(keep_rows), length(post$site_names)))
    write_csv(summarise_site_parameters(post, by = "site"), out, "site_summary.csv")

    # Country and all-site summaries weight sites by population at risk, read
    # from the site files (SITE_CACHE_DIR); skipped without site-file access.
    agg <- tryCatch(list(country = summarise_site_parameters(post, by = "country", weights = AGG_WEIGHTS),
                         all     = summarise_site_parameters(post, by = "all",     weights = AGG_WEIGHTS)),
                    error = function(e) { cat(sprintf("  country summaries skipped: %s\n", conditionMessage(e))); NULL })
    if (!is.null(agg)) {
      write_csv(agg$country, out, "country_summary.csv")
      write_csv(agg$all, out, "all_sites_summary.csv")
    }

    # Convergence across chains, on the sampled (log) scale
    conv_par <- c("mu_d", "sigma_d", grep("^log_(d|eir)_[0-9]+$", names(ps), value = TRUE))
    conv <- rbindlist(lapply(conv_par, function(p) {
      m <- sapply(mc$chains, function(ch) ch$samples[[p]])
      data.table(parameter = p, rhat = posterior::rhat(m), ess_bulk = posterior::ess_bulk(m),
                 ess_tail = posterior::ess_tail(m))
    }))
    idx <- suppressWarnings(as.integer(sub("^log_(d|eir)_", "", conv$parameter)))
    conv[, site := fifelse(is.na(idx), NA_character_, post$site_names[idx])]
    write_csv(conv, out, "convergence.csv")

    pr <- if (!is.null(mc$priors)) mc$priors else post$priors   # older fits keep priors in mcmc_inputs.RDS only
    pr <- pr[vapply(pr, function(x) is.atomic(x) && length(x) == 1, logical(1))]
    write_csv(data.table(prior = names(pr), value = vapply(pr, as.character, "")), out, "priors.csv")
    rm(mc, post, ps, site_draws); gc()
  }
}
cat(sprintf("\nPosterior exports in %s\n", POSTERIOR_DIR))
