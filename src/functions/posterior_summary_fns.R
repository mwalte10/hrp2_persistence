# =============================================================================
# src/functions/posterior_summary_fns.R
#
# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: posterior summaries)
#
# Helpers for pulling parameter summaries out of the primary hierarchical fit
# (06_hierarchical_mcmc.R -> outputs/06_hierarchical_fitting/mcmc_combined.RDS,
# shared_sigma_d model: mu_d -> sigma_d -> log_d_i, plus unpooled log_eir_i).
# Also works unchanged on the sigma_d prior-width sensitivity refits (09i/09j),
# which share the same schema. NOT for the retired country-spread fit
# (sensitivity_country_spread/), which has extra sigma_d_c[c] parameters.
#
# Scales (see hierarchical_mcmc_fns.R's own header):
#   - d   = malariasimulation hrp2_scale = MEAN duration of HRP2 positivity
#           (days); half-life = d * log(2) (cluster_fns.R::scale_to_halflife()).
#   - mu_d, log_d_i are on log(d); sigma_d is the SD on log(d).
#   - log_eir_i is log annual EIR. EIR is NOT pooled -- there is no
#     population-level EIR, so an unobserved site's EIR is just the prior.
#
# Country-level and all-site values: the model has no country-level
# parameters of its own, so a country's value is built per posterior draw as
# a WEIGHTED GEOMETRIC MEAN of its sites (weighted mean of log_d / log_eir,
# then exp()), and summarised across draws afterwards. Geometric because the
# model pools on log(d) -- an arithmetic mean is dominated by a country's few
# longest-tailed sites. Per-draw (not averaging per-site intervals) because
# sites share each draw's mu_d/sigma_d, so they are not independent.
#
# Weights have NO default anywhere -- every caller must choose "none" (equal
# weights), "par_pf" (population at risk, averaged over that site's own MIS
# survey years), or a named numeric vector over site_names.
#
# Typical use:
#   post <- load_hierarchical_posterior(MCMC_OUT_DIR, GRID_OUT_DIR)
#   summarise_site_parameters(post, by = "country", weights = "par_pf")
#   summarise_site_parameters(post, by = "all",     weights = "none")
#   summarise_site_parameters(post, by = "site")
#   summarise_site_spread(post, by = "country", weights = "par_pf")   # range ACROSS sites
#   get_population_parameters(post)
#   predict_unobserved_site(post)$summary
# =============================================================================

# Full country names for the 16 ISO3 codes that can appear in the fitted
# site set -- shared by 06 Section Q/R2 and 09_results_summary.R (previously
# duplicated in both). Extend here if a new country is added.
COUNTRY_NAME_LOOKUP <- c(
  BDI = "Burundi", KEN = "Kenya", MDG = "Madagascar", MWI = "Malawi",
  RWA = "Rwanda", UGA = "Uganda", BEN = "Benin", BFA = "Burkina Faso",
  GHA = "Ghana", GIN = "Guinea", LBR = "Liberia", MLI = "Mali",
  NGA = "Nigeria", SEN = "Senegal", TGO = "Togo", CIV = "Cote d'Ivoire"
)

source("./src/functions/paths.R")   # SITE_CACHE_DIR, DIR_STRUCTURE

# --- Loading -----------------------------------------------------------------

# load_hierarchical_posterior(): everything the functions below need, in one
# list. `mcmc_combined` can be passed in when it's already loaded (it's
# ~180MB) to avoid a second readRDS(). Site -> ISO mapping comes from
# adm1_combos_with_paths.RDS, the same authoritative source 06 Section Q uses,
# aligned 1:1 with site_names (index i in log_d_i <-> site_names[i]).
load_hierarchical_posterior <- function(mcmc_dir, grid_out_dir, mcmc_combined = NULL) {
  if (is.null(mcmc_combined)) {
    mcmc_combined <- readRDS(file.path(mcmc_dir, "mcmc_combined.RDS"))
  }
  mcmc_inputs <- readRDS(file.path(mcmc_dir, "mcmc_inputs.RDS"))
  site_names  <- mcmc_combined$chains[[1]]$site_names
  if (!"mu_d" %in% names(mcmc_combined$pooled_samples) ||
      any(grepl("^sigma_d_c", names(mcmc_combined$pooled_samples)))) {
    stop("load_hierarchical_posterior(): pooled_samples isn't the shared_sigma_d schema (country-spread fits aren't supported).")
  }

  adm1_combos <- data.table::as.data.table(readRDS(file.path(grid_out_dir, "adm1_combos_with_paths.RDS")))
  adm1_combos[, site_name := sprintf("%s_%s", iso, name_1)]
  adm1_combos <- adm1_combos[match(site_names, site_name)]
  if (anyNA(adm1_combos$site_name)) {
    stop(sprintf("load_hierarchical_posterior(): no adm1_combos entry for fitted site(s): %s",
                 paste(site_names[is.na(adm1_combos$site_name)], collapse = ", ")))
  }

  list(
    pooled_samples = mcmc_combined$pooled_samples,
    site_names     = site_names,
    site_iso       = adm1_combos$iso,
    site_name_1    = adm1_combos$name_1,
    log_d_bounds   = mcmc_inputs$log_d_bounds,
    priors         = mcmc_inputs$priors
  )
}

# --- Unit conversions --------------------------------------------------------

# Same conversion as cluster_fns.R::scale_to_halflife(), restated here so this
# file doesn't need cluster_fns.R (which loads malariasimulation).
d_to_halflife     <- function(d)     d * log(2)
log_d_to_halflife <- function(log_d) exp(log_d) * log(2)

# --- Site weights ------------------------------------------------------------

# get_adm1_population_weights(): par_pf (population at risk of P.
# falciparum) per (iso, name_1, year) from each country's site file. Reads
# each iso's site file ONCE (grouped, same efficiency reasoning as
# cluster_fns.R::get_adm1_eir_weighted_from_site()'s per-iso grouped reads).
# Moved here from 06 Sections K/S, which still use it for their per-year
# prevalence weighting.
get_adm1_population_weights <- function(isos, site_cache_dir = SITE_CACHE_DIR) {
  data.table::rbindlist(lapply(isos, function(i) {
    site <- readRDS(file.path(site_cache_dir, paste0(i, "_site.rds")))
    pop  <- data.table::data.table(site$population$population_total)
    pop[, .(iso = iso3c, name_1, year, weight = par_pf)]
  }))
}

# get_site_weights(): one non-negative weight per fitted site, named and
# ordered by post$site_names. `weights` deliberately has NO default:
#   "none"   -> equal weights
#   "par_pf" -> mean par_pf over the years that site has an MIS survey in
#               matched_site_mis_prevalence.RDS (ties each site's weight to
#               the period its own data come from)
#   numeric  -> a named vector covering every site in post$site_names
# Stops (listing the sites) rather than silently dropping any site without a
# weight.
get_site_weights <- function(post, weights,
                             site_cache_dir = SITE_CACHE_DIR,
                             mis_path = file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS")) {
  if (missing(weights)) {
    stop("get_site_weights(): `weights` must be given explicitly -- \"none\", \"par_pf\", or a named numeric vector over site_names.")
  }
  site_names <- post$site_names

  if (is.numeric(weights)) {
    if (is.null(names(weights)) || !all(site_names %in% names(weights))) {
      stop("get_site_weights(): a numeric `weights` vector must be named and cover every site in post$site_names.")
    }
    w <- weights[site_names]
  } else if (identical(weights, "none")) {
    w <- stats::setNames(rep(1, length(site_names)), site_names)
  } else if (identical(weights, "par_pf")) {
    survey_years <- unique(data.table::as.data.table(readRDS(mis_path))[
      !is.na(iso3c), .(site = sprintf("%s_%s", iso3c, name_1), year)])
    survey_years <- survey_years[site %in% site_names]
    pop <- get_adm1_population_weights(unique(post$site_iso), site_cache_dir)
    pop[, site := sprintf("%s_%s", iso, name_1)]
    site_w <- merge(survey_years, pop[, .(site, year, weight)], by = c("site", "year"), all.x = TRUE)
    site_w <- site_w[, .(weight = mean(weight)), by = site]
    w <- stats::setNames(site_w$weight, site_w$site)[site_names]
    names(w) <- site_names
  } else {
    stop("get_site_weights(): `weights` must be \"none\", \"par_pf\", or a named numeric vector.")
  }

  bad <- is.na(w) | w < 0
  if (any(bad)) {
    stop(sprintf("get_site_weights(): missing or negative weight for %d site(s): %s",
                 sum(bad), paste(site_names[bad], collapse = ", ")))
  }
  w
}

# --- Draw matrices and aggregation -------------------------------------------

# site_param_draws(): draws x sites matrix of log_d or log_eir (log scale),
# colnames = site_names.
site_param_draws <- function(post, param = c("log_d", "log_eir")) {
  param <- match.arg(param)
  cols  <- paste0(param, "_", seq_along(post$site_names))
  mat   <- as.matrix(post$pooled_samples[, ..cols])
  colnames(mat) <- post$site_names
  mat
}

# aggregate_site_draws(): per-draw weighted geometric mean of d (param =
# "log_d") or EIR (param = "log_eir") across the sites in each group, on the
# NATURAL scale. Returns a draws x groups matrix (groups = ISO3 codes for
# by = "country", or the single column "all"). Weights are normalised within
# each group, so a country's value only depends on its own sites' relative
# weights.
aggregate_site_draws <- function(post, param = c("log_d", "log_eir"), by = c("country", "all"), weights) {
  param <- match.arg(param)
  by    <- match.arg(by)
  w     <- get_site_weights(post, weights)

  groups   <- if (by == "country") unique(post$site_iso) else "all"
  group_of <- if (by == "country") post$site_iso else rep("all", length(post$site_names))
  W <- matrix(0, nrow = length(post$site_names), ncol = length(groups), dimnames = list(post$site_names, groups))
  for (g in groups) {
    idx <- group_of == g
    if (sum(w[idx]) <= 0) stop(sprintf("aggregate_site_draws(): weights sum to zero in group %s.", g))
    W[idx, g] <- w[idx] / sum(w[idx])
  }
  exp(site_param_draws(post, param) %*% W)
}

# summarise_draws(): mean / median / 95% equal-tailed interval of each column
# of a draws matrix, after an optional per-draw transform (e.g.
# d_to_halflife). Report medians alongside means: medians survive the
# log -> natural back-transform, means don't.
summarise_draws <- function(mat, transform = identity) {
  mat <- transform(mat)
  data.table::data.table(
    group  = colnames(mat),
    mean   = colMeans(mat),
    median = apply(mat, 2, stats::median),
    q025   = apply(mat, 2, stats::quantile, probs = 0.025, names = FALSE),
    q975   = apply(mat, 2, stats::quantile, probs = 0.975, names = FALSE)
  )
}

# --- User-facing summaries ---------------------------------------------------

# summarise_site_parameters(): one table of HRP2 half-life (days) and EIR,
# each with mean / median / 95% CrI, per group:
#   by = "country" -> one row per country (weighted geometric mean of its sites)
#   by = "all"     -> one row across every fitted site
#   by = "site"    -> one row per site (its own marginal posterior; `weights`
#                     isn't used and needn't be given)
summarise_site_parameters <- function(post, by = c("country", "all", "site"), weights) {
  by <- match.arg(by)

  if (by == "site") {
    d_draws   <- exp(site_param_draws(post, "log_d"))
    eir_draws <- exp(site_param_draws(post, "log_eir"))
  } else {
    w         <- get_site_weights(post, weights)   # resolve once, reuse for both params
    d_draws   <- aggregate_site_draws(post, "log_d",   by, w)
    eir_draws <- aggregate_site_draws(post, "log_eir", by, w)
  }

  hl  <- summarise_draws(d_draws, d_to_halflife)
  eir <- summarise_draws(eir_draws)
  data.table::setnames(hl,  c("mean", "median", "q025", "q975"), paste0("halflife_", c("mean", "median", "q025", "q975")))
  data.table::setnames(eir, c("mean", "median", "q025", "q975"), paste0("eir_",      c("mean", "median", "q025", "q975")))
  out <- merge(hl, eir, by = "group", sort = FALSE)

  if (by == "site") {
    out[, `:=`(iso = post$site_iso[match(group, post$site_names)],
               name_1 = post$site_name_1[match(group, post$site_names)])]
    data.table::setnames(out, "group", "site")
    data.table::setcolorder(out, c("site", "iso", "name_1"))
  } else if (by == "country") {
    out[, `:=`(country = unname(COUNTRY_NAME_LOOKUP[group]),
               n_sites = as.integer(table(post$site_iso)[group]))]
    data.table::setnames(out, "group", "iso")
    data.table::setcolorder(out, c("iso", "country", "n_sites"))
  } else {
    out[, n_sites := length(post$site_names)]
    data.table::setcolorder(out, c("group", "n_sites"))
  }
  out[]
}


# get_population_parameters(): posterior summaries of the population-level
# decay hyperparameters, on their own scale and as half-lives:
#   mu_d                  log(d) population centre
#   sigma_d               SD of log(d) across sites
#   halflife_pop_median   exp(mu_d) * log(2) -- the population MEDIAN site
#                         half-life (days), not a mean
#   halflife_spread_factor exp(1.96 * sigma_d) -- 95% of sites' half-lives
#                         lie within [median / factor, median * factor]
get_population_parameters <- function(post) {
  ps  <- post$pooled_samples
  mat <- cbind(
    mu_d                   = ps$mu_d,
    sigma_d                = ps$sigma_d,
    halflife_pop_median    = log_d_to_halflife(ps$mu_d),
    halflife_spread_factor = exp(1.96 * ps$sigma_d)
  )
  out <- summarise_draws(mat)
  data.table::setnames(out, "group", "parameter")
  out[]
}

# predict_unobserved_site(): posterior predictive HRP2 half-life for a new,
# unobserved site -- the expected value for a location with no data, from
# the population-level parameters alone. For each posterior draw, draws
# log_d_new ~ Normal(mu_d, sigma_d) (or the Student-t form when the fit used
# decay_family = "t"), truncated to the same log_d_bounds the site-level
# prior uses (rejection: out-of-bounds draws are redrawn from that same
# posterior draw's mu_d/sigma_d). This propagates the uncertainty in mu_d and
# sigma_d, unlike a plug-in mu_d +/- 1.96*sigma_d.
#
# Decay only: EIR isn't pooled, so there's no posterior information about an
# unobserved site's EIR beyond its prior.
#
# Returns list(draws = data.table(log_d, d, halflife), summary = 1-row table).
predict_unobserved_site <- function(post, n_per_draw = 1L, seed = 1L, max_tries = 1000L) {
  if (!is.null(seed)) set.seed(seed)
  mu  <- rep(post$pooled_samples$mu_d,    each = n_per_draw)
  sig <- rep(post$pooled_samples$sigma_d, each = n_per_draw)
  is_t <- identical(post$priors$decay_family, "t")
  draw <- function(idx) {
    if (is_t) mu[idx] + sig[idx] * stats::rt(length(idx), df = post$priors$decay_t_df)
    else      stats::rnorm(length(idx), mu[idx], sig[idx])
  }

  log_d <- numeric(length(mu))
  todo  <- seq_along(mu)
  for (i in seq_len(max_tries)) {
    log_d[todo] <- draw(todo)
    todo <- todo[log_d[todo] < post$log_d_bounds[1] | log_d[todo] > post$log_d_bounds[2]]
    if (length(todo) == 0L) break
  }
  if (length(todo) > 0L) {
    stop(sprintf("predict_unobserved_site(): %d draw(s) still outside log_d_bounds after %d tries.", length(todo), max_tries))
  }

  draws <- data.table::data.table(log_d = log_d, d = exp(log_d), halflife = log_d_to_halflife(log_d))
  summary <- data.table::data.table(
    quantity = "halflife_unobserved_site",
    mean     = mean(draws$halflife),
    median   = stats::median(draws$halflife),
    q025     = stats::quantile(draws$halflife, 0.025, names = FALSE),
    q975     = stats::quantile(draws$halflife, 0.975, names = FALSE)
  )
  list(draws = draws, summary = summary)
}
