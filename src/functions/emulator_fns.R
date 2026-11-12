# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: emulators)
#
# 16_fit_emulators_fns.R
# Pure function library for 16_fit_emulators.R's GP emulator fitting --
# fit site-specific Gaussian-process emulators to binomial log-likelihood
# surfaces from grid search over (d, EIR) -- d is the single HRP2 clearance
# scale malariasimulation 3.0.0's set_hrp2_parameters() now takes (the old
# separate treated/untreated decay rates, d_treated/d_untreated, were
# collapsed into one shared rate; see 15_launch_grid_persistence.R Section B).
#
# Each site's 2D likelihood array is reshaped to a design matrix, d/eir are
# log-transformed for a better-behaved GP input space, and a km() model is
# fitted. Emulators are saved to disk for downstream use.
#
# hrp2_deletion_pct is NOT a GP dimension: it's fixed per-country at
# data/hrp2_deletion.csv's central value in 15_launch_grid_persistence.R
# Section F, applied as a deterministic post-hoc scaling of simulated
# hrp2_prev rather than optimised or emulated. It WAS previously a 4th grid/
# GP dimension, profiled per-country and logit-transformed — that consistently
# produced degenerate covariance hyperparameters (every theta pinned at
# DiceKriging's default bound: 1e-10 lower or 2x-column-range upper, not a
# genuine interior optimum, confirmed identical to several decimal places
# across different sites/countries) regardless of optimizer settings tried
# (multistart, explicit wide bounds, biased subsampling) — see git history
# for the investigation. Dropping it as a free/fitted parameter entirely
# removes that failure mode, and matches the fact that pfhrp2/3 deletion
# prevalence is measured/reported at country level in the source data
# anyway (no admin1 resolution to profile against).
#
# Input expectations:
#   ll_arrays  - named list of 2D arrays, dim = [n_d x n_eir],
#                values are binomial log-likelihoods from the grid search.
#                Names must match site_names.
#   param_grids - EITHER a single list of two numeric vectors on natural
#                 scale, shared by every site:
#                   $d    (length = dim(array)[1])
#                   $eir  (length = dim(array)[2])
#                 OR a named list of such lists, keyed by site name, for
#                 per-site grids (see build_emulator_inputs_from_grid_search()
#                 below). eir is ABSOLUTE — 04_grid_simulation.R's EIR_GRID,
#                 the same shared range simulated identically for every adm1
#                 combo, not a per-adm1-baseline-relative multiplier.
#                 build_emulator_inputs_from_grid_search() still derives each
#                 site's grid from its own data (per-site form) rather than
#                 assuming this shared shape, but the shared-grid form is now
#                 the real, common case, not just illustrative/test data.
#   d/eir must be strictly positive (log-transform is applied).
#
# IMPORTANT: this file must contain ONLY function definitions (plus the
# library()/source() calls they need) -- no top-level submit/wait/collect
# side effects. It's sourced by hipercow_environment_create() (both the
# "emulator_fit" and "hierarchical_mcmc" environments), which means every
# worker task loading either environment re-executes this file's top level.
# The driver script (16_fit_emulators.R) sources this file and keeps all the
# if(TRUE)/if(FALSE)-gated submit/wait/collect sections for interactive use.

library(DiceKriging)
library(data.table)
library(parallel)
library(hipercow)

# validate_emulators()/loo_cv_emulators() below are used to sanity-check
# fits before they're handed to 17_run_hierarchical_mcmc.R (see the guarded
# blocks in 16_fit_emulators.R). Sourced here rather than left to callers
# (example_usage.R also sources it directly — redundant but harmless, since
# re-sourcing just redefines the same functions).
source("./src/functions/validate_emulators.R")

# plot_emulator_diagnostics()/save_emulator_diagnostics_png() render a raw-
# grid-vs-smooth-emulator heatmap comparison per site, for visual review
# alongside the numeric loo/validate checks above (see plot_emulators.R —
# a "clean smooth bump" vs. a striped pattern can reveal degenerate GP
# hyperparameters that grid-locked numeric checks can't detect at all).
source("./src/functions/plot_emulators.R")

# A "shared" param_grids has the three grid names directly at the top level
# (one grid for every site); a "per-site" param_grids is instead keyed by
# site name, each value being one of those three-grid lists. Detecting this
# lets fit_site_emulators()/validate_emulators()/flag_poor_emulator_regions()
# accept either shape without a separate code path per caller.
is_shared_param_grids <- function(param_grids) {
  all(c("d", "eir") %in% names(param_grids))
}
resolve_site_grid <- function(param_grids, site) {
  if (is_shared_param_grids(param_grids)) param_grids else param_grids[[site]]
}


# --- Core emulator fitting ---


fit_single_emulator <- function(ll_arr, param_grids, site_label = "",
                                 max_design_points = 2000L, seed = 123,
                                 cap_floor = TRUE, covtype = "matern5_2",
                                 multistart = 1, floor_below_max = 100) {
  # Reshape 2D FACTORIAL array to a tidy design, then delegate to
  # fit_single_emulator_design() below for the actual fit. Only valid when
  # ll_arr really is a dense rectangular array over param_grids' axes — an
  # irregular/space-filling design (e.g. 15_launch_grid_persistence.R's
  # GRID_METHOD = "lhs") has no such array to reshape (each row has its own
  # (d, eir) pair, not shared per-axis levels), so those go through
  # fit_single_emulator_design() directly instead; see
  # build_emulator_inputs_from_grid_search()'s rectangular-grid detection,
  # which routes each site to the right one automatically.
  # cap_floor: passed straight through — see fit_single_emulator_design().
  design <- as.data.table(expand.grid(
    d   = param_grids$d,
    eir = param_grids$eir
  ))
  design[, ll := as.vector(ll_arr)]

  fit_single_emulator_design(design, site_label, max_design_points, seed, cap_floor,
                              covtype = covtype, multistart = multistart,
                              floor_below_max = floor_below_max)
}


fit_single_emulator_design <- function(design, site_label = "",
                                        max_design_points = 2000L, seed = 123,
                                        cap_floor = TRUE, covtype = "matern5_2",
                                        multistart = 1, floor_below_max = 100) {
  # Core fit: design is a tidy data.table/data.frame with (at least) columns
  # d, eir, ll — one row per simulated (or grid) point. No array/expand.grid
  # involved, so this works identically for a factorial grid's flattened
  # design (via fit_single_emulator() above) AND an irregular design where
  # every row has its own unique combination.
  # max_design_points: km()'s covariance factorisation is O(n^3), so a real
  # 15_launch_grid_persistence.R-derived grid (thousands of finite points
  # per site) makes every extra training point expensive. A full factorial
  # grid is also heavily redundant for a smooth GP surface — neighbouring
  # grid points carry near-identical information — so a random subsample
  # below this cap keeps fit time bounded without materially hurting
  # emulator accuracy on a smooth binomial log-lik surface. NULL disables
  # subsampling (keep every finite point).
  # cap_floor: TRUE (default) applies the adaptive IQR floor below. FALSE
  # skips it entirely and trains on the raw log-lik values as-is — an
  # experiment worth comparing against the default, NOT a known-better
  # setting: the floor was added specifically because a deep, uninformative
  # tail dominating the training signal produced a degenerate fit (see the
  # floor's own comment below), and turning it off reintroduces exactly that
  # risk. Compare fit quality (loo_cv_emulators()/validate_emulators()) with
  # and without before trusting an uncapped fit.

  design <- as.data.table(design)

  # Log-transform d/eir: GP works on additive scale; log makes the input
  # space more uniform and prevents the kernel from being dominated by EIR
  # range.
  design[, `:=`(
    log_d   = log(d),
    log_eir = log(eir)
  )]

  # Drop grid points outside the likelihood support; -Inf arises where the
  # simulator produces prevalence of exactly 0 or 1 against non-zero data.
  # Also require the DESIGN columns themselves to be finite, not just ll —
  # a row can have a perfectly finite ll while one of its predictor columns
  # is NA/Inf (e.g. an NA grid value), and km() has no tolerance for that:
  # DiceKriging's internal covariance computation throws a low-level
  # "NA/NaN/Inf in foreign function call" error rather than a friendly
  # R-level message.
  n_total <- nrow(design)
  design   <- design[is.finite(ll) & is.finite(log_d) & is.finite(log_eir)]
  n_finite <- nrow(design)

  if (n_finite < 10) {
    stop(sprintf(
      "Site '%s': only %d finite likelihood points (of %d). Check grid coverage.",
      site_label, n_finite, n_total
    ))
  }

  if (n_finite < n_total) {
    message(sprintf(
      "  Site '%s': dropped %d non-finite grid points.",
      site_label, n_total - n_finite
    ))
  }

  if (!is.null(max_design_points) && n_finite > max_design_points) {
    set.seed(seed)
    design <- design[sample(n_finite, max_design_points)]
    message(sprintf(
      "  Site '%s': subsampled %d of %d finite grid points for GP fitting (max_design_points = %d).",
      site_label, max_design_points, n_finite, max_design_points
    ))
  }

  X <- as.matrix(design[, .(log_d, log_eir)])
  y <- design$ll

  # Cap log-likelihood at a per-site, scale-adaptive floor (Tukey's classic
  # 1.5*IQR fence, anchored at the max rather than Q3 since we care about
  # distance below the best point, not the median) rather than a fixed
  # constant. A flat cap (e.g. -50) doesn't work here: real per-site log-lik
  # ranges span anywhere from ~1e2 to ~2e4 units (verified against real
  # 15-derived data), so a small fixed cap clamps 80-99%+ of every site's
  # training points to one identical floor value -- destroying almost all
  # the GP's training signal outside a tiny sliver near the optimum (this is
  # what produced near-zero LOO-CV error but genuinely terrible held-out
  # validate_emulators() error: the "emulator" had collapsed to
  # near-constant almost everywhere). k=1.5 empirically clamps ~15-35% of
  # points across sites with wildly different scales (768 to 22,000-unit
  # ranges) -- discarding the genuinely posterior-irrelevant deep tail
  # without destroying the informative region near the optimum.
  #
  # 2026-10-09: replaced by a fixed distance below each site's own best,
  # floor_below_max = 100 log-lik units (points further down are <= e^-100 as
  # likely, so the MCMC never visits them). The 1.5*IQR fence sat hundreds of
  # units down at steep sites, so one stationary GP had to span the whole
  # drop: its fitted nugget reached SD ~50 log-lik units, larger than the
  # structure near the peak, and it invented spurious peaks there (GIN_Boké:
  # posterior chains stuck in separate modes, Rhat 1.53; GIN_Kindia: positive
  # predictions under matern5_2). Across all 248 sites (RDT sensitivity 0.9),
  # floor 100 with matern5_2 left near-peak LOO error unchanged at 240 sites,
  # improved it by > 1 unit at 8 and worsened none; overshoot > 5 fell from 5
  # sites to 1 (outputs/05_emulation/floor100_vs_current_all_sites.csv).
  # Whole-surface LOO/validation error is not comparable with earlier fits:
  # the deep tail now sits on the floor by design. Set floor_below_max = NULL
  # to restore the 1.5*IQR fence.
  if (cap_floor) {
    floor <- if (is.null(floor_below_max)) {
      max(y) - 1.5 * unname(diff(quantile(y, c(0.25, 0.75))))
    } else {
      max(y) - floor_below_max
    }
    y <- pmax(y, floor)
  }

  # Matern 5/2 is the standard, DEFAULT choice for smooth physical-model
  # outputs; nugget estimated to absorb numerical noise at the boundary of
  # the likelihood surface. covtype/multistart are per-call overrides (see
  # SITE_KERNEL_OVERRIDES in 05_emulation.R Section C) for the small number
  # of sites where matern5_2 has a confirmed extrapolation-spike artifact
  # (GIN_Kindia/GIN_Boké investigation: the smooth GP surface fabricates an
  # impossible POSITIVE log-likelihood between real grid points -- neither
  # theta bounds, more training data, nor multistart alone fixed it under
  # matern5_2, but the "exp" kernel does, at some LOO-CV cost). multistart
  # > 1 also matters on its own: a single BFGS start can converge to a
  # degenerate near-zero theta(log_eir) (confirmed for real, same
  # investigation) that happens to fit the training points fine but is
  # wildly unstable between them -- multistart=5+ reliably avoids that
  # local optimum. Both default to matern5_2/1, reproducing every existing
  # caller's behaviour exactly when neither is overridden.
  km(
    formula      = ~1,
    design       = data.frame(X),
    response     = y,
    covtype      = covtype,
    nugget       = 1e-6,
    nugget.estim = TRUE,
    control      = list(trace = FALSE),
    multistart   = multistart
  )
}


# --- Emulator prediction ---

predict_emulator_ll <- function(emulator, log_d, log_eir) {
  # Return the GP mean prediction of log-likelihood at a single point.
  # Uses simple kriging (SK) with the estimated nugget; prediction variance
  # is not returned here but can be retrieved for diagnostics if needed.
  newdata <- data.frame(
    log_d   = log_d,
    log_eir = log_eir
  )
  predict(emulator, newdata = newdata, type = "SK", checkNames = FALSE)$mean
}


# --- Convert 15_launch_grid_persistence.R's grid-search output into
# fit_site_emulators()'s (ll_arrays, param_grids) inputs ---

build_emulator_inputs_from_grid_search <- function(name1_grid_log_lik) {
  # name1_grid_log_lik: data.table from 04_grid_simulation.R Section F
  #   (readRDS("outputs/04_grid_simulation/hrp2_grid_log_lik.RDS")$by_name1_grid_point),
  #   one row per (iso, name_1, eir, hrp2_scale) grid point, with a summed
  #   log_lik. malariasimulation 3.0.0 collapsed the old separate treated/
  #   untreated HRP2 clearance rates into one shared hrp2_scale (see
  #   04_grid_simulation.R Section B), so there's only one HRP2 column here
  #   now, not two. hrp2_deletion_pct is NOT a column here either — it's
  #   fixed per-country at data/hrp2_deletion.csv's central value in 04
  #   Section F, applied as a deterministic post-hoc scaling of simulated
  #   hrp2_prev rather than gridded/emulated (see this file's header for why:
  #   it consistently produced degenerate covariance hyperparameters as an
  #   extra GP dimension).
  #
  # eir here is ABSOLUTE (04_grid_simulation.R's EIR_GRID, shared identically
  # across every adm1 combo — see that file's Section B header), trained on
  # directly as log(eir). This function used to also take an eir_in_by_site
  # argument to convert 04's absolute grid back to a per-site RELATIVE
  # multiplier (eir/eir_in) for training — removed when the pipeline moved
  # from relative to absolute EIR (relative EIR meant the emulator's own
  # axis, and therefore the MCMC's EIR posterior, was defined against an
  # external per-adm1 baseline that turned out to be unreliable for several
  # high-transmission sites; see hierarchical_mcmc_fns.R's header for the
  # full account).
  #
  # Grid axes (d/eir) are derived PER SITE from the data itself, not from
  # 04's original grid-construction globals — this makes no assumption that
  # they stay shared across sites (even though, for now, every site's grid
  # IS in fact the same EIR_GRID — see this file's earlier per-site-grid
  # comment for why per-site derivation is still the right default). A
  # missing (site, grid-point) combo — e.g. a cluster task that failed and
  # was dropped in 04 Section E — leaves that array cell NA, which
  # fit_single_emulator() already drops via its is.finite(ll) filter.
  #
  # PER SITE, this ALSO detects whether (d, eir) form a RECTANGULAR
  # (factorial) grid or an IRREGULAR one (04's GRID_METHOD = "lhs": each row
  # has its own unique pair, not shared per-axis levels). Building the usual
  # dense 2D array for an irregular design would try to allocate one cell
  # per (#unique d) x (#unique eir) combination — for a 500-point LHS design
  # that's up to ~500^2 = 250K cells PER SITE even though only ~500 rows of
  # real data exist, and per-site that adds up across a real
  # multi-hundred-site run. Irregular sites instead get a TIDY design table
  # (see fit_single_emulator_design()) — fit_site_emulators() and
  # validate_emulators() both dispatch on this shape automatically.
  #
  # Returns list(ll_arrays, param_grids, site_country) — the first two feed
  # fit_site_emulators()/validate_emulators() directly; site_country is a
  # site -> ISO map for any downstream per-country grouping. Site names are
  # "<iso>_<name_1>" since name_1 alone isn't unique across countries.

  dt <- data.table::as.data.table(name1_grid_log_lik)

  site_keys  <- unique(dt[, .(iso, name_1)])
  site_names <- sprintf("%s_%s", site_keys$iso, site_keys$name_1)

  ll_arrays   <- vector("list", length(site_names))
  param_grids <- vector("list", length(site_names))
  names(ll_arrays)   <- site_names
  names(param_grids) <- site_names

  for (k in seq_along(site_names)) {
    site_dt <- dt[iso == site_keys$iso[k] & name_1 == site_keys$name_1[k]]

    d_grid   <- sort(unique(site_dt$hrp2_scale))
    eir_grid <- sort(unique(site_dt$eir))

    # Tolerant of a handful of missing combos (failed/dropped cluster
    # tasks, 04 Section E) without misclassifying a near-complete factorial
    # grid as irregular — an LHS design's unique-pair count sits orders of
    # magnitude below the full axis product (~500 vs ~500^2 here), so a
    # generous 50% threshold cleanly separates the two regimes either way.
    n_pairs        <- nrow(unique(site_dt[, .(hrp2_scale, eir)]))
    is_rectangular <- n_pairs >= 0.5 * length(d_grid) * length(eir_grid)

    if (is_rectangular) {
      arr <- array(NA_real_, dim = c(length(d_grid), length(eir_grid)))
      arr[cbind(
        match(site_dt$hrp2_scale, d_grid),
        match(site_dt$eir,        eir_grid)
      )] <- site_dt$log_lik

      ll_arrays[[k]]   <- arr
      param_grids[[k]] <- list(
        d   = d_grid,
        eir = eir_grid
      )
    } else {
      ll_arrays[[k]] <- site_dt[, .(
        d   = hrp2_scale,
        eir = eir,
        ll  = log_lik
      )]
      # Ranges only (not a full grid) — enough for fit_site_emulators()'s
      # positivity/bounds checks and for anything that just needs bounds.
      param_grids[[k]] <- list(
        d   = range(site_dt$hrp2_scale),
        eir = range(site_dt$eir)
      )
    }
  }

  list(
    ll_arrays    = ll_arrays,
    param_grids  = param_grids,
    site_country = setNames(site_keys$iso, site_names)
  )
}


# --- Fallback: validate only sites LOO-CV flags as poor performers ---

identify_poor_loo_sites <- function(loo, near_peak_rmse_threshold = 2) {
  # loo: output of loo_cv_emulators() (needs rmse_near).
  # Flags sites whose LOO RMSE near the peak (within loo_cv_emulators()'s
  # near_peak_ll = 3 log-likelihood units of the best grid point) exceeds
  # near_peak_rmse_threshold, the same 2-unit convention as
  # flag_poor_emulator_regions()'s ll_error_threshold. Error near the peak is
  # what moves the posterior; whole-surface RMSE and coverage are dominated
  # by the steep tail and the 100-unit floor, which the MCMC never visits.
  # 2026-10-13: replaced the earlier rule (whole-surface LOO RMSE above
  # Q3 + 1.5 x IQR, or 95% coverage < 0.8; 33 sites at RDT sensitivity 0.9,
  # mostly in Uganda and Malawi with near-peak RMSE < 1).
  if (is.null(loo$rmse_near)) {
    stop("identify_poor_loo_sites(): loo has no rmse_near column -- rerun loo_cv_emulators().")
  }
  loo$site[!is.na(loo$rmse_near) & loo$rmse_near > near_peak_rmse_threshold]
}

validate_poor_loo_sites <- function(
  ll_arrays, param_grids, emulators,
  near_peak_rmse_threshold = 2,
  holdout_frac          = 0.2,
  seed                  = 123,
  max_design_points     = 2000L,
  loo                   = NULL
) {
  # Cheap-first fallback: loo_cv_emulators() is free (reuses each emulator's
  # own already-fitted training design, no refit) and always runs first.
  # validate_emulators() DOES refit per site — exactly the O(n^3) cost
  # max_design_points exists to bound in fit_single_emulator() — so it only
  # runs on sites LOO-CV actually flags as poor, not the full site list.
  # Returns list(loo, poor_sites, val) — val is NULL if nothing was flagged.
  #
  # ll_arrays/param_grids are passed through to validate_emulators()
  # UNSLICED: param_grids may be shared (one grid for every site) rather
  # than per-site (see resolve_site_grid() near the top of this file), so
  # slicing it here by poor_sites would break the shared case — only
  # `emulators` needs restricting, since validate_emulators() drives its
  # site loop off names(emulators).
  if (is.null(loo)) loo <- loo_cv_emulators(emulators)

  poor_sites <- identify_poor_loo_sites(loo, near_peak_rmse_threshold)

  if (length(poor_sites) == 0L) {
    message(sprintf(
      "loo_cv_emulators(): no poor performers (near-peak rmse > %.3g) among %d site(s); skipping validate_emulators().",
      near_peak_rmse_threshold, nrow(loo)
    ))
    return(list(loo = loo, poor_sites = character(0), val = NULL))
  }

  message(sprintf(
    "loo_cv_emulators(): %d/%d site(s) flagged as poor performers: %s — running validate_emulators() on them.",
    length(poor_sites), nrow(loo), paste(poor_sites, collapse = ", ")
  ))

  val <- validate_emulators(
    ll_arrays, param_grids, emulators[poor_sites],
    holdout_frac = holdout_frac, seed = seed, max_design_points = max_design_points
  )

  list(loo = loo, poor_sites = poor_sites, val = val)
}


# --- Save / load helpers ---

save_emulators <- function(emulators, path) {
  saveRDS(emulators, path)
  message("Saved ", length(emulators), " emulators to: ", path)
}

load_emulators <- function(path) {
  readRDS(path)
}
