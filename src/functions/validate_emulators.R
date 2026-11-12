# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: emulator validation)
#
# validate_emulators.R
# Held-out cross-validation of site-specific GP emulators.
# Assess whether emulators accurately reproduce log-likelihoods in the
# region of parameter space that matters.
#
# Two checks are implemented:
#   1. Leave-one-out CV on the training grid (fast, uses DiceKriging leaveOneOut)
#   2. Held-out validation on a random subset of grid points (recommended)
#
# Emulators are 2D: (log_d, log_eir) -- d is the single HRP2 clearance scale
# malariasimulation 3.0.0's set_hrp2_parameters() now takes (the old
# separate treated/untreated decay rates were collapsed into one shared
# rate; see 15_launch_grid_persistence.R Section B). deletion is NOT a GP
# dimension -- see 16_fit_emulators.R's header for why -- it's fixed
# per-country outside the emulator entirely.
#
# Call validate_emulators() after fit_site_emulators().

library(DiceKriging)
library(data.table)

# Must match 16_fit_emulators.R's is_shared_param_grids()/resolve_site_grid()
# exactly — param_grids here can be either one grid shared by every site, or
# a named list of per-site grids (see 16_fit_emulators.R's
# build_emulator_inputs_from_grid_search(), which produces the per-site form
# for 15_launch_grid_persistence.R's actual output).
is_shared_param_grids <- function(param_grids) {
  all(c("d", "eir") %in% names(param_grids))
}
resolve_site_grid <- function(param_grids, site) {
  if (is_shared_param_grids(param_grids)) param_grids else param_grids[[site]]
}


# --- Held-out validation ---

validate_emulators <- function(
  ll_arrays,
  param_grids,
  emulators,
  holdout_frac      = 0.2,
  seed              = 123,
  max_design_points = 2000L,
  cap_floor         = TRUE,
  floor_below_max   = 100,
  covtype           = "matern5_2",
  multistart        = 5
) {
  # Refits each emulator on a training subset and evaluates on held-out points.
  # floor_below_max/covtype/multistart: must match the production fit
  # (emulator_fns.R::fit_single_emulator_design(), called from
  # 3b_calibration_emulator/emulation.R with multistart = 5). Until
  # 2026-10-12 this refit still used the old 1.5*IQR floor and a single
  # optimiser start, so it validated a different emulator from the one the
  # MCMC used.
  # ll_arrays: named list of 2D arrays [n_d x n_eir].
  # param_grids: shared or per-site (see resolve_site_grid() above).
  # max_design_points: see validate_single_site() below — must mirror
  # 16_fit_emulators.R::fit_single_emulator()'s cap, or a real 15-derived
  # per-site grid (up to ~10,000 finite points) reintroduces exactly the
  # O(n^3) km() refit cost that cap exists to avoid, just here instead.
  # cap_floor: must mirror whatever fit_site_emulators()/fit_single_emulator*()
  # was actually called with (see 16_fit_emulators.R::fit_single_emulator_design()'s
  # comment) — validating against a differently-capped response than what
  # production trained on would make the reported error meaningless.
  # Returns list: $summary (data.table) and $detail (data.table, all sites combined).

  set.seed(seed)
  site_names <- names(emulators)

  result_list <- lapply(site_names, function(site) {
    message("Validating emulator: ", site)
    if (is.data.frame(ll_arrays[[site]])) {
      validate_single_site_design(
        ll_arrays[[site]], holdout_frac,
        max_design_points = max_design_points, cap_floor = cap_floor,
        floor_below_max = floor_below_max, covtype = covtype, multistart = multistart
      )
    } else {
      validate_single_site(
        ll_arrays[[site]], resolve_site_grid(param_grids, site), holdout_frac,
        max_design_points = max_design_points, cap_floor = cap_floor,
        floor_below_max = floor_below_max, covtype = covtype, multistart = multistart
      )
    }
  })

  detail <- rbindlist(result_list, idcol = "site")
  detail[, site := site_names[site]]

  # Summary per site: RMSE and mean absolute error
  summary_dt <- detail[, .(
    n_holdout = .N,
    rmse      = sqrt(mean(error^2)),
    mae       = mean(abs(error)),
    max_ae    = max(abs(error)),
    r2        = 1 - sum(error^2) / sum((ll_true - mean(ll_true))^2)
  ), by = site]

  message("\nEmulator validation summary:")
  print(summary_dt)

  list(summary = summary_dt, detail = detail)
}


validate_single_site <- function(ll_arr, param_grids, holdout_frac, max_design_points = 2000L,
                                  cap_floor = TRUE, floor_below_max = 100,
                                  covtype = "matern5_2", multistart = 5) {
  # Build 2D FACTORIAL design matrix matching fit_emulators.R conventions,
  # then delegate to validate_single_site_design() below. Only valid for a
  # rectangular grid — irregular designs (e.g. 15's GRID_METHOD = "lhs") go
  # through validate_single_site_design() directly (see validate_emulators()
  # above, which dispatches on ll_arrays[[site]]'s shape).
  design <- as.data.table(expand.grid(
    d   = param_grids$d,
    eir = param_grids$eir
  ))
  design[, ll := as.vector(ll_arr)]
  validate_single_site_design(design, holdout_frac, max_design_points = max_design_points,
                               cap_floor = cap_floor, floor_below_max = floor_below_max,
                               covtype = covtype, multistart = multistart)
}


validate_single_site_design <- function(design, holdout_frac, max_design_points = 2000L,
                                         cap_floor = TRUE, floor_below_max = 100,
                                         covtype = "matern5_2", multistart = 5) {
  # Core validation: design is a tidy data.table/data.frame with (at least)
  # columns d, eir, ll — see
  # 16_fit_emulators.R::fit_single_emulator_design()'s equivalent split.
  design <- as.data.table(design)
  design[, `:=`(
    log_d   = log(d),
    log_eir = log(eir)
  )]
  # Require the design columns themselves to be finite, not just ll — see
  # 16_fit_emulators.R for why a finite-ll/non-finite-predictor row still
  # breaks km()'s internal covariance computation.
  design <- design[is.finite(ll) & is.finite(log_d) & is.finite(log_eir)]

  # Mirrors 16_fit_emulators.R::fit_single_emulator()'s subsampling exactly —
  # this validates "how well would an emulator trained on ~max_design_points
  # points generalise", which is the actually-relevant question given
  # fit_site_emulators() trains on that same capped size, not the full grid.
  if (!is.null(max_design_points) && nrow(design) > max_design_points) {
    design <- design[sample(nrow(design), max_design_points)]
  }

  # Mirrors fit_single_emulator_design()'s floor exactly (floor_below_max
  # log-lik units below the max; NULL = the old 1.5*IQR fence) -- applied to the WHOLE design
  # before the train/holdout split, same as production applies it to the
  # whole (subsampled) design before fitting. This matters for TWO reasons:
  # y_train must match production's response shape for validation to mean
  # anything, and holdout points capped the same way keeps the reported
  # error focused on the posterior-relevant region validate_emulators() is
  # actually meant to check (per this file's header) rather than being
  # dominated by deep-tail extrapolation error that's irrelevant by design.
  # cap_floor = FALSE skips this (see fit_single_emulator_design()'s own
  # comment on that experiment) -- must match whatever production actually
  # trained with, or this validates against the wrong response entirely.
  if (cap_floor) {
    floor <- if (is.null(floor_below_max)) {
      max(design$ll) - 1.5 * unname(diff(quantile(design$ll, c(0.25, 0.75))))
    } else {
      max(design$ll) - floor_below_max
    }
    design[, ll := pmax(ll, floor)]
  }

  n       <- nrow(design)
  holdout <- sample(seq_len(n), size = floor(n * holdout_frac))
  train   <- setdiff(seq_len(n), holdout)

  X_train <- as.matrix(design[train,   .(log_d, log_eir)])
  y_train <- design$ll[train]
  X_test  <- as.matrix(design[holdout, .(log_d, log_eir)])
  y_test  <- design$ll[holdout]

  # Same km() call as fit_single_emulator_design()
  em_train <- km(
    formula      = ~1,
    design       = data.frame(X_train),
    response     = y_train,
    covtype      = covtype,
    nugget       = 1e-6,
    nugget.estim = TRUE,
    control      = list(trace = FALSE),
    multistart   = multistart
  )

  pred <- predict(
    em_train,
    newdata    = data.frame(X_test),
    type       = "SK",
    checkNames = FALSE
  )

  data.table(
    ll_true   = y_test,
    ll_pred   = pred$mean,
    ll_sd     = pred$sd,
    error     = pred$mean - y_test,
    # Standardised residual: large values suggest emulator uncertainty is underestimated
    std_error = (pred$mean - y_test) / (pred$sd + 1e-10)
  )
}


# --- Leave-one-out CV using DiceKriging built-in ---

loo_cv_emulators <- function(emulators, near_peak_ll = 3) {
  # Fast LOO-CV using DiceKriging's leaveOneOut.predict().
  # Less stringent than held-out validation but requires no refit.
  # Returns per-site RMSE and coverage of 95% prediction intervals over the
  # whole training grid, plus rmse_near / n_near over the training points
  # within near_peak_ll log-likelihood units of the site's best point. A drop
  # of 3 units is ~ the 95% joint region for two parameters (chi-squared(2)
  # 95% critical value 5.99 / 2), i.e. where the posterior sits; the
  # emulator flag (identify_poor_loo_sites()) uses rmse_near.

  site_names <- names(emulators)

  rbindlist(lapply(site_names, function(site) {
    em   <- emulators[[site]]
    loo  <- leaveOneOut.km(em, type = "SK", trend.reestim = FALSE)

    y_true <- em@y
    y_pred <- loo$mean
    y_sd   <- loo$sd

    within_95 <- abs(y_true - y_pred) < 1.96 * y_sd
    near      <- y_true >= max(y_true) - near_peak_ll

    data.table(
      site       = site,
      rmse       = sqrt(mean((y_pred - y_true)^2)),
      mae        = mean(abs(y_pred - y_true)),
      coverage95 = mean(within_95),
      n_near     = sum(near),
      rmse_near  = sqrt(mean((y_pred[near] - y_true[near])^2))
    )
  }))
}
