# =============================================================================
# src/3b_calibration_emulator/emulation.R              (was src/05_emulation.R)
#
# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (stage 2: Gaussian process emulators)
#
# Step 3b, calibration stage 2: Gaussian-process emulators. For each site, a
# GP (DiceKriging::km, Matern 5/2 kernel, multistart = 5) is fitted to the
# grid log-likelihood surface from step 3a over (log d, log EIR), so the MCMC
# (step 3c) can evaluate the likelihood anywhere in the grid's range without
# re-running malariasimulation. One hipercow task per site.
#
# Sections
#   A  hipercow environment (once per hipercow root)
#   C  submit one emulator fit per site
#   D  collect -> emulators.rds; LOO-CV / hold-out validation -> emulator_validation.rds
# Sections C-D run once per value in RDT_SENSITIVITY_RUN (main analysis 0.9;
# sensitivity/rdt_sensitivity_emulation.R runs 0.8 and 1.0). Each reads its
# own grid log-likelihood (cluster_fns.R::log_lik_file()) and writes to its
# own folder (DIR_EMULATOR for 0.9, DIR_EMULATOR/sensitivity_rdt_se_NNN/
# otherwise).
#
# Emulator fitting details (emulator_fns.R::fit_single_emulator*()): up to
# 2,000 design points per site; points more than 100 log-likelihood units
# below the site's best are floored (CAP_FLOOR = TRUE), so the GP spends its
# flexibility near the peak rather than on the steep, irrelevant tails.
#
# The emulator vetting figure (raw grid vs emulator vs 1:1 LOO scatter, per
# RDT sensitivity) is drawn in src/plots.R. The SIR posterior-predictive
# check (sample each emulator surface and re-run malariasimulation), the
# per-site diagnostic PNGs and the 2-10y EIR-baseline runs are in
# archive/src/05_emulation.R.
# =============================================================================
library(data.table)
library(hipercow)

source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/emulator_fns.R")
source("./src/functions/pipeline_fns.R")

# Settings ======================================================================
OVR <- take_overrides()
# Sensitivity analyses (defaults = main analysis): GRID_OUT_DIR is where the
# grid log-likelihood is read; EMULATOR_BASE_DIR where emulators are written.
# REFIT_ISOS (NULL: every site) refits only those countries' sites; the other
# sites' emulators are copied from BASE_EMULATOR_DIR's emulators.rds.
GRID_OUT_DIR      <- setting(OVR, "GRID_OUT_DIR", DIR_GRID)
EMULATOR_BASE_DIR <- setting(OVR, "EMULATOR_BASE_DIR", DIR_EMULATOR)
REFIT_ISOS        <- setting(OVR, "REFIT_ISOS", NULL)
BASE_EMULATOR_DIR <- setting(OVR, "BASE_EMULATOR_DIR", DIR_EMULATOR)
# Bundle names carry the output folder for non-main runs, so they stay unique
BUNDLE_PREFIX     <- if (identical(EMULATOR_BASE_DIR, DIR_EMULATOR)) "emulator_fit" else paste0("emulator_fit_", basename(EMULATOR_BASE_DIR))
dir.create(EMULATOR_BASE_DIR, showWarnings = FALSE, recursive = TRUE)
SUBMIT_EMULATOR <- setting(OVR, "SUBMIT_EMULATOR", FALSE)   # Section C
AGG_EMULATOR    <- setting(OVR, "AGG_EMULATOR", FALSE)      # Section D
RDT_SENSITIVITY_RUN <- setting(OVR, "RDT_SENSITIVITY_RUN", RDT_SENSITIVITY_MAIN)

CAP_FLOOR <- TRUE

# Kernel per site. All sites use matern5_2 with multistart = 5. (GIN_Kindia and
# GIN_Boké once had an "exp" override; their spikes came from the old 1.5*IQR
# training floor, which floor_below_max = 100 removed.) Add a site here to
# override its kernel on every refit.
SITE_KERNEL_OVERRIDES <- list()
site_kernel_override <- function(site_label) {
  ov <- SITE_KERNEL_OVERRIDES[[site_label]]
  if (is.null(ov)) list(covtype = "matern5_2", multistart = 5) else ov
}

# Task wrapper: a site's ll_arr is a dense array on a rectangular grid (every
# site here) or a data.table for an irregular design; the irregular case goes
# straight to fit_single_emulator_design().
fit_emulator_task <- function(task_args) {
  if (is.data.frame(task_args$ll_arr)) {
    fit_single_emulator_design(
      task_args$ll_arr, task_args$site_label,
      max_design_points = task_args$max_design_points, seed = task_args$seed,
      cap_floor = task_args$cap_floor,
      covtype = task_args$covtype, multistart = task_args$multistart
    )
  } else {
    do.call(fit_single_emulator, task_args)
  }
}

# Hipercow setup ===============================================================
# SECTION A -- run once per hipercow root (interactively). Sources the
# function library only, never this driver script.
if (FALSE) {
  hipercow_environment_create(name = "emulator_fit", sources = "src/functions/emulator_fns.R")
}

# Loop over RDT sensitivity variants ===========================================
for (RDT_SENSITIVITY in RDT_SENSITIVITY_RUN) {
LOG_LIK_FILE     <- log_lik_file(RDT_SENSITIVITY)
EMULATOR_OUT_DIR <- rdt_sensitivity_dir(EMULATOR_BASE_DIR, RDT_SENSITIVITY)
BUNDLE_TAG       <- rdt_sensitivity_tag(RDT_SENSITIVITY)
cat(sprintf("\n########## 3b: RDT sensitivity %.2f -> %s ##########\n", RDT_SENSITIVITY, EMULATOR_OUT_DIR))

# Submit emulator fitting ======================================================
# SECTION C -- one task per site
if (SUBMIT_EMULATOR) {
  lik_obj     <- readRDS(file.path(GRID_OUT_DIR, LOG_LIK_FILE))

  # One row per (site, grid point); mean() is a no-op unless a grid point
  # appears twice.
  grid_point_dt <- lik_obj$by_name1_grid_point[
    , .(log_lik = mean(log_lik), n_obs = n_obs[1]),
    by = .(iso, name_1, eir, hrp2_scale)
  ]

  emulator_inputs <- build_emulator_inputs_from_grid_search(grid_point_dt)
  site_names <- names(emulator_inputs$ll_arrays)
  saveRDS(site_names, file.path(EMULATOR_OUT_DIR, "site_names.RDS"))
  # Seeds stay tied to each site's position in the full site list, so a
  # refitted site gets the same seed it would in a full run
  fit_idx <- if (is.null(REFIT_ISOS)) seq_along(site_names) else which(sub("_.*", "", site_names) %in% REFIT_ISOS)
  saveRDS(site_names[fit_idx], file.path(EMULATOR_OUT_DIR, "site_names_fitted.RDS"))

  # One list per site, named after fit_single_emulator()'s arguments.
  tasks <- lapply(fit_idx, function(i) {
    ov <- site_kernel_override(site_names[i])
    list(
      ll_arr            = emulator_inputs$ll_arrays[[i]],
      param_grids       = emulator_inputs$param_grids[[i]],
      site_label        = site_names[i],
      max_design_points = 2000L,
      seed              = 123L + i,
      cap_floor         = CAP_FLOOR,
      covtype           = ov$covtype,
      multistart        = ov$multistart
    )
  })

  bundle_fit <- task_create_bulk_call(
    fit_emulator_task,
    data        = tasks,
    environment = "emulator_fit",
    bundle_name = sprintf("%s_%s", BUNDLE_PREFIX, BUNDLE_TAG),
    # Single-threaded fits: 1 core/task lets the scheduler pack many per node.
    resources   = hipercow_resources(cores = 1)
  )
  saveRDS(bundle_fit, file.path(EMULATOR_OUT_DIR, "bundle_fit_new.RDS"))
  cat(sprintf("Submitted %d emulator-fit tasks | bundle id: %s\n",
              length(tasks), bundle_fit$name))
}

# Aggregate emulators ==========================================================
# SECTION D -- wait, retry transient failures, collect, validate
if (AGG_EMULATOR) {
  bundle_fit <- readRDS(file.path(EMULATOR_OUT_DIR, "bundle_fit_new.RDS"))
  site_names <- readRDS(file.path(EMULATOR_OUT_DIR, "site_names.RDS"))
  fitted_path <- file.path(EMULATOR_OUT_DIR, "site_names_fitted.RDS")
  # Older bundles fitted every site and have no site_names_fitted.RDS
  if (file.exists(fitted_path)) site_names <- readRDS(fitted_path)
  coll        <- collect_bundle(bundle_fit, label = "emulator-fit task")
  raw_results <- coll$raw
  ok          <- coll$ok

  emulators <- setNames(raw_results[ok], site_names[ok])
  if (!is.null(REFIT_ISOS)) {
    base_emulators <- load_emulators(file.path(rdt_sensitivity_dir(BASE_EMULATOR_DIR, RDT_SENSITIVITY), "emulators.rds"))
    kept <- base_emulators[!names(base_emulators) %in% site_names]
    cat(sprintf("Refitted %d sites (%s); %d other sites' emulators copied from %s.\n", sum(ok),
                paste(REFIT_ISOS, collapse = ", "), length(kept), BASE_EMULATOR_DIR))
    emulators <- c(kept, emulators)
  }
  save_emulators(emulators, file.path(EMULATOR_OUT_DIR, "emulators.rds"))

  # Failed sites are written to emulator_fit_failures.csv and raised as a
  # warning, so a dropped site can't go unnoticed.
  if (any(!ok)) {
    fail_msgs <- vapply(raw_results[!ok], conditionMessage, character(1))
    failure_log <- data.table::data.table(site = site_names[!ok], error = fail_msgs)
    cat("\n--- Emulator-fit FAILURES (dropped from emulators.rds) ---\n")
    print(failure_log)
    fwrite(failure_log, file.path(EMULATOR_OUT_DIR, "emulator_fit_failures.csv"))
    warning(sprintf(
      "%d/%d sites have NO emulator in emulators.rds (see emulator_fit_failures.csv in %s): %s",
      sum(!ok), length(ok), EMULATOR_OUT_DIR, paste(site_names[!ok], collapse = ", ")
    ), call. = FALSE)
  } else if (file.exists(file.path(EMULATOR_OUT_DIR, "emulator_fit_failures.csv"))) {
    # Every site succeeded this run -- remove any stale failure log from a
    # previous run so it can't be mistaken for a still-current problem.
    file.remove(file.path(EMULATOR_OUT_DIR, "emulator_fit_failures.csv"))
  }
  cat(sprintf("\nemulators.rds: %d/%d expected sites present.\n", length(emulators),
              length(readRDS(file.path(EMULATOR_OUT_DIR, "site_names.RDS")))))

  # Validation: leave-one-out CV on every site (cheap), then a hold-out refit
  # only for the sites LOO flags (validate_poor_loo_sites()); refitting all
  # ~250 sites would take hours.
  lik_obj     <- readRDS(file.path(GRID_OUT_DIR, LOG_LIK_FILE))

  grid_point_dt <- lik_obj$by_name1_grid_point[
    , .(log_lik = mean(log_lik), n_obs = n_obs[1]),
    by = .(iso, name_1, eir, hrp2_scale)
  ]

  emulator_inputs <- build_emulator_inputs_from_grid_search(grid_point_dt)

  fallback <- validate_poor_loo_sites(
    emulator_inputs$ll_arrays[names(emulators)],
    emulator_inputs$param_grids[names(emulators)],
    emulators
  )
  print(fallback$loo)
  saveRDS(fallback, file.path(EMULATOR_OUT_DIR, "emulator_validation.rds"))
}
}   # end RDT sensitivity loop
