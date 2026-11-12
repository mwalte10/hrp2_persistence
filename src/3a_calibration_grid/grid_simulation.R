# =============================================================================
# src/3a_calibration_grid/grid_simulation.R           (was src/04_grid_simulation.R)
#
# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (stage 1: likelihood grid)
#
# Step 3a, calibration stage 1: grid search. For every admin-1 site with MIS
# data, malariasimulation is run over a factorial grid of HRP2 half-life x
# absolute EIR on the DIDE cluster (hipercow), and the observation
# log-likelihood of each site's MIS surveys is computed at every grid point.
#
# Sections
#   A  hipercow set-up (once per hipercow root)
#   B  grid design and the site list (always runs; cheap)
#   C  build and cache each site's base parameters (site files, ITN, SMC, ...)
#   D  submit the grid simulations (only points not already simulated)
#   E  collect the grid simulations -> hrp2_grid_results_<SETUP_TAG>.RDS
#   F  log-likelihood per site and grid point, for each RDT sensitivity in
#      RDT_SENSITIVITY_RUN -> hrp2_grid_log_lik[_rdt_se_NNN].RDS
# The RDT sensitivity analysis (0.8, 1.0) only repeats Section F; run it with
# sensitivity/rdt_sensitivity_log_lik.R. The simulations do not depend on RDT
# sensitivity.
#
# Observation model (Section F, LIKELIHOOD = "multinomial_deff"): one
# multinomial over the joint RDT x microscopy 2x2 (TP/FP/FN/TN) for children
# tested by both, each survey down-weighted by its design effect and by the
# number of fitted sites sharing its DHS-region counts; surveys with only one
# test use that test's weighted binomial. Deletion and RDT sensitivity act on
# RDT only: an HRP2-positive child reads RDT-positive with probability
# (1 - deletion) x sensitivity.
#
# Grid (Section B): 20 log-spaced half-lives (GRID_HALFLIFE_RANGE, 5-365 d)
# x 24 log-spaced absolute EIRs (grid_eir_values(), ~0.5-2,069), shared by
# every site, plus extra EIR values for the Rwanda sites whose likelihood peak
# is narrower than the EIR spacing (cluster_fns.R::GRID_SITE_EIR_REFINEMENT).
# Both axes are log-spaced to concentrate points at short half-lives and low
# EIR, where the likelihood surface curves most.
#
# Simulation setup (SETUP_TAG "mapitn_burnin5"): ITN usage converted with
# MAP's net-loss curve and a 5-year burn-in, as in the site-file EIR
# calibration (cluster_fns.R::site_itn_input_dist(), build_adm1_base_pars()).
#
# Outputs (DIR_GRID): adm1_combos_with_paths.RDS, base_pars_cache_<tag>/,
#   hrp2_grid_results_<tag>.RDS, hrp2_grid_log_lik[_rdt_se_NNN].RDS
#
# Archived (archive/src/04_grid_simulation.R): the Latin-hypercube design
# option, the edge-pinned site refinement rounds (Sections G0/G/H) and the
# comparison against the 2026-09-11 original-setup grid.
# =============================================================================
library(malariasimulation)
library(data.table)
library(hipercow)

source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

OUT_DIR             <- DIR_GRID
SETUP_TAG           <- GRID_SETUP_TAG   # cluster_fns.R, shared with steps 3b-4
BASE_PARS_CACHE_DIR <- file.path(OUT_DIR, sprintf("base_pars_cache_%s", SETUP_TAG))
BASE_PARS_PATH_COL  <- sprintf("base_pars_path_%s", SETUP_TAG)
GRID_SETUP          <- SETUP_TAG
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Settings ======================================================================
# Every cluster section is off by default; turn on one at a time, or run e.g.
#   Rscript run_sections.R src/3a_calibration_grid/grid_simulation.R SUBMIT_GRID=TRUE
OVR <- take_overrides()
BUILD_CACHE1    <- setting(OVR, "BUILD_CACHE1", FALSE)      # Section C: submit base_pars builds
BUILD_CACHE2    <- setting(OVR, "BUILD_CACHE2", FALSE)      # Section C: collect base_pars paths
SUBMIT_GRID     <- setting(OVR, "SUBMIT_GRID", FALSE)       # Section D
AGG_GRID        <- setting(OVR, "AGG_GRID", FALSE)          # Section E
CALC_LIKELIHOOD <- setting(OVR, "CALC_LIKELIHOOD", FALSE)   # Section F
RDT_SENSITIVITY_RUN <- setting(OVR, "RDT_SENSITIVITY_RUN", RDT_SENSITIVITY_MAIN)
# Section F inputs/outputs for the sensitivity analyses (defaults = main analysis):
#   GRID_RESULTS_PATH  grid results to read (NULL: this setup's results file in OUT_DIR)
#   LOG_LIK_OUT_DIR    folder the log-likelihood files are written to
#   S_U_RATIO          s_U / s, the RDT detection probability of HRP2-flagged
#                      children in U relative to other flagged children (NULL:
#                      no split, as the main analysis; needs p_fp_u in the results)
#   S_U_ISOS           countries the s_U split applies to (NULL: all); elsewhere
#                      fp_U = 0, which leaves the main cells unchanged
#   GRID_POINTS_FROM   a log-likelihood file whose grid points to keep (NULL: the current design)
GRID_RESULTS_PATH <- setting(OVR, "GRID_RESULTS_PATH", NULL)
GRID_POINTS_FROM  <- setting(OVR, "GRID_POINTS_FROM", NULL)
LOG_LIK_OUT_DIR   <- setting(OVR, "LOG_LIK_OUT_DIR", OUT_DIR)
S_U_RATIO         <- setting(OVR, "S_U_RATIO", NULL)
S_U_ISOS          <- setting(OVR, "S_U_ISOS", NULL)
# Observation model for Section F:
#   "binomial"         independent RDT and LM binomials on raw counts (v1)
#   "binomial_deff"    the same, down-weighted by design effect and shared regions (v2)
#   "multinomial_deff" joint RDT x LM multinomial, same down-weighting (v3, main analysis)
LIKELIHOOD <- setting(OVR, "LIKELIHOOD", "multinomial_deff")

# Hipercow setup ===============================================================
# SECTION A -- run once per hipercow root (interactively): creates the
# "calibration" environment the simulation tasks run in and installs packages.
if (FALSE) {
  hipercow_init(root = PROJECT_ROOT, driver = "dide-windows")
  hipercow_environment_create(name = "calibration", sources = "src/functions/cluster_fns.R")
  hipercow::hipercow_provision("pkgdepends")
}

# Build grid ===================================================================
# SECTION B -- GRID DESIGN AND SITE LIST
# hrp2_deletion_pct is not a grid axis: Section F applies each country's fixed
# deletion prevalence to the simulated RDT outcomes, so it costs no simulation.
# Bounds come from cluster_fns.R (GRID_HALFLIFE_RANGE, grid_eir_values()) so
# steps 3b-4 see the range the emulators are trained on.
hrp2_halflife_grid <- exp(seq(log(GRID_HALFLIFE_RANGE[1]), log(GRID_HALFLIFE_RANGE[2]), length.out = 20))
hrp2_scale_grid     <- halflife_to_scale(hrp2_halflife_grid)
EIR_GRID            <- grid_eir_values()   # 24 points, see cluster_fns.R (original 20 + 2 below + 2 above)

# The current grid's (hrp2_scale, eir) design points, shared by every site.
current_grid_points <- function() data.table::CJ(hrp2_scale = hrp2_scale_grid, eir = EIR_GRID)
# Join key robust to floating-point noise between recomputed and stored values.
grid_key <- function(hrp2_scale, eir) sprintf("%.6f_%.6f", log(hrp2_scale), log(eir))

# Every design point for one site: the shared grid plus that site's extra EIR
# values (cluster_fns.R::GRID_SITE_EIR_REFINEMENT) at every grid half-life.
# Sections D, E and F all use this, so refined sites are simulated, merged and
# kept in the likelihood like any other grid point.
site_grid_points <- function(iso, name_1) {
  base  <- current_grid_points()
  extra <- site_refinement_eir_values(paste(iso, name_1, sep = "_"))
  if (length(extra) == 0) return(base)
  unique(rbind(base, data.table::CJ(hrp2_scale = sort(unique(base$hrp2_scale)), eir = extra)))
}
# All sites' design points, keyed by site and grid point.
all_site_grid_points <- function(sites) {
  data.table::rbindlist(lapply(seq_len(nrow(sites)), function(i)
    site_grid_points(sites$iso[i], sites$name_1[i])[, `:=`(iso = sites$iso[i], name_1 = sites$name_1[i])]))
}

# get_adm1_eir_weighted() only has an EIR baseline for the isos in the OLD
# cache (EIR_SITE_CACHE_DIR) — restricts scope to those 16 isos even though
# the new primary cache (SITE_CACHE_DIR) covers 67.
new_cache_isos <- unlist(lapply(
  list.files(SITE_CACHE_DIR, pattern = "_site\\.rds$"),
  gsub, pattern = "_site.rds", replacement = ""
))
eir_cache_isos <- unlist(lapply(
  list.files(EIR_SITE_CACHE_DIR, pattern = "_site\\.RDS$"),
  gsub, pattern = "_site.RDS", replacement = ""
))
isos <- intersect(new_cache_isos, eir_cache_isos)
cat(sprintf("%d isos in the new site cache, %d with an EIR baseline in the old cache -> %d usable isos.\n",
            length(new_cache_isos), length(eir_cache_isos), length(isos)))

# One row per (iso, name_1) — the unit of the expensive base_pars build.
# Reading each cached site file locally is cheap, so this stays off the
# cluster. eir_in is looked up here too, kept only as an independent external
# reference value for reference only (not used to scale the grid or by the MCMC).
adm1_combos <- data.table::rbindlist(lapply(isos, list_adm1_combos))
cat(sprintf("%d iso x adm1 combinations.\n", nrow(adm1_combos)))

# Grouped by iso so each old-cache site file is read once (16 reads) rather
# than once per name_1 (304 reads) — get_adm1_eir_weighted_from_site() takes
# the already-loaded object.
adm1_combos[, eir_in := {
  old_site <- readRDS(file.path(EIR_SITE_CACHE_DIR, paste0(.BY$iso, "_site.RDS")))
  mapply(get_adm1_eir_weighted_from_site, name_1 = name_1,
         MoreArgs = list(old_site = old_site, iso = .BY$iso))
}, by = iso]
if (anyNA(adm1_combos$eir_in) || any(lengths(adm1_combos$eir_in) == 0)) {
  stop("Missing site eir_in for some adm1 combos — check old-cache $eir coverage.")
}
adm1_combos[, eir_in := as.numeric(eir_in)]

# Only simulate sites with MIS data
mis <- readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS"))
adm1_combos <- adm1_combos[name_1 %in% unique(mis[!is.na(pfpr_lm_num) | !is.na(pfpr_rdt_num),name_1]),]

cat(sprintf("%d hrp2_scale grid points per adm1.\n", length(hrp2_scale_grid)))
cat(sprintf("%d EIR grid points per adm1 (absolute, %.2f-%.2f, shared identically across all adm1s).\n",
            length(EIR_GRID), min(EIR_GRID), max(EIR_GRID)))

n_grid_points_per_adm1 <- length(hrp2_scale_grid) * length(EIR_GRID)
n_total_sims <- nrow(adm1_combos) * n_grid_points_per_adm1
cat(sprintf("%s grid points per adm1 (hrp2_scale x eir); %s total simulations across the full grid.\n",
            formatC(n_grid_points_per_adm1, format = "d", big.mark = ","),
            formatC(n_total_sims, format = "d", big.mark = ",")))

# Build cache ==================================================================
# SECTION C -- build and cache base_pars once per site (the expensive,
# grid-independent part of each simulation), then record the cache paths in
# column BASE_PARS_PATH_COL of adm1_combos_with_paths.RDS.
if (BUILD_CACHE1) {

  bundle_build <- task_create_bulk_call(
    build_adm1_base_pars,
    # Only the columns build_adm1_base_pars() actually takes as arguments —
    # task_create_bulk_call() passes every data.frame column as a named arg,
    # and the function has no ... to absorb extras like eir_in.
    data        = as.data.frame(adm1_combos[, .(iso, name_1)]),
    args        = list(base_pars_cache_dir = BASE_PARS_CACHE_DIR),
    environment = "calibration",
    bundle_name = sprintf("base_pars_build_%s", SETUP_TAG),
    # Each task is single-threaded (confirmed: run_simulation()'s user/elapsed
    # time ratio is ~1 regardless of RcppParallel thread settings), so
    # requesting 1 core/task lets the scheduler pack ~32 tasks per node on a
    # 32-core-per-node cluster instead of reserving whole nodes per task.
    resources   = hipercow_resources(cores = 1)
  )
  saveRDS(bundle_build, file.path(OUT_DIR, sprintf("bundle_build_%s.RDS", SETUP_TAG)))
  cat(sprintf("Submitted %d base_pars build tasks | bundle id: %s\n",
              nrow(adm1_combos), bundle_build$name))
}

# --- wait for Phase 1, collect cached base_pars paths ---
if (BUILD_CACHE2) {

  bundle_build <- readRDS(file.path(OUT_DIR, sprintf("bundle_build_%s.RDS", SETUP_TAG)))
  hipercow_bundle_wait(bundle_build, timeout = Inf)

  statuses <- hipercow_bundle_status(bundle_build)
  if (any(statuses != "success")) {
    failed <- which(statuses != "success")
    stop(sprintf("%d base_pars build task(s) did not succeed: indices %s",
                 length(failed), paste(failed, collapse = ", ")))
  }

  # Write this setup's paths into their own column of the EXISTING table
  # (keeping base_pars_path, the original setup's, and any other columns).
  new_paths <- data.table::data.table(adm1_combos[, .(iso, name_1)],
                                      path = unlist(hipercow_bundle_result(bundle_build)))
  paths_file <- file.path(OUT_DIR, "adm1_combos_with_paths.RDS")
  combos_out <- data.table::data.table(readRDS(paths_file))
  if (BASE_PARS_PATH_COL %in% names(combos_out)) combos_out[, (BASE_PARS_PATH_COL) := NULL]
  combos_out <- merge(combos_out, data.table::setnames(new_paths, "path", BASE_PARS_PATH_COL),
                      by = c("iso", "name_1"), all.x = TRUE, sort = FALSE)
  saveRDS(combos_out, paths_file)
  cat(sprintf("Cached base_pars (%s) for %d adm1 combos -> column %s of %s.\n",
              SETUP_TAG, nrow(new_paths), BASE_PARS_PATH_COL, paths_file))
}

# Submit grid ==================================================================
# SECTION D -- run the eir x hrp2 grid, chunked per site
CHUNK_SIZE <- 10L
bundle_grid_path      <- file.path(OUT_DIR, sprintf("bundle_grid_%s.RDS", SETUP_TAG))
setup_results_path    <- file.path(OUT_DIR, grid_results_file(SETUP_TAG))
bundle_collected_path <- file.path(OUT_DIR, sprintf("bundle_grid_%s_collected.txt", SETUP_TAG))
bundle_is_collected <- function(bundle_name) {
  file.exists(bundle_collected_path) && bundle_name %in% readLines(bundle_collected_path)
}
site_point_key <- function(iso, name_1, hrp2_scale, eir) paste(iso, name_1, grid_key(hrp2_scale, eir), sep = "|")

if (SUBMIT_GRID) {

  if (file.exists(bundle_grid_path)) {
    prev_bundle <- readRDS(bundle_grid_path)
    if (!bundle_is_collected(prev_bundle$name)) {
      # Bundles collected before this log existed: the results file being newer
      # than the bundle file means Section E already ran on it.
      if (file.exists(setup_results_path) && file.mtime(setup_results_path) > file.mtime(bundle_grid_path)) {
        cat(prev_bundle$name, "\n", file = bundle_collected_path, append = TRUE, sep = "")
      } else {
        stop(sprintf("Section D: previous bundle '%s' hasn't been collected -- run Section E (AGG_GRID) first.", prev_bundle$name))
      }
    }
  }

  adm1_paths <- data.table::data.table(readRDS(file.path(OUT_DIR, "adm1_combos_with_paths.RDS")))
  if (!BASE_PARS_PATH_COL %in% names(adm1_paths)) {
    stop(sprintf("Section D: no %s column -- build this setup's base_pars cache first (BUILD_CACHE1, then BUILD_CACHE2).", BASE_PARS_PATH_COL))
  }
  grid_sites <- merge(adm1_combos[, .(iso, name_1)],
                      adm1_paths[, .(iso, name_1, base_pars_path = get(BASE_PARS_PATH_COL))],
                      by = c("iso", "name_1"))
  if (anyNA(grid_sites$base_pars_path)) stop("Section D: sites without a cached base_pars_path for this setup.")

  # Full design (every site's own design points: shared grid + any EIR
  # refinement), minus what's already simulated
  design <- merge(all_site_grid_points(grid_sites[, .(iso, name_1)]), grid_sites, by = c("iso", "name_1"))
  n_design <- nrow(design)
  if (file.exists(setup_results_path)) {
    done <- unique(data.table::data.table(readRDS(setup_results_path))[, .(iso, name_1, hrp2_scale, eir)])
    done_keys <- site_point_key(done$iso, done$name_1, done$hrp2_scale, done$eir)
    design <- design[!site_point_key(iso, name_1, hrp2_scale, eir) %in% done_keys]
  }
  cat(sprintf("Section D: %d sites, %s design points (%d shared grid points per site + %d refinement points); %s already simulated; %s to submit.\n",
              nrow(grid_sites), formatC(n_design, format = "d", big.mark = ","), nrow(current_grid_points()),
              n_design - nrow(grid_sites) * nrow(current_grid_points()),
              formatC(n_design - nrow(design), format = "d", big.mark = ","),
              formatC(nrow(design), format = "d", big.mark = ",")))
  if (nrow(design) > 0) {
    cat("New points by site:\n"); print(design[, .(points = .N, eir_values = uniqueN(eir)), by = .(iso, name_1)])
  }

  if (nrow(design) == 0) {
    cat("Section D: nothing to submit -- every site already has the full current grid.\n")
  } else {
    tasks <- unlist(lapply(split(design, by = c("iso", "name_1")), function(d) {
      lapply(split(d, ceiling(seq_len(nrow(d)) / CHUNK_SIZE)), function(chunk) list(
        base_pars_path = chunk$base_pars_path[1],
        hrp2_grid      = chunk[, .(hrp2_scale, eir)],
        iso            = chunk$iso[1],
        name_1         = chunk$name_1[1],
        joint          = TRUE
      ))
    }), recursive = FALSE)
    cat(sprintf("Submitting %s grid tasks (chunks of up to %d; setup %s).\n",
                formatC(length(tasks), format = "d", big.mark = ","), CHUNK_SIZE, SETUP_TAG))

    run_hrp2_grid_task <- function(task_args) do.call(run_hrp2_grid_adm1, task_args)

    bundle_grid <- task_create_bulk_call(
      run_hrp2_grid_task,
      data        = tasks,
      environment = "calibration",
      bundle_name = sprintf("hrp2_grid_%s_%s", SETUP_TAG, format(Sys.time(), "%Y%m%d_%H%M")),
      # See Phase 1: single-threaded per task, so 1 core/task lets ~32 tasks
      # pack onto each node of a 32-core-per-node cluster.
      resources   = hipercow_resources(cores = 1)
    )
    # Keep the previous (collected) bundle's file under its own name
    if (file.exists(bundle_grid_path)) {
      file.copy(bundle_grid_path, sub("\\.RDS$", sprintf("__%s.RDS", readRDS(bundle_grid_path)$name), bundle_grid_path))
    }
    saveRDS(bundle_grid, bundle_grid_path)
    cat(sprintf("Submitted %d grid tasks | bundle id: %s\n", length(tasks), bundle_grid$name))
  }
}

# Aggregate grid ===============================================================
# SECTION E -- collect, merge into the setup's results file
if (AGG_GRID) {
  bundle_grid <- readRDS(bundle_grid_path)
  coll <- collect_bundle(bundle_grid, label = "grid task")
  cat(sprintf("Compiling results from %d/%d successful tasks.\n", sum(coll$ok), length(coll$ok)))
  results <- data.table::rbindlist(coll$results)
  missing_joint <- setdiff(c("p_tp", "p_fp", "p_fn", "p_tn"), names(results))
  if (length(missing_joint) > 0) {
    stop(sprintf("Section E: results lack joint column(s) %s -- was this bundle submitted with joint = TRUE?",
                 paste(missing_joint, collapse = ", ")))
  }
  cat(sprintf("Joint cells sum to 1 within %.2g.\n", max(abs(results[, p_tp + p_fp + p_fn + p_tn] - 1))))
  # Merge into the existing results file (Section D submits only missing
  # points): new rows replace any existing rows for the same site/grid point.
  if (file.exists(setup_results_path)) {
    existing  <- data.table::data.table(readRDS(setup_results_path))
    new_keys  <- unique(site_point_key(results$iso, results$name_1, results$hrp2_scale, results$eir))
    existing  <- existing[!site_point_key(iso, name_1, hrp2_scale, eir) %in% new_keys]
    cat(sprintf("Merging %s new rows into %s existing rows.\n",
                formatC(nrow(results), format = "d", big.mark = ","), formatC(nrow(existing), format = "d", big.mark = ",")))
    results <- data.table::rbindlist(list(existing, results), use.names = TRUE, fill = TRUE)
  }
  pts <- unique(results[, .(iso, name_1, hrp2_scale, eir)])
  want <- all_site_grid_points(unique(pts[, .(iso, name_1)]))
  pts[, on_design := site_point_key(iso, name_1, hrp2_scale, eir) %in% site_point_key(want$iso, want$name_1, want$hrp2_scale, want$eir)]
  cover <- merge(pts[on_design == TRUE, .(have = .N), by = .(iso, name_1)], want[, .(want = .N), by = .(iso, name_1)], by = c("iso", "name_1"))
  cat(sprintf("Design coverage: %d sites; %d complete, %d incomplete (shared grid %d points per site; refined sites have more).\n",
              nrow(cover), cover[have == want, .N], cover[have < want, .N], nrow(current_grid_points())))

  saveRDS(results, setup_results_path)
  cat(bundle_grid$name, "\n", file = bundle_collected_path, append = TRUE, sep = "")
  cat(sprintf("Saved %s rows to %s\n",
              formatC(nrow(results), format = "d", big.mark = ","), setup_results_path))
}

# Calculate likelihood =========================================================
# SECTION F -- log-likelihood of the observed MIS data at every grid point
if (CALC_LIKELIHOOD) {
  include <- readRDS(file.path(DIR_STRUCTURE, "countries_included.RDS"))

  stopifnot(LIKELIHOOD %in% c("binomial", "binomial_deff", "multinomial_deff"))
  cat(sprintf("Section F: LIKELIHOOD = '%s' on the '%s' grid.\n", LIKELIHOOD, GRID_SETUP))
  # Always read from disk: an in-memory `results` left over from Section E or
  # an earlier session could be a different grid or setup.
  results_path <- if (is.null(GRID_RESULTS_PATH)) file.path(OUT_DIR, grid_results_file(GRID_SETUP)) else GRID_RESULTS_PATH
  results_file <- basename(results_path)
  if (!file.exists(results_path)) {
    stop(sprintf("Section F: %s not found -- run Sections D and E (SUBMIT_GRID, then AGG_GRID) for this setup first.", results_path))
  }
  results <- data.table::data.table(readRDS(results_path))

  # Current design only: rows from older or wider grids in the results file
  # (e.g. the archived 2026-09 refinement rounds) are dropped.
  # Each site's own design (shared grid + any EIR refinement, site_grid_points())
  want    <- all_site_grid_points(unique(results[, .(iso, name_1)]))
  # GRID_POINTS_FROM: keep only the grid points of an existing log-likelihood
  # file, so a sensitivity analysis uses exactly the design its reference fit did
  if (!is.null(GRID_POINTS_FROM)) {
    ref_pts <- unique(data.table::as.data.table(readRDS(GRID_POINTS_FROM)$by_name1_grid_point)[, .(iso, name_1, eir, hrp2_scale)])
    ref_keys <- site_point_key(ref_pts$iso, ref_pts$name_1, ref_pts$hrp2_scale, ref_pts$eir)
    n_want <- nrow(want)
    want <- want[site_point_key(iso, name_1, hrp2_scale, eir) %in% ref_keys]
    cat(sprintf("Section F: design restricted to the %s grid points of %s (%s of the current design's %s).\n",
                formatC(nrow(ref_pts), format = "d", big.mark = ","), basename(GRID_POINTS_FROM),
                formatC(nrow(want), format = "d", big.mark = ","), formatC(n_want, format = "d", big.mark = ",")))
  }
  want_keys <- site_point_key(want$iso, want$name_1, want$hrp2_scale, want$eir)
  on_grid <- site_point_key(results$iso, results$name_1, results$hrp2_scale, results$eir) %in% want_keys
  if (any(!on_grid)) {
    cat(sprintf("Section F: dropped %s rows not in any site's current design.\n",
                formatC(sum(!on_grid), format = "d", big.mark = ",")))
  }
  results <- results[on_grid]
  pts_per_site <- unique(results[, .(iso, name_1, hrp2_scale, eir)])[, .(have = .N), by = .(iso, name_1)]
  pts_per_site <- merge(pts_per_site, want[, .(want = .N), by = .(iso, name_1)], by = c("iso", "name_1"))
  incomplete   <- pts_per_site[have < want]
  if (nrow(incomplete) > 0) {
    print(incomplete)
    stop(sprintf("Section F: %d site(s) are missing design points (failed tasks, or refinement not yet simulated?).", nrow(incomplete)))
  }
  cat(sprintf("Section F: %d sites (%d with EIR refinement), %d-%d design points per site, from %s.\n",
              nrow(pts_per_site), pts_per_site[want > nrow(current_grid_points()), .N],
              min(pts_per_site$want), max(pts_per_site$want), basename(results_path)))

  # Not restricted to a single iso: name_1 alone isn't a unique key across
  # countries (two isos can share an admin1 name), so every join/grouping below
  # keys on (iso, name_1), not name_1 alone.
  mis <- readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS"))
  mis <- mis[country %in% c(include$country),]
  mis_dt <- unique(data.table::data.table(mis)[
    (!is.na(pfpr_rdt_num) & !is.na(pfpr_rdt_denom)) |
      (!is.na(pfpr_lm_num)  & !is.na(pfpr_lm_denom)),
    .(iso = iso3c, name_1, year, pfpr_rdt_num, pfpr_rdt_denom, pfpr_lm_num, pfpr_lm_denom,
      joint_n, joint_tp, joint_fp, joint_fn, joint_tn,
      src_rdt, src_lm, src_joint, deff_rdt, deff_lm, deff_joint)
  ])

  sim_cols <- c("iso", "name_1", "year", "eir", "hrp2_scale", "hrp2_prev", "lm_prev",
                if (LIKELIHOOD == "multinomial_deff") c("p_tp", "p_fp", "p_fn", "p_tn"),
                if (!is.null(S_U_RATIO)) "p_fp_u")
  if (!is.null(S_U_RATIO) && LIKELIHOOD != "multinomial_deff") stop("Section F: S_U_RATIO needs LIKELIHOOD = 'multinomial_deff'.")
  missing_sim <- setdiff(sim_cols, names(results))
  if (length(missing_sim) > 0) {
    stop(sprintf("Section F (LIKELIHOOD = '%s'): grid results lack column(s) %s.", LIKELIHOOD, paste(missing_sim, collapse = ", ")))
  }
  admin1_prev <- data.table::data.table(results)[, ..sim_cols]
  if (!is.null(S_U_RATIO)) {
    # Only the replicate runs (RWA, SEN) rendered U; other sites keep the main cells
    if (!is.null(S_U_ISOS)) admin1_prev[!iso %in% S_U_ISOS, p_fp_u := 0]
    if (anyNA(admin1_prev$p_fp_u)) {
      stop(sprintf("Section F: p_fp_u missing for %s -- set S_U_ISOS to the countries whose runs rendered U.",
                   paste(unique(admin1_prev[is.na(p_fp_u), iso]), collapse = ", ")))
    }
  }

  lik <- merge(admin1_prev, mis_dt, by = c("iso", "name_1", "year"))

  # --- Observation weights: w = 1 / (deff x n_share), see this section's header ---
  # n_share counts the fitted sites (those actually in lik, i.e. with
  # simulations) drawing on each (survey, DHS region) source row. A survey
  # with no design-effect estimate (no usable regions in 02b) falls back to
  # the median across surveys.
  obs_keys <- unique(lik[, .(iso, name_1, year, src_rdt, src_lm, src_joint)])
  for (t in c("rdt", "lm", "joint")) {
    src <- paste0("src_", t)
    obs_keys[, (paste0("n_share_", t)) := data.table::fifelse(is.na(get(src)), 1L, .N), by = src]
  }
  lik <- merge(lik, obs_keys[, .(iso, name_1, year, n_share_rdt, n_share_lm, n_share_joint)],
               by = c("iso", "name_1", "year"))
  for (t in c("rdt", "lm", "joint")) {
    dc <- paste0("deff_", t)
    per_source <- unique(data.table::data.table(src = lik[[paste0("src_", t)]], d = lik[[dc]]))
    fallback <- stats::median(per_source$d, na.rm = TRUE)
    n_fill <- uniqueN(lik[is.na(get(dc)), .(iso, name_1, year)])
    if (n_fill > 0) cat(sprintf("deff_%s: %d site-year(s) with no survey design effect -> median %.2f.\n", t, n_fill, fallback))
    lik[is.na(get(dc)), (dc) := fallback]
  }
  if (LIKELIHOOD == "binomial") {
    lik[, `:=`(w_rdt = 1, w_lm = 1, w_joint = 1)]
  } else {
    lik[, `:=`(w_rdt   = 1 / (deff_rdt   * n_share_rdt),
               w_lm    = 1 / (deff_lm    * n_share_lm),
               w_joint = 1 / (deff_joint * n_share_joint))]
  }
  wsum <- unique(lik[, .(iso, name_1, year, w_rdt, w_lm, w_joint, n_share_rdt, deff_rdt)])
  cat(sprintf("LIKELIHOOD = '%s': median obs weight RDT %.3f, LM %.3f, joint %.3f; %d of %d site-years share their RDT source row.\n",
              LIKELIHOOD, median(wsum$w_rdt), median(wsum$w_lm), median(wsum$w_joint),
              sum(wsum$n_share_rdt > 1), nrow(wsum)))
  
  # --- HRP2 (pfhrp2/3) deletion post-processing, BEFORE the likelihood ---
  # Deletion-positive parasites don't express HRP2, so they're false negatives
  # on an HRP2-based RDT — the observable HRP2-positive population is the
  # non-deleted fraction of the true simulated positives.
  #
  # hrp2_deletion_pct is FIXED at each country's central value in
  # data/hrp2_deletion.csv (0% if the country has no row;
  # cluster_fns.R::get_hrp2_deletion_mean()). It is reported at country level
  # only, and as a 4th emulator dimension it gave degenerate GP
  # hyperparameters, so it is not a grid axis.
  deletion_dt <- data.table::fread(HRP2_DELETION_CSV)
  lik_isos <- unique(lik$iso)
  deletion_info_by_iso <- data.table::rbindlist(lapply(lik_isos, function(i) {
    m <- get_hrp2_deletion_mean(i, deletion_dt)
    data.table::data.table(iso = i, hrp2_deletion_pct = m$mean_pct, source = m$source)
  }))
  cat(sprintf("HRP2 deletion: %d isos from file, %d using the %g%% fallback (no row).\n",
              sum(deletion_info_by_iso$source == "file"),
              sum(deletion_info_by_iso$source == "fallback"),
              formals(get_hrp2_deletion_mean)$fallback_mean_pct))
  
  lik_base <- merge(lik, deletion_info_by_iso, by = "iso")

  # One likelihood per value in RDT_SENSITIVITY_RUN (main analysis: 0.9; the
  # sensitivity script runs 0.8 and 1.0).
  # Sensitivity acts on HRP2-positive children after deletion, so an
  # HRP2-positive child reads RDT-positive with probability (1 - del) x se.
  # Each value is saved to cluster_fns.R::log_lik_file(se); the main value
  # (RDT_SENSITIVITY_MAIN) keeps the unsuffixed hrp2_grid_log_lik.RDS.
  for (rdt_se in RDT_SENSITIVITY_RUN) {
  cat(sprintf("\n=== RDT sensitivity %.2f%s ===\n", rdt_se, if (is_main_rdt_sensitivity(rdt_se)) " (main)" else ""))
  lik <- data.table::copy(lik_base)
  lik[, rdt_sensitivity := rdt_se]
  lik[, hrp2_prev_adj := rdt_prev_observed(hrp2_prev, hrp2_deletion_pct, rdt_se)]
  # s_U: flagged U children read RDT-positive with r_U = (1 - del) x s x S_U_RATIO.
  # For RDT-only surveys the U share is taken from the joint (6-59 month) band.
  if (!is.null(S_U_RATIO)) {
    lik[, hrp2_prev_adj := hrp2_prev_adj - (1 - hrp2_deletion_pct / 100) * rdt_se * (1 - S_U_RATIO) * p_fp_u]
  }

  # --- joint binomial log-likelihood vs observed MIS RDT + LM prevalence ---
  lik[, prob_rdt := pmin(pmax(hrp2_prev_adj, 1e-6), 1 - 1e-6)]
  lik[, prob_lm  := pmin(pmax(lm_prev,       1e-6), 1 - 1e-6)]
  
  lik[, log_lik_rdt := data.table::fifelse(
    !is.na(pfpr_rdt_num) & !is.na(pfpr_rdt_denom),
    w_rdt * stats::dbinom(x = pfpr_rdt_num, size = pfpr_rdt_denom, prob = prob_rdt, log = TRUE),
    0
  )]
  lik[, log_lik_lm := data.table::fifelse(
    !is.na(pfpr_lm_num) & !is.na(pfpr_lm_denom),
    w_lm * stats::dbinom(x = pfpr_lm_num, size = pfpr_lm_denom, prob = prob_lm, log = TRUE),
    0
  )]

  if (LIKELIHOOD == "multinomial_deff") {
    # Joint 2x2 cell probabilities, deletion and RDT sensitivity acting on
    # RDT only: an HRP2-positive child reads RDT-positive with probability
    # r = (1-d) x se, so tp' = r tp, fp' = r fp, fn' = fn + (1-r) tp,
    # tn' = tn + (1-r) fp. The direct generalisation of hrp2_prev_adj
    # above -- tp' + fp' = r (tp + fp).
    lik[, r_detect := (1 - hrp2_deletion_pct / 100) * rdt_sensitivity]
    if (is.null(S_U_RATIO)) {
      lik[, `:=`(q_tp = r_detect * p_tp, q_fp = r_detect * p_fp,
                 q_fn = p_fn + (1 - r_detect) * p_tp, q_tn = p_tn + (1 - r_detect) * p_fp)]
    } else {
      # U children are LM-negative, so they sit only in FP (flagged) and TN:
      #   fp' = r fp_nonU + r_U fp_U,  tn' = tn + (1-r) fp_nonU + (1-r_U) fp_U
      lik[, r_detect_u := r_detect * S_U_RATIO]
      lik[, `:=`(q_tp = r_detect * p_tp,
                 q_fp = r_detect * (p_fp - p_fp_u) + r_detect_u * p_fp_u,
                 q_fn = p_fn + (1 - r_detect) * p_tp,
                 q_tn = p_tn + (1 - r_detect) * (p_fp - p_fp_u) + (1 - r_detect_u) * p_fp_u)]
      cat(sprintf("S_U_RATIO = %.2f: cells sum to 1 within %.2g before flooring; p_fp_u <= p_fp everywhere: %s.\n",
                  S_U_RATIO, max(abs(lik[, q_tp + q_fp + q_fn + q_tn] - 1)), all(lik$p_fp_u <= lik$p_fp + 1e-12)))
    }
    q_cols <- c("q_tp", "q_fp", "q_fn", "q_tn")
    lik[, (q_cols) := lapply(.SD, function(q) pmax(q, 1e-6)), .SDcols = q_cols]
    lik[, q_sum := q_tp + q_fp + q_fn + q_tn]
    lik[, (q_cols) := lapply(.SD, function(q) q / q_sum), .SDcols = q_cols]
    # Multinomial log-pmf, vectorised (dmultinom() isn't).
    lik[, log_lik_joint := data.table::fifelse(
      !is.na(joint_n) & joint_n > 0,
      w_joint * (lgamma(joint_n + 1) - lgamma(joint_tp + 1) - lgamma(joint_fp + 1) - lgamma(joint_fn + 1) - lgamma(joint_tn + 1) +
                 joint_tp * log(q_tp) + joint_fp * log(q_fp) + joint_fn * log(q_fn) + joint_tn * log(q_tn)),
      NA_real_
    )]
    # Rows with joint counts use the multinomial in place of both binomials;
    # rows with only one test keep that test's (weighted) binomial.
    lik[, log_lik := data.table::fifelse(!is.na(log_lik_joint), log_lik_joint, log_lik_rdt + log_lik_lm)]
  } else {
    lik[, log_lik := log_lik_rdt + log_lik_lm]
  }
  
  # total joint log-likelihood per (iso, name_1, eir, hrp2_scale), summed
  # across that admin1's survey years only. eir/hrp2_scale are still chosen
  # freely PER (iso, name_1) — treatment-seeking/RDT behavior and transmission
  # intensity plausibly vary spatially within a country. hrp2_deletion_pct is
  # NOT part of this grid at all now (fixed per-country above, see the
  # comment on deletion_info_by_iso), so there's no profiling step needed for
  # it — best_by_name1 is a direct argmax per admin1.
  name1_grid_log_lik <- lik[, .(
    log_lik = sum(log_lik),
    n_obs   = .N
  ), by = .(iso, name_1, eir, hrp2_scale)]
  
  best_by_name1 <- name1_grid_log_lik[
    , .SD[which.max(log_lik)], by = .(iso, name_1)
  ][, .(iso, name_1, eir, hrp2_scale, log_lik, n_obs)]
  best_by_name1 <- merge(best_by_name1, deletion_info_by_iso[, .(iso, hrp2_deletion_pct)], by = "iso")

  dir.create(LOG_LIK_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(LOG_LIK_OUT_DIR, log_lik_file(rdt_se))
  saveRDS(list(by_obs = lik, by_name1_grid_point = name1_grid_log_lik,
               best_by_name1 = best_by_name1,
               deletion_info_by_iso = deletion_info_by_iso,
               likelihood = LIKELIHOOD, grid_setup = GRID_SETUP, results_file = results_file,
               rdt_sensitivity = rdt_se, s_u_ratio = S_U_RATIO, results_path = results_path),
          out_path)
  cat(sprintf("Saved log-likelihood for %d iso x name_1 x grid-point combos to %s\n",
              nrow(name1_grid_log_lik), out_path))
  cat(sprintf("Best-point half-life: median %.1fd, %d of %d sites at the %gd grid floor.\n",
              median(best_by_name1$hrp2_scale * log(2)),
              sum(best_by_name1$hrp2_scale * log(2) < min(GRID_HALFLIFE_RANGE) + 0.01), nrow(best_by_name1),
              min(GRID_HALFLIFE_RANGE)))
  }   # end RDT sensitivity loop

  cat("\n--- HRP2 deletion_pct used per country (fixed at the file's central value, not optimised) ---\n")
  print(deletion_info_by_iso)
}
