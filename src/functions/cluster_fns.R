# =============================================================================
# src/functions/cluster_fns.R
#
# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: simulation, likelihood, RDT sensitivity settings)
#
# Simulation functions run on the DIDE cluster (hipercow environment
# "calibration") and shared settings read by steps 3a-4:
#   - build_adm1_base_pars(): per-site malariasimulation parameters from the
#     site files (interventions, ITN conversion, seasonality), cached once per site
#   - run_hrp2_grid_adm1() / run_simulation_joint(): run a set of (HRP2 scale,
#     EIR) points for one site; returns annual RDT (HRP2), LM and joint
#     RDT x LM prevalence in 6-59-month-olds
#   - grid settings (GRID_SETUP_TAG, GRID_HALFLIFE_RANGE, grid_eir_values(),
#     GRID_SITE_EIR_REFINEMENT) and RDT sensitivity settings and helpers
#     (RDT_SENSITIVITY_MAIN/OPTIONS, rdt_prev_observed(), joint_cells_observed(),
#     log_lik_file(), rdt_sensitivity_dir())
# Must stay free of side effects beyond library()/source() and definitions:
# every cluster task re-sources it. Functions used only by archived analyses
# are in archive/src/functions/cluster_fns_removed.R.
# =============================================================================

library(data.table)
library(malariasimulation)
source("./src/functions/paths.R")   # SITE_CACHE_DIR, EIR_SITE_CACHE_DIR, DATA_DIR

# Shared RDT/LM series colours -- dataviz skill categorical slots 1/2
# (blue/orange). Previously redefined independently in 04_grid_simulation.R,
# 05_emulation.R, 06_hierarchical_mcmc.R (three times), and
# 07_visualisation.R; centralised here since cluster_fns.R is already
# sourced by all of them.
COLOR_HRP2_RDT <- "#2a78d6"
COLOR_LM       <- "#eb6834"

# "<iso>_<name_1>" is this pipeline's site key from 04_grid_simulation.R
# onward (build_emulator_inputs_from_grid_search()'s naming convention) --
# used as a plain sprintf() at ~19 call sites across 04/05/06/
# check_site_coverage.R. Added here for new call sites; existing ones were
# left as-is (a one-line sprintf() vs a one-line function call has no real
# duplication cost, so retrofitting every existing site wasn't worth the
# re-verification risk on already-tested scripts).
make_site_name <- function(iso, name_1) sprintf("%s_%s", iso, name_1)

# Standard ggplot2 theme tokens repeated (as an inline theme() block, not a
# shared call) across most plots in 04_grid_simulation.R, 05_emulation.R,
# 06_hierarchical_mcmc.R, and 07_visualisation.R -- muted grid, dark bold
# title, muted subtitle/axis text. Available for new plots via
# `+ project_theme()`; existing plots' own inline theme() calls were left
# as-is (same reasoning as make_site_name() above -- retrofitting ~10
# already-tested plots for a cosmetic-only change wasn't worth the
# re-verification cost, but new plots should use this instead of pasting
# the block again).
project_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.3),
      plot.title       = ggplot2::element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "#52514e", size = 9),
      axis.text        = ggplot2::element_text(colour = "#898781"),
      axis.title       = ggplot2::element_text(colour = "#52514e"),
      legend.position  = "bottom"
    )
}

# Convert an HRP2 half-life (days) to the "scale" parameter malariasimulation
# actually takes (hrp2_scale in set_hrp2_parameters() — one shared clearance
# rate for treated and untreated infections alike as of malariasimulation
# 3.0.0; there is no separate treated-specific rate any more).
# Clearance there is a constant-hazard exponential process with p_clear =
# 1 - exp(-1/scale) per day, i.e. scale IS the distribution's MEAN duration
# (confirmed from the fork's docs/source), not its half-life — for an
# exponential with mean `scale`, half_life = scale * log(2). Grids are more
# intuitive to reason about in half-life terms, so build them there and
# convert at the point they're fed to set_hrp2_parameters()/run_hrp2_grid_adm1().
halflife_to_scale <- function(halflife_days) {
  halflife_days / log(2)
}

# Inverse of halflife_to_scale() — for display: results/likelihood tables
# carry hrp2_scale in scale units (what the simulation actually used), but
# half-life reads more intuitively on a plot.
scale_to_halflife <- function(scale) {
  scale * log(2)
}


# Central pfhrp2/3 deletion prevalence (%) for one iso -- the file's mean if
# a row exists, else fallback_mean_pct (assume no deletion where none is
# reported, rather than guessing a nonzero rate). Companion to
# build_hrp2_deletion_grid() above (same fallback rule) but returns just the
# point estimate, not a whole grid: deletion is no longer treated as a free/
# optimised parameter in 15_launch_grid_persistence.R Section F or fit as a
# GP dimension in 16_fit_emulators.R -- fitting it as a 4th emulator
# dimension consistently produced degenerate covariance hyperparameters
# (every theta pinned at DiceKriging's default bound, not a genuine interior
# optimum), so it's fixed at this file-reported central value instead.
get_hrp2_deletion_mean <- function(iso, deletion_dt, fallback_mean_pct = 0) {
  row <- deletion_dt[ISO == iso]
  if (nrow(row) == 0L) {
    list(mean_pct = fallback_mean_pct, source = "fallback")
  } else {
    list(mean_pct = row$mean[1], source = "file")
  }
}


# =============================================================================
# EIR emulator: build once from a grid of run_one_eir() results; call
# interpolate_eir() in place of run_one_eir() during MCMC.
# =============================================================================


# =============================================================================
# Worker function for district-level 2000-2024 simulations (script 09)
#
# Runs one malariasimulation for a given district × urban_rural combo and
# returns state counts averaged over each MIS survey year present for that combo.
# params_built must already have the full time-varying treatment schedule applied
# (built locally via build_params in script 09) but EIR not yet set —
# set_equilibrium is called here on the cluster.
# =============================================================================


# Primary site cache: no urban/rural split (Z:/Maggie/2026_08_07_Maggie/,
# <ISO>_site.rds). Used for everything except the EIR baseline.
# SITE_CACHE_DIR is set in paths.R.

# Old site cache: STILL urban/rural split, and the only source with real
# (non-NA) baseline EIR values ($eir$eir is entirely NA in the new cache).
# Used only by get_adm1_eir_weighted() below, and only covers 16 isos —
# adm1 combos from SITE_CACHE_DIR outside that set have no EIR baseline to
# build a relative grid from.
# EIR_SITE_CACHE_DIR is set in paths.R.

# Current grid simulation setup, read by 04-10 so they all agree.
# GRID_SETUP_TAG names the base_pars cache (base_pars_cache_<tag>/, column
# base_pars_path_<tag> of adm1_combos_with_paths.RDS), the grid bundle and
# the grid results file. GRID_*_RANGE are the factorial grid's bounds
# (04 Section B), i.e. the range the emulators are trained on.
GRID_SETUP_TAG      <- "mapitn_burnin5"   # MAP ITN net loss + 5-year burn-in
GRID_HALFLIFE_RANGE <- c(5, 365)          # HRP2 half-life, days

# EIR axis: the original 20 log-spaced points from 1 to 1000, extended by
# GRID_EIR_EXTEND further steps of the same log spacing below and above
# (2026-10-08: sites were pinned at EIR 1 and 1000). Built as
# exp(k * step) from the same step as the original seq(), so the original 20
# values are bit-identical and their simulations are reused -- 04 Section D
# only submits the points missing from the results file.
GRID_EIR_BASE_RANGE <- c(1, 1000)
GRID_EIR_BASE_N     <- 20L
GRID_EIR_EXTEND     <- c(below = 2L, above = 2L)
grid_eir_values <- function() {
  step <- (log(GRID_EIR_BASE_RANGE[2]) - log(GRID_EIR_BASE_RANGE[1])) / (GRID_EIR_BASE_N - 1L)
  k    <- seq(-GRID_EIR_EXTEND[["below"]], GRID_EIR_BASE_N - 1L + GRID_EIR_EXTEND[["above"]])
  exp(log(GRID_EIR_BASE_RANGE[1]) + k * step)
}
GRID_EIR_RANGE <- range(grid_eir_values())   # absolute EIR, ~0.48-2069

# Site-specific EIR refinement (2026-10-09): extra EIR values, simulated at
# every grid half-life, for sites whose likelihood peak is narrower in EIR than
# the grid spacing (x1.44 between EIR values) -- the Rwanda sites, where the
# emulator otherwise guesses the peak from one or two EIR columns. Each range is
# the union across RDT sensitivity variants of the site's EIR values within 50
# log-lik units of its best, one grid step either side. Simulations don't depend
# on RDT sensitivity, so one refinement serves every variant. Each site's grid
# stays rectangular (20 half-lives x shared + extra EIRs), so 05's emulator
# inputs need no change.
GRID_SITE_EIR_REFINEMENT <- list(
  "RWA_Amajyaruguru"     = c(0.6952, 4.281),
  "RWA_Amajyepfo"        = c(6.158, 18.33),
  "RWA_Iburasirazuba"    = c(12.74, 37.93),
  "RWA_Iburengerazuba"   = c(1, 6.158),
  "RWA_Umujyi wa Kigali" = c(4.281, 18.33)
)
GRID_REFINE_EIR_STEP <- 1.1   # spacing of the extra EIR values (vs x1.44 on the shared grid)

# Extra EIR values for a site ("<iso>_<name_1>"), skipping any within 3% of a
# shared grid EIR; numeric(0) for sites without refinement.
site_refinement_eir_values <- function(site) {
  r <- GRID_SITE_EIR_REFINEMENT[[site]]
  if (is.null(r)) return(numeric(0))
  v  <- exp(seq(log(r[1]), log(r[2]), by = log(GRID_REFINE_EIR_STEP)))
  es <- grid_eir_values()
  v[vapply(v, function(x) min(abs(log(x) - log(es))) > log(1.03), logical(1))]
}

# Grid results file for a setup: "original" is the 2026-09-11 netz-ITN /
# 20-year burn-in grid (no joint RDT x LM outputs).
grid_results_file <- function(setup = GRID_SETUP_TAG) {
  if (setup == "original") "hrp2_grid_results.RDS" else sprintf("hrp2_grid_results_%s.RDS", setup)
}

# RDT sensitivity among HRP2-expressing infections, separate from deletion:
# an HRP2-positive child reads RDT-positive with probability
# (1 - deletion) x sensitivity, so missed RDTs (lot quality, storage, reading)
# give RDT-/LM+ results without forcing them through the decay parameter.
# 04 Section F computes the likelihood for every value in
# RDT_SENSITIVITY_OPTIONS; RDT_SENSITIVITY_MAIN is the main analysis and keeps
# the unsuffixed file and folder names, the others are sensitivity analyses
# that 05/06 pick up by setting their own RDT_SENSITIVITY.
RDT_SENSITIVITY_MAIN    <- 0.90
RDT_SENSITIVITY_OPTIONS <- c(0.80, 0.90, 1.00)

# Observed RDT prevalence from simulated HRP2 prevalence.
rdt_prev_observed <- function(hrp2_prev, deletion_pct, rdt_sensitivity = RDT_SENSITIVITY_MAIN) {
  hrp2_prev * (1 - deletion_pct / 100) * rdt_sensitivity
}

# Observed joint RDT x LM cells from simulated ones -- the same transform as
# 04 Section F's multinomial: an HRP2-positive child reads RDT-positive with
# probability r = (1 - deletion) x sensitivity, so a fraction (1 - r) of each
# RDT-positive cell moves to the RDT-negative cell with the same LM result.
# Returns a list (q_tp, q_fp, q_fn, q_tn), for use with data.table `:=`.
joint_cells_observed <- function(p_tp, p_fp, p_fn, p_tn, deletion_pct, rdt_sensitivity = RDT_SENSITIVITY_MAIN) {
  r <- (1 - deletion_pct / 100) * rdt_sensitivity
  list(q_tp = r * p_tp, q_fp = r * p_fp,
       q_fn = p_fn + (1 - r) * p_tp, q_tn = p_tn + (1 - r) * p_fp)
}

# Wu et al. (2015) comparator, used by step 4 Section I and src/plots.R.
# Wu et al. 2015 (Nature 528:S86) Table 1, RDT:microscopy, under 15 years:
#   logit(RDT prevalence) = 0.109 + 0.908 x logit(microscopy prevalence).
# Table 1 gives point estimates only (no credible intervals or between-survey
# variance), so predictive uncertainty is built from: (1) our posterior LM
# draws, (2) a survey-level deviation from the line, N(0, wu_sd_c), and
# (3) binomial resampling at the survey's effective size, as for our own
# intervals. wu_sd_c is the root mean square of observed logit(RDT) minus the
# Wu prediction in every OTHER country (leave-one-country-out; joint-tested
# children). This favours the Wu comparator: the residuals already contain
# survey sampling noise, which layer (3) adds again, and the RMS includes any
# systematic offset from the Wu line in our data.
WU_INTERCEPT <- 0.109
WU_SLOPE     <- 0.908
# Which uncertainty the Wu comparator carries:
#   "lm_only"  -- only our posterior uncertainty in LM prevalence, passed
#                 through the Wu line (no deviation from the line, no survey
#                 sampling noise). Current choice.
#   "loco_sd"  -- also the N(0, wu_sd_c) survey-level deviation and binomial
#                 resampling described above.
WU_UNCERTAINTY <- "lm_only"

# RDT prevalence from LM prevalence on the Wu line, plus the survey-level
# deviation when WU_UNCERTAINTY is "loco_sd" (sd: one value per element).
wu_rdt_from_lm <- function(lm, sd = 0, uncertainty = WU_UNCERTAINTY) {
  eta <- WU_INTERCEPT + WU_SLOPE * stats::qlogis(pmin(pmax(lm, 1e-6), 1 - 1e-6))
  if (identical(uncertainty, "loco_sd")) eta <- eta + stats::rnorm(length(lm), 0, sd)
  stats::plogis(eta)
}

# Leave-one-country-out RMS of logit(RDT) residuals around the Wu line, from
# observed joint-tested prevalence (obs: iso, rdt, lm). Survey-years with RDT
# or LM prevalence of 0 or 1 have no logit and are left out.
wu_residual_sd_loco <- function(obs, a = WU_INTERCEPT, b = WU_SLOPE) {
  ok <- obs[rdt > 0 & rdt < 1 & lm > 0 & lm < 1]
  ok[, resid := qlogis(rdt) - (a + b * qlogis(lm))]
  isos <- sort(unique(obs$iso))
  data.table(iso = isos,
             wu_sd = vapply(isos, function(i) sqrt(mean(ok[iso != i, resid]^2)), 0),
             n_resid = vapply(isos, function(i) ok[iso != i, .N], 0L))
}

is_main_rdt_sensitivity <- function(rdt_sensitivity) isTRUE(all.equal(rdt_sensitivity, RDT_SENSITIVITY_MAIN))
rdt_sensitivity_tag <- function(rdt_sensitivity) sprintf("rdt_se_%03d", round(100 * rdt_sensitivity))

# 04 Section F's log-likelihood file for a sensitivity value.
log_lik_file <- function(rdt_sensitivity = RDT_SENSITIVITY_MAIN) {
  if (is_main_rdt_sensitivity(rdt_sensitivity)) "hrp2_grid_log_lik.RDS"
  else sprintf("hrp2_grid_log_lik_%s.RDS", rdt_sensitivity_tag(rdt_sensitivity))
}

# Output folder for a sensitivity value: the main analysis writes to base_dir,
# others to base_dir/sensitivity_rdt_se_<NNN> (same pattern as 06's other
# sensitivity_* folders).
rdt_sensitivity_dir <- function(base_dir, rdt_sensitivity = RDT_SENSITIVITY_MAIN) {
  if (is_main_rdt_sensitivity(rdt_sensitivity)) return(base_dir)
  out <- file.path(base_dir, paste0("sensitivity_", rdt_sensitivity_tag(rdt_sensitivity)))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  out
}

# Point adm1_combos$base_pars_path at a setup's cache. adm1_combos_with_paths.RDS
# keeps the original cache in base_pars_path, so every reader that simulates
# must call this, or it silently runs the old setup.
use_setup_base_pars <- function(adm1_combos, setup = GRID_SETUP_TAG) {
  col <- sprintf("base_pars_path_%s", setup)
  if (!col %in% names(adm1_combos)) {
    stop(sprintf("adm1_combos has no %s column -- build that setup's base_pars cache (04 Section C).", col))
  }
  adm1_combos <- data.table::copy(data.table::as.data.table(adm1_combos))
  adm1_combos[, base_pars_path := get(col)]
  if (anyNA(adm1_combos$base_pars_path)) stop(sprintf("adm1_combos: NA in %s.", col))
  adm1_combos
}

# List the name_1 combinations present for one iso's cached site object.
# Cheap and local — used to enumerate the build-phase task list. No
# urban_rural column in this cache (unlike EIR_SITE_CACHE_DIR).
list_adm1_combos <- function(iso, site_cache_dir = SITE_CACHE_DIR) {
  site <- readRDS(file.path(site_cache_dir, paste0(iso, "_site.rds")))
  combos <- unique(data.table::data.table(site$sites)[, .(name_1)])
  combos[, iso := iso]
  combos[]
}

# Look up the OLD cache's baseline "pf" EIR for (name_1, urban_rural) — used
# only by get_adm1_eir_weighted() below. Falls back to the OTHER urban_rural
# category's eir when this one has none recorded (some adm1s only have a
# usable eir estimate for one of urban/rural). Verified equivalent to
# filtering site::subset_site()'s $eir for the same adm1.
get_site_eir <- function(site_raw, name_1, urban_rural) {
  target_name_1      <- name_1
  target_urban_rural <- urban_rural
  other_urban_rural   <- if (target_urban_rural == "urban") "rural" else "urban"

  eir_dt <- data.table::data.table(site_raw$eir)
  eir_in <- eir_dt[eir_dt$sp == "pf" & eir_dt$name_1 == target_name_1 &
                     eir_dt$urban_rural == target_urban_rural, eir]
  if (length(eir_in) == 0) {
    eir_in <- eir_dt[eir_dt$sp == "pf" & eir_dt$name_1 == target_name_1 &
                       eir_dt$urban_rural == other_urban_rural, eir]
  }
  eir_in
}

# Population-weighted admin1-level baseline EIR for (name_1), combining the
# OLD cache's urban + rural eir estimates (get_site_eir()) since the new,
# primary site cache (SITE_CACHE_DIR) has no urban/rural split and no usable
# eir at all. Weights are each stratum's mean population share over
# 2000-2025 (OLD cache's population$population_total) — a single baseline
# EIR isn't time-varying, so an across-year average weight is used rather
# than picking one arbitrary year. Takes an already-loaded OLD site object
# (see get_adm1_eir_weighted() below for the file-reading wrapper) so callers
# working through many name_1s for the same iso (e.g. Section B) can read
# that iso's file once instead of once per name_1.
get_adm1_eir_weighted_from_site <- function(old_site, iso, name_1) {
  target_name_1 <- name_1

  eir_urban <- get_site_eir(old_site, target_name_1, "urban")
  eir_rural <- get_site_eir(old_site, target_name_1, "rural")
  if (length(eir_urban) == 0 || length(eir_rural) == 0) {
    stop(sprintf("No old-cache EIR found for %s / %s (urban or rural).", iso, target_name_1))
  }

  pop_dt <- data.table::data.table(old_site$population$population_total)
  pop_dt <- pop_dt[pop_dt$name_1 == target_name_1, ]
  pop_by_ur <- pop_dt[, .(mean_pop = mean(pop)), by = urban_rural]

  if (nrow(pop_by_ur) == 1) {
    # Only one of urban/rural has population data for this adm1 (e.g.
    # BDI/Cankuzo) — no weighting possible/needed, use that category's eir.
    only_ur <- pop_by_ur$urban_rural[1]
    return(get_site_eir(old_site, target_name_1, only_ur))
  }
  if (nrow(pop_by_ur) != 2) {
    stop(sprintf("Expected 1 or 2 urban/rural population rows for %s / %s, found %d.",
                 iso, target_name_1, nrow(pop_by_ur)))
  }
  w_urban <- pop_by_ur[urban_rural == "urban", mean_pop] / sum(pop_by_ur$mean_pop)

  w_urban * eir_urban + (1 - w_urban) * eir_rural
}

# File-reading wrapper around get_adm1_eir_weighted_from_site() for one-off
# lookups (e.g. build_adm1_base_pars(), which only needs one iso's file per
# cluster task anyway). Only covers the 16 isos in EIR_SITE_CACHE_DIR; errors
# clearly if the iso or name_1 isn't found there.
get_adm1_eir_weighted <- function(iso, name_1, eir_site_cache_dir = EIR_SITE_CACHE_DIR) {
  old_site_path <- file.path(eir_site_cache_dir, paste0(iso, "_site.RDS"))
  if (!file.exists(old_site_path)) {
    stop(sprintf("No old-cache EIR source for iso '%s' (name_1 '%s') — %s not found.",
                 iso, name_1, old_site_path))
  }
  get_adm1_eir_weighted_from_site(readRDS(old_site_path), iso, name_1)
}


# site::site_parameters()'s intervention-schedule fields (SMC/MDA/TBV/mass
# PEV/EPI PEV campaign days) are meant to be whole-number simulation-day
# indices, but can come back with ~1e-14-scale floating-point noise (e.g.
# 175.99999999999997 instead of exactly 176) -- confirmed directly against
# real cached base_pars objects, and traced to a genuine crash: malaria-
# simulation's own internal event scheduler rounds a *_timesteps value to an
# integer tick when it actually SCHEDULES/fires that event, but several of
# its listeners (create_mda_listeners() -- covers both MDA and SMC, since
# attach_event_listeners() dispatches both through it -- and, by the same
# pattern, create_tbv_listener()/create_mass_pev_listener()/
# attach_pev_dose_listeners()) then re-derive "which scheduled round is
# this" via an EXACT `which(timesteps == timestep)` equality check against
# the ORIGINAL, un-rounded value. That check fails on the noisy entries,
# `which()` returns integer(0), and code of the shape
# `if (time_index == 0) { ... }` (create_mda_listeners()'s own pattern)
# then throws "argument is of length zero" -- confirmed by direct
# reproduction: 20 real sites all showed this exact crash AND non-integer
# smc_timesteps, while >10 other real SMC-enabled sites with clean integer
# smc_timesteps all simulated successfully; rounding smc_timesteps and
# re-running the SAME site end-to-end (full 25-year horizon) resolved it
# with no other error. The likely upstream source is site::site_parameters()
# (or the site cache's own campaign-date arithmetic) introducing this noise
# during date-to-timestep conversion for specific dates -- not this
# project's own code, and not something to patch inside the installed
# package.
#
# Applied to ALL of the *_timesteps fields malariasimulation's own
# attach_event_listeners() dispatches through this same "exact match against
# a stored timestep vector" pattern, not just smc_timesteps -- the failure
# was only OBSERVED on smc_timesteps in this project's data, but mda/tbv/PEV
# campaigns are built by the same site::site_parameters() code path and
# would be silently vulnerable to the identical bug the moment a site
# actually uses one of them. Rounding an already near-integer numeric
# vector to the nearest integer changes nothing else about its meaning --
# every one of these fields is a day-index by construction.
round_intervention_timesteps <- function(pars) {
  ts_fields <- c(
    "mda_timesteps", "smc_timesteps", "tbv_timesteps",
    "mass_pev_timesteps", "pev_epi_timesteps"
  )
  for (f in ts_fields) {
    if (!is.null(pars[[f]]) && length(pars[[f]]) > 0) {
      pars[[f]] <- round(pars[[f]])
    }
  }
  pars
}

# Build and cache the malariasimulation parameters for one (iso, name_1)
# combination. Everything here — site subsetting, the ITN input
# distribution, and site::site_parameters() — is independent of hrp2_scale
# AND of the EIR used to initialise equilibrium, so it is done ONCE per adm1
# combo and cached to disk, rather than being rebuilt for every point of the
# grid (run_hrp2_grid_adm1() below re-uses the cached result across the full
# grid for that adm1).
#
# set_equilibrium() is NOT called here (unlike the earlier version of this
# function): EIR is now a free grid parameter alongside hrp2_scale (see
# run_hrp2_grid_adm1()), and set_equilibrium()'s
# init_EIR determines the cached object's initial immunity/age structure, so
# it can't be fixed until EIR is known. The `eir` used in site::site_parameters()
# below is a separate, fixed site property (the seasonal transmission/vector
# calibration baseline, population-weighted from the OLD cache via
# get_adm1_eir_weighted() since SITE_CACHE_DIR has none) — unrelated to the
# init_EIR grid.
# ITN usage -> model net distribution, EXACTLY as the site-file EIR
# calibration builds it (Z:/Pete/foresite-orderly/src/calibration/
# calibration_utils.R::calibrate_site()): site::site_usage_to_model_distribution()
# with MAP's net-loss curve (netz::net_loss_map) and the site file's
# retention_half_life passed as a HALF-life. Replaces the earlier
# netz::usage_to_model_distribution() call, which used netz's default
# exponential loss (net_loss_exp) and passed the half-life as a MEAN retention
# time -- a different net-retention curve, hence different modelled net use
# for the same target usage, so the site files' calibrated EIRs didn't carry
# over (08_description_of_data.R sections 12-15).
site_itn_input_dist <- function(site_sub) {
  itn <- site_sub$interventions$itn
  site::site_usage_to_model_distribution(
    usage                    = itn$use$itn_use,
    usage_year               = itn$use$year,
    usage_day_of_year        = itn$use$usage_day_of_year,
    distribution_year        = itn$implementation$year,
    distribution_day_of_year = itn$implementation$distribution_day_of_year,
    distribution_lower       = itn$implementation$distribution_lower,
    distribution_upper       = itn$implementation$distribution_upper,
    net_loss_function        = netz::net_loss_map,
    half_life                = itn$retention_half_life
  )
}

# burn_in_years: 5, matching the calibration's own EIR search
# (calibrate_site()'s calibration_burnin = 5, start_year = 2000 - 5) -- the
# EIRs in the site files were fitted under a 5-year burn-in. Stored on the
# cached object as base_pars$burn_in_years so run_hrp2_grid_adm1() and friends
# map timesteps -> calendar years correctly whichever cache they're given
# (caches built before this field existed used 20 years).
build_adm1_base_pars <- function(iso, name_1, base_pars_cache_dir,
                                  site_cache_dir = SITE_CACHE_DIR,
                                  sim_start = 2000L, sim_end = 2025L,
                                  human_pop = 50000L,
                                  burn_in_years = 5L,
                                  age_min = 183L, age_max = 1825L) {
  site_raw <- readRDS(file.path(site_cache_dir, paste0(iso, "_site.rds")))

  # Local copy with a non-colliding name: this installed data.table version's
  # "..var" NSE prefix doesn't reliably resolve function-argument names that
  # match column names (name_1), so avoid relying on it.
  target_name_1 <- name_1
  sites_dt <- data.table::data.table(site_raw$sites)
  site_filter <- data.frame(sites_dt[sites_dt$name_1 == target_name_1, ])
  site_adm2 <- site::subset_site(site = site_raw, site_filter = site_filter)
  site_adm2$interventions$treatment$implementation$day_of_year <- 1

  site_adm2$interventions$itn$implementation$itn_input_dist <- site_itn_input_dist(site_adm2)

  eir_in <- get_adm1_eir_weighted(iso, target_name_1)
  # burn_in_years extra years are prepended by shifting start_year back --
  # same convention as run_baseline_lm_2_10_urban_rural()'s own burn-in (see
  # that function's header for the rationale, matching Z:/Pete/foresite-
  # orderly/src/calibration/calibration_utils.R::calibrate_site()'s own
  # diagnostic run). run_hrp2_grid_adm1() below must be told the SAME
  # burn_in_years used here (it only reads base_pars_path, not this
  # function's own args) to correctly map timestep -> calendar year and trim
  # the burn-in back off -- both default to 20L so this stays consistent by
  # default, but a caller overriding one MUST override the other identically.
  base_pars <- site::site_parameters(
    interventions = site_adm2$interventions,
    demography    = site_adm2$demography,
    vectors       = site_adm2$vectors,
    seasonality   = site_adm2$seasonality,
    eir           = eir_in,
    overrides     = list(human_population = human_pop),
    start_year    = sim_start - burn_in_years,
    end_year      = sim_end
  )
  base_pars <- round_intervention_timesteps(base_pars)
  # Two simultaneous prevalence age bands -- the pipeline's own established
  # 6m-5y (183-1825d) band PLUS a 2-10y (730-3650d) one, rendered in the SAME
  # run rather than a separate resimulation (formerly only ever computed in
  # 05_emulation.R's own dedicated urban/rural EIR-baseline function -- see
  # that function's header for why the URBAN/RURAL-split version still lives
  # there: this combined, non-split grid has no urban/rural distinction to
  # draw on at all). malariasimulation renders one n_detect_lm_<tag> column
  # per (min,max) pair when these are given as parallel vectors -- see
  # to_annual_hrp2_lm()'s own lm_age_tag_2 handling below.
  base_pars$prevalence_rendering_min_ages  <- c(age_min, 730L)
  base_pars$prevalence_rendering_max_ages  <- c(age_max, 3650L)
  base_pars$state_count_rendering_min_ages <- age_min
  base_pars$state_count_rendering_max_ages <- age_max

  # HRP2 age band (0-5y) and the fixed (non-grid) HRP2 parameters are set
  # once here; hrp2_scale is varied later per grid point in
  # run_hrp2_grid_adm1(). set_epi_outputs()'s age bands are unioned into
  # the n_age_<band> renderer automatically, so n_age_0_1824 is available
  # alongside n_hrp2_positive_0_1824 with no extra rendering setup.
  base_pars <- malariasimulation::set_epi_outputs(
    base_pars,
    hrp2 = c(0, 5 * 365, 100 * 365)
  )
  base_pars <- malariasimulation::set_hrp2_parameters(
    base_pars,
    ##assume that asymptomatic infections are as likely to be
    ##dxd by rdt as lm
    hrp2_asymptomatic_prob = 1
  )

  base_pars$burn_in_years <- burn_in_years   # read back by run_hrp2_grid_adm1() etc.

  dir.create(base_pars_cache_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(base_pars_cache_dir, sprintf("%s_%s.RDS", iso, name_1))
  saveRDS(base_pars, out_path)
  out_path
}

# Burn-in a cached base_pars object was built with: its own burn_in_years
# field (build_adm1_base_pars(), current), else 20 -- every cache built before
# that field existed used build_adm1_base_pars()'s then-default of 20 years.
cached_burn_in_years <- function(base_pars) {
  if (is.null(base_pars$burn_in_years)) 20L else as.integer(base_pars$burn_in_years)
}

# Annual mean HRP2 and light-microscopy prevalence for their respective
# rendered age bands (these differ: hrp2_age_tag is the 0-5y HRP2 band set in
# build_adm1_base_pars() via set_epi_outputs(hrp2 = ...); lm_age_tag is the
# prevalence_rendering_min_ages/max_ages band, e.g. "183_1825" for 6-59mo).
# lm_age_tag_2: a SECOND, optional LM age band -- build_adm1_base_pars() now
# renders 2-10y (730-3650d, tag "730_3650") alongside the primary lm_age_tag
# band in every grid simulation (see that function's own header), so this is
# computed here too, as lm_prev_2_10, whenever those columns are present.
# Guarded (not required) so this still works unchanged against OLDER cached
# base_pars/sim objects that only ever rendered the one band.
# by_month = TRUE: the same means by calendar year AND month (adds a month
# column), for the survey-month matching pilot
# (archived: archive/src/3a_calibration_grid/sensitivity/survey_month_pilot.R); the default is the annual summary.
to_annual_hrp2_lm <- function(sim, hrp2_age_tag, lm_age_tag, lm_age_tag_2 = "730_3650",
                               baseline_year = 2000L, by_month = FALSE) {
  hrp2_pos_col <- paste0("n_hrp2_positive_", hrp2_age_tag)
  hrp2_age_col <- paste0("n_age_",          hrp2_age_tag)
  lm_col       <- paste0("n_detect_lm_",    lm_age_tag)
  lm_age_col   <- paste0("n_age_",          lm_age_tag)
  lm_col_2     <- paste0("n_detect_lm_",    lm_age_tag_2)
  lm_age_col_2 <- paste0("n_age_",          lm_age_tag_2)
  origin       <- as.Date(paste0(baseline_year, "-01-01"))
  dt           <- data.table::as.data.table(sim)
  dt[, year := as.integer(format(origin + timestep - 1L, "%Y"))]
  if (by_month) dt[, month := as.integer(format(origin + timestep - 1L, "%m"))]
  grp <- if (by_month) c("year", "month") else "year"
  out <- dt[, .(
    hrp2_prev = mean(get(hrp2_pos_col) / get(hrp2_age_col), na.rm = TRUE),
    lm_prev   = mean(get(lm_col)       / get(lm_age_col),   na.rm = TRUE)
  ), by = grp]
  if (all(c(lm_col_2, lm_age_col_2) %in% names(dt))) {
    out_2 <- dt[, .(lm_prev_2_10 = mean(get(lm_col_2) / get(lm_age_col_2), na.rm = TRUE)), by = grp]
    out <- merge(out, out_2, by = grp)
  }
  # Joint RDT x LM cell probabilities (run_simulation_joint() only), annual
  # mean of the daily expected proportions in the joint band (6-59mo):
  # p_tp = RDT+/LM+, p_fp = RDT+/LM-, p_fn = RDT-/LM+, p_tn = RDT-/LM-.
  jt <- function(x) paste0("joint_", x, "_", lm_age_tag)
  if (all(vapply(c("n", "rdt", "lm", "tp"), jt, character(1)) %in% names(dt))) {
    out_j <- dt[, .(
      p_tp = mean(get(jt("tp")) / get(jt("n")), na.rm = TRUE),
      p_fp = mean((get(jt("rdt")) - get(jt("tp"))) / get(jt("n")), na.rm = TRUE),
      p_fn = mean((get(jt("lm"))  - get(jt("tp"))) / get(jt("n")), na.rm = TRUE)
    ), by = grp]
    out_j[, p_tn := 1 - p_tp - p_fp - p_fn]
    out <- merge(out, out_j, by = grp)
  }
  # p_fp_u: the part of p_fp from HRP2-flagged children in U (render_u runs only)
  if (jt("rdt_u") %in% names(dt)) {
    out <- merge(out, dt[, .(p_fp_u = mean(get(jt("rdt_u")) / get(jt("n")), na.rm = TRUE)), by = grp], by = grp)
  }
  out
}

# --- Joint RDT (HRP2) x light-microscopy outcomes ----------------------------
# malariasimulation renders HRP2 positivity and LM detection only as separate
# marginals, but the MIS likelihood (04_grid_simulation.R Section F,
# LIKELIHOOD = "multinomial_deff") needs the joint 2x2 for the SAME children.
# Per individual the fork (mwalte10/malariasimulation@age_spec_states,
# fa6a541) stores HRP2 status as a 0/1 variable (variables$hrp2), and LM
# detection is 1 for D/Tr, probability_of_detection(age, id) for A, 0
# otherwise -- drawn independently of HRP2 given state -- so the EXPECTED
# joint counts are exact sums over individuals:
#   joint_rdt = #(hrp2 = 1)                 joint_lm = #(D/Tr) + sum_A p_lm
#   joint_tp  = #(D/Tr & hrp2) + sum_A p_lm * hrp2
# (fp = rdt - tp, fn = lm - tp, tn = n - tp - fp - fn). Expected rather than
# Bernoulli-drawn counts, so this renderer consumes NO random numbers: a
# run_simulation_joint() run replays exactly the same trajectory (and the
# same hrp2_prev/lm_prev) as run_simulation() under the same seed.
# render_u = TRUE also renders joint_rdt_u: HRP2-flagged children in U. LM
# never detects U, so these are exactly the U share of the false positives
# (the s_U sensitivity analysis, 3a sensitivity/su_grid.R). Off by default,
# leaving the main grid's output unchanged.
create_joint_rdt_lm_renderer <- function(renderer, variables, parameters,
                                         min_age = 183L, max_age = 1825L, render_u = FALSE) {
  get_age                  <- malariasimulation:::get_age
  probability_of_detection <- malariasimulation:::probability_of_detection
  in_age_range             <- malariasimulation:::in_age_range
  tag <- paste0(min_age, "_", max_age)
  function(timestep) {
    in_age <- in_age_range(variables$birth, timestep, min_age, max_age)
    hrp2   <- variables$hrp2$get_index_of(1)$and(in_age)
    clin   <- variables$state$get_index_of(c("Tr", "D"))$and(in_age)
    asym   <- variables$state$get_index_of("A")$and(in_age)
    p_lm   <- probability_of_detection(get_age(variables$birth$get_values(asym), timestep),
                                       variables$id$get_values(asym), parameters)
    h_asym <- variables$hrp2$get_values(asym)
    renderer$render(paste0("joint_n_",   tag), in_age$size(), timestep)
    renderer$render(paste0("joint_rdt_", tag), hrp2$size(), timestep)
    renderer$render(paste0("joint_lm_",  tag), clin$size() + sum(p_lm), timestep)
    renderer$render(paste0("joint_tp_",  tag), clin$copy()$and(hrp2)$size() + sum(p_lm * h_asym), timestep)
    if (render_u) {
      renderer$render(paste0("joint_rdt_u_", tag), variables$state$get_index_of("U")$and(hrp2)$size(), timestep)
    }
  }
}

# run_simulation_joint(): malariasimulation::run_simulation() plus the joint
# renderer above. A copy of malariasimulation:::run_resumable_simulation()
# (fork commit fa6a541; no initial_state/resume branch) evaluated in the
# package namespace so its internals resolve, with one extra process appended
# to create_processes()'s list. All processes see the same pre-update state
# within a timestep (individual applies queued updates afterwards), so the
# joint counts are consistent with the built-in renderers. Re-check this
# against run_resumable_simulation() if the fork is ever updated.
run_simulation_joint <- local({
  f <- function(timesteps, parameters, make_joint_renderer,
                joint_min_age = 183L, joint_max_age = 1825L) {
    random_seed(ceiling(runif(1) * .Machine$integer.max))
    correlations <- get_correlation_parameters(parameters)
    variables <- create_variables(parameters)
    events <- create_events(parameters)
    initialise_events(events, variables, parameters)
    renderer <- individual::Render$new(timesteps)
    populate_incidence_rendering_columns(renderer, parameters)
    attach_event_listeners(events, variables, parameters, correlations, renderer)
    vector_models <- parameterise_mosquito_models(parameters, timesteps)
    solvers <- parameterise_solvers(vector_models, parameters)
    lagged_eir <- create_lagged_eir(variables, solvers, parameters)
    lagged_infectivity <- create_lagged_infectivity(variables, parameters)
    stateful_objects <- list(RandomState$new(FALSE), correlations, vector_models,
                             solvers, lagged_eir, lagged_infectivity)
    processes <- create_processes(renderer, variables, events, parameters, vector_models,
                                  solvers, correlations, lagged_eir, lagged_infectivity, timesteps)
    processes <- c(processes, make_joint_renderer(renderer, variables, parameters,
                                                  joint_min_age, joint_max_age))
    individual::simulation_loop(processes = processes, variables = variables, events = events,
                                timesteps = timesteps, state = NULL, restore_random_state = FALSE)
    renderer$to_dataframe()
  }
  environment(f) <- asNamespace("malariasimulation")
  function(timesteps, parameters, joint_min_age = 183L, joint_max_age = 1825L, render_u = FALSE) {
    f(timesteps, parameters,
      function(renderer, variables, parameters, min_age, max_age)
        create_joint_rdt_lm_renderer(renderer, variables, parameters, min_age, max_age, render_u = render_u),
      joint_min_age, joint_max_age)
  }
})

# Run one chunk of the hrp2_scale x eir grid for a single adm1 combo,
# re-using its cached base_pars (see build_adm1_base_pars) instead of
# rebuilding it per grid point. hrp2_grid must have hrp2_scale AND eir
# columns — EIR is a free/optimisable parameter here (set_equilibrium() is
# called per row, since build_adm1_base_pars() no longer fixes it). Returns
# a compact annual summary (hrp2_prev, lm_prev, and now lm_prev_2_10) rather
# than the raw sim object, since millions of these results get collected
# back from the cluster.
# burn_in_years MUST match whatever build_adm1_base_pars() actually used to
# build base_pars_path (both default 20L) -- base_pars_path's own internal
# calendar already starts sim_start - burn_in_years (see that function's own
# header), so this is only used here to correctly map timestep -> calendar
# year and drop the burn-in years back off before returning, not to change
# base_pars itself.
run_hrp2_grid_adm1 <- function(base_pars_path, hrp2_grid, iso, name_1,
                                sim_start = 2000L, sim_end = 2025L,
                                burn_in_years = NULL,
                                hrp2_age_tag = "0_1824", lm_age_tag = "183_1825",
                                seeds = NULL, human_population = NULL, joint = FALSE,
                                monthly = FALSE, render_u = FALSE) {
  # render_u = TRUE (joint runs only): also return p_fp_u, the false positives
  # from HRP2-flagged children in U (see create_joint_rdt_lm_renderer())
  # joint = TRUE: simulate via run_simulation_joint() and also return the
  # joint RDT x LM cell probabilities p_tp/p_fp/p_fn/p_tn (lm_age_tag band;
  # see to_annual_hrp2_lm()). Same trajectory as joint = FALSE under the same
  # seed -- the joint renderer draws no random numbers.
  # monthly = TRUE (survey-month pilot only): returns list(annual, monthly),
  # the monthly table holding the same means by calendar month. The default
  # returns the annual table alone, as every main-pipeline caller expects.
  base_pars <- readRDS(base_pars_path)
  if (is.null(burn_in_years)) burn_in_years <- cached_burn_in_years(base_pars)   # NULL = the cache's own burn-in
  # human_population: optional override of the CACHED base_pars' own
  # population (set at build_adm1_base_pars() time, default 50,000) --
  # applied before set_equilibrium() so the equilibrium state itself is
  # computed at the smaller population, not just downstream rendering.
  # Rebuilding base_pars_cache at a different population is the "real" way
  # to change this permanently; this override exists for a caller (e.g. a
  # cheaper posterior-predictive check) that wants a smaller population for
  # ONE run without touching the shared cache every other caller reuses.
  # NULL (default) leaves base_pars completely unchanged -- every existing
  # caller keeps whatever population its own base_pars_path was built with.
  if (!is.null(human_population)) base_pars$human_population <- human_population
  n_days    <- 365L * (sim_end - sim_start + burn_in_years)

  # seeds: optional per-row RNG seed vector (length nrow(hrp2_grid)) -- lets
  # a caller vary malariasimulation's own intrinsic stochastic noise across
  # grid rows, rather than every row replaying the IDENTICAL random stream.
  # Defaults to the fixed 123 for every row -- UNCHANGED behaviour for 04's
  # own grid search (its own callers never pass this), which deliberately
  # wants every grid point evaluated under the same simulated "weather" so
  # the likelihood surface reflects genuine parameter differences, not
  # confounded by between-point simulation noise. 06_hierarchical_mcmc.R's
  # own posterior-predictive tasks (Sections F/L) DO pass a real per-draw
  # seed vector -- a posterior predictive interval should reflect both
  # parameter uncertainty AND the simulator's own aleatoric noise, and
  # replaying one fixed seed across every draw was silently dropping the
  # second of those (confirmed directly: every draw's simulated trajectory
  # previously differed ONLY by its posterior parameter values).
  if (is.null(seeds)) seeds <- rep(123L, nrow(hrp2_grid))
  stopifnot(length(seeds) == nrow(hrp2_grid))

  results <- vector("list", nrow(hrp2_grid))
  results_monthly <- vector("list", nrow(hrp2_grid))
  for (i in seq_len(nrow(hrp2_grid))) {
    pars_i <- malariasimulation::set_equilibrium(base_pars, init_EIR = hrp2_grid$eir[i])
    pars_i <- malariasimulation::set_hrp2_parameters(
      pars_i,
      hrp2_scale = hrp2_grid$hrp2_scale[i]
    )
    set.seed(seeds[i])
    sim_i  <- if (joint) {
      run_simulation_joint(timesteps = n_days, parameters = pars_i,
                           joint_min_age = as.integer(sub("_.*", "", lm_age_tag)),
                           joint_max_age = as.integer(sub(".*_", "", lm_age_tag)),
                           render_u = render_u)
    } else {
      malariasimulation::run_simulation(timesteps = n_days, parameters = pars_i)
    }
    annual <- to_annual_hrp2_lm(sim_i, hrp2_age_tag = hrp2_age_tag, lm_age_tag = lm_age_tag,
                                 baseline_year = sim_start - burn_in_years)
    annual <- annual[year >= sim_start]
    annual[, `:=`(
      iso        = iso,
      name_1     = name_1,
      eir        = hrp2_grid$eir[i],
      hrp2_scale = hrp2_grid$hrp2_scale[i]
    )]
    results[[i]] <- annual
    if (monthly) {
      mon <- to_annual_hrp2_lm(sim_i, hrp2_age_tag = hrp2_age_tag, lm_age_tag = lm_age_tag,
                               baseline_year = sim_start - burn_in_years, by_month = TRUE)
      mon <- mon[year >= sim_start]
      mon[, `:=`(iso = iso, name_1 = name_1, eir = hrp2_grid$eir[i], hrp2_scale = hrp2_grid$hrp2_scale[i])]
      results_monthly[[i]] <- mon
    }
  }
  if (monthly) return(list(annual = data.table::rbindlist(results),
                           monthly = data.table::rbindlist(results_monthly)))
  data.table::rbindlist(results)
}


