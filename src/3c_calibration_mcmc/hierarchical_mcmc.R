# =============================================================================
# src/3c_calibration_mcmc/hierarchical_mcmc.R
#   (was src/06_hierarchical_mcmc.R Sections A-E, and 09i/09j for the prior-SD fits)
#
# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (stage 3: Bayesian hierarchical model)
#
# Step 3c, calibration stage 3: hierarchical Bayesian MCMC across all sites.
# Blocked Metropolis-Hastings (hierarchical_mcmc_fns.R::run_mcmc()) with each
# site's likelihood given by its GP emulator from step 3b.
#
# Model (priors from pipeline_fns.R::build_mcmc_priors()):
#   log d_i   ~ Normal(mu_d, sigma_d), truncated to the emulators' trained range
#   mu_d      ~ Normal(MU_D_PRIOR_MEAN, 0.75), or its logit-normal match on the
#               emulators' log d range (PRIOR_FAMILY); main MU_D_PRIOR_MEAN = log(15 / log 2)
#   sigma_d   ~ Normal(0.3, SIGMA_D_PRIOR_SD), truncated at 0
#   log EIR_i ~ Normal(eir_mu, eir_sd), truncated to the grid, or its logit-normal
#               match on the grid (PRIOR_FAMILY); no pooling
# PRIOR_FAMILY: main "normal"; "logit_normal" (no truncation) is planned to
# become the main analysis -- see pipeline_fns.R "Prior family".
# d is the mean HRP2 persistence time (half-life = d x log 2). Partial pooling
# acts on decay only: EIR is site-specific.
#
# Sections
#   A    hipercow environment (once per hipercow root)
#   B/C  submit the prep task (build inputs + tune proposal SDs; blocks, ~2 h,
#        unless WAIT_FOR_PREP = FALSE), then N_CHAINS chain tasks (~8.5 h). Resumable: an existing
#        task_prep.RDS / bundle_chains.RDS in the output folder is reused, so a
#        new fit must start from a folder without them.
#   D    collect chains -> chains.RDS (every chain must succeed)
#   E    R-hat, pool chains -> mcmc_combined.RDS (summaries and shareable
#        draws are written by src/export_posteriors.R)
#
# One fit per run, chosen by RDT_SENSITIVITY, SIGMA_D_PRIOR_SD, MU_D_PRIOR_MEAN
# and PRIOR_FAMILY. The main analysis (0.9, 0.01, NULL, "normal") writes to DIR_MCMC.
# Sensitivity fits:
#   sensitivity/rdt_sensitivity_mcmc.R   RDT sensitivity 0.8, 1.0 -> DIR_MCMC/sensitivity_rdt_se_NNN/
#   sensitivity/sigma_prior_sd_mcmc.R    sigma_d prior SD 0.05, 0.1 (RDT 0.9)
#                                        -> DIR_MCMC/sensitivity_sigma_d_sd_<sd>_rdt_se_090/
# The mu_d prior-centre analysis (not in the manuscript) is archived:
#   archive/src/3c_calibration_mcmc/sensitivity/mu_d_prior_mcmc.R
# Posterior draws and summaries for every fit are exported by src/export_posteriors.R.
#
# Plots are drawn in src/plots.R; the posterior predictive check is step 4.
# =============================================================================
library(data.table)
library(hipercow)

source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/hierarchical_mcmc_fns.R")
source("./src/functions/pipeline_fns.R")

# Settings ======================================================================
OVR <- take_overrides()
RUN_MCMC     <- setting(OVR, "RUN_MCMC", FALSE)       # Section B/C
COLLECT_MCMC <- setting(OVR, "COLLECT_MCMC", FALSE)   # Section D
POOL_MCMC    <- setting(OVR, "POOL_MCMC", FALSE)      # Section E
RDT_SENSITIVITY  <- setting(OVR, "RDT_SENSITIVITY", RDT_SENSITIVITY_MAIN)
SIGMA_D_PRIOR_SD <- setting(OVR, "SIGMA_D_PRIOR_SD", SIGMA_D_PRIOR_SD_MAIN)
MU_D_PRIOR_MEAN  <- setting(OVR, "MU_D_PRIOR_MEAN", NULL)   # log(d); NULL = main, log(15 / log 2)
PRIOR_FAMILY     <- setting(OVR, "PRIOR_FAMILY", PRIOR_FAMILY_MAIN)   # mu_d and EIR priors: "normal" or "logit_normal"
WAIT_FOR_PREP    <- setting(OVR, "WAIT_FOR_PREP", TRUE)   # FALSE: submit the prep task and return (Section B/C)

# Names this fit's cluster bundles: observation-likelihood version (3a's
# LIKELIHOOD) plus the variant tag.
MCMC_RUN_TAG <- "v3_multinomial_deff"

N_CHAINS        <- 4L
CHAIN_N_ITER    <- 10000L
CHAIN_N_WARMUP  <- 5000L
CHAIN_THIN      <- 1L
CHAIN_BASE_SEED <- 123L
MAX_RETRIES     <- 3L

# Upper bound on site half-life (days); NULL = the emulators' full trained range.
MAX_DECAY_HALFLIFE_DAYS <- NULL

PRIORS <- if (is.null(MU_D_PRIOR_MEAN)) build_mcmc_priors(SIGMA_D_PRIOR_SD, prior_family = PRIOR_FAMILY) else
  build_mcmc_priors(SIGMA_D_PRIOR_SD, mu_d_mean = MU_D_PRIOR_MEAN, prior_family = PRIOR_FAMILY)
# GRID_OUT_DIR / EMULATOR_OUT_DIR / MCMC_OUT_DIR can be overridden by the
# model-structure sensitivity analyses (s_U, replicate runs), which read their
# own likelihood and emulators and write to their own folder; BUNDLE_SUFFIX
# keeps their cluster bundle names distinct. Defaults = this variant's folders.
GRID_OUT_DIR     <- setting(OVR, "GRID_OUT_DIR", DIR_GRID)
EMULATOR_OUT_DIR <- setting(OVR, "EMULATOR_OUT_DIR", rdt_sensitivity_dir(DIR_EMULATOR, RDT_SENSITIVITY))
MCMC_OUT_DIR     <- setting(OVR, "MCMC_OUT_DIR", mcmc_variant_dir(RDT_SENSITIVITY, SIGMA_D_PRIOR_SD,
                                                                  mu_d_mean = PRIORS$mu_d_mean, prior_family = PRIOR_FAMILY))
BUNDLE_TAG       <- sprintf("%s_%s%s", MCMC_RUN_TAG,
                            mcmc_variant_tag(RDT_SENSITIVITY, SIGMA_D_PRIOR_SD, PRIORS$mu_d_mean, PRIOR_FAMILY),
                            setting(OVR, "BUNDLE_SUFFIX", ""))
dir.create(MCMC_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("\n########## 3c: RDT sensitivity %.2f, sigma_d prior SD %s, priors %s, mu_d median %.3f (half-life %.1f days) -> %s ##########\n",
            RDT_SENSITIVITY, format(SIGMA_D_PRIOR_SD), PRIOR_FAMILY, PRIORS$mu_d_mean,
            exp(PRIORS$mu_d_mean) * log(2), MCMC_OUT_DIR))

# Hipercow setup ===============================================================
# SECTION A -- run once per hipercow root (interactively). Sources the
# function libraries only, never this driver script.
if (FALSE) {
  hipercow_environment_create(
    name    = "hierarchical_mcmc",
    sources = c("src/functions/hierarchical_mcmc_fns.R", "src/functions/emulator_fns.R")
  )
}

# Task wrapper (top level, so hipercow doesn't serialise a local environment)
run_one_chain_task <- function(task_args) do.call(run_one_chain, task_args)

# Run MCMC =====================================================================
# SECTION B/C -- prep task (blocks until done), then N_CHAINS chain tasks
if (RUN_MCMC) {
  prep_path   <- file.path(MCMC_OUT_DIR, "task_prep.RDS")
  chains_path <- file.path(MCMC_OUT_DIR, "bundle_chains.RDS")
  if (file.exists(prep_path)) {
    task_prep <- readRDS(prep_path)
    cat(sprintf("Reusing prep task | id: %s\n", task_prep))
  } else {
    task_prep <- task_create_call(
      build_and_tune_mcmc_inputs,
      args = list(
        grid_out_dir            = GRID_OUT_DIR,
        emulator_out_dir        = EMULATOR_OUT_DIR,
        n_tune                  = 2000L,
        target                  = 0.30,
        seed                    = 123L,
        priors                  = PRIORS,
        max_decay_halflife_days = MAX_DECAY_HALFLIFE_DAYS,
        log_lik_file            = log_lik_file(RDT_SENSITIVITY)
      ),
      environment = "hierarchical_mcmc",
      resources   = hipercow_resources(cores = 1)
    )
    saveRDS(task_prep, prep_path)
    cat(sprintf("Submitted prep task | id: %s\n", task_prep))
  }
  # WAIT_FOR_PREP = FALSE: stop here unless the prep task has already finished,
  # so several fits' prep tasks can be submitted at once; run RUN_MCMC again
  # (it reuses task_prep.RDS) to submit the chains.
  submit_chains <- WAIT_FOR_PREP || task_status(task_prep) %in% c("success", "failure", "cancelled")
  if (!submit_chains) {
    cat(sprintf("Prep task %s is %s -- not waiting (WAIT_FOR_PREP = FALSE); run RUN_MCMC again to submit the chains.\n",
                task_prep, task_status(task_prep)))
  }
}

if (RUN_MCMC && submit_chains) {
  task_wait(task_prep, timeout = Inf)
  status <- task_status(task_prep)
  if (status != "success") {
    stop(sprintf("Prep task did not succeed (status: %s) -- check task_log_show('%s').", status, task_prep))
  }
  mcmc_inputs <- task_result(task_prep)
  mcmc_inputs$rdt_sensitivity <- RDT_SENSITIVITY   # the likelihood variant these emulators came from
  saveRDS(mcmc_inputs, file.path(MCMC_OUT_DIR, "mcmc_inputs.RDS"))
  cat(sprintf("Prep done: %d sites, tuned mu_d proposal_sd = %.3g, sigma_d proposal_sd = %.3g\n",
              length(mcmc_inputs$site_names), mcmc_inputs$proposal_sds["mu_d"], mcmc_inputs$proposal_sds["sigma_d"]))

  chain_tasks <- lapply(seq_len(N_CHAINS), function(k) list(
    emulator_out_dir = EMULATOR_OUT_DIR,
    log_eir_bounds   = mcmc_inputs$log_eir_bounds,
    log_d_bounds     = mcmc_inputs$log_d_bounds,
    priors           = mcmc_inputs$priors,
    init_theta       = mcmc_inputs$init_theta,
    proposal_sds     = mcmc_inputs$proposal_sds,
    n_iter           = CHAIN_N_ITER,
    n_warmup         = CHAIN_N_WARMUP,
    thin             = CHAIN_THIN,
    seed             = CHAIN_BASE_SEED + k
  ))

  if (file.exists(chains_path)) {
    bundle_chains <- readRDS(chains_path)
    cat(sprintf("Reusing chain bundle | id: %s\n", bundle_chains$name))
  } else {
    bundle_chains <- task_create_bulk_call(
      run_one_chain_task,
      data        = chain_tasks,
      environment = "hierarchical_mcmc",
      bundle_name = safe_bundle_name(sprintf("hierarchical_mcmc_chains_%s", BUNDLE_TAG)),
      resources   = hipercow_resources(cores = 1)
    )
    saveRDS(bundle_chains, chains_path)
    cat(sprintf("Submitted %d chain tasks | bundle id: %s\n", N_CHAINS, bundle_chains$name))
  }
}

# Collect MCMC =================================================================
# SECTION D -- every chain must succeed (R-hat needs all of them)
if (COLLECT_MCMC) {
  bundle_chains <- readRDS(file.path(MCMC_OUT_DIR, "bundle_chains.RDS"))
  chains <- collect_bundle(bundle_chains, MAX_RETRIES, allow_failures = FALSE, label = "chain")$results
  saveRDS(chains, file.path(MCMC_OUT_DIR, "chains.RDS"))
  cat(sprintf("Collected %d chains.\n", length(chains)))
}

# Pool MCMC posterior ==========================================================
# SECTION E -- R-hat on mu_d, sigma_d and each site's log d_i and log EIR_i;
# investigate values above ~1.1 before using the pooled posterior.
if (POOL_MCMC) {
  if (!exists("chains") || !isTRUE(COLLECT_MCMC)) chains <- readRDS(file.path(MCMC_OUT_DIR, "chains.RDS"))
  site_names <- chains[[1]]$site_names
  n_sites    <- length(site_names)

  rhat_mu_d    <- compute_rhat(lapply(chains, function(ch) ch$samples$mu_d))
  rhat_sigma_d <- compute_rhat(lapply(chains, function(ch) ch$samples$sigma_d))
  rhat_site <- function(prefix) setNames(vapply(seq_len(n_sites), function(i) {
    compute_rhat(lapply(chains, function(ch) ch$samples[[paste0(prefix, i)]]))
  }, numeric(1)), site_names)
  rhat_decay <- rhat_site("log_d_")
  rhat_eir   <- rhat_site("log_eir_")

  cat(sprintf("R-hat mu_d: %.3f | sigma_d: %.3f\n", rhat_mu_d, rhat_sigma_d))
  cat(sprintf("R-hat site decay: median %.3f, max %.3f (%d/%d > 1.1)\n",
              median(rhat_decay), max(rhat_decay), sum(rhat_decay > 1.1), n_sites))
  cat(sprintf("R-hat site EIR:   median %.3f, max %.3f (%d/%d > 1.1)\n",
              median(rhat_eir), max(rhat_eir), sum(rhat_eir > 1.1), n_sites))

  pooled_samples <- data.table::rbindlist(lapply(seq_along(chains), function(k) {
    dt <- data.table::copy(chains[[k]]$samples)
    dt[, chain := k]
    dt
  }))

  saveRDS(list(
    chains             = chains,
    pooled_samples     = pooled_samples,
    rhat_mu_d          = rhat_mu_d,
    rhat_sigma_d       = rhat_sigma_d,
    rhat_decay         = rhat_decay,
    rhat_eir           = rhat_eir,
    sigma_d_prior_mean = PRIORS$sigma_d_prior_mean,
    sigma_d_prior_sd   = PRIORS$sigma_d_prior_sd,
    mu_d_mean          = PRIORS$mu_d_mean,
    mu_d_sd            = PRIORS$mu_d_sd,
    prior_family       = PRIOR_FAMILY,
    priors             = PRIORS                   # full prior list, incl. logit-normal parameters
  ), file.path(MCMC_OUT_DIR, "mcmc_combined.RDS"))

  cat(sprintf("Saved mcmc_combined.RDS to %s (posterior summaries: src/export_posteriors.R)\n", MCMC_OUT_DIR))
}
