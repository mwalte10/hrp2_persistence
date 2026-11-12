# =============================================================================
# src/4_sample_posterior/sample_posterior.R      (was src/06_hierarchical_mcmc.R Sections F-I)
#
# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (posterior predictive check)
#
# Step 4. Posterior predictive check: draw N_POSTERIOR_PRED_DRAWS rows from a
# fit's pooled joint posterior (each row gives every site's own d_i and EIR_i),
# run malariasimulation at each site's (d_i, EIR_i) on the cluster, and
# summarise the simulated RDT, LM and joint RDT x LM (TP/FP/FN/TN) prevalence
# per site and year as a mean and 95% posterior predictive interval.
#
# Survey-year intervals include survey sampling noise: each draw's prevalence
# is resampled as a binomial on the survey's effective size (denominator /
# design effect) before taking quantiles; years without a survey use the
# quantiles of the draws. Deletion and RDT sensitivity are applied as in the
# step 3a likelihood, using the RDT sensitivity the fit was run with.
#
# Sections (each runs once per variant in VARIANTS)
#   F/G  sample posterior draws, build one task per (site, chunk of draws), submit
#   H    collect -> posterior_predictive/posterior_predictive_results.RDS
#   I    aggregate -> posterior_predictive/posterior_predictive_summary.RDS
#
# Section I also builds the comparisons for the PPC paragraph: marginal RDT
# (TP + FP) and LM (TP + FN) prevalence from the joint cells, and RDT
# prevalence predicted from the same LM draws with the Wu et al. (2015)
# RDT:microscopy relationship (see "Wu et al. comparator" below).
#
# Runs for the main fit only (RDT sensitivity 0.9, sigma_P prior SD 0.01):
# the manuscript reports posterior predictive coverage for the main fit
# (appendix fig:scatter, \nm{ppc.*} keys). RDT_SENSITIVITY_RUN can add an
# RDT sensitivity fit if needed.
# =============================================================================
library(data.table)
library(hipercow)

source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

# Settings ======================================================================
OVR <- take_overrides()
SUBMIT_POST  <- setting(OVR, "SUBMIT_POST", FALSE)    # Sections F/G: sample and submit
COLLECT_POST <- setting(OVR, "COLLECT_POST", FALSE)   # Section H
AGG_POST     <- setting(OVR, "AGG_POST", FALSE)       # Section I
RDT_SENSITIVITY_RUN <- setting(OVR, "RDT_SENSITIVITY_RUN", RDT_SENSITIVITY_MAIN)

VARIANTS <- data.table::data.table(rdt_sensitivity = RDT_SENSITIVITY_RUN, sigma_d_prior_sd = SIGMA_D_PRIOR_SD_MAIN)

N_POSTERIOR_PRED_DRAWS <- 100L
CHUNK_SIZE             <- 30L
DRAW_SEED              <- 123L
MAX_RETRIES            <- 3L
# Simulated population per run (the grid used 50,000). 10,000 keeps the
# runs affordable; the survey-noise resampling dominates the interval width.
POSTERIOR_PRED_HUMAN_POPULATION <- 10000L

# Wu et al. comparator: WU_INTERCEPT, WU_SLOPE and wu_residual_sd_loco() are in
# src/functions/cluster_fns.R (shared with src/plots.R).
WU_SEED      <- 2015L
# Section I summary prefix -> per-draw column
MARGINAL_COLS <- c(rdtj = "rdt_joint", lmj = "lm_joint", rdtwu = "rdt_wu")

# Helpers ======================================================================

# One row per (draw, site): draws WITHOUT replacement from the pooled chains.
sample_joint_posterior_draws <- function(pooled_samples, site_names, n_draws, seed = 123) {
  n_sites <- length(site_names)
  stopifnot(n_draws <= nrow(pooled_samples))
  set.seed(seed)
  draws <- pooled_samples[sample.int(nrow(pooled_samples), n_draws)]
  draws[, draw_id := seq_len(.N)]
  long_d <- data.table::melt(draws[, c("draw_id", paste0("d_", seq_len(n_sites))), with = FALSE],
                             id.vars = "draw_id", variable.name = "site_col", value.name = "d")
  long_d[, site_idx := as.integer(sub("^d_", "", site_col))]
  long_eir <- data.table::melt(draws[, c("draw_id", paste0("eir_", seq_len(n_sites))), with = FALSE],
                               id.vars = "draw_id", variable.name = "site_col", value.name = "eir")
  long_eir[, site_idx := as.integer(sub("^eir_", "", site_col))]
  long <- merge(long_d[, .(draw_id, site_idx, d)], long_eir[, .(draw_id, site_idx, eir)], by = c("draw_id", "site_idx"))
  long[, site_name := site_names[site_idx]]
  long[, .(draw_id, site_name, d, eir)]
}

# 95% interval of prevalence draws as a survey of `denom` children with design
# effect `deff` would observe them; plain quantiles when there is no survey.
binom_resampled_quantiles <- function(prevs, denom, deff = 1) {
  if (is.na(denom)) return(stats::quantile(prevs, c(0.025, 0.975)))
  if (is.na(deff)) deff <- 1
  n_eff <- max(1, round(denom / deff))
  resampled <- rbinom(length(prevs), size = n_eff, prob = prevs) / n_eff
  stats::quantile(resampled, c(0.025, 0.975))
}

run_hrp2_grid_task <- function(task_args) do.call(run_hrp2_grid_adm1, task_args)

# Loop over variants ===========================================================
for (v in seq_len(nrow(VARIANTS))) {
RDT_SENSITIVITY  <- VARIANTS$rdt_sensitivity[v]
SIGMA_D_PRIOR_SD <- VARIANTS$sigma_d_prior_sd[v]
MCMC_OUT_DIR     <- mcmc_variant_dir(RDT_SENSITIVITY, SIGMA_D_PRIOR_SD)
POSTERIOR_PRED_OUT_DIR <- file.path(MCMC_OUT_DIR, "posterior_predictive")
BUNDLE_TAG       <- sprintf("v3_multinomial_deff_%s", mcmc_variant_tag(RDT_SENSITIVITY, SIGMA_D_PRIOR_SD))
cat(sprintf("\n########## 4: RDT sensitivity %.2f, sigma_d prior SD %s -> %s ##########\n",
            RDT_SENSITIVITY, format(SIGMA_D_PRIOR_SD), POSTERIOR_PRED_OUT_DIR))
dir.create(POSTERIOR_PRED_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# The fit's own RDT sensitivity must match (fits before 2026-10-08 recorded none: 1)
if (SUBMIT_POST || AGG_POST) {
  mcmc_inputs <- readRDS(file.path(MCMC_OUT_DIR, "mcmc_inputs.RDS"))
  se_fit <- if (is.null(mcmc_inputs$rdt_sensitivity)) 1 else mcmc_inputs$rdt_sensitivity
  if (!isTRUE(all.equal(se_fit, RDT_SENSITIVITY))) {
    stop(sprintf("%s was fitted with RDT sensitivity %.2f, not %.2f.", MCMC_OUT_DIR, se_fit, RDT_SENSITIVITY))
  }
}

# Sample + submit ==============================================================
# SECTION F/G -- duplicated (d, EIR) pairs within a site (rejected MH proposals)
# are simulated once; draw_lookup keeps every draw for the merge in Section I.
# Each unique combination gets its own seed so simulator noise is included.
if (SUBMIT_POST) {
  mcmc_combined  <- readRDS(file.path(MCMC_OUT_DIR, "mcmc_combined.RDS"))
  pooled_samples <- mcmc_combined$pooled_samples
  site_names     <- mcmc_combined$chains[[1]]$site_names
  rm(mcmc_combined)

  adm1_combos <- use_setup_base_pars(readRDS(file.path(DIR_GRID, "adm1_combos_with_paths.RDS")))   # current grid setup's cache
  adm1_combos[, site_name := sprintf("%s_%s", iso, name_1)]
  adm1_combos <- adm1_combos[site_name %in% site_names]
  if (nrow(adm1_combos) != length(site_names)) {
    stop(sprintf("adm1_combos has %d of the fit's %d sites.", nrow(adm1_combos), length(site_names)))
  }
  adm1_combos <- merge(adm1_combos, deletion_by_iso(unique(adm1_combos$iso))[, .(iso, deletion_pct)], by = "iso")

  draws <- sample_joint_posterior_draws(pooled_samples, site_names, N_POSTERIOR_PRED_DRAWS, seed = DRAW_SEED)
  draws <- merge(draws, adm1_combos[, .(site_name, iso, name_1, base_pars_path, deletion_pct)], by = "site_name")
  draws[, hrp2_scale := d]
  draw_lookup <- draws[, .(iso, name_1, eir, hrp2_scale, draw_id, deletion_pct)]

  tasks <- list()
  seed_counter <- 0L
  for (site in site_names) {
    row          <- adm1_combos[site_name == site]
    unique_draws <- unique(draws[site_name == site], by = c("eir", "hrp2_scale"))
    chunk_id     <- ceiling(seq_len(nrow(unique_draws)) / CHUNK_SIZE)
    for (chunk in split(unique_draws, chunk_id)) {
      chunk_seeds  <- DRAW_SEED + seed_counter + seq_len(nrow(chunk))
      seed_counter <- seed_counter + nrow(chunk)
      tasks[[length(tasks) + 1L]] <- list(
        base_pars_path   = row$base_pars_path,
        hrp2_grid        = chunk[, .(eir, hrp2_scale)],
        iso              = row$iso,
        name_1           = row$name_1,
        seeds            = chunk_seeds,
        human_population = POSTERIOR_PRED_HUMAN_POPULATION,
        joint            = TRUE   # TP/FP/FN/TN cells for Section I
      )
    }
  }
  cat(sprintf("Built %d tasks: %d sites, %d unique (d, EIR) combinations of %d draws.\n",
              length(tasks), length(site_names),
              nrow(unique(draw_lookup, by = c("iso", "name_1", "eir", "hrp2_scale"))), nrow(draw_lookup)))

  bundle_posterior_pred <- task_create_bulk_call(
    run_hrp2_grid_task,
    data        = tasks,
    environment = "calibration",
    bundle_name = safe_bundle_name(sprintf("joint_posterior_predictive_%s", BUNDLE_TAG)),
    resources   = hipercow_resources(cores = 1)
  )
  saveRDS(bundle_posterior_pred, file.path(POSTERIOR_PRED_OUT_DIR, "bundle_posterior_pred.RDS"))
  saveRDS(draw_lookup, file.path(POSTERIOR_PRED_OUT_DIR, "draw_lookup.RDS"))
  cat(sprintf("Submitted %d tasks | bundle id: %s\n", length(tasks), bundle_posterior_pred$name))
}

# Collect ======================================================================
# SECTION H
if (COLLECT_POST) {
  bundle_posterior_pred <- readRDS(file.path(POSTERIOR_PRED_OUT_DIR, "bundle_posterior_pred.RDS"))
  coll    <- collect_bundle(bundle_posterior_pred, MAX_RETRIES, label = "posterior-predictive task")
  results <- data.table::rbindlist(coll$results, use.names = TRUE, fill = TRUE)
  saveRDS(results, file.path(POSTERIOR_PRED_OUT_DIR, "posterior_predictive_results.RDS"))
  cat(sprintf("Compiled %s rows from %d/%d successful tasks.\n",
              formatC(nrow(results), format = "d", big.mark = ","), sum(coll$ok), length(coll$ok)))
}

# Aggregate ====================================================================
# SECTION I -- mean + 95% PPI per (site, year)
if (AGG_POST) {
  if (!COLLECT_POST) results <- readRDS(file.path(POSTERIOR_PRED_OUT_DIR, "posterior_predictive_results.RDS"))
  draw_lookup <- readRDS(file.path(POSTERIOR_PRED_OUT_DIR, "draw_lookup.RDS"))

  # One row per (site, combination, year) in results; the merge fans out to
  # every draw that shares a combination.
  results <- merge(results, draw_lookup, by = c("iso", "name_1", "eir", "hrp2_scale"), allow.cartesian = TRUE)
  results[, hrp2_prev_adj := rdt_prev_observed(hrp2_prev, deletion_pct, RDT_SENSITIVITY)]
  has_joint <- all(c("p_tp", "p_fp", "p_fn", "p_tn") %in% names(results))
  if (has_joint) {
    results[, c("q_tp", "q_fp", "q_fn", "q_tn") := joint_cells_observed(p_tp, p_fp, p_fn, p_tn, deletion_pct, RDT_SENSITIVITY)]
  } else {
    cat("Section I: no joint cells in the results (runs submitted without joint = TRUE) -- TP/FP/FN/TN skipped.\n")
  }

  mis <- data.table::as.data.table(readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS")))
  mis_obs_denoms <- unique(mis[, .(iso = iso3c, name_1, year, pfpr_rdt_denom, pfpr_lm_denom, deff_rdt, deff_lm, joint_n, deff_joint)])
  results <- merge(results, mis_obs_denoms, by = c("iso", "name_1", "year"), all.x = TRUE)

  # Marginal RDT/LM from the joint cells, and the Wu et al. comparator, on the
  # same children (joint_n) as the TP/FP/FN/TN check. LM is unaffected by
  # deletion and RDT sensitivity, so it comes from the unadjusted cells.
  if (has_joint) {
    results[, `:=`(rdt_joint = q_tp + q_fp, lm_joint = p_tp + p_fn)]
    obs_joint <- unique(mis[!is.na(joint_n) & joint_n > 0,
                            .(iso = iso3c, name_1, year, rdt = (joint_tp + joint_fp) / joint_n, lm = (joint_tp + joint_fn) / joint_n)])
    wu_sd <- wu_residual_sd_loco(obs_joint)
    saveRDS(wu_sd, file.path(POSTERIOR_PRED_OUT_DIR, "wu_residual_sd.RDS"))
    results <- merge(results, wu_sd[, .(iso, wu_sd)], by = "iso", all.x = TRUE)
    # One deviation per (draw, site, year): the Wu line is a survey-level relationship
    set.seed(WU_SEED)
    results[, rdt_wu := wu_rdt_from_lm(lm_joint, wu_sd)]
    cat(sprintf("Wu comparator (%s): leave-one-country-out residual RMS %.2f-%.2f on the logit scale%s.\n",
                WU_UNCERTAINTY, min(wu_sd$wu_sd), max(wu_sd$wu_sd),
                if (identical(WU_UNCERTAINTY, "lm_only")) " (not used)" else ""))
  }

  set.seed(123)
  posterior_pred_agg <- results[, {
    hrp2_q <- binom_resampled_quantiles(hrp2_prev_adj, pfpr_rdt_denom[1], deff_rdt[1])
    lm_q   <- binom_resampled_quantiles(lm_prev,       pfpr_lm_denom[1],  deff_lm[1])
    out <- list(hrp2_mean = mean(hrp2_prev_adj), hrp2_lo = hrp2_q[[1]], hrp2_hi = hrp2_q[[2]],
                lm_mean   = mean(lm_prev),       lm_lo   = lm_q[[1]],   lm_hi   = lm_q[[2]],
                n_draws   = .N)
    if (has_joint) {
      for (cl in c("tp", "fp", "fn", "tn")) {
        q  <- get(paste0("q_", cl))
        qq <- binom_resampled_quantiles(q, joint_n[1], deff_joint[1])
        out[[paste0(cl, "_mean")]] <- mean(q)
        out[[paste0(cl, "_lo")]]   <- qq[[1]]
        out[[paste0(cl, "_hi")]]   <- qq[[2]]
      }
    }
    out
  }, by = .(iso, name_1, year)]

  # Marginal RDT (rdtj = TP + FP) and LM (lmj = TP + FN) from our model, and
  # RDT from the Wu et al. line (rdtwu). A separate pass with its own seed, so
  # adding these leaves the TP/FP/FN/TN intervals above unchanged.
  if (has_joint) {
    set.seed(WU_SEED + 1L)
    marg_agg <- results[, {
      out <- list()
      for (k in names(MARGINAL_COLS)) {
        x  <- get(MARGINAL_COLS[[k]])
        # "lm_only": the Wu interval is the spread of the LM draws through the line, no survey noise
        n_k <- if (k == "rdtwu" && identical(WU_UNCERTAINTY, "lm_only")) NA_real_ else joint_n[1]
        qq <- binom_resampled_quantiles(x, n_k, deff_joint[1])
        out[[paste0(k, "_mean")]] <- mean(x)
        out[[paste0(k, "_lo")]]   <- qq[[1]]
        out[[paste0(k, "_hi")]]   <- qq[[2]]
      }
      out
    }, by = .(iso, name_1, year)]
    posterior_pred_agg <- merge(posterior_pred_agg, marg_agg, by = c("iso", "name_1", "year"))
  }

  saveRDS(posterior_pred_agg, file.path(POSTERIOR_PRED_OUT_DIR, "posterior_predictive_summary.RDS"))
  cat(sprintf("Aggregated to %d (site, year) rows across %d sites.\n",
              nrow(posterior_pred_agg), uniqueN(posterior_pred_agg[, .(iso, name_1)])))
}
}   # end variant loop
