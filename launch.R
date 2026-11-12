# =============================================================================
# launch.R
#
# ANALYSIS COMPONENT: PIPELINE CONTROL (runs any stage of any analysis)
#
# Central launcher: runs any stage of the pipeline, for the main analysis and
# (optionally) the sensitivity analyses, without editing the step scripts.
#
# From a terminal (project root):
#   Rscript launch.R list                         # show every action
#   Rscript launch.R 3c_collect 3c_pool           # run actions in order
#   Rscript launch.R --dry-run 3a_likelihood 3b_submit
#   Rscript launch.R --main-only 4_submit         # skip the sensitivity fits
# Interactively: RUN <- "export_posteriors"; source("launch.R")
#
# Cluster actions come in submit / collect pairs. A submit returns as soon as
# the tasks are queued (except 3c_run, which waits ~2 h per fit for the prep
# task before submitting chains); a collect waits until the tasks finish.
# Typical order:
#   1_data_prep -> 2_regression_fit -> 2_regression_selection
#   3a_cache_submit -> 3a_cache_collect -> 3a_grid_submit -> 3a_grid_collect -> 3a_likelihood
#   3b_submit -> 3b_collect
#   3c_run -> 3c_collect -> 3c_pool
#   4_submit -> 4_collect -> 4_aggregate
#   plots -> results_numbers -> export_posteriors -> manuscript
# The sections behind each action (and what they write) are described in each
# step script's header.
# =============================================================================
source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

# Settings ======================================================================
# Actions to run when sourced interactively. Set RUN before source("launch.R"),
# e.g. RUN <- c("3c_collect", "3c_pool"); a RUN left over from an earlier
# source() is cleared at the end of each run.
if (!exists("RUN", envir = globalenv(), inherits = FALSE)) RUN <- character(0)
INCLUDE_RDT_SENSITIVITY <- TRUE     # SENSITIVITY ANALYSIS: also run RDT sensitivity 0.8 and 1.0 (steps 3a-3c)
INCLUDE_SIGMA_PRIOR_SD  <- TRUE     # SENSITIVITY ANALYSIS: also run sigma_d prior SD 0.05 and 0.1 (step 3c)
DRY_RUN <- FALSE                    # print what would run, without running it

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if ("--dry-run" %in% args)   DRY_RUN <- TRUE
  if ("--main-only" %in% args) INCLUDE_RDT_SENSITIVITY <- INCLUDE_SIGMA_PRIOR_SD <- FALSE
  RUN <- setdiff(args, c("--dry-run", "--main-only"))
}

# Which fits each step covers ===================================================
rdt_values <- function() if (INCLUDE_RDT_SENSITIVITY) RDT_SENSITIVITY_OPTIONS else RDT_SENSITIVITY_MAIN
mcmc_fits <- function() {
  v <- mcmc_variants()
  v[(is_main_sigma_sd(sigma_d_prior_sd) & rdt_sensitivity %in% rdt_values()) |
      (INCLUDE_SIGMA_PRIOR_SD & !is_main_sigma_sd(sigma_d_prior_sd))]
}
mcmc_overrides <- function(flags) {
  v <- mcmc_fits()
  lapply(seq_len(nrow(v)), function(i) c(flags, list(RDT_SENSITIVITY = v$rdt_sensitivity[i],
                                                      SIGMA_D_PRIOR_SD = v$sigma_d_prior_sd[i])))
}

S1  <- c("src/1_data_prep/01_format_dhs_mis_data.R",
         "src/1_data_prep/02_survey_design_effects.R",
         "src/1_data_prep/03_structure_mis_data.R")
S2  <- "src/2_regression/fit_fdr_fpr_regression.R"
S2_SEL <- "src/2_regression/fixed_effect_selection.R"
S3A <- "src/3a_calibration_grid/grid_simulation.R"
S3B <- "src/3b_calibration_emulator/emulation.R"
S3C <- "src/3c_calibration_mcmc/hierarchical_mcmc.R"
S4  <- "src/4_sample_posterior/sample_posterior.R"

# Actions =======================================================================
# Each action is a script plus a list of override sets (one run of the script
# per set), or "source" for scripts without settings.
ACTIONS <- list(
  "1_data_prep"          = list(desc = "Format DHS/MIS data, design effects, site-level prevalence (step 1)",
                                fn = function() list(script = S1, overrides = "source")),
  "2_regression_fit"     = list(desc = "Fit the brms regression locally (hours) and compute FPR/FDR (step 2)",
                                fn = function() list(script = S2, overrides = list(list(FIT_MODEL = TRUE, COVARIATE_FPFN = TRUE)))),
  "2_regression_fpfn"    = list(desc = "FPR/FDR summaries from the existing fit (step 2)",
                                fn = function() list(script = S2, overrides = list(list(FIT_MODEL = FALSE, COVARIATE_FPFN = TRUE)))),
  "2_regression_selection" = list(desc = "Fixed-effect selection: LOO refits (hours per model) and ELPD table (step 2)",
                                fn = function() list(script = S2_SEL, overrides = list(list(FIT = TRUE, COMPARE = TRUE)))),
  "3a_cache_submit"      = list(desc = "Submit per-site base_pars builds (step 3a C)",
                                fn = function() list(script = S3A, overrides = list(list(BUILD_CACHE1 = TRUE)))),
  "3a_cache_collect"     = list(desc = "Collect base_pars cache paths (step 3a C)",
                                fn = function() list(script = S3A, overrides = list(list(BUILD_CACHE2 = TRUE)))),
  "3a_grid_submit"       = list(desc = "Submit the grid simulations, missing points only (step 3a D; ~3 h to queue)",
                                fn = function() list(script = S3A, overrides = list(list(SUBMIT_GRID = TRUE)))),
  "3a_grid_collect"      = list(desc = "Collect the grid simulations (step 3a E)",
                                fn = function() list(script = S3A, overrides = list(list(AGG_GRID = TRUE)))),
  "3a_likelihood"        = list(desc = "Grid log-likelihood for each RDT sensitivity (step 3a F)",
                                fn = function() list(script = S3A, overrides = list(list(CALC_LIKELIHOOD = TRUE, RDT_SENSITIVITY_RUN = rdt_values())))),
  "3b_submit"            = list(desc = "Submit GP emulator fits for each RDT sensitivity (step 3b C)",
                                fn = function() list(script = S3B, overrides = list(list(SUBMIT_EMULATOR = TRUE, RDT_SENSITIVITY_RUN = rdt_values())))),
  "3b_collect"           = list(desc = "Collect and validate the emulators (step 3b D)",
                                fn = function() list(script = S3B, overrides = list(list(AGG_EMULATOR = TRUE, RDT_SENSITIVITY_RUN = rdt_values())))),
  "3c_run"               = list(desc = "MCMC prep task (waits ~2 h each) then chains, for every fit (step 3c B/C)",
                                fn = function() list(script = S3C, overrides = mcmc_overrides(list(RUN_MCMC = TRUE)))),
  "3c_collect"           = list(desc = "Collect MCMC chains for every fit (step 3c D)",
                                fn = function() list(script = S3C, overrides = mcmc_overrides(list(COLLECT_MCMC = TRUE)))),
  "3c_pool"              = list(desc = "R-hat, pool chains, posterior summaries for every fit (step 3c E)",
                                fn = function() list(script = S3C, overrides = mcmc_overrides(list(POOL_MCMC = TRUE)))),
  "4_submit"             = list(desc = "Posterior predictive runs for the main fit: sample and submit (step 4 F/G)",
                                fn = function() list(script = S4, overrides = list(list(SUBMIT_POST = TRUE)))),
  "4_collect"            = list(desc = "Collect posterior predictive runs (step 4 H)",
                                fn = function() list(script = S4, overrides = list(list(COLLECT_POST = TRUE)))),
  "4_aggregate"          = list(desc = "Posterior predictive intervals per site and year (step 4 I)",
                                fn = function() list(script = S4, overrides = list(list(AGG_POST = TRUE)))),
  "plots"                = list(desc = "Manuscript figures (src/plots.R -> figures/)",
                                fn = function() list(script = "src/plots.R", overrides = "source")),
  "results_numbers"      = list(desc = "Knit src/results_numbers.Rmd (RDT sensitivity 0.9) -> results/",
                                fn = function() list(script = "src/results_numbers.Rmd", overrides = "render")),
  "export_posteriors"    = list(desc = "Posterior draws and summaries, regression and every calibration fit -> results/posteriors/",
                                fn = function() list(script = "src/export_posteriors.R", overrides = "source")),
  "manuscript"           = list(desc = "Copy figures and numbers.tex into HRP2_MANUSCRIPT_DIR (src/sync_manuscript.R)",
                                fn = function() list(script = "src/sync_manuscript.R", overrides = "source"))
)

list_actions <- function() {
  cat("Actions (run in the order given):\n")
  for (nm in names(ACTIONS)) cat(sprintf("  %-18s %s\n", nm, ACTIONS[[nm]]$desc))
  cat(sprintf("\nRDT sensitivity fits included: %s | sigma_d prior SD fits included: %s\n",
              INCLUDE_RDT_SENSITIVITY, INCLUDE_SIGMA_PRIOR_SD))
}

run_action <- function(name) {
  if (!name %in% names(ACTIONS)) stop(sprintf("Unknown action '%s' -- run `Rscript launch.R list`.", name))
  a <- ACTIONS[[name]]$fn()
  cat(sprintf("\n==================== %s: %s ====================\n", name, ACTIONS[[name]]$desc))
  describe <- function(o) paste(sprintf("%s = %s", names(o), vapply(o, function(x) paste(format(x), collapse = ","), "")), collapse = "; ")
  if (identical(a$overrides, "source")) {
    for (s in a$script) {
      cat(sprintf("  source(%s)\n", s))
      if (!DRY_RUN) source(s, local = globalenv())
    }
  } else if (identical(a$overrides, "render")) {
    cat(sprintf("  render(%s)\n", a$script))
    if (!DRY_RUN) {
      dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
      # The knitted HTML is a by-product: write it to a temporary folder so
      # only results/numbers.tex lands in the repository
      rmarkdown::render(a$script, output_dir = tempdir(), intermediates_dir = tempdir(), envir = new.env())
    }
  } else if (DRY_RUN) {
    for (o in a$overrides) cat(sprintf("  %s | %s\n", a$script, describe(o)))
  } else {
    run_step(a$script, a$overrides)
  }
  invisible(NULL)
}

# Run ===========================================================================
if (length(RUN) == 0 || identical(RUN, "list")) {
  list_actions()
} else {
  bad <- setdiff(RUN, names(ACTIONS))
  if (length(bad)) stop(sprintf("Unknown action(s): %s -- run `Rscript launch.R list`.", paste(bad, collapse = ", ")))
  for (nm in RUN) run_action(nm)
  cat(sprintf("\nlaunch.R finished: %s\n", paste(RUN, collapse = ", ")))
  if (!length(commandArgs(trailingOnly = TRUE))) rm(RUN, envir = globalenv())   # next interactive source() lists actions unless RUN is set again
}
