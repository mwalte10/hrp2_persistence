# =============================================================================
# src/3c_calibration_mcmc/sensitivity/sigma_prior_sd_mcmc.R
#   (was src/09i_sensitivity_sigma_d_sd_0.05.R and 09j_sensitivity_sigma_d_sd_0.1.R)
#
# ANALYSIS COMPONENT: SENSITIVITY ANALYSIS (calibration: sigma_P prior SD 0.05 and 0.1, stage 3)
#
# sigma_d prior sensitivity analysis, step 3c. The main fit's prior,
# sigma_d ~ Normal(0.3, 0.01), is deliberately tight; these fits loosen its SD
# to 0.05 (5x) and 0.1 (10x, the width used by the retired country-spread
# model) to show how much of the fitted between-site spread comes from the
# prior. Everything else matches the main fit, at RDT sensitivity 0.9 only.
# Writes DIR_MCMC/sensitivity_sigma_d_sd_0.05_rdt_se_090/ and
#        DIR_MCMC/sensitivity_sigma_d_sd_0.1_rdt_se_090/.
# =============================================================================
source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

RUN_MCMC     <- FALSE
COLLECT_MCMC <- FALSE
POOL_MCMC    <- FALSE

run_step("src/3c_calibration_mcmc/hierarchical_mcmc.R", lapply(SIGMA_D_PRIOR_SD_OPTIONS, function(sd) list(
  RUN_MCMC = RUN_MCMC, COLLECT_MCMC = COLLECT_MCMC, POOL_MCMC = POOL_MCMC,
  RDT_SENSITIVITY = RDT_SENSITIVITY_MAIN, SIGMA_D_PRIOR_SD = sd
)))
