# =============================================================================
# src/3c_calibration_mcmc/sensitivity/rdt_sensitivity_mcmc.R
#
# ANALYSIS COMPONENT: SENSITIVITY ANALYSIS (calibration: RDT sensitivity 0.8 and 1.0, stage 3)
#
# RDT sensitivity analysis, step 3c: hierarchical MCMC on the RDT sensitivity
# 0.8 and 1.0 emulators (sigma_d prior SD 0.01, as the main fit).
# Writes DIR_MCMC/sensitivity_rdt_se_080/ and .../sensitivity_rdt_se_100/.
#
# RUN_MCMC blocks while each variant's prep task runs (~2 h each), then
# submits its chains. Turn on COLLECT_MCMC and POOL_MCMC once the chains finish.
# =============================================================================
source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

RUN_MCMC     <- FALSE
COLLECT_MCMC <- FALSE
POOL_MCMC    <- FALSE

RDT_SENSITIVITY_SENS <- setdiff(RDT_SENSITIVITY_OPTIONS, RDT_SENSITIVITY_MAIN)   # 0.8, 1.0

run_step("src/3c_calibration_mcmc/hierarchical_mcmc.R", lapply(RDT_SENSITIVITY_SENS, function(se) list(
  RUN_MCMC = RUN_MCMC, COLLECT_MCMC = COLLECT_MCMC, POOL_MCMC = POOL_MCMC,
  RDT_SENSITIVITY = se, SIGMA_D_PRIOR_SD = SIGMA_D_PRIOR_SD_MAIN
)))
