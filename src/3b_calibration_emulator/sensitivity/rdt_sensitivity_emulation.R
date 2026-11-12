# =============================================================================
# src/3b_calibration_emulator/sensitivity/rdt_sensitivity_emulation.R
#
# ANALYSIS COMPONENT: SENSITIVITY ANALYSIS (calibration: RDT sensitivity 0.8 and 1.0, stage 2)
#
# RDT sensitivity analysis, step 3b: fit the GP emulators to the RDT
# sensitivity 0.8 and 1.0 grid log-likelihoods (from
# 3a_calibration_grid/sensitivity/rdt_sensitivity_log_lik.R).
# Writes DIR_EMULATOR/sensitivity_rdt_se_080/ and .../sensitivity_rdt_se_100/.
#
# Run SUBMIT first, then COLLECT once the cluster tasks have finished.
# =============================================================================
source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

SUBMIT  <- FALSE   # submit one emulator fit per site, per variant
COLLECT <- FALSE   # collect + validate (blocks until the tasks finish)

RDT_SENSITIVITY_SENS <- setdiff(RDT_SENSITIVITY_OPTIONS, RDT_SENSITIVITY_MAIN)   # 0.8, 1.0

run_step("src/3b_calibration_emulator/emulation.R",
         list(list(SUBMIT_EMULATOR = SUBMIT, AGG_EMULATOR = COLLECT, RDT_SENSITIVITY_RUN = RDT_SENSITIVITY_SENS)))
