# =============================================================================
# src/3a_calibration_grid/sensitivity/rdt_sensitivity_log_lik.R
#
# ANALYSIS COMPONENT: SENSITIVITY ANALYSIS (calibration: RDT sensitivity 0.8 and 1.0, stage 1)
#
# RDT sensitivity analysis, step 3a: recompute the grid log-likelihood with
# RDT sensitivity 0.8 and 1.0 (main analysis: 0.9). Only Section F of
# grid_simulation.R depends on RDT sensitivity, so the simulations are reused.
#
# Writes DIR_GRID/hrp2_grid_log_lik_rdt_se_080.RDS and _rdt_se_100.RDS.
# =============================================================================
source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")
source("./src/functions/pipeline_fns.R")

RDT_SENSITIVITY_SENS <- setdiff(RDT_SENSITIVITY_OPTIONS, RDT_SENSITIVITY_MAIN)   # 0.8, 1.0

run_step("src/3a_calibration_grid/grid_simulation.R",
         list(list(CALC_LIKELIHOOD = TRUE, RDT_SENSITIVITY_RUN = RDT_SENSITIVITY_SENS)))
