# =============================================================================
# run_sections.R
#
# ANALYSIS COMPONENT: PIPELINE CONTROL (run one script with changed settings)
#
# Run a step script from the command line with some of its settings changed,
# without editing it:
#   Rscript run_sections.R src/3a_calibration_grid/grid_simulation.R SUBMIT_GRID=TRUE
#   Rscript run_sections.R src/3c_calibration_mcmc/hierarchical_mcmc.R RUN_MCMC=TRUE RDT_SENSITIVITY=0.8
# Each NAME=VALUE becomes STEP_OVERRIDES$NAME, which the script reads through
# take_overrides()/setting() (src/functions/pipeline_fns.R). Values are parsed
# as R expressions, so c(0.8, 1) and "binomial_deff" work (quote for the shell).
#
# Run from the project root. Long cluster submissions (e.g. the ~11,000-task
# grid) are best run this way in your own terminal.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_sections.R <script.R> NAME=VALUE ...")
script <- args[1]
STEP_OVERRIDES <- list()
for (kv in args[-1]) {
  i <- regexpr("=", kv, fixed = TRUE)
  if (i < 0) stop(sprintf("Expected NAME=VALUE, got '%s'", kv))
  STEP_OVERRIDES[[substr(kv, 1, i - 1)]] <- eval(parse(text = substring(kv, i + 1)))
}
source(script, echo = FALSE)
cat(sprintf("%s DONE: %s\n", script, paste(args[-1], collapse = " ")))
