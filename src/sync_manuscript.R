# =============================================================================
# src/sync_manuscript.R
#
# ANALYSIS COMPONENT: OUTPUTS (manuscript build inputs)
#
# Copies every code-generated file the manuscript uses from figures/ and
# results/ into the manuscript folder, under the names main.tex and
# appendix.tex reference, so the manuscript never holds a hand-copied or stale
# figure or number. The manuscript is kept outside this repository; set
# HRP2_MANUSCRIPT_DIR in .Renviron to its local folder (e.g. an Overleaf git
# clone). Run after src/plots.R and src/results_numbers.Rmd
# (launch.R action "manuscript"). Stops if a source file is missing.
#
# Hand-written manuscript sources are not touched: main.tex, appendix.tex,
# figures/nmcp_guidance*.tex (TikZ; numbers via \nm{} keys), *.bib, lancet.bst.
# =============================================================================
source("./src/functions/paths.R")


# Without a local manuscript folder there is nothing to copy: the figures stay
# in figures/ and the numbers in results/numbers.tex.
if (is.na(MANUSCRIPT_DIR) || !dir.exists(MANUSCRIPT_DIR)) {
  message("HRP2_MANUSCRIPT_DIR is not set (or does not exist), so nothing was copied. ",
          "Figures are in figures/ and numbers in results/numbers.tex; set HRP2_MANUSCRIPT_DIR ",
          "in .Renviron to copy them into a local manuscript folder (see README).")
} else {
  # source (relative to the project root)          -> destination (relative to MANUSCRIPT_DIR)
  SYNC <- c(
    "results/numbers.tex"                          = "numbers.tex",
    "figures/schematic.png"                        = "figures/schematic.png",
    "figures/fpr_fdr_maps_ssa.png"                 = "figures/fpr_fdr_maps_ssa.png",
    "figures/combined_ratio_density_and_forest_plot.png" = "figures/combined_ratio_density_and_forest_plot.png",
    "figures/data_stream.png"                      = "figures/data stream.png",
    "figures/fpfn_ame_forest.png"                  = "figures/fpfn_ame_forest.png",
    "figures/static/model_framework.png"           = "figures/model framework.png",
    "figures/prevalence_rdt_lm_wu_scatter.png"     = "figures/prevalence_rdt_lm_wu_scatter.png"
  )

  src <- file.path(PROJECT_ROOT, names(SYNC))
  missing <- src[!file.exists(src)]
  if (length(missing)) stop("Not yet produced (run src/plots.R and src/results_numbers.Rmd first):\n  ",
                            paste(missing, collapse = "\n  "))
  dst <- file.path(MANUSCRIPT_DIR, SYNC)
  for (d in unique(dirname(dst))) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ok <- file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
  if (!all(ok)) stop("Could not copy: ", paste(src[!ok], collapse = ", "))
  cat(sprintf("Copied %d files into %s\n", length(dst), MANUSCRIPT_DIR))
}
