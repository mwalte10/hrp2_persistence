# =============================================================================
# src/plots.R
#
# ANALYSIS COMPONENT: OUTPUTS (manuscript figures)
#
# Draws every figure in the manuscript and appendix except the model
# framework diagram (figures/static/, see README). Reads the analysis outputs
# of steps 1-4 only (no cluster jobs, no refitting) and writes to FIGURE_DIR.
# src/sync_manuscript.R then copies them into the manuscript (HRP2_MANUSCRIPT_DIR).
#
#   Flag                        File (FIGURE_DIR/)                          Manuscript
#   PLOT_SCHEMATIC              schematic.png                               Fig 1 (fig:dynamics), illustrative
#   PLOT_FPR_FDR_MAPS           fpr_fdr_maps_ssa.png                        Fig 2 (fig:regression)
#   PLOT_DENSITY_FOREST         combined_ratio_density_and_forest_plot.png  Fig 3 (fig:halflife)
#   PLOT_DATA_STREAM            data_stream.png                             Appendix fig:data_stream
#   PLOT_AME                    fpfn_ame_forest.png                         Appendix fig:app_fpr_fdr
#   PLOT_PREVALENCE_WU_SCATTER  prevalence_rdt_lm_wu_scatter.png            Appendix fig:scatter
#
# PLOT_FPR_FDR_MAPS also writes DIR_REGRESSION/fpr_fdr_by_subregion.csv
# (sub-regional rates read by results_numbers.Rmd; derived from DHS/MIS data,
# so it stays in the git-ignored OUTPUT_ROOT).
# Settings can be overridden through run_step() / run_sections.R.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
  library(ggrepel)
  library(ltc)
})
source("./src/functions/paths.R")
source("./src/functions/cluster_fns.R")             # RDT sensitivity settings, colours
source("./src/functions/posterior_summary_fns.R")   # country aggregation, unobserved-site prediction
source("./src/functions/plot_mcmc_diagnostics.R")   # simulate_hierarchical_priors()
source("./src/functions/pipeline_fns.R")            # mcmc_variant_dir()
source("./src/functions/schematic_fns.R")           # Figure 1
source("./src/functions/data_stream_fns.R")         # appendix data stream

# Settings ======================================================================
OVR <- take_overrides()
PLOT_SCHEMATIC             <- setting(OVR, "PLOT_SCHEMATIC", TRUE)
PLOT_DATA_STREAM           <- setting(OVR, "PLOT_DATA_STREAM", TRUE)
PLOT_AME                   <- setting(OVR, "PLOT_AME", TRUE)
PLOT_FPR_FDR_MAPS          <- setting(OVR, "PLOT_FPR_FDR_MAPS", TRUE)
# Age for the FPR/FDR maps: the oldest group in the data (4-year-olds,
# 48-59 months); age is continuous in the model, so evaluated at the midpoint.
FPR_FDR_MAP_AGE_MONTHS     <- setting(OVR, "FPR_FDR_MAP_AGE_MONTHS", 54)
FPR_FDR_MAP_AGE_LABEL      <- setting(OVR, "FPR_FDR_MAP_AGE_LABEL",
                                      "aged 4 years (48-59 months, evaluated at 54 months)")
PLOT_DENSITY_FOREST        <- setting(OVR, "PLOT_DENSITY_FOREST", TRUE)
PLOT_PREVALENCE_WU_SCATTER <- setting(OVR, "PLOT_PREVALENCE_WU_SCATTER", TRUE)

# Calibration figures use the main fit: RDT sensitivity 0.9, sigma_d prior SD 0.01, main priors.
fit_dir <- function(rdt_sensitivity, sigma_d_prior_sd) {
  mcmc_variant_dir(rdt_sensitivity, sigma_d_prior_sd, mu_d_mean = MU_D_PRIOR_MEAN_MAIN, prior_family = PRIOR_FAMILY_MAIN)
}

# Site weights for country half-life summaries: "par_pf" = population at
# risk averaged over each site's MIS survey years; "none" = equal weights.
AGG_WEIGHTS <- "par_pf"

dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 0. FIGURE 1: ILLUSTRATIVE SCHEMATIC (no data; src/functions/schematic_fns.R)
# =============================================================================
if (PLOT_SCHEMATIC) {
  ggsave(file.path(FIGURE_DIR, "schematic.png"), plot_schematic(simulate_schematic_cohort()),
         width = 9.4, height = 10.5, dpi = 600, bg = "white")
  cat("Saved schematic.png\n")
}

# =============================================================================
# 0b. APPENDIX DATA STREAM (steps 1a, 1c, 2, 3c; src/functions/data_stream_fns.R)
# =============================================================================
if (PLOT_DATA_STREAM) {
  ds <- data_stream_counts(DIR_PULL, DIR_STRUCTURE,
                           readRDS(file.path(DIR_REGRESSION, "fit_brms_int_full_local.RDS"))$data,
                           readRDS(file.path(fit_dir(RDT_SENSITIVITY_MAIN, SIGMA_D_PRIOR_SD_MAIN), "mcmc_inputs.RDS"))$site_names)
  ggsave(file.path(FIGURE_DIR, "data_stream.png"), plot_data_stream(ds), width = 9, height = 8, dpi = 600, bg = "white")
  cat("Saved data_stream.png\n")
}

# =============================================================================
# 1. AVERAGE MARGINAL EFFECTS FOREST PLOT (step 2: fpfn_ame.RDS)
# =============================================================================
# Complements the maps (levels) with the magnitude of each covariate's
# influence, in percentage points averaged over the surveyed children.
# Colour follows the appendix tables: blue = increase, red = decrease,
# lighter where the 95% CrI crosses zero.
if (PLOT_AME) {
  ame <- as.data.table(readRDS(file.path(DIR_REGRESSION, "fpfn_ame.RDS")))
  ame[, outcome := factor(fifelse(outcome == "P(RDT+ | microscopy-)", "False positive rate, P(RDT+ | LM-)",
                                  "False discovery rate, P(LM- | RDT+)"),
                          levels = c("False positive rate, P(RDT+ | LM-)", "False discovery rate, P(LM- | RDT+)"))]
  AME_LABELS <- c("Prevalence: 10-20% vs <10%" = "10-20%",
                  "Prevalence: >20% vs <10%"   = ">20%",
                  "PfHRP2: Medium vs Slight"   = "Medium",
                  "PfHRP2: High vs Slight"     = "High",
                  "Fever"                    = "Fever",
                  "Formal care/ACT"          = "Formal care or ACT\n(febrile children)",
                  "SMC"                      = "SMC",
                  "During/after peak season" = "During/after peak season",
                  "Age (+1 year)"            = "Age (+1 year)")
  if (!all(ame$covariate %in% names(AME_LABELS))) stop("PLOT_AME: unrecognised covariate label")
  ame[, cov_lab := factor(AME_LABELS[covariate], levels = rev(AME_LABELS))]
  # Row groups: prevalence and deletion risk are each compared with their
  # reference level, every child keeping their own level of the other factor
  ame[, grp := factor(fcase(startsWith(covariate, "Prevalence"), "Sub-regional\nprevalence\n(vs <10%)",
                            startsWith(covariate, "PfHRP2"),     "Deletion risk\n(vs slight)",
                            default = "Other\ncovariates"),
                      levels = c("Sub-regional\nprevalence\n(vs <10%)", "Deletion risk\n(vs slight)", "Other\ncovariates"))]
  ame[, shade := fcase(mean < 0 & hi < 0, "dec", mean < 0, "dec_ns", mean > 0 & lo > 0, "inc", default = "inc_ns")]
  ame[, lab := sprintf("%s%.1f (%.1f to %.1f)", fifelse(round(mean, 1) > 0, "+", ""), mean, lo, hi)]
  AME_COLS <- c(dec = "#9C0006", dec_ns = "#C0504D", inc = "#1F4E79", inc_ns = "#5B7DA8")
  p_ame <- ggplot(ame, aes(x = mean, y = cov_lab, colour = shade)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1) +
    geom_point(shape = 15, size = 3) +
    geom_text(aes(label = lab), vjust = -1, size = 3.2, show.legend = FALSE) +
    facet_grid(grp ~ outcome, scales = "free", space = "free_y") +
    # "Insignificant" = 95% CrI includes 0; increases on the first legend row
    scale_colour_manual(values = AME_COLS, limits = c("inc", "inc_ns", "dec", "dec_ns"), name = NULL,
                        labels = c(inc = "Increase", inc_ns = "Increase, insignificant effect",
                                   dec = "Decrease", dec_ns = "Decrease, insignificant effect")) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "0", sprintf("%+g", x)), breaks = scales::breaks_pretty(n = 6),
                       expand = expansion(mult = 0.22)) +
    scale_y_discrete(expand = expansion(add = c(0.5, 0.8))) +
    labs(x = "Average marginal effect (percentage points)", y = NULL) +
    theme_bw(base_size = 13) +
    theme(strip.text = element_text(face = "bold"), strip.text.y = element_text(angle = 0, hjust = 0),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")   # TODO: confirm legend position with Maggie
  ggsave(file.path(FIGURE_DIR, "fpfn_ame_forest.png"), p_ame, width = 10.5, height = 7.4, dpi = 600, bg = "white")
  cat("Saved fpfn_ame_forest.png\n")
}

# =============================================================================
# 1b. FPR / FDR MAPS (step 2: fit_brms_int_full_local.RDS)
# =============================================================================
# Every country's most recent survey on one sub-Saharan Africa page. Method,
# panels and boundary sources: src/functions/fpr_fdr_map_fns.R.
if (PLOT_FPR_FDR_MAPS) {
  source("./src/functions/brms_int_cluster_fns.R")   # prepare_brms_full_data()
  source("./src/functions/fpr_fdr_map_fns.R")
  fit_int_full <- readRDS(file.path(DIR_REGRESSION, "fit_brms_int_full_local.RDS"))
  map_dt <- fpr_fdr_subregion_rates(fit_int_full, FPR_FDR_MAP_AGE_MONTHS)
  # tryCatch: the CSV is often open in Excel, which locks it on Windows
  tryCatch(fwrite(map_dt, file.path(DIR_REGRESSION, "fpr_fdr_by_subregion.csv")),
           error = function(e) warning("fpr_fdr_by_subregion.csv not written (open in another program?): ", conditionMessage(e)))

  surveys <- unique(map_dt[, .(country, survey_year)])[order(country)]
  ssa_list <- list()
  for (i in seq_len(nrow(surveys))) {
    ctry <- surveys$country[i]
    yr   <- surveys$survey_year[i]
    m    <- country_map_sf(ctry, yr, map_dt)
    if (length(m$unmatched) > 0) {
      message(sprintf("%s %d: no boundary for %s", ctry, yr, paste(m$unmatched, collapse = ", ")))
    }
    ssa_list[[ctry]] <- m$sf[, c("country", "prev_hrp2", "smc_lab", "prev_fever",
                                 "prop_care", "fpr", "fdr", "abs_fp")]
  }

  ssa_note <- paste0(
    "\nB: SMC coverage > 1%. D: formal care or ACT, among children with reported fever in the past two weeks (only asked after recent fever).",
    "\nSurveys: ", paste(strwrap(paste(sprintf("%s %d", surveys$country, surveys$survey_year),
                                       collapse = ", "), 160), collapse = "\n"), ".",
    "\nBoundaries: Natural Earth admin-1, except Burundi, Cote d'Ivoire, Ghana and Togo (geoBoundaries gbOpen)."
  )
  # border_lwd = 0: country outlines only, no sub-region borders
  p_ssa <- build_fpr_fdr_page(do.call(rbind, ssa_list), "Sub-Saharan Africa",
                              "Each country's most recent survey (predictions include that survey's random effect)",
                              ssa_note, FPR_FDR_MAP_AGE_LABEL, bg_sf = ssa_background(), border_lwd = 0,
                              publication = TRUE)
  ggsave(file.path(FIGURE_DIR, "fpr_fdr_maps_ssa.png"), p_ssa, width = 14, height = 12.5, dpi = 300, bg = "white")
  cat(sprintf("Saved %s\n", file.path(FIGURE_DIR, "fpr_fdr_maps_ssa.png")))
}

# =============================================================================
# 2. HRP2 HALF-LIFE FOREST WITH OBSERVED RDT:LM RATIO DENSITY (steps 1c, 3c, 4)
# =============================================================================

# --- Observed MIS RDT vs LM, one row per distinct (iso, year, RDT, LM) value --
build_mis_scatter <- function() {
  # RDT vs LM scatter: OBSERVED MIS survey data, NOT the simulator's
  # posterior mean (deliberately different data source from every other
  # plot in this section, which all show malariasimulation's own fitted
  # estimate) -- one point per (site, survey year) with BOTH RDT and LM
  # actually measured in that same round (a site-year missing either
  # channel can't contribute a scatter point at all). LM on x / RDT on y
  # (same convention as the ratio map above, RDT/LM -- numerator on y),
  # coloured by country (14 countries is a legible number of legend
  # entries; individual site/year labels on the plot itself would be
  # unreadable clutter). 1:1 line (dashed) marks where RDT would exactly
  # equal LM.
  SCATTER_COUNTRY_NAME_LOOKUP <- c(
    BDI = "Burundi", KEN = "Kenya", MDG = "Madagascar", MWI = "Malawi",
    RWA = "Rwanda", UGA = "Uganda", BEN = "Benin", BFA = "Burkina Faso",
    GHA = "Ghana", GIN = "Guinea", LBR = "Liberia", MLI = "Mali",
    NGA = "Nigeria", SEN = "Senegal", TGO = "Togo", CIV = "Cote d'Ivoire"
  )
  # Deduplicated on the (iso, year, pfpr_rdt, pfpr_lm) VALUE combination,
  # deliberately NOT including name_1 -- the DHS-location-to-site-file
  # mapping upstream (02_structure_mis_data.R) legitimately assigns the SAME
  # source DHS survey region's prevalence to several different site-file
  # admin1 units (a real, intentional many-to-one mapping, not a bug there),
  # but that means the SAME single real observation would otherwise appear
  # as several near-identical scatter points/regression-fit rows here --
  # confirmed directly: MWI 2012 alone repeats one (rdt, lm) pair across 13
  # different name_1 rows. Deduplicating here (858 -> 435 rows) treats it as
  # the ONE real measurement it actually is.
  mis_scatter <- unique(data.table::as.data.table(
    readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS"))
  )[, .(iso = iso3c, year,
        pfpr_rdt = pfpr_rdt_num / pfpr_rdt_denom,
        pfpr_lm  = pfpr_lm_num  / pfpr_lm_denom)])
  mis_scatter <- mis_scatter[!is.na(pfpr_rdt) & !is.na(pfpr_lm)]
  mis_scatter[, country := SCATTER_COUNTRY_NAME_LOOKUP[iso]]
  if (anyNA(mis_scatter$country)) {
    stop(sprintf("No country name mapping for: %s -- extend SCATTER_COUNTRY_NAME_LOOKUP above.",
                 paste(unique(mis_scatter[is.na(country), iso]), collapse = ", ")))
  }

  # Per-country best-fit line, FORCED THROUGH THE ORIGIN (y ~ 0 + x, not the
  # usual y ~ x with its own intercept): 0 LM prevalence should mean 0 RDT
  # prevalence for either channel, so a single slope per country is a
  # complete, sufficient summary of that country's own RDT/LM relationship
  # (an unconstrained intercept+slope fit would need BOTH numbers reported,
  # not just "the slope").
  country_slopes <- mis_scatter[, .(slope = coef(lm(pfpr_rdt ~ 0 + pfpr_lm))[[1]]), by = country]
  mis_scatter <- merge(mis_scatter, country_slopes, by = "country")
  mis_scatter[, country_label := sprintf("%s (%.1f)", country, slope)]
  # Ordered by slope (descending) so the legend itself reads as a ranked
  # list of which countries' RDT reads furthest above LM, not alphabetical.
  country_label_levels <- unique(mis_scatter[order(-slope), country_label])
  mis_scatter[, country_label := factor(country_label, levels = country_label_levels)]
  mis_scatter
}

# --- Country half-life summaries, prior / posterior predictive bands and
# posterior predictive coverage (was 06 Section Q) ---------------------------
# Country values: per posterior draw, the AGG_WEIGHTS-weighted geometric mean of
# the country's site d_i, then the median and 95% CrI across draws.
build_forest_inputs <- function(mcmc_dir) {
  mcmc_combined  <- readRDS(file.path(mcmc_dir, "mcmc_combined.RDS"))
  mcmc_inputs    <- readRDS(file.path(mcmc_dir, "mcmc_inputs.RDS"))
  pooled_samples <- mcmc_combined$pooled_samples
  site_names     <- mcmc_combined$chains[[1]]$site_names
  post           <- load_hierarchical_posterior(mcmc_dir, DIR_GRID, mcmc_combined = mcmc_combined)
  site_iso       <- post$site_iso
  d_mat <- as.matrix(pooled_samples[, paste0("d_", seq_along(site_names)), with = FALSE])
  colnames(d_mat) <- site_names

  country_summary <- summarise_site_parameters(post, by = "country", weights = AGG_WEIGHTS)
  country_summary[, n_sites_in_country := n_sites]

  # Posterior predictive coverage per country: share of observed site-years
  # inside their 95% PPI -- TP/FP/FN/TN when the summary has joint cells, else
  # RDT/LM. A summary older than the fit belongs to an earlier fit and is ignored.
  ppi_path <- file.path(mcmc_dir, "posterior_predictive", "posterior_predictive_summary.RDS")
  have_ppi <- file.exists(ppi_path) && file.mtime(ppi_path) >= file.mtime(file.path(mcmc_dir, "mcmc_combined.RDS"))
  ppi <- if (have_ppi) data.table::as.data.table(readRDS(ppi_path)) else NULL
  mis <- data.table::as.data.table(readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS")))
  country_summary[, ppc_text := NA_character_]
  if (have_ppi) {
    cells <- c("tp", "fp", "fn", "tn")
    if (all(paste0(cells, "_lo") %in% names(ppi))) {
      obs <- unique(mis[!is.na(joint_n) & joint_n > 0, .(iso = iso3c, name_1, year,
              obs_tp = joint_tp / joint_n, obs_fp = joint_fp / joint_n, obs_fn = joint_fn / joint_n, obs_tn = joint_tn / joint_n)])
      chk <- merge(obs, ppi, by = c("iso", "name_1", "year"))
      for (cl in cells) chk[, (paste0(cl, "_in")) := get(paste0("obs_", cl)) >= get(paste0(cl, "_lo")) & get(paste0("obs_", cl)) <= get(paste0(cl, "_hi"))]
      cov <- chk[, lapply(.SD, function(x) 100 * mean(x)), by = iso, .SDcols = paste0(cells, "_in")]
      for (cl in cells) country_summary[, (paste0("ppi_coverage_", cl, "_pct")) := cov[[paste0(cl, "_in")]][match(iso, cov$iso)]]
      country_summary[, ppc_text := sprintf("TP %.0f%%, FP %.0f%%<br>FN %.0f%%, TN %.0f%%",
                                            ppi_coverage_tp_pct, ppi_coverage_fp_pct, ppi_coverage_fn_pct, ppi_coverage_tn_pct)]
    } else {
      obs <- unique(mis[, .(iso = iso3c, name_1, year, pfpr_rdt = pfpr_rdt_num / pfpr_rdt_denom, pfpr_lm = pfpr_lm_num / pfpr_lm_denom)])
      chk <- merge(obs, ppi[, .(iso, name_1, year, hrp2_lo, hrp2_hi, lm_lo, lm_hi)], by = c("iso", "name_1", "year"))
      cov <- chk[, .(rdt = 100 * mean(pfpr_rdt >= hrp2_lo & pfpr_rdt <= hrp2_hi, na.rm = TRUE),
                     lm  = 100 * mean(pfpr_lm  >= lm_lo   & pfpr_lm  <= lm_hi,   na.rm = TRUE)), by = iso]
      country_summary[, ppc_text := sprintf("RDT %.0f%%, LM %.0f%%", cov$rdt[match(iso, cov$iso)], cov$lm[match(iso, cov$iso)])]
    }
  } else {
    cat(sprintf("  No current posterior_predictive_summary.RDS in %s -- forest plot without coverage.\n", mcmc_dir))
  }
  country_summary[, country_full_name := country]
  country_summary <- country_summary[order(halflife_median)]

  # Prior predictive interval for a new site's half-life (marginal prior on d)
  priors <- simulate_hierarchical_priors(mcmc_inputs$log_eir_bounds, mcmc_inputs$log_d_bounds, priors = mcmc_inputs$priors)
  prior_halflife <- priors$d * log(2)

  list(mcmc_inputs = mcmc_inputs, pooled_samples = pooled_samples, chains = mcmc_combined$chains,
       site_names = site_names, site_iso = site_iso,
       post = post, d_mat = d_mat, country_summary = country_summary, priors = priors, ppi = ppi,
       prior_pred_lo = stats::quantile(prior_halflife, 0.025), prior_pred_hi = stats::quantile(prior_halflife, 0.975))
}

# --- Observed RDT vs LM scatter (left) + country half-life forest (right) ----
# (was 06 Section R2)
# Text size multipliers on top of SCALE_A4 (2026-10-12, for readability at
# print size): COMBINED_TEXT_SCALE for titles, axes, country detail lines,
# legend and caption; COMBINED_LABEL_SCALE for panel A's direct country line
# labels. Point sizes and line widths are unchanged. The figure is taller
# than the original 4.98 in so panel B's four-line country labels still fit.
COMBINED_TEXT_SCALE  <- 1.6
COMBINED_LABEL_SCALE <- 2.2
COMBINED_WIDTH_IN    <- 10.5
COMBINED_HEIGHT_IN   <- 5.6
COMBINED_DENSITY_HEIGHT_IN <- 5.9   # density version: two square panels;
COMBINED_DENSITY_WIDTH_IN  <- 9.6   # width sized to fit them with no blank sides
COMBINED_DENSITY_PANEL_WIDTHS <- c(1, 1)   # patchwork widths are panel widths, so A and B are drawn the same size
COMBINED_DENSITY_GAP_PT <- 14              # space between A and B
COMBINED_DENSITY_TEXT_SCALE <- 1.25        # density version: text this much larger again than COMBINED_TEXT_SCALE
COMBINED_PANEL_WIDTHS <- c(31, 26.6)   # panel A : panel B (was 24 : 26.6); A is square, so its width sets its size
COMBINED_FOREST_ASPECT <- 1            # panel B height / width; fixed so A and B are drawn the same height
# RDT:LM ratio map (ratio_map != NULL replaces panel A): log scale centred on
# 1, teal where LM is higher and red where RDT is higher (the RDT and LM
# colours of Figure 1, src/functions/schematic_fns.R); values outside the
# limits are clamped to the end colours.
RATIO_MAP_LIMITS <- c(0.25, 4)
RATIO_MAP_COL_LM  <- "#00798C"
RATIO_MAP_COL_RDT <- "#D1495B"
RATIO_MAP_COL_NO_DATA <- "#E6DFCF"   # countries / sub-regions outside the analysis
# ratio_map: NULL for the scatter version, or rdt_lm_ratio_sf()'s output plus
# $bg (ssa_background()) for the map version; panel B is the same either way.
# ratio_forest = TRUE replaces panel A with an RDT:LM ratio forest in panel B's
# country order (takes precedence over ratio_map). Site-year ratios beyond
# RATIO_FOREST_LIMITS are clamped to the axis ends.
RATIO_FOREST_LIMITS <- c(0.125, 16)
# ratio_density = TRUE replaces panel A with the density of the site-year
# RDT:LM ratios in panel B's rows, with site points sized by population at risk.
plot_combined_scatter_forest <- function(fi, mis_scatter, ratio_map = NULL, ratio_forest = FALSE,
                                         ratio_density = FALSE) {
  list2env(fi, environment())
  # Whole figure uniformly rescaled to an A4-appropriate physical print
  # size (50.6x24in -> 10.5x4.98in, a scale factor of ~0.2075) -- SCALE_A4
  # multiplies every absolute font/point-size/linewidth/margin value below
  # by that SAME factor as the width/height change, so every relative
  # proportion already tuned (text-to-panel ratio, line weight vs panel
  # border weight, etc.) stays exactly as-is; only the physical print size
  # changes. dpi is left unchanged (350) so print resolution/crispness is
  # unaffected -- only the pixel COUNT shrinks accordingly, which is the
  # correct/expected consequence of a smaller physical size at a fixed dpi.
  SCALE_A4 <- 10.5 / 50.6
  if (ratio_density) COMBINED_TEXT_SCALE <- COMBINED_TEXT_SCALE * COMBINED_DENSITY_TEXT_SCALE

  # Publication spelling (circumflex) for the combined figure only --
  # relabels the factor LEVELS, so the `country` column the label-placement
  # rules below key on ("Cote d'Ivoire") is untouched.
  levels(mis_scatter$country_label) <- gsub("Cote d'Ivoire", "C\u00f4te d'Ivoire", levels(mis_scatter$country_label), fixed = TRUE)

  # 14-colour extension of the ltc package's own "minou" palette (6 anchor
  # hues) via colorRampPalette -- every generated colour sits ON the same
  # hue journey between minou's own actual colours, so all 14 are "close to
  # or on" the palette rather than an unrelated auto-generated rainbow.
  # Named to line_label_levels (same descending-slope order the legend/
  # labels already use) so both geom_point/geom_smooth and the direct
  # labels below share one consistent mapping.
  minou_country_colors <- setNames(
    grDevices::colorRampPalette(ltc::palettes$minou)(nlevels(mis_scatter$country_label)),
    levels(mis_scatter$country_label)
  )

  # --- LHS: scatter, alpha'd down so overplotted points don't compete with
  # the 14 direct line-end labels. Only a SUBSET of countries use ggrepel's
  # force-directed placement + leader-line segment (REPEL_COUNTRIES) --
  # everyone else gets a plain geom_text placed right at their own line's
  # end point with a small fixed nudge and NO segment. Restricting repel to
  # 5 of 14 countries (rather than all 14, as before) means far less mutual
  # crowding for ggrepel to resolve, which was the root cause of the prior
  # version's leader lines visibly drifting away from their true anchor --
  # with fewer competing labels the force-directed solver settles close to
  # a straight line from anchor to text, so the leader line reads as
  # actually connected to its own country's line tip instead of appearing
  # to terminate in open space.
  # Rwanda/Senegal/Malawi's own lines are short (low max observed
  # prevalence) and sit right in the crowded low-prevalence corner where
  # most countries' lines start, so their labels need to travel much
  # further than Ghana/Liberia's own already fairly well-separated line
  # ends -- an EXTRA nudge (nudge_x/nudge_y, a one-time starting push
  # before ggrepel's own iterative repulsion even begins) is layered on
  # top for JUST those three countries specifically.
  # scatter_axis_max_combined widened (1.15x, was 1.03x) + coord_equal
  # (clip="off") together give those long leader lines real panel space to
  # extend into rather than being squeezed against the data's own edge.
  line_label_df <- mis_scatter[, .(x_end = max(pfpr_lm), slope = slope[1]), by = .(country, country_label)]
  line_label_df[, y_end := slope * x_end]
  REPEL_COUNTRIES <- c("Rwanda", "Senegal", "Malawi", "Ghana", "Liberia")
  LONG_LABEL_COUNTRIES <- c("Rwanda", "Senegal", "Malawi")
  repel_label_df <- line_label_df[country %in% REPEL_COUNTRIES]
  plain_label_df <- line_label_df[!country %in% REPEL_COUNTRIES]
  # Senegal's own anchor point (x_end=0.178, y_end=0.173) sits just barely
  # below the 1:1 line already, but the shared RWA/SEN/MWI long-nudge
  # (up and to the left) was pushing its LABEL up into y > x territory.
  # Senegal gets its own nudge instead -- down and to the right, keeping a
  # similar total displacement magnitude (~0.27) so its leader line stays
  # just as long, but now lands below the diagonal.
  repel_label_df[, nudge_x := data.table::fcase(
    country == "Senegal", 0.21,
    country %in% LONG_LABEL_COUNTRIES, -0.22,
    default = 0
  )]
  repel_label_df[, nudge_y := data.table::fcase(
    country == "Senegal", -0.17,
    country %in% LONG_LABEL_COUNTRIES, 0.16,
    default = 0
  )]
  # Uganda and Burkina Faso's line ends sit almost exactly on top of each
  # other (y_end differs by <0.001) -- with the default small nudge both
  # plain labels landed on the same spot and overlapped. Flip Uganda's
  # label to sit to the LEFT of its point (hjust=1), now with only a
  # SMALL nudge (just enough to clear Burkina Faso) so it reads as
  # attached to its own line rather than floating off it. Cote d'Ivoire
  # sits in the same crowded Nigeria/Guinea/Burundi/Togo cluster and gets
  # the same left-attached treatment, both to read as attached to its own
  # line and to stop competing for space with Nigeria's label immediately
  # to its right.
  plain_label_df[, hjust   := data.table::fifelse(country %in% c("Uganda", "Cote d'Ivoire"), 1, 0)]
  plain_label_df[, nudge_x := data.table::fifelse(country %in% c("Uganda", "Cote d'Ivoire"), -0.02, 0.015)]
  plain_label_df[, nudge_y := 0.01]
  # At COMBINED_LABEL_SCALE the Nigeria/Guinea, Burundi/Togo and Togo/Mali
  # labels touch: lift Nigeria above its line end, drop Togo below its own,
  # and drop Mali further.
  plain_label_df[country == "Nigeria", nudge_y := 0.04]
  plain_label_df[country == "Togo",    nudge_y := -0.03]
  plain_label_df[country == "Mali",    nudge_y := -0.04]
  scatter_axis_max_combined <- max(c(mis_scatter$pfpr_rdt, mis_scatter$pfpr_lm), na.rm = TRUE) * 1.15

  p_scatter_combined <- ggplot(mis_scatter, aes(x = pfpr_lm, y = pfpr_rdt, colour = country_label)) +
    geom_abline(slope = 1, intercept = 0, colour = "#898781", linewidth = 1.1 * SCALE_A4, linetype = "dashed") +
    geom_point(size = 3.7 * SCALE_A4, alpha = 0.2) +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, linewidth = 1.25 * SCALE_A4) +
    geom_text(
      data = plain_label_df, aes(x = x_end, y = y_end, label = country_label, colour = country_label,
                                  hjust = hjust),
      inherit.aes = FALSE, size = 5.9 * SCALE_A4 * COMBINED_LABEL_SCALE, fontface = "bold", show.legend = FALSE,
      nudge_x = plain_label_df$nudge_x, nudge_y = plain_label_df$nudge_y
    ) +
    ggrepel::geom_text_repel(
      data = repel_label_df, aes(x = x_end, y = y_end, label = country_label, colour = country_label),
      inherit.aes = FALSE, size = 5.9 * SCALE_A4 * COMBINED_LABEL_SCALE, fontface = "bold", show.legend = FALSE,
      nudge_x = repel_label_df$nudge_x, nudge_y = repel_label_df$nudge_y,
      force = 10, force_pull = 0.3, box.padding = 1, point.padding = 0.6,
      segment.size = 0.5 * SCALE_A4, segment.alpha = 0.9,
      max.iter = 50000, max.time = 5, min.segment.length = 0, max.overlaps = Inf, seed = 1L
    ) +
    scale_colour_manual(values = minou_country_colors) +
    coord_equal(xlim = c(0, scatter_axis_max_combined), ylim = c(0, scatter_axis_max_combined), clip = "off") +
    labs(title = "Underlying MIS data",
         subtitle = "Number in brackets: RDT:LM scalar, the slope of the country's\nbest-fit line through the origin.",
         x = "Observed MIS LM prevalence", y = "Observed MIS RDT (HRP2) prevalence") +
    theme_minimal(base_size = 23 * SCALE_A4 * COMBINED_TEXT_SCALE) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "#0b0b0b", fill = NA, linewidth = 0.85 * SCALE_A4),
      legend.position  = "none",
      plot.title        = element_text(face = "bold", size = 30 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b"),
      plot.subtitle     = element_text(size = 21 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e", lineheight = 1.1),
      axis.text.x      = element_text(size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE),
      axis.text.y      = element_text(size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE),
      axis.title.x     = element_text(size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b"),
      axis.title.y     = element_text(size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b"),
      plot.margin      = margin(t = 4 * SCALE_A4, r = 4 * SCALE_A4, b = 4 * SCALE_A4, l = 4 * SCALE_A4)
    )

  # --- RHS: same country_summary/band_df/etc. Section Q already computed,
  # rebuilt at LHS-matching text size (base_size = 11, same as the scatter
  # above) -- Section Q's own p_country_forest_base is deliberately sized
  # small (base_size 5.25) to fit its own compact standalone page, which
  # would look mismatched/inconsistent next to the scatter at full size, so
  # this composite gets its own text-sized copy rather than reusing that
  # object's baked-in small theme. country_label_big is built fresh here
  # (not derived from country_label) with two deliberate differences: (1)
  # the "n=X"/coverage detail lines are their OWN smaller, dark-grey
  # <span> -- reduces their size RELATIVE to the (bigger, default-ink)
  # country name above them, and stops them competing visually with it; (2)
  # the "Posterior predictive coverage:" explanatory phrase is only written
  # out on Cote d'Ivoire's own row (this plot's top row) -- every other
  # country just shows the bare RDT/LM numbers, since the explanation only
  # needs stating once, not 14 times.
  # --- RHS data, rebuilt here from Section Q's own d_mat/site_iso/
  # pooled_samples/priors rather than reusing country_summary$halflife_mean
  # and band_df directly (2026-09-24 revision, to match the methods text):
  #
  # (1) Country summary = per-draw GEOMETRIC mean of site d (exp of the mean
  # of log d), not the arithmetic mean. The model lives on log(d), so an
  # arithmetic mean is dominated by a country's few longest-tailed sites
  # (e.g. a few edge-of-emulator sites can set CIV's whole country value).
  # Reported as posterior median + 95% equal-tailed CrI, since medians
  # survive the log -> natural back-transform and means do not.
  #
  # (2) Each site's own posterior median is drawn behind its country's
  # estimate -- makes the site-level spread within a country visible
  # directly (no country-specific spread parameter under this model -- see
  # hierarchical_mcmc_fns.R's own header; every site shares one population
  # sigma_d, so a country's own site-to-site scatter here is purely
  # descriptive, not a fitted quantity).
  #
  # (3) The posterior predictive band propagates the SAME hierarchy the
  # prior band does (simulate_hierarchical_priors(): mu_d -> sigma_d ->
  # log d_new, truncated to log_d_bounds), draw by draw -- posterior_summary_
  # fns.R::predict_unobserved_site(), the same function Section Q now uses.
  forest_dt <- data.table::copy(country_summary)
  forest_dt[, country_full_name := gsub("Cote d'Ivoire", "C\u00f4te d'Ivoire", country_full_name, fixed = TRUE)]

  # Section Q's country_summary is already the AGG_WEIGHTS-weighted
  # geometric mean (posterior_summary_fns.R::summarise_site_parameters()).
  forest_dt[, `:=`(
    halflife_gm_med  = halflife_median,
    halflife_gm_q025 = halflife_q025,
    halflife_gm_q975 = halflife_q975
  )]
  forest_dt <- forest_dt[order(halflife_gm_med)]

  # Country name only (2026-10-12: the site-count and posterior predictive
  # coverage lines, and the PPC footnote, were removed; coverage is reported
  # in the site table in results_numbers.Rmd instead).
  NAME_SPAN_SIZE <- 25 * SCALE_A4 * COMBINED_TEXT_SCALE
  forest_dt[, country_label_big := sprintf("<span style='font-size:%.2fpt'>%s</span>", NAME_SPAN_SIZE, country_full_name)]
  forest_dt[, country_label_big := factor(country_label_big, levels = country_label_big)]

  # Site-level posterior medians, coloured with the SAME minou mapping the
  # LHS scatter uses (keyed by iso) so a country reads as one colour across
  # both panels.
  scatter_iso_labels <- unique(mis_scatter[, .(iso, country_label)])
  country_colors_by_iso <- setNames(
    unname(minou_country_colors[as.character(scatter_iso_labels$country_label)]),
    scatter_iso_labels$iso
  )
  # Point size = the site's own AGG_WEIGHTS weight (par_pf: population at
  # risk, mean over its MIS survey years) -- the same weight that site
  # carries in its country's geometric mean, so bigger dots pull harder.
  site_forest_dt <- data.table::data.table(
    iso          = site_iso,
    halflife_med = apply(d_mat, 2, median) * log(2),
    site_weight  = unname(get_site_weights(post, AGG_WEIGHTS))
  )
  site_forest_dt[, country_label_big := forest_dt$country_label_big[match(iso, forest_dt$iso)]]

  # Country point size = that country's total weight (sum of its sites'
  # par_pf), area-proportional. Set directly rather than mapped: the plot's
  # one size scale belongs to the site dots (no ggnewscale for a second), and
  # country totals on that scale would dwarf every site.
  COUNTRY_POINT_SIZE_RANGE <- c(2.5, 9) * SCALE_A4
  country_weight <- site_forest_dt[, .(country_weight = sum(site_weight)), by = iso]
  forest_dt[, country_weight := country_weight$country_weight[match(iso, country_weight$iso)]]
  forest_dt[, country_point_size := pmax(COUNTRY_POINT_SIZE_RANGE[1],
                                         COUNTRY_POINT_SIZE_RANGE[2] * sqrt(country_weight / max(country_weight)))]
  site_forest_dt[, site_colour := data.table::fifelse(is.na(country_colors_by_iso[iso]), "#898781",
                                                      country_colors_by_iso[iso])]

  # Posterior predictive for an unobserved site -- one draw of log d_new per
  # retained posterior draw, propagating mu_d/sigma_d's own joint
  # uncertainty directly (no country layer under this model -- see
  # hierarchical_mcmc_fns.R's own header), same structure/truncation as the
  # prior simulator so the two bands are like-for-like.
  post_pred_halflife <- predict_unobserved_site(post, seed = 1L)$draws$halflife

  COMBINED_COLOR_PRIOR_BAND <- "#2e4057"
  COMBINED_COLOR_POST_BAND  <- "#edae49"
  # Prior drawn as an OUTLINE only and posterior as the only fill: two
  # overlapping semi-transparent fills blended into a third, unexplained
  # colour in the previous version.
  # Row order = drawing order: posterior fill first, prior outline on top
  # so the fill never hides the outline's edge (factor levels keep the
  # legend reading Prior, Posterior).
  band_combined_dt <- data.table::data.table(
    band = factor(c("Posterior", "Prior"), levels = c("Prior", "Posterior")),
    xmin = c(quantile(post_pred_halflife, 0.025), prior_pred_lo),
    xmax = c(quantile(post_pred_halflife, 0.975), prior_pred_hi)
  )
  # Population median half-life, exp(mu_d) * log(2) -- median, not mean,
  # to match the country/site summaries above.
  median_line_dt <- data.table::data.table(
    what = factor(c("Prior median", "Posterior median"), levels = c("Prior median", "Posterior median")),
    x    = c(median(priors$mu_d_natural), median(pooled_samples$mu_d_natural)) * log(2)
  )
  cat(sprintf(
    "Posterior predictive half-life, new site: median %.1fd, 95%% PI [%.1f, %.1f]d (prior [%.1f, %.1f]d)\n",
    median(post_pred_halflife), band_combined_dt[band == "Posterior", xmin], band_combined_dt[band == "Posterior", xmax],
    prior_pred_lo, prior_pred_hi
  ))

  # Log x-axis: decay is modelled on log(d), and on a linear 0-300 axis 11
  # of 14 countries were squeezed into 5-40 days.
  HALFLIFE_BREAKS <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)

  p_forest_natural_combined <- ggplot() +
    # Prior interval band and Prior median line both dropped (per request)
    # -- only the Posterior predictive band + Posterior median remain.
    # Legend removed entirely (legend.position = "none" below), so these
    # are plotted with fixed aesthetics rather than data-driven aes()/scale
    # mappings -- there is no longer a legend key to build.
    geom_rect(data = band_combined_dt[band == "Posterior"],
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = COMBINED_COLOR_POST_BAND, colour = NA, linewidth = 1.1 * SCALE_A4, alpha = 0.3) +
    geom_vline(xintercept = median_line_dt[what == "Posterior median", x],
               colour = COMBINED_COLOR_POST_BAND, linewidth = 1.1 * SCALE_A4, linetype = "dashed") +
    geom_point(data = site_forest_dt, aes(x = halflife_med, y = country_label_big, colour = site_colour, size = site_weight),
               position = position_jitter(width = 0, height = 0.18, seed = 1L), alpha = 0.5) +
    scale_colour_identity() +
    # Area proportional to weight; the 1.2 floor keeps the smallest admin1s
    # (~13k at risk vs ~12M for the largest) visible.
    scale_size(range = c(1.2, 9) * SCALE_A4,
               name = if (identical(AGG_WEIGHTS, "par_pf")) "Site population at risk" else "Site weight",
               labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    geom_segment(data = forest_dt,
                 aes(x = halflife_gm_q025, xend = halflife_gm_q975, y = country_label_big, yend = country_label_big),
                 colour = "black", linewidth = 1.5 * SCALE_A4) +
    geom_point(data = forest_dt, aes(x = halflife_gm_med, y = country_label_big), colour = "black",
               size = forest_dt$country_point_size) +
    scale_x_log10(breaks = HALFLIFE_BREAKS, labels = HALFLIFE_BREAKS) +
    # expand = expansion(add = 0.6): trims the near-full blank "row" a
    # discrete y-scale pads above the top and below the bottom country.
    scale_y_discrete(expand = expansion(add = 0.6)) +
    labs(
      title    = "Fitted HRP2 decay half-life",
      subtitle = "Black: country population-weighted geometric mean (median, 95% CrI).\nPoints: site posterior medians. Point size: population at risk.",
      x = "HRP2 half-life (days, log scale)", y = NULL
    ) +
    theme_bw(base_size = 23 * SCALE_A4 * COMBINED_TEXT_SCALE) +
    theme(
      # linewidth scaled from the LHS panel's 0.85 by the RHS:LHS column
      # width ratio (26.6/24 ~= 1.11) so the two boxes read equally bold --
      # RHS's column is only modestly wider now (it was proportionally much
      # wider in an earlier version of this figure, when 1.2 was tuned;
      # that ratio has since narrowed and 1.2 started reading heavier than
      # the LHS box, not matched to it).
      panel.border       = element_rect(colour = "#0b0b0b", fill = NA, linewidth = 0.95 * SCALE_A4),
      panel.background   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "#e1e0d9", linewidth = 0.42 * SCALE_A4),
      plot.title          = element_text(face = "bold", size = 30 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b"),
      plot.subtitle       = element_text(size = 21 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e", lineheight = 1.1),
      # plot.caption.position = "plot" (not the default "panel"): spans the
      # caption across the FULL plot width, including the y-axis area the
      # country labels sit in -- with the default "panel" position, hjust=0
      # alone still starts the caption at the PANEL's own left edge, which
      # sits well right of the country labels and doesn't read as aligned
      # with them.
      plot.caption        = element_text(size = 23 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e", hjust = 1),
      plot.caption.position = "plot",
      axis.text.y        = ggtext::element_markdown(size = 22 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b", lineheight = 1.0),
      axis.text.x        = element_text(colour = "#898781", size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE),
      axis.title.x       = element_text(colour = "#0b0b0b", size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE),
      axis.ticks         = element_blank(),
      # Only the size key shows (colour is scale_colour_identity(), no guide).
      legend.position    = "bottom",
      legend.title       = element_text(size = 28 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e"),
      legend.text        = element_text(size = 26 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e"),
      legend.key         = element_blank(),
      legend.key.size    = grid::unit(1.4, "lines"),
      legend.margin      = margin(t = 0, b = 0),
      plot.margin        = margin(t = 4 * SCALE_A4, r = 4 * SCALE_A4, b = 4 * SCALE_A4, l = 4 * SCALE_A4),
      # Fixed panel shape, like panel A's coord_equal: two fixed-aspect panels
      # in one patchwork row share the row height, so both panels are drawn
      # exactly the same height.
      aspect.ratio       = COMBINED_FOREST_ASPECT
    )

  # --- Map version of panel A: observed RDT:LM ratio by sub-region, each
  # country's most recent survey. Same title/subtitle/text sizes as the
  # scatter, and a square frame (fixed_frame) so the layout widths still fit.
  if (!is.null(ratio_map) && !ratio_forest) {
    # Boundaries with no survey row (e.g. northern Mali, north-eastern Kenya,
    # not sampled) fall back to the background; grey is kept for LM = 0
    # Only the countries in panel B (the fitted sites' countries)
    map_sf <- ratio_map$sf[!is.na(ratio_map$sf$prev_lm) &
                           ratio_map$sf$country %in% country_summary$country, ]
    missing_b <- setdiff(country_summary$country, map_sf$country)
    if (length(missing_b) > 0) warning("Panel B countries with no ratio map data: ", paste(missing_b, collapse = ", "))
    map_sf$ratio_plot <- pmin(pmax(map_sf$ratio, RATIO_MAP_LIMITS[1]), RATIO_MAP_LIMITS[2])
    ratio_breaks <- c(0.25, 0.5, 1, 2, 4)
    p_scatter_combined <- ggplot(map_sf) +
      # No-data fill is a warm sand: a light grey read as the pale teal end
      # of the colour bar
      geom_sf(data = ratio_map$bg, inherit.aes = FALSE, fill = RATIO_MAP_COL_NO_DATA, colour = NA) +
      geom_sf(aes(fill = ratio_plot), colour = NA) +
      geom_sf(data = ratio_map$bg, inherit.aes = FALSE, fill = NA, colour = "grey40",
              linewidth = 0.85 * SCALE_A4) +
      scale_fill_gradient2(
        low = RATIO_MAP_COL_LM, mid = "white", high = RATIO_MAP_COL_RDT, midpoint = 1,
        transform = "log2", limits = RATIO_MAP_LIMITS, breaks = ratio_breaks,
        labels = c("≤0.25", "0.5", "1", "2", "≥4"), na.value = "grey60",
        name = "RDT:LM prevalence ratio",
        guide = guide_colourbar(title.position = "top", title.hjust = 0.5)
      ) +
      fixed_frame(ratio_map$bg, aspect = 1) +
      labs(title = "Observed RDT:LM prevalence",
           subtitle = sprintf("Sub-region ratio, each country's most recent survey.%s",
                              if (any(map_sf$prev_lm == 0)) "\nGrey: LM prevalence 0%." else "")) +
      theme_void(base_size = 23 * SCALE_A4 * COMBINED_TEXT_SCALE) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        # Left margin clears the panel tag (theme_void has no axis area to do it)
        plot.title       = element_text(face = "bold", size = 30 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b",
                                        margin = margin(l = 16, b = 2)),
        plot.subtitle    = element_text(size = 21 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e", lineheight = 1.1,
                                        margin = margin(l = 16)),
        legend.position  = "bottom",
        legend.title     = element_text(size = 26 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e"),
        legend.text      = element_text(size = 24 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e"),
        legend.key.width = grid::unit(2, "cm"),
        legend.key.height = grid::unit(0.3, "cm"),
        legend.margin    = margin(t = 0, b = 0),
        plot.margin      = margin(t = 4 * SCALE_A4, r = 4 * SCALE_A4, b = 4 * SCALE_A4, l = 4 * SCALE_A4)
      )
  }

  # --- Ratio-forest version of panel A: observed RDT:LM by country, laid out
  # like panel B (same rows, order, colours, theme and country point sizes).
  # Points: each site-year's ratio (mis_scatter: all surveys in the fit,
  # deduplicated as for the scatter). Black: the country RDT:LM scalar, the
  # slope of RDT on LM through the origin (the scatter version's bracketed
  # number), with its 95% CI.
  if (ratio_forest) {
    ctry_levels <- levels(forest_dt$country_label_big)
    rf_pts <- mis_scatter[pfpr_lm > 0, .(iso, ratio = pfpr_rdt / pfpr_lm)]
    rf_pts[, ratio_plot := pmin(pmax(ratio, RATIO_FOREST_LIMITS[1]), RATIO_FOREST_LIMITS[2])]
    rf_pts[, country_label_big := factor(forest_dt$country_label_big[match(iso, forest_dt$iso)], levels = ctry_levels)]
    rf_pts[, site_colour := data.table::fifelse(is.na(country_colors_by_iso[iso]), "#898781", country_colors_by_iso[iso])]
    n_clamped <- rf_pts[ratio != ratio_plot, .N]
    rf_ctry <- mis_scatter[, {
      f <- lm(pfpr_rdt ~ 0 + pfpr_lm)
      ci <- suppressWarnings(confint(f))
      .(slope = coef(f)[[1]], lo = ci[1, 1], hi = ci[1, 2])
    }, by = iso]
    rf_ctry[, country_label_big := factor(forest_dt$country_label_big[match(iso, forest_dt$iso)], levels = ctry_levels)]
    rf_ctry[, point_size := forest_dt$country_point_size[match(iso, forest_dt$iso)]]
    rf_pts  <- rf_pts[!is.na(country_label_big)]
    rf_ctry <- rf_ctry[!is.na(country_label_big)]
    cat(sprintf("Ratio forest: %d site-years (%d with LM = 0 dropped, %d clamped to the axis)\n",
                nrow(rf_pts), mis_scatter[pfpr_lm == 0, .N], n_clamped))

    p_scatter_combined <- ggplot() +
      geom_vline(xintercept = 1, colour = "#898781", linewidth = 1.1 * SCALE_A4, linetype = "dashed") +
      geom_point(data = rf_pts, aes(x = ratio_plot, y = country_label_big, colour = site_colour),
                 position = position_jitter(width = 0, height = 0.18, seed = 1L), alpha = 0.5,
                 size = 2.6 * SCALE_A4) +
      scale_colour_identity() +
      geom_segment(data = rf_ctry, aes(x = lo, xend = hi, y = country_label_big, yend = country_label_big),
                   colour = "black", linewidth = 1.5 * SCALE_A4) +
      geom_point(data = rf_ctry, aes(x = slope, y = country_label_big), colour = "black",
                 size = rf_ctry$point_size) +
      scale_x_log10(limits = RATIO_FOREST_LIMITS, breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
                    labels = c("≤1/8", "1/4", "1/2", "1", "2", "4", "8", "≥16")) +
      scale_y_discrete(drop = FALSE, expand = expansion(add = 0.6)) +
      labs(title = "Observed RDT:LM prevalence",
           subtitle = "Black: country RDT:LM scalar (slope through the origin, 95% CI).\nPoints: site-year ratios, all surveys.",
           x = "RDT:LM prevalence ratio (log scale)", y = NULL) +
      p_forest_natural_combined$theme +
      theme(legend.position = "none")
  }

  # --- Ratio density panel (ratio_density = TRUE): replaces panel A, so the
  # figure is density (A) | half-life forest (B). Density of the site-year
  # RDT:LM ratios (as in the ratio forest) on the log scale, one ridge per
  # country in panel B's rows. Drawn by hand on a numeric y axis (row i spans
  # i - 0.5 to i + 0.5, the discrete axis positions) so each ridge stays
  # inside its own row and the rows line up with panel B's; ggridges' ridges
  # rise above their row and would overflow the top one.
  # Points: each fitted site's observed ratio (mean RDT over mean LM
  # prevalence across its survey years), sized by population at risk on the
  # same scale as panel B's site points. The site-year ratios behind the
  # density are deduplicated on value (build_mis_scatter()), so they can't be
  # tied to one site's population; the points are per site instead.
  if (ratio_density) {
    ctry_levels <- levels(forest_dt$country_label_big)
    dens_pts <- mis_scatter[pfpr_lm > 0, .(iso, ratio = pfpr_rdt / pfpr_lm)]
    dens_pts[, row := match(iso, forest_dt$iso)]
    dens_pts <- dens_pts[!is.na(row)]
    n_clamped <- dens_pts[ratio < RATIO_FOREST_LIMITS[1] | ratio > RATIO_FOREST_LIMITS[2], .N]
    dens_pts[, log_ratio := log10(pmin(pmax(ratio, RATIO_FOREST_LIMITS[1]), RATIO_FOREST_LIMITS[2]))]
    # Densities cut at the axis ends; heights scaled to the tallest ridge so
    # countries compare in shape, the tallest filling 90% of its row
    dens_dt <- dens_pts[, {
      d <- stats::density(log_ratio, from = log10(RATIO_FOREST_LIMITS[1]), to = log10(RATIO_FOREST_LIMITS[2]), n = 256)
      .(x = 10^d$x, y = d$y, n = .N)
    }, by = .(iso, row)]
    dens_dt[, y := 0.9 * y / max(y)]
    dens_dt[, fill := data.table::fifelse(is.na(country_colors_by_iso[iso]), "#898781", country_colors_by_iso[iso])]
    dens_med <- dens_pts[, .(x = 10^median(log_ratio)), by = .(iso, row)]

    site_obs <- data.table::as.data.table(readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS")))[
      !is.na(iso3c), .(site = sprintf("%s_%s", iso3c, name_1), year,
                       pfpr_rdt = pfpr_rdt_num / pfpr_rdt_denom, pfpr_lm = pfpr_lm_num / pfpr_lm_denom)]
    site_obs <- unique(site_obs[!is.na(pfpr_rdt) & !is.na(pfpr_lm)])
    site_ratio <- site_obs[, .(ratio = mean(pfpr_rdt) / mean(pfpr_lm)), by = site]
    site_dens_dt <- data.table::copy(site_forest_dt)[, site := site_names]
    site_dens_dt <- merge(site_dens_dt, site_ratio, by = "site")[is.finite(ratio)]
    site_dens_dt[, ratio_plot := pmin(pmax(ratio, RATIO_FOREST_LIMITS[1]), RATIO_FOREST_LIMITS[2])]
    site_dens_dt[, row := match(iso, forest_dt$iso)]
    cat(sprintf(paste0("Ratio density: %d site-years in %d countries (%d with LM = 0 dropped, %d clamped to the axis); ",
                       "fewest per country %d. Site points: %d of %d sites (%d clamped)\n"),
                nrow(dens_pts), uniqueN(dens_pts$iso), mis_scatter[pfpr_lm == 0, .N], n_clamped,
                min(dens_pts[, .N, by = iso]$N), nrow(site_dens_dt), length(site_names),
                site_dens_dt[ratio != ratio_plot, .N]))

    p_scatter_combined <- ggplot() +
      geom_vline(xintercept = 1, colour = "#898781", linewidth = 1.1 * SCALE_A4, linetype = "dashed") +
      geom_ribbon(data = dens_dt, aes(x = x, ymin = row - 0.45, ymax = row - 0.45 + y, group = iso, fill = fill),
                  colour = "#0b0b0b", linewidth = 0.6 * SCALE_A4, alpha = 0.3) +
      geom_point(data = site_dens_dt, aes(x = ratio_plot, y = row - 0.3, colour = site_colour, size = site_weight),
                 position = position_jitter(width = 0, height = 0.08, seed = 1L), alpha = 0.5) +   # near the row bottom; alpha as panel B, so the size legends merge
      geom_segment(data = dens_med, aes(x = x, xend = x, y = row - 0.45, yend = row + 0.45),
                   colour = "black", linewidth = 1.5 * SCALE_A4) +
      scale_fill_identity() +
      scale_colour_identity() +
      # Same size scale as panel B, so patchwork collects one legend
      scale_size(range = c(1.2, 9) * SCALE_A4, limits = range(site_forest_dt$site_weight),
                 name = if (identical(AGG_WEIGHTS, "par_pf")) "Population at risk" else "Weight",
                 labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      scale_x_log10(limits = RATIO_FOREST_LIMITS, breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
                    labels = c("≤1/8", "1/4", "1/2", "1", "2", "4", "8", "≥16")) +
      # Same span as panel B's discrete axis (rows 1..n plus 0.6 either side),
      # labelled with panel B's country names
      scale_y_continuous(limits = c(0.4, length(ctry_levels) + 0.6), expand = c(0, 0),
                         breaks = seq_along(ctry_levels), labels = ctry_levels) +
      labs(title = "DHS/MIS observed RDT:LM ratio",
           subtitle = "Ridges: density of site-year ratios, all surveys. Black line: country median.\nPoints: site ratios across survey years.",
           x = "RDT:LM prevalence ratio (log scale)", y = NULL) +
      p_forest_natural_combined$theme +
      theme(plot.margin = margin(t = 4 * SCALE_A4, r = COMBINED_DENSITY_GAP_PT, b = 4 * SCALE_A4, l = 4 * SCALE_A4))
    # Country names are on panel A; B drops its own, so its title needs a left
    # margin to clear the panel tag
    p_forest_natural_combined <- p_forest_natural_combined +
      labs(subtitle = "Black: country population-weighted geometric mean (median, 95% CrI).\nPoints: site posterior medians.") +
      theme(
      axis.text.y   = element_blank(),
      plot.title    = element_text(face = "bold", size = 30 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#0b0b0b",
                                   margin = margin(l = 16, b = 2)),
      plot.subtitle = element_text(size = 21 * SCALE_A4 * COMBINED_TEXT_SCALE, colour = "#52514e", lineheight = 1.1,
                                   margin = margin(l = 16)))
    # Both sizes keyed to the full weight range, so the two legends merge
    p_forest_natural_combined <- p_forest_natural_combined +
      scale_size(range = c(1.2, 9) * SCALE_A4, limits = range(site_forest_dt$site_weight),
                 name = if (identical(AGG_WEIGHTS, "par_pf")) "Population at risk" else "Weight",
                 labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  }

  # RHS (the paper's main result) gets the larger width share -- LHS's
  # column is sized to just match the figure height (its coord_equal panel
  # is forced square), everything left over goes to RHS. Panel tags A/B for
  # referencing from the text.
  # Ratio forest: no fixed panel shape and B's size legend collected below
  # both panels, so patchwork aligns the two panels and their country rows
  # line up; equal widths since both carry the country labels.
  # Density version: both panels square and the same size (patchwork widths
  # are panel widths); the figure is sized to fit them so the fixed shape
  # leaves no blank sides.
  if (ratio_forest) {
    p_scatter_combined        <- p_scatter_combined + theme(aspect.ratio = NULL)
    p_forest_natural_combined <- p_forest_natural_combined + theme(aspect.ratio = NULL)
  }
  if (ratio_density) {
    p_scatter_combined        <- p_scatter_combined + theme(aspect.ratio = 1)
    p_forest_natural_combined <- p_forest_natural_combined + theme(aspect.ratio = 1)
  }
  p_combined_scatter_forest <- (p_scatter_combined | p_forest_natural_combined) +
    patchwork::plot_layout(widths = if (ratio_density) COMBINED_DENSITY_PANEL_WIDTHS else if (ratio_forest) c(1, 1) else COMBINED_PANEL_WIDTHS,
                           guides = if (ratio_forest || ratio_density) "collect" else "auto")
  p_combined_scatter_forest <- p_combined_scatter_forest +
    patchwork::plot_annotation(tag_levels = "A") &
    # plot.tag.location = "plot": the tag sits in each subplot's top-left
    # corner, level with its title. (Anchored to the panel corner, as before
    # 2026-10-12, it overlapped panel A's subtitle.)
    theme(plot.tag = element_text(size = 36 * SCALE_A4 * COMBINED_TEXT_SCALE, face = "bold"),
          plot.tag.location = "plot",
          plot.tag.position = c(0, 1))
  if (ratio_forest || ratio_density) p_combined_scatter_forest <- p_combined_scatter_forest & theme(legend.position = "bottom")
  # Density version: legend tucked up under the axis titles
  if (ratio_density) p_combined_scatter_forest <- p_combined_scatter_forest & theme(legend.box.spacing = grid::unit(2, "pt"))
  p_combined_scatter_forest
}

# --- Driver -------------------------------------------------------------------
if (PLOT_DENSITY_FOREST) {
  mis_scatter <- build_mis_scatter()
  mcmc_dir <- fit_dir(RDT_SENSITIVITY_MAIN, SIGMA_D_PRIOR_SD_MAIN)
  if (!file.exists(file.path(mcmc_dir, "mcmc_combined.RDS"))) {
    cat(sprintf("No mcmc_combined.RDS in %s -- combined_ratio_density_and_forest_plot.png skipped.\n", mcmc_dir))
  } else {
    fi <- build_forest_inputs(mcmc_dir)
    ggsave(file.path(FIGURE_DIR, "combined_ratio_density_and_forest_plot.png"),
           plot_combined_scatter_forest(fi, mis_scatter, ratio_density = TRUE),
           width = COMBINED_DENSITY_WIDTH_IN, height = COMBINED_DENSITY_HEIGHT_IN, dpi = 350, limitsize = FALSE, bg = "white")
    cat(sprintf("Saved %s\n", file.path(FIGURE_DIR, "combined_ratio_density_and_forest_plot.png")))
  }
}

# =============================================================================
# 3. POSTERIOR PREDICTIVE COVERAGE: OBSERVED VS PREDICTED PREVALENCE, THIS
#    MODEL VS WU ET AL. (steps 1c, 4; main fit)
# =============================================================================
# Observed: per site and survey year, RDT (TP + FP) and LM (TP + FN)
# prevalence among children tested by both (sites sharing a DHS region source
# row show the same values).
# This model: per posterior draw, site RDT (q_tp + q_fp, deletion and RDT
# sensitivity applied) and LM (p_tp + p_fn); 95% posterior predictive
# intervals, each draw resampled as a binomial on the site's effective size
# (joint_n / design effect; the site's mean in years without a survey), as in
# step 4.
# Wu et al.: per draw, RDT predicted from the same site LM draws with the Wu
# et al. line (wu_rdt_from_lm() in cluster_fns.R). With WU_UNCERTAINTY =
# "lm_only" its band is the spread of our LM draws through the line only;
# "loco_sd" adds a survey-level deviation and binomial resampling.
PREV_MEASURE_COLOURS <- c(RDT = "#0F86B6", LM = "#E1AF00")
PREV_MODEL_LINETYPES <- c("This analysis" = "solid", "Wu et al" = "dashed")
PREV_MODEL_ALPHA     <- c("This analysis" = 0.12, "Wu et al" = 0.08)

# Site-level observed and predicted prevalence (used for the coverage intervals
# when step 4's summary lacks the Section I marginal columns).
prevalence_wu_site_data <- function(mcmc_dir, rdt_sensitivity, mis, seed = 2015L) {
  pp_dir  <- file.path(mcmc_dir, "posterior_predictive")
  results <- readRDS(file.path(pp_dir, "posterior_predictive_results.RDS"))
  results <- merge(results, readRDS(file.path(pp_dir, "draw_lookup.RDS")),
                   by = c("iso", "name_1", "eir", "hrp2_scale"), allow.cartesian = TRUE)
  results[, c("q_tp", "q_fp", "q_fn", "q_tn") := joint_cells_observed(p_tp, p_fp, p_fn, p_tn, deletion_pct, rdt_sensitivity)]
  results[, `:=`(rdt = q_tp + q_fp, lm = p_tp + p_fn)]

  mis <- data.table::as.data.table(mis)
  obs <- unique(mis[!is.na(joint_n) & joint_n > 0,
    .(iso = iso3c, country, name_1, year, n = joint_n, RDT = joint_tp + joint_fp, LM = joint_tp + joint_fn,
      n_eff = joint_n / data.table::fifelse(is.na(deff_joint), 1, deff_joint))])
  obs <- obs[unique(results[, .(iso, name_1)]), on = c("iso", "name_1"), nomatch = NULL]
  site_n_eff <- obs[, .(iso, name_1, year, n_eff)]
  obs <- data.table::melt(obs[, !"n_eff"], id.vars = c("iso", "country", "name_1", "year", "n"),
                          measure.vars = c("RDT", "LM"), variable.name = "measure", value.name = "k")
  obs[, `:=`(prev = k / n, lo = stats::qbeta(0.025, k, n - k + 1), hi = stats::qbeta(0.975, k + 1, n - k))]

  wu_sd <- wu_residual_sd_loco(unique(mis[!is.na(joint_n) & joint_n > 0,
    .(iso = iso3c, name_1, year, rdt = (joint_tp + joint_fp) / joint_n, lm = (joint_tp + joint_fn) / joint_n)]))
  pred <- merge(results[, .(iso, name_1, year, draw_id, rdt, lm)], wu_sd[, .(iso, wu_sd)], by = "iso")
  set.seed(seed)
  pred[, rdt_wu := wu_rdt_from_lm(lm, wu_sd)]
  pred <- data.table::rbindlist(list(
    pred[, .(iso, name_1, year, draw_id, measure = "RDT", model = "This analysis", prev = rdt)],
    pred[, .(iso, name_1, year, draw_id, measure = "LM",  model = "This analysis", prev = lm)],
    pred[, .(iso, name_1, year, draw_id, measure = "RDT", model = "Wu et al",  prev = rdt_wu)]))
  pred <- merge(pred, site_n_eff, by = c("iso", "name_1", "year"), all.x = TRUE)
  pred <- merge(pred, site_n_eff[, .(n_eff_avg = mean(n_eff)), by = .(iso, name_1)], by = c("iso", "name_1"))
  pred[, n_eff := max(1, round(data.table::fifelse(is.na(n_eff), n_eff_avg, n_eff))), by = .(iso, name_1, year)]
  pred <- pred[, {
    resampled <- if (model[1] == "Wu et al" && identical(WU_UNCERTAINTY, "lm_only")) prev else
      stats::rbinom(.N, size = n_eff[1], prob = prev) / n_eff[1]
    .(mean = mean(prev), lo = stats::quantile(resampled, 0.025), hi = stats::quantile(resampled, 0.975))
  }, by = .(iso, name_1, year, measure, model)]

  yr_range <- c(min(obs$year) - 2L, max(obs$year) + 2L)
  pred <- pred[year >= yr_range[1] & year <= yr_range[2]]
  country_names <- unique(obs[, .(iso, country)])
  country_names[, country := gsub("Cote d'Ivoire", "Côte d'Ivoire", country, fixed = TRUE)]
  for (dt in list(obs, pred)) {
    dt[, measure := factor(measure, levels = names(PREV_MEASURE_COLOURS))]
    dt[country_names, country := i.country, on = "iso"]
  }
  pred[, model := factor(model, levels = names(PREV_MODEL_LINETYPES))]
  list(obs = obs, pred = pred, yr_range = yr_range,
       y_max = max(c(pred$hi, obs$hi), na.rm = TRUE))
}


# Coverage of the 95% interval per (site, survey year) for RDT (this model),
# LM (this model) and RDT (Wu et al.). Uses step 4's
# posterior_predictive_summary.RDS when it has the Section I marginal columns,
# so the counts match results_numbers.Rmd; otherwise the site intervals built
# for the site figure (same construction, different random draws).
prevalence_wu_coverage <- function(mcmc_dir, mis, site_data) {
  obs <- unique(data.table::as.data.table(mis)[!is.na(joint_n) & joint_n > 0,
    .(iso = iso3c, name_1, year, rdt_obs = (joint_tp + joint_fp) / joint_n, lm_obs = (joint_tp + joint_fn) / joint_n)])
  s_path <- file.path(mcmc_dir, "posterior_predictive", "posterior_predictive_summary.RDS")
  s <- if (file.exists(s_path)) data.table::as.data.table(readRDS(s_path)) else NULL
  if (!is.null(s) && all(c("rdtj_lo", "lmj_lo", "rdtwu_lo") %in% names(s))) {
    iv <- s[, .(iso, name_1, year, rdt_mean = rdtj_mean, rdt_lo = rdtj_lo, rdt_hi = rdtj_hi,
                lm_mean = lmj_mean, lm_lo = lmj_lo, lm_hi = lmj_hi,
                wu_mean = rdtwu_mean, wu_lo = rdtwu_lo, wu_hi = rdtwu_hi)]
    src <- "step 4 summary"
  } else {
    p <- site_data$pred
    one <- function(m, mo, nm) {
      x <- p[measure == m & model == mo, .(iso, name_1, year, mean, lo, hi)]
      data.table::setnames(x, c("mean", "lo", "hi"), paste0(nm, c("_mean", "_lo", "_hi")))
    }
    iv <- Reduce(function(a, b) merge(a, b, by = c("iso", "name_1", "year")),
                 list(one("RDT", "This analysis", "rdt"), one("LM", "This analysis", "lm"), one("RDT", "Wu et al", "wu")))
    src <- "site figure intervals -- rerun step 4 Section I (AGG_POST = TRUE) for counts matching results_numbers.Rmd"
  }
  cv <- merge(obs, iv, by = c("iso", "name_1", "year"))
  cv[, `:=`(rdt_in = rdt_obs >= rdt_lo & rdt_obs <= rdt_hi,
            lm_in  = lm_obs  >= lm_lo  & lm_obs  <= lm_hi,
            wu_in  = rdt_obs >= wu_lo  & rdt_obs <= wu_hi)]
  cat(sprintf("Coverage from %s: %d site-years.\n", src, nrow(cv)))
  cv
}
ppc_kn <- function(x) sprintf("%s (%d/%d)", scales::percent(mean(x), accuracy = 1), sum(x), length(x))
ppc_panel_labels <- function(cv, by) {
  cv[, .(label = sprintf("LM: %s\nRDT: %s\nRDT (Wu et al): %s", ppc_kn(lm_in), ppc_kn(rdt_in), ppc_kn(wu_in))), by = by]
}
ppc_subtitle <- function(cv, where) {
  sprintf("RDT coverage%s: this analysis %s; Wu et al %s", where, ppc_kn(cv$rdt_in), ppc_kn(cv$wu_in))
}
# Top-right corner of each panel, on a translucent white box
ppc_label_layer <- function(labs_cov, size) {
  geom_label(data = labs_cov, aes(x = Inf, y = Inf, label = label), hjust = 1, vjust = 1, size = size,
             lineheight = 0.9, label.size = 0, fill = scales::alpha("white", 0.75), inherit.aes = FALSE)
}

# Observed MIS prevalence (x) vs predicted mean and 95% interval (y) per
# site-year, from prevalence_wu_coverage(): one panel each for LM (this
# analysis), RDT (this analysis) and RDT (Wu et al). Styled as
# plot_mcmc_diagnostics.R::plot_joint_ppi_all_sites(): points inside their
# interval blue, outside gold, dashed y = x, coverage per panel.
# Every panel shares the same x and y range.
plot_prevalence_wu_scatter <- function(cv) {
  panels <- c(lm = "LM prevalence (this analysis)", rdt = "RDT prevalence (this analysis)",
              wu = "RDT prevalence (Wu et al.)")
  obs_col <- c(lm = "lm_obs", rdt = "rdt_obs", wu = "rdt_obs")
  d <- data.table::rbindlist(lapply(names(panels), function(k) cv[, .(
    panel = panels[[k]], observed = get(obs_col[[k]]), mean = get(paste0(k, "_mean")),
    lo = get(paste0(k, "_lo")), hi = get(paste0(k, "_hi")), inside = get(paste0(k, "_in")))]))
  d[, panel := factor(panel, levels = panels)]
  d[, status := factor(data.table::fifelse(inside, "Inside 95% interval", "Outside 95% interval"),
                       levels = c("Inside 95% interval", "Outside 95% interval"))]
  lim <- c(0, max(c(d$observed, d$hi), na.rm = TRUE))
  # Per panel: coverage of the 95% interval
  cov <- d[, .(label = sprintf(" Coverage %s%%", formatC(100 * mean(inside), format = "f", digits = 0))), by = panel]
  ggplot(d, aes(x = observed, y = mean)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = lo, ymax = hi, colour = status), alpha = 0.35, linewidth = 0.3) +
    geom_point(aes(colour = status), size = 0.9, alpha = 0.7) +
    geom_text(data = cov, aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1.15,
              size = 3.6, lineheight = 0.95, inherit.aes = FALSE) +
    facet_wrap(~ panel, nrow = 1) +
    coord_equal(xlim = lim, ylim = lim) +
    scale_colour_manual(values = c("Inside 95% interval" = "#156082", "Outside 95% interval" = "#E1AF00"), name = NULL) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1, suffix = "")) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1, suffix = "")) +
    labs(x = "Observed prevalence, MIS/DHS (6mo-5y, %)",
         y = "Predicted prevalence (6mo-5y, %)\n(mean and 95% interval)") +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(size = rel(0.9)), legend.position = "bottom",
          panel.grid.minor = element_blank())
}


if (PLOT_PREVALENCE_WU_SCATTER) {
  mis       <- readRDS(file.path(DIR_STRUCTURE, "matched_site_mis_prevalence.RDS"))
  mcmc_dir  <- fit_dir(RDT_SENSITIVITY_MAIN, SIGMA_D_PRIOR_SD_MAIN)
  fit_path  <- file.path(mcmc_dir, "mcmc_combined.RDS")
  pp_path   <- file.path(mcmc_dir, "posterior_predictive", "posterior_predictive_results.RDS")
  if (!file.exists(fit_path) || !file.exists(pp_path) || file.mtime(pp_path) < file.mtime(fit_path)) {
    cat(sprintf("prevalence_rdt_lm_wu_scatter: no current posterior predictive results in %s -- skipped.\n", mcmc_dir))
  } else {
    site_data <- prevalence_wu_site_data(mcmc_dir, RDT_SENSITIVITY_MAIN, mis)
    coverage  <- prevalence_wu_coverage(mcmc_dir, mis, site_data)
    ggsave(file.path(FIGURE_DIR, "prevalence_rdt_lm_wu_scatter.png"), plot_prevalence_wu_scatter(coverage),
           width = 9.4, height = 4, dpi = 600, bg = "white")
    cat(sprintf("Saved %s\n", file.path(FIGURE_DIR, "prevalence_rdt_lm_wu_scatter.png")))
  }
}
