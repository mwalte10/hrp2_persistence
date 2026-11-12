# =============================================================================
# src/functions/data_stream_fns.R
#
# ANALYSIS COMPONENT: OUTPUTS (function library: appendix data stream figure
# and its numbers)
#
# Counts the DHS/MIS surveys at each step from all surveys measuring malaria
# prevalence to the regression and calibration data sets, using the same
# rules as the pipeline:
#   - survey_test_years(): the survey-year table step 1c uses to choose the
#     countries with more than one paired survey (countries_included.RDS)
#   - regression: the surveys in the fitted step 2 model
#   - calibration: the surveys with paired site-level RDT and LM prevalence in
#     the countries of the fitted step 3c sites
# Used by src/plots.R (figure) and src/results_numbers.Rmd (data.stream.* keys),
# so the figure and the text always agree.
# =============================================================================

# Country-years with RDT and/or LM prevalence, one row each, with the same
# exclusions as step 1c section 1.1 (single-child rows; the 2008 Senegal
# survey, which has no paired testing in the household survey file).
survey_test_years <- function(comb) {
  comb <- data.table::as.data.table(comb)[n_denominator != 1]
  comb <- comb[!(survey_year == 2008 & country == "Senegal")]
  yi <- unique(comb[indicator %in% c("rdt", "microscopy"), .(country, survey_year, indicator)])
  yi[, .(has_rdt = "rdt" %in% indicator, has_micro = "microscopy" %in% indicator), by = .(country, survey_year)]
}

# Every count in the data stream figure.
#   pull_dir, structure_dir: step 1 output folders (DIR_PULL, DIR_STRUCTURE)
#   reg_data:   the data the step 2 model was fitted to (fit$data)
#   site_names: fitted calibration sites, "<ISO>_<name_1>" (mcmc_inputs.RDS)
data_stream_counts <- function(pull_dir, structure_dir, reg_data, site_names) {
  sy <- survey_test_years(data.table::fread(file.path(pull_dir, "dhs_mis_combined.csv")))
  paired <- sy[has_rdt & has_micro]
  per_country <- paired[, .N, by = country]
  keep <- per_country[N > 1, country]
  multi <- paired[country %in% keep]

  # Regression: surveys in the fitted model
  reg_surveys <- unique(as.character(reg_data$survid))
  reg_years <- as.integer(sub("^.*_", "", reg_surveys))
  multi_id <- paste0(multi$country, "_", multi$survey_year)

  # Calibration: surveys with paired site-level RDT and LM prevalence in the
  # fitted sites' countries
  mis <- data.table::as.data.table(readRDS(file.path(structure_dir, "matched_site_mis_prevalence.RDS")))
  cal_isos <- unique(sub("_.*$", "", site_names))
  cal <- unique(mis[iso3c %in% cal_isos & !is.na(pfpr_rdt_num) & !is.na(pfpr_lm_num), .(country, survey_year = year)])
  cal_id <- paste0(cal$country, "_", cal$survey_year)
  cal_excl <- multi[!paste0(country, "_", survey_year) %in% cal_id]

  # Children tested, from the national survey-weighted estimates (6-59 months)
  nat_n <- function(file) {
    d <- data.table::fread(file.path(pull_dir, file))[tolower(sub_region) == "national"]
    sum(d[cal, on = .(country, survey_year), nomatch = 0]$n_tested)
  }

  list(
    n_any            = nrow(sy),
    n_rdt_only       = sy[has_rdt & !has_micro, .N],
    n_lm_only        = sy[!has_rdt & has_micro, .N],
    n_paired         = nrow(paired),
    n_paired_countries = data.table::uniqueN(paired$country),
    n_single_surveys = per_country[N == 1, sum(N)],
    n_single_countries = per_country[N == 1, .N],
    n_multi          = nrow(multi),
    n_multi_countries = length(keep),
    n_reg            = length(reg_surveys),
    n_reg_excluded   = sum(!multi_id %in% reg_surveys),
    n_reg_countries  = data.table::uniqueN(sub("_[0-9]+$", "", reg_surveys)),
    reg_years        = range(reg_years),
    n_reg_children   = nrow(reg_data),
    n_cal            = nrow(cal),
    n_cal_countries  = data.table::uniqueN(cal$country),
    cal_years        = range(cal$survey_year),
    n_cal_excluded   = nrow(cal_excl),
    cal_excl_countries = sort(unique(cal_excl$country)),
    n_cal_sites      = length(site_names),
    n_cal_lm         = nat_n("dhs_mis_microscopy.csv"),
    n_cal_rdt        = nat_n("dhs_mis_rdt.csv")
  )
}

# Flow diagram of the counts (appendix fig:data_stream).
plot_data_stream <- function(k) {
  int <- function(x) formatC(x, format = "d", big.mark = ",")
  and_list <- function(x) if (length(x) <= 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
  W <- 3.4; WS <- 3.0
  boxes <- data.table::data.table(
    id = c("any", "rdt_only", "paired", "single", "multi", "reg_excl", "cal_excl", "reg", "cal"),
    x  = c(0, 4.2, 0, 4.2, 0, -4.2, 4.2, -2.2, 2.2),
    y  = c(10, 8.75, 7.5, 6.25, 5, 3.4, 3.4, 0.9, 0.9),
    w  = c(W, WS, W, WS, W, WS, WS, 3.8, 3.8),
    h  = c(1, 0.8, 1, 0.8, 1, 0.8, 0.8, 2.2, 2.2),
    kind = c("main", "excl", "main", "excl", "main", "excl", "excl", "reg", "cal"),
    label = c(
      sprintf("%s DHS/MIS surveys measuring malaria prevalence by RDT or LM", int(k$n_any)),
      sprintf("%s surveys with RDT only%s", int(k$n_rdt_only),
              if (k$n_lm_only > 0) sprintf(" and %s with LM only", int(k$n_lm_only)) else ""),
      sprintf("%s surveys with paired RDT and LM results (%s countries)", int(k$n_paired), int(k$n_paired_countries)),
      sprintf("%s surveys from %s countries with only one paired survey", int(k$n_single_surveys), int(k$n_single_countries)),
      sprintf("%s surveys from %s countries with more than one paired survey", int(k$n_multi), int(k$n_multi_countries)),
      sprintf("%s surveys without all individual-level determinants", int(k$n_reg_excluded)),
      sprintf("%s surveys from %s, which report results by ecological zone", int(k$n_cal_excluded), and_list(k$cal_excl_countries)),
      sprintf("Regression analysis\n%s surveys from %s countries, %d-%d\n%s children with RDT and LM results and all individual-level determinants",
              int(k$n_reg), int(k$n_reg_countries), k$reg_years[1], k$reg_years[2], int(k$n_reg_children)),
      sprintf("Calibration analysis\n%s surveys from %s countries, %d-%d, %s admin-1 regions\n%s children tested by LM\n%s children tested by RDT",
              int(k$n_cal), int(k$n_cal_countries), k$cal_years[1], k$cal_years[2], int(k$n_cal_sites),
              int(k$n_cal_lm), int(k$n_cal_rdt)))
  )
  boxes[, label := vapply(label, function(s) paste(vapply(strsplit(s, "\n")[[1]],
                    function(l) paste(strwrap(l, width = round(w * 11)), collapse = "\n"), ""), collapse = "\n"), ""),
        by = id]
  boxes <- boxes[!(id == "reg_excl" & k$n_reg_excluded == 0)]
  fills <- c(main = "#DDE2E8", excl = "white", reg = "#E0ECE5", cal = "#F6DBDF")
  edges <- c(main = "#2E4057", excl = "#8D96A3", reg = "#66A182", cal = "#D1495B")
  # Main line down the centre, side exclusions at right angles, and a split
  # into the regression (left) and calibration (right) branches
  lines <- data.table::data.table(x = c(0, -2.2, 0), y = c(4.5, 4.1, 4.5), xend = c(0, 2.2, 0), yend = c(4.1, 4.1, 4.1))
  arrows <- data.table::data.table(
    x    = c(0, 0, 0, 0, -2.2, 2.2, -2.2, 2.2),
    y    = c(9.5, 8.75, 7, 6.25, 4.1, 4.1, 3.4, 3.4),
    xend = c(0, 4.2 - WS / 2, 0, 4.2 - WS / 2, -2.2, 2.2, -4.2 + WS / 2, 4.2 - WS / 2),
    yend = c(8, 8.75, 5.5, 6.25, 2, 2, 3.4, 3.4))
  if (k$n_reg_excluded == 0) arrows <- arrows[-7]
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = lines, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = "#2E4057", linewidth = 0.5) +
    ggplot2::geom_segment(data = arrows, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.18, "cm"), type = "closed"),
                          colour = "#2E4057", linewidth = 0.5) +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                                                  fill = kind, colour = kind), linewidth = 0.7) +
    ggplot2::geom_text(data = boxes, ggplot2::aes(x = x, y = y, label = label), size = 3.3, lineheight = 0.95,
                       colour = "#232333") +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::scale_colour_manual(values = edges, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-6, 6), ylim = c(-0.3, 10.6), expand = FALSE) +
    ggplot2::theme_void()
}
