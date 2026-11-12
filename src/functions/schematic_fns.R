# =============================================================================
# src/functions/schematic_fns.R
#
# ANALYSIS COMPONENT: OUTPUTS (function library: manuscript Figure 1)
#
# Illustrative schematic of detectability by light microscopy (LM) and rapid
# diagnostic test (RDT) after infection (main text fig:dynamics). The figure
# is conceptual: no survey data or model fits are used. All three panels are
# drawn from ONE artificial data set, simulate_schematic_cohort(), so they
# are internally consistent:
#   A  the illustrative individual is the cohort member with no timing offset;
#      its parasitaemia and circulating HRP2 curves, the LM and RDT detection
#      limits and the event times all come from the same parameters
#   B  the share of the cohort detectable by LM and by RDT at each day
#   C  the share of the cohort in each joint LM x RDT outcome at each day
#
# HRP2 clearance follows the model structure used in the calibration
# (malariasimulation fork): after clinical clearance an individual stays
# HRP2-positive for an exponentially distributed time with mean delta. The
# illustrative individual is the average cohort member, so its circulating
# HRP2 is scaled to fall below the RDT detection limit exactly delta days
# after clinical clearance (the mean time to a negative RDT).
# RDT sensitivity among HRP2-expressing infections is SCHEMATIC_PARS$rdt_sens.
# =============================================================================

SCHEMATIC_PARS <- list(
  t_max        = 70,              # days since infection shown
  dt           = 0.05,            # time step (days)
  t_treat      = 10.5,            # treatment, at peak parasitaemia (days)
  t_clear      = 15.5,            # clinical clearance: parasitaemia falls below the LM limit (days)
  lm_limit     = 0.35,            # LM detection limit (relative parasitaemia)
  rdt_limit    = 0.06,            # RDT detection limit (relative HRP2); RDTs detect lower antigen levels
  symptom_frac = 0.5,             # symptoms start when parasitaemia reaches this share of its peak
  delta        = 15 / log(2),     # mean HRP2 clearance time (days): 15-day median, main prior centre
  hrp2_peak    = 0.75,            # peak circulating HRP2 (relative units)
  rdt_sens     = 0.9,             # RDT sensitivity among HRP2-expressing infections (main analysis)
  n_cohort     = 5000,            # cohort size
  offset_sd    = 0.75,            # SD of each individual's timing offset (days)
  seed         = 20261112
)

# Relative parasitaemia of the illustrative individual: gamma-shaped growth to
# a peak of 1 at treatment, then exponential decline that reaches the LM limit
# at clinical clearance.
schematic_parasitaemia <- function(t, p = SCHEMATIC_PARS) {
  a    <- 6                                               # growth shape
  tau  <- (p$t_clear - p$t_treat) / log(1 / p$lm_limit)   # post-treatment decline time constant
  rise <- ifelse(t > 0, (t / p$t_treat)^a * exp(a * (1 - t / p$t_treat)), 0)
  ifelse(t <= p$t_treat, rise, exp(-(t - p$t_treat) / tau))
}

# Circulating HRP2 of the illustrative individual: produced in proportion to
# parasitaemia and removed with time constant tau, scaled to hrp2_peak. tau is
# solved so that HRP2 falls below the RDT detection limit exactly delta days
# after clinical clearance: the average individual turns RDT-negative at
# t_clear + delta, the cohort mean.
schematic_hrp2 <- function(p = SCHEMATIC_PARS) {
  t   <- seq(0, p$t_max + 30, by = p$dt)
  par <- schematic_parasitaemia(t, p)
  curve <- function(tau) {
    h <- numeric(length(t))
    for (i in seq_along(t)[-1]) h[i] <- h[i - 1] + p$dt * (par[i - 1] - h[i - 1] / tau)
    h / max(h) * p$hrp2_peak
  }
  gap <- function(tau) stats::approx(t, curve(tau), xout = p$t_clear + p$delta)$y - p$rdt_limit
  tau <- stats::uniroot(gap, c(1, 60))$root
  list(t = t, h = curve(tau), tau = tau, rdt_limit = p$rdt_limit)
}

# Event times of the illustrative individual (timeline under the figure).
schematic_events <- function(p = SCHEMATIC_PARS) {
  t <- seq(0, p$t_treat, by = p$dt)
  data.table::data.table(
    event = c("Infection", "Symptom\nonset", "Treatment", "Clinical\nclearance",
              "Mean time to\nnegative rapid\ndiagnostic test"),
    t = c(0, t[which(schematic_parasitaemia(t, p) >= p$symptom_frac)[1]], p$t_treat,
          p$t_clear, p$t_clear + p$delta))
}

# The artificial cohort: every individual follows the illustrative curves
# shifted by a timing offset, keeps (with probability rdt_sens) an RDT that
# detects HRP2, and stays HRP2-positive for Exp(mean delta) days after its own
# clinical clearance. Returns the individual curves (panel A), the cohort
# detectability (panel B) and the joint outcome shares (panel C).
simulate_schematic_cohort <- function(p = SCHEMATIC_PARS) {
  set.seed(p$seed)
  hr      <- schematic_hrp2(p)
  times   <- seq(0, p$t_max, by = 0.25)
  offset  <- stats::rnorm(p$n_cohort, 0, p$offset_sd)
  rdt_ok  <- stats::runif(p$n_cohort) < p$rdt_sens
  hrp2_end <- p$t_clear + offset + stats::rexp(p$n_cohort, 1 / p$delta)   # HRP2-negative after this

  lm_pos <- outer(offset, times, function(o, t) schematic_parasitaemia(t - o, p) > p$lm_limit)
  hrp2_rising <- outer(offset, times, function(o, t)
    stats::approx(hr$t, hr$h, xout = pmax(t - o, 0), rule = 2)$y >= hr$rdt_limit & (t - o) <= p$t_clear)
  after_clear <- outer(offset, times, function(o, t) (t - o) > p$t_clear)
  hrp2_pos <- hrp2_rising | (after_clear & outer(hrp2_end, times, ">"))
  rdt_pos  <- hrp2_pos & rdt_ok

  share <- function(m) colMeans(m)
  outcomes <- data.table::data.table(
    t = times,
    `True positive`  = share(lm_pos & rdt_pos),
    `False positive` = share(!lm_pos & rdt_pos),
    `False negative` = share(lm_pos & !rdt_pos),
    `True negative`  = share(!lm_pos & !rdt_pos))

  individual <- data.table::data.table(t = times,
    parasitaemia = schematic_parasitaemia(times, p),
    hrp2 = stats::approx(hr$t, hr$h, xout = times)$y)

  list(individual = individual,
       detect = data.table::data.table(t = times, LM = share(lm_pos), RDT = share(rdt_pos)),
       outcomes = outcomes, events = schematic_events(p),
       lm_limit = p$lm_limit, rdt_limit = hr$rdt_limit, pars = p)
}

# Three panels plus the event timeline.
plot_schematic <- function(sim = simulate_schematic_cohort()) {
  COL_LM <- "#00798C"; COL_RDT <- "#D1495B"
  OUTCOME_COLS <- c(`True positive` = "#D1495B", `False positive` = "#E8A3AD",
                    `False negative` = "#8CC5CE", `True negative` = "#00798C")
  OUTCOME_LABS <- c(`True positive`  = "True positive\n(light microscopy +,\nrapid diagnostic test +)",
                    `False positive` = "False positive\n(light microscopy -,\nrapid diagnostic test +)",
                    `False negative` = "False negative\n(light microscopy +,\nrapid diagnostic test -)",
                    `True negative`  = "True negative\n(light microscopy -,\nrapid diagnostic test -)")
  ev <- sim$events
  t_max <- sim$pars$t_max
  vlines <- ggplot2::geom_vline(xintercept = ev$t[-1], linetype = "dotted", colour = "grey30", linewidth = 0.4)
  x_scale <- ggplot2::scale_x_continuous(limits = c(0, t_max), breaks = seq(0, t_max, 10), expand = c(0, 0))
  theme_s <- ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.key.height = ggplot2::unit(1.6, "lines"),
                   strip.text = ggplot2::element_text(size = ggplot2::rel(1.05)))

  # A: illustrative individual
  ind <- data.table::melt(sim$individual, id.vars = "t", variable.name = "series", value.name = "level")
  ind[, series := factor(series, c("parasitaemia", "hrp2"), c("Parasitaemia", "Circulating HRP2"))]
  delta_y <- 1.0
  pA <- ggplot2::ggplot(ind, ggplot2::aes(t, level)) +
    ggplot2::geom_area(ggplot2::aes(fill = series), alpha = 0.25, position = "identity") +
    ggplot2::geom_line(ggplot2::aes(colour = series), linewidth = 1.2) +
    ggplot2::geom_hline(yintercept = sim$lm_limit, linetype = "dashed", colour = COL_LM) +
    ggplot2::geom_hline(yintercept = sim$rdt_limit, linetype = "dashed", colour = COL_RDT) +
    ggplot2::annotate("text", x = t_max - 1, y = sim$lm_limit, label = "Light microscopy detection limit",
                      hjust = 1, vjust = -0.5, size = 3.6, colour = "grey20") +
    ggplot2::annotate("text", x = t_max - 1, y = sim$rdt_limit, label = "Rapid diagnostic test detection limit",
                      hjust = 1, vjust = -0.5, size = 3.6, colour = "grey20") +
    ggplot2::annotate("segment", x = ev$t[4], xend = ev$t[5], y = delta_y * 0.8, yend = delta_y * 0.8,
                      colour = "#66A182", linewidth = 0.8,
                      arrow = ggplot2::arrow(ends = "both", length = ggplot2::unit(0.15, "cm"), type = "closed")) +
    ggplot2::annotate("label", x = (ev$t[4] + ev$t[5]) / 2, y = delta_y * 0.94,
                      label = 'atop(bold(delta*":"~mean~HRP2), bold(clearance~time))', parse = TRUE,
                      fill = "#66A182", colour = "white", size = 3.6, label.size = 0) +
    vlines + x_scale +
    ggplot2::scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = c(Parasitaemia = COL_LM, `Circulating HRP2` = COL_RDT), name = NULL) +
    ggplot2::scale_fill_manual(values = c(Parasitaemia = COL_LM, `Circulating HRP2` = COL_RDT), name = NULL) +
    ggplot2::facet_wrap(~ "Illustrative individual") +
    ggplot2::labs(x = NULL, y = "Relative level\n(arbitrary units)") + theme_s

  # B: cohort detectability
  det <- data.table::melt(sim$detect, id.vars = "t", variable.name = "test", value.name = "share")
  det[, test := factor(test, c("LM", "RDT"), c("Light\nmicroscopy", "Rapid\ndiagnostic test"))]
  plateau_t <- (ev$t[3] + ev$t[4]) / 2 - 2
  pB <- ggplot2::ggplot(det, ggplot2::aes(t, share, colour = test)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::annotate("segment", x = plateau_t, xend = plateau_t, y = 1, yend = sim$pars$rdt_sens,
                      arrow = ggplot2::arrow(ends = "both", length = ggplot2::unit(0.12, "cm")), colour = "grey20") +
    ggplot2::annotate("segment", x = plateau_t, xend = plateau_t + 14, y = (1 + sim$pars$rdt_sens) / 2,
                      yend = (1 + sim$pars$rdt_sens) / 2, colour = "grey20", linewidth = 0.3) +
    ggplot2::annotate("text", x = plateau_t + 14.5, y = (1 + sim$pars$rdt_sens) / 2, hjust = 0, size = 3.6,
                      label = "Rapid diagnostic test sensitivity", colour = "grey20") +
    vlines + x_scale +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1.05), expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(COL_LM, COL_RDT), levels(det$test)), name = "Diagnostic") +
    ggplot2::facet_wrap(~ "Cohort") +
    ggplot2::labs(x = NULL, y = "Proportion of cohort\ndetectable (%)") + theme_s

  # C: cohort joint outcomes
  out <- data.table::melt(sim$outcomes, id.vars = "t", variable.name = "outcome", value.name = "share")
  out[, outcome := factor(outcome, names(OUTCOME_COLS))]
  pC <- ggplot2::ggplot(out, ggplot2::aes(t, share, fill = outcome)) +
    ggplot2::geom_area(colour = "white", linewidth = 0.2) +   # first level (true positive) stacks on top
    vlines + x_scale +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = OUTCOME_COLS, labels = OUTCOME_LABS, name = "Outcome") +
    ggplot2::facet_wrap(~ "Cohort") +
    ggplot2::labs(x = "Days since infection", y = "Proportion of\ncohort (%)") + theme_s

  # Event timeline (alternating label rows so close events do not overlap)
  ev[, row := c(1, 2, 1, 2, 1)]
  pT <- ggplot2::ggplot(ev, ggplot2::aes(t, 0)) +
    ggplot2::annotate("segment", x = 0, xend = t_max, y = 0, yend = 0, colour = "#2E4057", linewidth = 1,
                      arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"), type = "closed")) +
    ggplot2::geom_point(size = 3, colour = "#2E4057") +
    ggplot2::geom_text(ggplot2::aes(y = ifelse(row == 1, -0.35, -1.0), label = event), vjust = 1,
                       fontface = "bold", size = 3.7, colour = "#2E4057", lineheight = 0.9) +
    ggplot2::scale_x_continuous(limits = c(0, t_max), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-2.2, 0.2)) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme_void()

  patchwork::wrap_plots(pA, pB, pC, pT, ncol = 1, heights = c(1, 1, 1.15, 0.45)) +
    patchwork::plot_annotation(tag_levels = list(c("A", "B", "C", ""))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 20))
}
