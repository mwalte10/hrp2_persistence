# ANALYSIS COMPONENT: CALIBRATION ANALYSIS (function library: emulator surface helpers)
#
# plot_emulators.R
# Visual diagnostic: for a given site, compare the SMOOTH GP emulator's
# predicted log-likelihood surface against the RAW gridded log-likelihood it
# was fit to, as a pair of heatmaps over the (d, eir) input space -- d is
# the single HRP2 clearance scale malariasimulation 3.0.0's
# set_hrp2_parameters() now takes (the old separate treated/untreated decay
# rates were collapsed into one shared rate; see
# 15_launch_grid_persistence.R Section B), so (d, eir) is the WHOLE input
# space now -- there's no third dimension left to fix/profile the way an
# older 3-axis version of this file did.
#
# The raw panel renders an exact grid for a rectangular/factorial design; for
# an irregular (LHS) design, where the points aren't on a discrete grid at
# all, it instead bins d/eir into an n_bins x n_bins pseudo-grid and takes
# the best ll per bin. See plot_emulator_slice()'s is_irregular branch.
# hrp2_deletion_pct is NOT one of these dimensions -- it's fixed per-country
# outside the emulator entirely (see 16_fit_emulators_fns.R's header for why).
#
# Requires 16_fit_emulators_fns.R already sourced (predict_emulator_ll()) --
# call after fit_site_emulators().

library(ggplot2)
library(data.table)
library(patchwork)

# Must match 16_fit_emulators.R's is_shared_param_grids()/resolve_site_grid()
# exactly.
is_shared_param_grids <- function(param_grids) {
  all(c("d", "eir") %in% names(param_grids))
}
resolve_site_grid <- function(param_grids, site) {
  if (is_shared_param_grids(param_grids)) param_grids else param_grids[[site]]
}

# Sequential blue, light->dark -- same values as
# 15_launch_grid_persistence.R Section G/H's own heatmaps, for visual
# consistency across this project's diagnostic plots.
SEQ_BLUE          <- c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b")
REL_LOG_LIK_FLOOR <- -20

# Warm orange, deliberately a DIFFERENT hue from SEQ_BLUE's fill -- the
# validation-point overlay (plot_emulator_slice()'s validation_points arg)
# is a second, orthogonal quantity (|true - predicted| log-lik error at an
# external held-out point) layered on top of the tile fill (relative
# log-lik), not a variant of it, so it gets its own colour channel
# (geom_point's outline colour, not fill) and its own scale.
VALIDATION_POINT_COLOUR <- c("#fdf1e9", "#eb6834", "#8a3616")

AXES        <- c("d", "eir")
AXIS_LABELS <- c(
  d   = "hrp2_scale",
  eir = "eir"
)
LOG_AXES <- c("d", "eir")   # both are log-spaced

# Natural-scale axis value -> the value the GP was actually trained on (log
# for both axes) — matches fit_single_emulator()'s own transform exactly.
to_gp_scale <- function(axis, value) log(value)

GP_INPUT_NAME <- c(
  d   = "log_d",
  eir = "log_eir"
)


# --- Single axis-pair comparison: raw grid vs. smooth emulator ---

# Builds the raw/smooth panel PAIR without combining or annotating them --
# factored out of plot_emulator_slice() so a caller assembling a bigger
# multi-panel page (e.g. a per-site vetting PDF combining this with other,
# unrelated diagnostics) can lay these two panels out itself instead of
# inheriting plot_emulator_slice()'s own fixed side-by-side + title/subtitle
# composition. plot_emulator_slice() below is now a thin wrapper around this.
build_emulator_slice_panels <- function(
  ll_arr, param_grids, emulator,
  n_fine = 60,
  n_bins = 12,
  validation_points = NULL,
  relative_to_best = TRUE
) {
  # relative_to_best: TRUE (default, unchanged behaviour for every existing
  # caller) rescales fill to log-lik RELATIVE to this site's own best grid
  # point (rel_ll = ll - global_best), with a FIXED colour scale (SEQ_BLUE,
  # limits c(REL_LOG_LIK_FLOOR, 0)) -- the right choice when panels from
  # DIFFERENT sites (wildly different absolute log-lik magnitudes/sample
  # sizes) need to sit on one shared, comparable scale (Sections G/H of 04,
  # the batch per-site diagnostic PNGs here). FALSE instead fills with the
  # RAW log-lik value and lets the colour scale auto-range to this site's own
  # data -- useful when a caller wants the heatmap directly comparable, ON
  # THE SAME NUMBERS, to something else plotted in raw log-lik units (e.g.
  # 05_emulation.R's vetting-page LOO scatter, whose axes are raw log-lik,
  # not relative-to-best -- a reader comparing the two side by side otherwise
  # has to mentally subtract global_best from one of them).
  # validation_points: optional data.table, one row per EXTERNAL held-out
  # point for this site (04_grid_simulation.R's Section I -- real simulated
  # log-lik at NEW (hrp2_scale, eir) draws NOT on this emulator's own
  # training grid, vs. this emulator's own prediction at the same points).
  # Columns: d (hrp2_scale), eir (absolute), abs_error (|log_lik_true -
  # log_lik_pred|). Overlaid as outlined markers on BOTH panels so a point
  # sitting in a region where raw and smooth visibly disagree can be read
  # directly against its own error, rather than trusting the RMSE summary
  # number alone. NULL (default) draws neither panel's overlay -- fully
  # backward compatible with every existing caller.
  # ll_arr/param_grids: one site's 2D array and grid (see 16_fit_emulators.R).
  # (d, eir) is the WHOLE input space now (malariasimulation 3.0.0 collapsed
  # the old separate treated/untreated HRP2 clearance rates into one shared
  # d), so there's no third axis left to fix or profile the way an older
  # 3-axis version of this function did (axis_x/axis_y/fix_at args all gone
  # -- always plots AXES = c("d", "eir")). eir is ABSOLUTE (04_grid_
  # simulation.R's shared EIR_GRID), so param_grids$eir is already on the
  # natural scale a reader wants -- no per-site baseline rescaling needed
  # (this function used to take an eir_in arg to annotate the axis label with
  # the absolute range implied by a RELATIVE eir grid; removed now that eir
  # itself is absolute).
  # n_fine: resolution of the smooth emulator panel's prediction grid.
  # n_bins: irregular (LHS) designs only -- number of bins per axis for the
  # raw panel's binned heatmap (see is_irregular branch below).
  axis_x <- AXES[1]
  axis_y <- AXES[2]

  # ll_arr is EITHER a dense 2D array (rectangular/factorial grid) or a
  # tidy data.frame/data.table (irregular design, e.g. 15's GRID_METHOD =
  # "lhs" — see 16_fit_emulators.R::build_emulator_inputs_from_grid_search()'s
  # rectangular-grid detection). Only the RAW panel differs between the two
  # (an exact grid vs. a binned heatmap) — the SMOOTH panel is identical
  # either way, since predict_emulator_ll() only needs the already-fitted
  # emulator, not the training data's shape.
  is_irregular <- is.data.frame(ll_arr)

  if (is_irregular) {
    design <- as.data.table(ll_arr)
  } else {
    # Reuses fit_single_emulator()'s own expand.grid ordering exactly, so
    # this is guaranteed to reshape ll_arr consistently with how it was built.
    design <- as.data.table(expand.grid(
      d   = param_grids$d,
      eir = param_grids$eir
    ))
    design[, ll := as.vector(ll_arr)]
  }
  finite_design <- design[is.finite(ll)]
  global_best   <- max(finite_design$ll)

  if (is_irregular) {
    # No discrete grid to slice -- an LHS design's points are ~continuous in
    # both d and eir, so there's no small discrete set of values to group by.
    # Instead: bin d/eir into an n_bins x n_bins grid (evenly spaced on the
    # GP's own log scale, matching how the factorial grids themselves are
    # built) and take the BEST (max) ll within each bin -- this uses ALL of
    # the site's simulated points rather than a lossy nearest-neighbour
    # subset near one arbitrary reference value.
    x_t   <- to_gp_scale(axis_x, finite_design[[axis_x]])
    y_t   <- to_gp_scale(axis_y, finite_design[[axis_y]])
    x_edges <- seq(min(x_t), max(x_t), length.out = n_bins + 1L)
    y_edges <- seq(min(y_t), max(y_t), length.out = n_bins + 1L)
    x_mids  <- (head(x_edges, -1) + tail(x_edges, -1)) / 2
    y_mids  <- (head(y_edges, -1) + tail(y_edges, -1)) / 2

    raw_slice <- copy(finite_design)
    raw_slice[, x_bin := pmin(n_bins, pmax(1L, findInterval(x_t, x_edges, all.inside = TRUE)))]
    raw_slice[, y_bin := pmin(n_bins, pmax(1L, findInterval(y_t, y_edges, all.inside = TRUE)))]
    raw_slice <- raw_slice[, .(ll = max(ll)), by = .(x_bin, y_bin)]
    raw_slice[, (axis_x) := exp(x_mids[x_bin])]
    raw_slice[, (axis_y) := exp(y_mids[y_bin])]
    raw_panel_title <- sprintf("Raw grid search (binned, %d bins/axis)", n_bins)
  } else {
    # No slicing needed: (d, eir) is already the whole input space, so
    # finite_design IS the full 2D grid.
    raw_slice <- finite_design
    raw_panel_title <- "Raw grid search"
  }
  raw_slice[, rel_ll := if (relative_to_best) ll - global_best else ll]

  # --- SMOOTH panel: emulator prediction on a fine continuous grid over
  # the full (d, eir) input space ---
  fine_seq <- function(axis) {
    rng <- range(param_grids[[axis]])
    if (axis %in% LOG_AXES) exp(seq(log(rng[1]), log(rng[2]), length.out = n_fine))
    else seq(rng[1], rng[2], length.out = n_fine)
  }
  fine_dt <- as.data.table(expand.grid(setNames(
    list(fine_seq(axis_x), fine_seq(axis_y)), c(axis_x, axis_y)
  )))
  fine_dt[, ll_pred := predict_emulator_ll(
    emulator,
    to_gp_scale(axis_x, fine_dt[[axis_x]]),
    to_gp_scale(axis_y, fine_dt[[axis_y]])
  )]
  fine_dt[, rel_ll := if (relative_to_best) ll_pred - global_best else ll_pred]

  scale_x <- if (axis_x %in% LOG_AXES) scale_x_log10() else scale_x_continuous()
  scale_y <- if (axis_y %in% LOG_AXES) scale_y_log10() else scale_y_continuous()

  # fill_range: the SHARED colour-scale range both panels are locked to --
  # c(REL_LOG_LIK_FLOOR, 0) for the relative case (unchanged), or the union
  # of the raw grid's and the smooth GP surface's own log-lik values for the
  # raw case (computed here, not left to each geom_tile() to auto-scale
  # independently -- the raw grid's own worst points and the smooth surface
  # rarely span the same range, so independent auto-scaling would put raw
  # and smooth back on two DIFFERENT scales even though both are nominally
  # "raw log-lik" -- exactly the "which scale am I looking at" confusion
  # this relative_to_best = FALSE path exists to remove). Returned below so
  # a caller (e.g. a companion scatter plot of the same raw/predicted
  # values) can match its own axis limits to this SAME range.
  fill_range <- if (relative_to_best) c(REL_LOG_LIK_FLOOR, 0) else range(c(raw_slice$rel_ll, fine_dt$rel_ll), finite = TRUE)
  fill_scale <- function(name) scale_fill_gradientn(
    colours = SEQ_BLUE, limits = fill_range, oob = scales::squish, name = name
  )
  fill_legend_name <- if (relative_to_best) "log-lik\n(rel. to best)" else "log-lik"

  axis_labels <- AXIS_LABELS
  base_theme <- theme_minimal(base_size = 9) + theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "#e1e0d9", fill = NA, linewidth = 0.3),
    plot.title      = element_text(colour = "#0b0b0b", size = 9),
    axis.text       = element_text(colour = "#898781", size = 6.5),
    axis.title      = element_text(colour = "#52514e", size = 8),
    legend.position = "bottom"
  )

  # Both branches above now produce the same shape (one row per grid cell/
  # bin, with axis_x/axis_y/rel_ll) -- factorial's cells sit on the exact
  # grid, the irregular case's on evenly log-spaced bin midpoints -- so one
  # geom_tile() rendering covers both. Not the same "explicit width in raw
  # units" pitfall 15_launch_grid_persistence.R's heatmaps hit: no width/
  # height aesthetic is mapped here, so geom_tile()'s default resolution()-
  # based sizing is computed on the scale_x_log10()/scale_y_log10()-
  # transformed positions and comes out uniform (verified against
  # ggplot_build() output) -- correct as-is.
  # Fixed shape (open diamond) so the overlay reads as "external point", not
  # as another tile -- colour (not fill) carries |true - predicted| error,
  # leaving the tile's own fill legend (rel_ll) untouched.
  validation_layer <- if (!is.null(validation_points) && nrow(validation_points) > 0) {
    list(
      geom_point(
        data = validation_points, aes(x = .data[[axis_x]], y = .data[[axis_y]], colour = abs_error),
        inherit.aes = FALSE, shape = 5, size = 2.4, stroke = 1
      ),
      scale_colour_gradientn(colours = VALIDATION_POINT_COLOUR, name = "|error|\n(held-out)")
    )
  } else {
    NULL
  }

  p_raw <- ggplot(raw_slice, aes(x = .data[[axis_x]], y = .data[[axis_y]], fill = rel_ll)) +
    geom_tile() +
    scale_x + scale_y +
    fill_scale(fill_legend_name) +
    validation_layer +
    labs(title = raw_panel_title, x = axis_labels[[axis_x]], y = axis_labels[[axis_y]]) +
    base_theme

  p_smooth <- ggplot(fine_dt, aes(x = .data[[axis_x]], y = .data[[axis_y]], fill = rel_ll)) +
    geom_tile() +
    scale_x + scale_y +
    fill_scale(fill_legend_name) +
    validation_layer +
    labs(title = "GP emulator (smooth)", x = axis_labels[[axis_x]], y = axis_labels[[axis_y]]) +
    base_theme

  caption <- if (is_irregular) {
    sprintf(
      "Binned into %d x %d cells (best ll per bin) -- no exact grid slice for an irregular design.",
      n_bins, n_bins
    )
  } else {
    "Full 2D input space (d x eir) -- no other axes to fix or profile."
  }
  if (!is.null(validation_points) && nrow(validation_points) > 0) {
    caption <- paste0(caption, sprintf(
      " Diamonds: %d external held-out points (04_grid_simulation.R Section I), coloured by |true - predicted| log-lik.",
      nrow(validation_points)
    ))
  }

  # design: the (d, eir, ll) points actually used for the raw panel (finite
  # values only) -- exposed so a caller can e.g. compute predict_emulator_ll()
  # at these SAME points for an LOO-style raw-vs-predicted scatter, without
  # re-deriving the rectangular-vs-irregular grid extraction logic above.
  list(raw = p_raw, smooth = p_smooth, caption = caption, axis_x = axis_x, axis_y = axis_y,
       design = finite_design, global_best = global_best, fill_range = fill_range)
}

# --- Single axis-pair comparison: raw grid vs. smooth emulator, combined +
# annotated onto one row (the original standalone-figure entry point;
# build_emulator_slice_panels() above is now the reusable building block) ---


# --- Single (d x eir) row, for one site ---
#
# There's only ONE axis pair now (d x eir is the whole input space --
# malariasimulation 3.0.0 collapsed the old separate treated/untreated
# decay rates into one shared d, see plot_emulator_slice()), so this is a
# single row rather than the old three-orthogonal-pair stack. Kept as its
# own wrapper (rather than inlining plot_emulator_slice() at call sites)
# for API stability and in case a future dimension gets added back.


# --- Batch: one PNG per site, in a folder ---

# Site names ("<iso>_<name_1>", see build_emulator_inputs_from_grid_search())
# can contain characters that are awkward or unsafe in filenames (spaces,
# slashes in some admin1 names) -- sanitised to a safe filename stem, kept
# unique via a numeric suffix on collision (distinct site names that
# sanitise to the same stem, e.g. differing only by punctuation).
sanitise_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)


