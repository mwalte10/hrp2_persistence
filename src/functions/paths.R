# =============================================================================
# src/functions/paths.R
#
# ANALYSIS COMPONENT: SHARED CONFIGURATION (used by every step)
#
# Every file-system location the pipeline reads from or writes to. Sourced by
# every script and by the function libraries that run on the cluster, so it
# must stay free of side effects (assignments, plus reading the project
# .Renviron).
#
# All scripts assume the working directory is PROJECT_ROOT (open
# hrp2_persistence.Rproj, or setwd(PROJECT_ROOT) before sourcing).
# =============================================================================

# Locations come from environment variables, set in a .Renviron file in the
# project root (git-ignored; copy .Renviron.example), read below, so every
# script and cluster task sees the same paths. Defaults make the repository
# self-contained: outputs under outputs/, site files under data/site_files/.
#   HRP2_PROJECT_ROOT       project root (default: the working directory)
#   HRP2_OUTPUT_ROOT        pipeline outputs (default: <root>/outputs)
#   HRP2_SITE_CACHE_DIR     site files, <ISO>_site.rds (default: <root>/data/site_files)
#   HRP2_UGA_SNT_SITE_FILE  Uganda subnational site file (default: <root>/data/site_files/UGA_snt_site.RDS)
#   HRP2_MANUSCRIPT_DIR     optional: local copy of the manuscript LaTeX sources, kept
#                           outside this repository (default: unset)
# R reads .Renviron only when a session starts in the project folder, so load
# it here as well: this covers a session started elsewhere and then setwd(),
# and cluster tasks.
if (file.exists(".Renviron")) readRenviron(".Renviron")
env_or <- function(var, default) { v <- Sys.getenv(var); if (nzchar(v)) v else default }

PROJECT_ROOT <- normalizePath(env_or("HRP2_PROJECT_ROOT", getwd()), winslash = "/", mustWork = FALSE)

# Pipeline outputs. Everything here is derived from DHS/MIS microdata or site
# files and is git-ignored (outputs/). adm1_combos_with_paths.RDS stores
# absolute paths to the base_pars cache, so moving OUTPUT_ROOT needs step 3a
# Section C re-run.
OUTPUT_ROOT <- env_or("HRP2_OUTPUT_ROOT", file.path(PROJECT_ROOT, "outputs"))

# Inputs ----------------------------------------------------------------------
DATA_DIR          <- file.path(PROJECT_ROOT, "data")            # small inputs kept in the repo
HRP2_DELETION_CSV <- file.path(DATA_DIR, "hrp2_deletion.csv")    # country pfhrp2 deletion prevalence
HRP2_RISK_CSV     <- file.path(DATA_DIR, "pfhrp2_deletion_risk.csv")
UGA_NAME1_MAP_CSV <- file.path(DATA_DIR, "name_1_map.csv")
GEOBOUNDARIES_DIR <- file.path(DATA_DIR, "geoboundaries")       # admin-1 boundaries for the FPR/FDR maps
# Site files: restricted, available on request (README "Data access").
# data/site_files/ is git-ignored.
SITE_FILE_DIR     <- file.path(DATA_DIR, "site_files")
UGA_SNT_SITE_FILE <- env_or("HRP2_UGA_SNT_SITE_FILE", file.path(SITE_FILE_DIR, "UGA_snt_site.RDS"))   # Uganda subnational site file (~40 MB)
SITE_CACHE_DIR    <- env_or("HRP2_SITE_CACHE_DIR", SITE_FILE_DIR)                                   # site files, <ISO>_site.rds, no urban/rural split

# Step output folders (names kept from mal_proj so existing outputs resolve) ---
DIR_PULL       <- file.path(OUTPUT_ROOT, "01_pull_mis_data")            # step 1: 01
DIR_STRUCTURE  <- file.path(OUTPUT_ROOT, "02_structure_mis_data")       # step 1: 02, 03
DIR_REGRESSION <- file.path(OUTPUT_ROOT, "03b_brms_int_cluster_fit")    # step 2
DIR_GRID       <- file.path(OUTPUT_ROOT, "04_grid_simulation")          # step 3a
DIR_EMULATOR   <- file.path(OUTPUT_ROOT, "05_emulation")                # step 3b
DIR_MCMC       <- file.path(OUTPUT_ROOT, "06_hierarchical_fitting")     # step 3c and 4

# Old urban/rural site cache, built by step 1 (03_structure_mis_data.R); the
# only source of baseline EIR values.
EIR_SITE_CACHE_DIR <- file.path(DIR_STRUCTURE, "site_cache")

# Figures, manuscript numbers/tables and posterior estimates live in the
# repo, not in OUTPUT_ROOT. None of them contain survey microdata.
FIGURE_DIR    <- file.path(PROJECT_ROOT, "figures")      # src/plots.R
RESULTS_DIR   <- file.path(PROJECT_ROOT, "results")      # src/results_numbers.Rmd, fixed_effect_selection.R
POSTERIOR_DIR <- file.path(RESULTS_DIR, "posteriors")    # src/export_posteriors.R

# The manuscript is kept outside this repository. When HRP2_MANUSCRIPT_DIR
# points to it, results_numbers.Rmd writes only the \nm{} keys it cites and
# src/sync_manuscript.R copies the figures and numbers.tex into it.
MANUSCRIPT_DIR <- env_or("HRP2_MANUSCRIPT_DIR", NA_character_)
