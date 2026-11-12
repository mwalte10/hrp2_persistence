# =============================================================================
# src/1_data_prep/01_format_dhs_mis_data.R            (was src/01_pull_mis_data.R)
#
# ANALYSIS COMPONENT: DATA CLEANING AND TRANSFORMATION (DHS/MIS pull)
#
# Step 1a. Pulls MIS/DHS household survey data from the DHS API (rdhs) and
# formats it for the later steps: survey-weighted RDT, microscopy and joint
# RDT x microscopy (concordance) prevalence at national and sub-national level,
# fever, care-seeking and ACT indicators, subnational boundary polygons and
# cluster GPS points, and the individual-level child file used by the
# regression (step 2).
#
# REPULL_DATA = FALSE (default) reloads the cached CSVs in DIR_PULL and only
# rewrites the outputs; TRUE re-downloads everything (needs DHS credentials,
# see set_rdhs_config() below).
#
# Outputs (DIR_PULL): dhs_mis_{rdt,microscopy,concordance,fever,...}.csv,
#   dhs_mis_individual.csv, dhs_subnational_boundaries.RDS, and the data-check
#   files written by check_fns.R (run_manifest.csv, check_*.csv,
#   survey_summary*.csv). PCR, cluster GPS, fever-stratified concordance,
#   age-specific and non-malarial fever indicators from the original pull are
#   not needed for the manuscript (archive/src/01_pull_mis_data.R).
# The descriptive plots the original script drew are in archive/src/01_pull_mis_data.R.
# =============================================================================
source("./src/functions/paths.R")
source("./src/functions/check_fns.R")
library(data.table)
library(tidyverse)
library(rdhs)
library(haven)
library(survey)
library(ggrepel)
library(sf)

# =============================================================================
# DHS LOGIN CREDENTIALS
# The DHS/MIS microdata are not distributed with this repository. Each user
# must register with the DHS Program (https://dhsprogram.com/data/new-user-registration.cfm),
# request access to the surveys listed in the appendix for
# their own project, and supply their own credentials. Set them in
# .Renviron in the project root (copy .Renviron.example; never in this file):
#   DHS_EMAIL=you@institution.org
#   DHS_PROJECT=Your DHS project title, exactly as registered
# The password is prompted for at run time. rdhs.json (written by
# set_rdhs_config() in the project root) and the raw-data cache
# (OUTPUT_ROOT/01_pull_mis_data/dhs_raw_cache/) are git-ignored.
# =============================================================================
DHS_EMAIL   <- Sys.getenv("DHS_EMAIL")
DHS_PROJECT <- Sys.getenv("DHS_PROJECT")

dir.create(DIR_PULL, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# CACHE CONTROL
# REPULL_DATA <- TRUE  : re-download everything from DHS API (requires auth).
# REPULL_DATA <- FALSE : load existing CSVs and skip all DHS API calls.
# =============================================================================
REPULL_DATA <- FALSE   # TRUE re-downloads every survey from the DHS API (hours; needs DHS credentials)

KEY_CACHE <- c(
  file.path(DIR_PULL, "dhs_mis_rdt.csv"),
  file.path(DIR_PULL, "dhs_mis_microscopy.csv"),
  file.path(DIR_PULL, "dhs_mis_concordance.csv"),
  file.path(DIR_PULL, "dhs_mis_fever.csv"),
  file.path(DIR_PULL, "dhs_subnational_boundaries.RDS")
)

# Whether this run pulls from the DHS API, so the individual file at the end
# is rewritten whenever extract_kr is freshly built, not only when
# REPULL_DATA = TRUE
PULL_FRESH <- REPULL_DATA || !all(file.exists(KEY_CACHE))
RUN_START  <- Sys.time()   # files written after this go in run_manifest.csv

if (PULL_FRESH) {

# Authenticate with the user's own DHS credentials (see DHS LOGIN CREDENTIALS
# above). An existing rdhs configuration is reused.
if (is.null(rdhs::get_rdhs_config())) {
  if (!nzchar(DHS_EMAIL) || !nzchar(DHS_PROJECT))
    stop("Set DHS_EMAIL and DHS_PROJECT in .Renviron in the project root (your own DHS Program login; see .Renviron.example) before pulling DHS/MIS data.")
  set_rdhs_config(
    email           = DHS_EMAIL,
    project         = DHS_PROJECT,
    password_prompt = TRUE,
    global          = FALSE,
    cache_path      = file.path(DIR_PULL, "dhs_raw_cache")
  )
}

# =============================================================================
# SURVEYS
# =============================================================================
surveys <- dhs_surveys() %>%
  filter(SurveyType %in% c("MIS", "DHS") & RegionName ==  "Sub-Saharan Africa" & SurveyYear > 2000) %>%
  data.table()

if (nrow(surveys) == 0)
  stop("No MIS surveys found — check credentials and downloads.")

# Isolate which survey ID(s) trigger the DHS API 403
# data.table syntax throughout, per project convention

library(rdhs)
library(data.table)

check_survey_ids <- function(survey_ids) {

  results <- vector("list", length(survey_ids))

  for (i in seq_along(survey_ids)) {

    sid <- survey_ids[i]

    res <- tryCatch({
      dt <- data.table(dhs_datasets(fileType = "PR",
                                    fileFormat = "flat",
                                    surveyIds = sid))
      list(SurveyId = sid, status = "ok", n_files = nrow(dt))
    }, error = function(e) {
      list(SurveyId = sid, status = "FAILED", n_files = NA_integer_,
           error_msg = conditionMessage(e))
    })

    results[[i]] <- res
    # print as we go so a slow/hanging call is visible immediately, not just at the end
    cat(sprintf("[%d/%d] %s: %s\n", i, length(survey_ids), sid, res$status))

  }

  rbindlist(results, fill = TRUE)

}

survey_check <- check_survey_ids(surveys$SurveyId)

# surveys that failed
surveys <- surveys %>%
  filter(!(SurveyId %in% survey_check[status == "FAILED", SurveyId]))


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
to_result <- function(x) {
  v <- as.integer(zap_labels(x))
  fifelse(v %in% 0:1, v, NA_integer_)
}

.make_design <- function(dat) {
  options(survey.lonely.psu = "adjust")
  ids_var    <- if ("hv001" %in% names(dat)) "hv001" else "v001"
  strata_var <- if ("hv022" %in% names(dat)) "hv022" else "v022"
  ids_f    <- as.formula(paste0("~", ids_var))
  strata_f <- as.formula(paste0("~", strata_var))
  if (strata_var %in% names(dat) && !all(is.na(dat[[strata_var]]))) {
    svydesign(ids = ids_f, strata = strata_f, weights = ~sample_weight,
              data = dat, nest = TRUE)
  } else {
    svydesign(ids = ids_f, weights = ~sample_weight, data = dat)
  }
}

# =============================================================================
# PR DATA
#
# Variables:
#   hv001/hv002/hv005/hv022/hv024 — design variables + region
#   hvidx                           — individual line number (for KR merge)
#   hml35/sb115/sh214              — RDT result (standard then country-specific)
#   hml32                           — microscopy result
# =============================================================================
chunk_size <- 20
id_chunks <- split(surveys$SurveyId,
                   ceiling(seq_along(surveys$SurveyId) / chunk_size))

datasets_pr_list <- lapply(id_chunks, function(ids) {
  data.table(dhs_datasets(fileType = "PR", fileFormat = "flat", surveyIds = ids))
})

datasets_pr <- rbindlist(datasets_pr_list, fill = TRUE)
get_datasets(datasets_pr$FileName)

# Variable audit — review output, then ensure RDT_VARS / MICRO_VARS below are complete
questions <- search_variable_labels(
  datasets_pr$FileName,
  search_terms = c("rapid test", "rapid diagnostic", "parachec", "rdt",
                   "malaria test", "malaria result",
                   "blood smear", "microscopy", "blood slide", "thick film",
                   "pcr", "polymerase", "molecular",
                   "brand", "test kit", "type of test", "month")
)
print(unique(data.table(questions)[, .(variable, description)]))

RDT_VARS   <- c("hml35", "sb115", "sh214")
MICRO_VARS <- c("hml32")
PCR_VARS   <- c("hml36", "spcr", "sspec")

vars_pr <- search_variables(
  datasets_pr$FileName,
  variables = c(
    "hv001", "hv002", "hv005", "hv006", "hv007", "hv022", "hv024", "hvidx",
    "hv103", "hv105",   # de-facto resident + age (needed for PR-based fever fallback)
    "hc1",              # age in months (children under 5; ages PR-only children in the KR merge)
    "h22",              # fever in last 2 weeks (present in MIS PR files without KR)
    RDT_VARS, MICRO_VARS, PCR_VARS,
    "hml37"
  )
)

extract_pr <- extract_dhs(vars_pr, add_geo = FALSE)
extract_pr <- rbindlist(lapply(extract_pr, function(x) {
  dt <- data.table(x)
  dt[, adm_1 := tolower(as.character(as_factor(hv024)))]
  if ("hml37" %in% names(dt))
    dt[, rdt_type := tolower(as.character(as_factor(hml37)))]
  dt
}), fill = TRUE, ignore.attr = TRUE)

extract_pr[, sample_weight := as.numeric(hv005) / 1e6]


# PR-based fever flag — present in MIS PR files that lack a separate KR file
if ("h22" %in% names(extract_pr))
  extract_pr[, fever_pr := to_result(h22)]

rdt_cols   <- intersect(RDT_VARS,   names(extract_pr))
micro_cols <- intersect(MICRO_VARS, names(extract_pr))
pcr_cols   <- intersect(PCR_VARS,   names(extract_pr))

if (length(rdt_cols) > 0) {
  extract_pr[, rdt_result := do.call(fcoalesce, lapply(rdt_cols,   function(v) to_result(.SD[[v]])))]
} else {
  warning("No RDT variables found in extracted PR data.")
  extract_pr[, rdt_result := NA_integer_]
}

if (length(micro_cols) > 0) {
  extract_pr[, micro_result := do.call(fcoalesce, lapply(micro_cols, function(v) to_result(.SD[[v]])))]
} else {
  warning("No microscopy variables found in extracted PR data.")
  extract_pr[, micro_result := NA_integer_]
}

if (length(pcr_cols) > 0) {
  extract_pr[, pcr_result := do.call(fcoalesce, lapply(pcr_cols, function(v) to_result(.SD[[v]])))]
} else {
  message("hml36 (PCR) not present in these surveys — pcr_result set to NA.")
  extract_pr[, pcr_result := NA_integer_]
}

# Test variables are optional (many DHS surveys did no malaria testing)
check_var_coverage(extract_pr,
                   required = c("hv001", "hv002", "hv005", "hv024", "hvidx"),
                   optional = c("rdt_result", "micro_result", "hc1", "hv105", "hv006", "hv007"),
                   label = "PR", out_csv = file.path(DIR_PULL, "check_pr_var_coverage.csv"))


# Household-level interview date from PR (hv006/hv007), used as fallback for surveys
# where v006/v007 is absent from the KR file.  One row per SurveyId + cluster + HH.
pr_date_lookup <- NULL
if (all(c("hv006", "hv007") %in% names(extract_pr))) {
  pr_date_lookup <- unique(extract_pr[, .(
    SurveyId,
    hv001              = as.integer(zap_labels(hv001)),
    hv002              = as.integer(zap_labels(hv002)),
    interview_month_pr = as.integer(zap_labels(hv006)),
    interview_year_pr  = as.integer(zap_labels(hv007))
  )])
}

# Slim slice used to merge test results onto KR children via household + line number
# rdt_type included when hml37 was decoded (propagates per-child brand info).
# Region, weight, stratum and age ride along (suffixed _pr) so tested children
# with no KR record -- mother not interviewed, or not her de jure child -- keep
# them after the merge instead of ending up with only test results. Suffixed
# so they don't collide with the KR columns or switch .make_design() to hv022.
slim_cols <- c("SurveyId", "hv001", "hv002", "hvidx", "hv006", "hv007", "rdt_result", "micro_result")
if ("rdt_type" %in% names(extract_pr)) slim_cols <- c(slim_cols, "rdt_type")
extract_pr_slim <- extract_pr[!is.na(rdt_result) | !is.na(micro_result), ..slim_cols]
extract_pr_slim[, `:=`(
  adm_1_pr         = extract_pr[!is.na(rdt_result) | !is.na(micro_result), adm_1],
  sample_weight_pr = extract_pr[!is.na(rdt_result) | !is.na(micro_result), sample_weight],
  stratum_pr       = extract_pr[!is.na(rdt_result) | !is.na(micro_result), as.integer(zap_labels(hv022))],
  # hc1 (months) where the survey has it, else the midpoint of hv105 (years)
  age_months_pr    = extract_pr[!is.na(rdt_result) | !is.na(micro_result),
                                fcoalesce(if ("hc1" %in% names(extract_pr)) as.integer(zap_labels(hc1)) else NA_integer_,
                                          fifelse(as.integer(zap_labels(hv105)) < 95L,
                                                  as.integer(zap_labels(hv105)) * 12L + 6L, NA_integer_))]
)]
rm(slim_cols)

# =============================================================================
# SUBNATIONAL BOUNDARY POLYGONS
#
# adm_1 (sub_region) above is just the DHS region label/factor — it carries no
# geometry. Pulling the DHS Spatial Data Repository's own region polygons here
# means later scripts can spatially match sub_region against the GADM
# shapefiles bundled in the site package (site::fetch_site()), instead of
# relying solely on the fuzzy name-string crosswalk used in
# 01c_site_seasonality.R. No spatial join happens in this file — this just
# pulls and caches the raw polygons for that downstream match.
# =============================================================================
cat(sprintf("\nDownloading DHS subnational boundary shapefiles for %d surveys (~5s pause/survey — this takes a while)...\n",
            nrow(surveys)))

BOUNDARY_TIMEOUT_S <- 30   # hard cap per survey so one bad request can't block the whole loop

get_survey_boundary <- function(survid) {
  res <- tryCatch({
    setTimeLimit(elapsed = BOUNDARY_TIMEOUT_S, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
    download_boundaries(surveyId = survid, method = "sf", quiet_download = TRUE)
  }, error = function(e) {
    warning(sprintf("No boundary shapefile for survey %s: %s", survid, conditionMessage(e)))
    NULL
  })
  if (is.null(res) || length(res) == 0) return(NULL)
  bnd <- res[[1]]   # first (only) layer — typically named sdr_subnational_boundaries
  bnd$SurveyId    <- survid
  bnd$country     <- surveys[SurveyId == survid, CountryName]
  bnd$survey_year <- surveys[SurveyId == survid, SurveyYear]
  bnd
}

n_surv       <- length(surveys$SurveyId)
boundary_list <- vector("list", n_surv)
for (i in seq_len(n_surv)) {
  sid <- surveys$SurveyId[i]
  t0  <- Sys.time()
  boundary_list[[i]] <- get_survey_boundary(sid)
  cat(sprintf("  [%d/%d] %s: %s (%.1fs)\n", i, n_surv, sid,
              if (is.null(boundary_list[[i]])) "FAILED" else "ok",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
boundary_list <- Filter(Negate(is.null), boundary_list)

if (length(boundary_list) > 0) {
  dhs_boundaries <- dplyr::bind_rows(boundary_list)
  cat(sprintf("  Boundary columns available (review to pick the region-name field for downstream matching): %s\n",
              paste(setdiff(names(dhs_boundaries), "geometry"), collapse = ", ")))
  saveRDS(dhs_boundaries, file.path(DIR_PULL, "dhs_subnational_boundaries.RDS"))
  cat(sprintf("  Collected boundaries for %d/%d surveys.\n",
              length(unique(dhs_boundaries$SurveyId)), nrow(surveys)))
} else {
  warning("No DHS boundary shapefiles could be downloaded for any survey.")
  dhs_boundaries <- NULL
}

# =============================================================================
# KR DATA
#
# Variables:
#   v001/v002/v005/v022/v024 — design variables + region
#   b5/b16/b19               — child alive, line no. in PR roster, age in months
#   v008/b3                  — interview and birth date (CMC), for age where b19 is absent
#   h22                      — fever in last two weeks
#   h32*                     — care-seeking provider indicators
# =============================================================================
datasets_kr <- dhs_datasets(fileType = "KR", fileFormat = "flat",
                             surveyIds = surveys$SurveyId) %>% data.table()
get_datasets(datasets_kr$FileName)

vars_kr <- search_variables(
  datasets_kr$FileName,
  variables = c(
    "v001", "v002", "v005", "v006", "v007", "v008", "v022", "v024",
    "b3",   "b5",   "b16",  "b19",
    "h22",
    "h32a", "h32b", "h32c", "h32k",
    "h32d", "h32e", "h32f", "h32i", "h32j",
    "h32r", "h32s",
    "h37j", "h37g"
  )
)

extract_kr <- extract_dhs(vars_kr, add_geo = FALSE)
extract_kr <- rbindlist(lapply(extract_kr, function(x) {
  dt <- data.table(x)
  dt[, adm_1 := tolower(as.character(as_factor(v024)))]
  dt
}), fill = TRUE, ignore.attr = TRUE)


extract_kr[, sample_weight   := as.numeric(v005) / 1e6]
extract_kr[, interview_month := as.integer(zap_labels(v006))]
extract_kr[, interview_year  := as.integer(zap_labels(v007))]

# Fallback: if v006/v007 is absent from a survey's KR dictionary, every row in that
# survey will be NA.  Detect those surveys and overwrite from the PR-derived lookup.
if (!is.null(pr_date_lookup)) {
  kr_date_na  <- extract_kr[, .(all_na = all(is.na(interview_month))), by = SurveyId]
  flagged     <- kr_date_na[all_na == TRUE, SurveyId]
  if (length(flagged) > 0L) {
    cat(sprintf("  interview_month fallback (v006 absent): %d survey(s): %s\n",
                length(flagged), paste(flagged, collapse = ", ")))
    extract_kr <- merge(
      extract_kr, pr_date_lookup,
      by.x  = c("SurveyId", "v001", "v002"),
      by.y  = c("SurveyId", "hv001", "hv002"),
      all.x = TRUE
    )
    extract_kr[SurveyId %in% flagged,
               `:=`(interview_month = interview_month_pr,
                    interview_year  = interview_year_pr)]
    extract_kr[, c("interview_month_pr", "interview_year_pr") := NULL]
    # Warn if PR also had no date for any flagged survey
    still_na <- extract_kr[SurveyId %in% flagged,
                            .(any_filled = any(!is.na(interview_month))), by = SurveyId]
    bad <- still_na[any_filled == FALSE, SurveyId]
    if (length(bad) > 0L)
      warning(sprintf("hv006/hv007 also absent from PR for survey(s): %s — interview_month remains NA",
                      paste(bad, collapse = ", ")))
  }
}

# Restrict to alive children under 5. b19 (age in months) only exists from
# DHS-7 recodes; earlier KR files have no b19, so filtering on it alone dropped
# every child in those surveys (most surveys before 2016), taking their fever,
# care seeking and sub_region with them. Age there is interview date minus
# birth date, both in century-month codes.
extract_kr[, age_months := fcoalesce(
  if ("b19" %in% names(extract_kr)) as.integer(zap_labels(b19)) else NA_integer_,
  as.integer(zap_labels(v008)) - as.integer(zap_labels(b3))
)]
check_var_coverage(extract_kr,
                   required = c("v001", "v002", "v005", "v024", "b5", "age_months"),
                   optional = c("b16", "b19", "h22", "v006", "v007"),
                   label = "KR", out_csv = file.path(DIR_PULL, "check_kr_var_coverage.csv"))
kr_n_before <- extract_kr[, .N, by = SurveyId]
extract_kr <- extract_kr[
  as.integer(zap_labels(b5)) == 1L &
  age_months < 60L
]
check_filter_loss(kr_n_before, extract_kr[, .N, by = SurveyId], label = "KR under-5 alive filter")
rm(kr_n_before)

# fever
extract_kr[, fever := to_result(h22)]

# formal care-seeking: any public or private provider (excludes traditional healers)
to_01 <- function(x) as.integer(!is.na(x) & as.integer(zap_labels(x)) == 1L)

pub_cols  <- intersect(c("h32a", "h32b", "h32c", "h32k"),              names(extract_kr))
priv_cols <- intersect(c("h32d", "h32e", "h32f", "h32i", "h32j", "h32r", "h32s"), names(extract_kr))

pub_flag  <- rowSums(sapply(pub_cols,  function(v) to_01(extract_kr[[v]])), na.rm = TRUE) > 0
priv_flag <- rowSums(sapply(priv_cols, function(v) to_01(extract_kr[[v]])), na.rm = TRUE) > 0
extract_kr[, formal_care := as.integer(pub_flag | priv_flag)]

# ACT treatment in last two weeks (h37j = ACT; h37g = other artesunate combination)
act_cols <- intersect(c("h37j", "h37g"), names(extract_kr))
if (length(act_cols) > 0) {
  extract_kr[, act_treatment := as.integer(
    rowSums(sapply(act_cols, function(v) to_01(.SD[[v]])), na.rm = TRUE) > 0
  )]
} else {
  warning("No ACT treatment variables found in KR data.")
  extract_kr[, act_treatment := NA_integer_]
}

# Merge PR test results onto KR using household + individual line number
# KR key: SurveyId + v001 (cluster) + v002 (HH) + b16 (child's line in PR roster)
# PR key: SurveyId + hv001           + hv002       + hvidx
extract_kr <- merge(
  extract_kr,
  extract_pr_slim,
  by.x = c("SurveyId", "v001", "v002", "b16"),
  by.y = c("SurveyId", "hv001", "hv002", "hvidx"),

  all.x = TRUE, all.y = T
)
extract_kr[,interview_month := ifelse(is.na(hv006), v006, hv006)]
extract_kr[,interview_year := ifelse(is.na(hv007), v007, hv007)]

# Tested children with no KR record (all.y): take region, weight, stratum and
# age from PR. Their fever, care seeking and ACT stay NA (KR-only questions),
# so the KR-based indicators below still use KR children only. KR children keep
# their KR values (women's weight v005).
extract_kr[, `:=`(
  adm_1         = fcoalesce(adm_1, adm_1_pr),
  sample_weight = fcoalesce(sample_weight, sample_weight_pr),
  v022          = fcoalesce(as.integer(zap_labels(v022)), stratum_pr),
  age_months    = fcoalesce(age_months, age_months_pr)
)]
extract_kr[, c("adm_1_pr", "sample_weight_pr", "stratum_pr", "age_months_pr") := NULL]
# PR-only rows include older children and adults where a survey tests all ages
# (e.g. Rwanda); keep under 5s, plus the few with no recorded age
extract_kr <- extract_kr[is.na(age_months) | age_months < 60L]
check_kr_pr_match(extract_kr, kr_ids = datasets_kr$SurveyId,
                  out_csv = file.path(DIR_PULL, "check_kr_pr_match_by_survey.csv"))

cat(sprintf("KR rows: %d  (%d with all three results non-missing)\n",
            nrow(extract_kr),
            extract_kr[!is.na(fever) & !is.na(rdt_result) & !is.na(micro_result), .N]))

# =============================================================================
# INDICATOR FUNCTIONS — PR-BASED
# get_rdt(), get_microscopy(), get_concordance()
# Each returns sub-regional rows + one "national" row using svyby()/svymean().
# =============================================================================
get_rdt <- function(survid) {
  dat <- copy(extract_pr[SurveyId == survid & !is.na(rdt_result) & !is.na(adm_1)])
  if (nrow(dat) == 0) { warning(sprintf("No RDT observations for survey %s", survid)); return(NULL) }

  des    <- .make_design(dat)
  counts <- dat[, .(n_tested   = .N,
                    n_positive = sum(rdt_result == 1L),
                    n_negative = sum(rdt_result == 0L)), by = adm_1]

  sub_est <- as.data.table(svyby(~rdt_result, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, c("rdt_result", "se"), c("prev_rdt", "se_prev_rdt"))

  nat_est <- svymean(~rdt_result, des, na.rm = TRUE)
  nat_row <- dat[, .(adm_1 = "national", n_tested = .N,
                     n_positive = sum(rdt_result == 1L), n_negative = sum(rdt_result == 0L))]
  nat_row[, prev_rdt    := coef(nat_est)[["rdt_result"]]]
  nat_row[, se_prev_rdt := SE(nat_est)]

  out <- rbindlist(list(merge(counts, sub_est[, .(adm_1, prev_rdt, se_prev_rdt)], by = "adm_1"),
                        nat_row), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}

get_microscopy <- function(survid) {
  dat <- copy(extract_pr[SurveyId == survid & !is.na(micro_result) & !is.na(adm_1)])
  if (nrow(dat) == 0) { warning(sprintf("No microscopy observations for survey %s", survid)); return(NULL) }

  des    <- .make_design(dat)
  counts <- dat[, .(n_tested   = .N,
                    n_positive = sum(micro_result == 1L),
                    n_negative = sum(micro_result == 0L)), by = adm_1]

  sub_est <- as.data.table(svyby(~micro_result, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, c("micro_result", "se"), c("prev_micro", "se_prev_micro"))

  nat_est <- svymean(~micro_result, des, na.rm = TRUE)
  nat_row <- dat[, .(adm_1 = "national", n_tested = .N,
                     n_positive = sum(micro_result == 1L), n_negative = sum(micro_result == 0L))]
  nat_row[, prev_micro    := coef(nat_est)[["micro_result"]]]
  nat_row[, se_prev_micro := SE(nat_est)]

  out <- rbindlist(list(merge(counts, sub_est[, .(adm_1, prev_micro, se_prev_micro)], by = "adm_1"),
                        nat_row), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}

# General concordance helper — computes weighted 2×2 table for any two binary
# PR-level test variables. label1/label2 set the output column name prefixes.
.get_concordance_pair <- function(survid, var1, var2, label1, label2) {
  dat <- copy(extract_pr[SurveyId == survid &
    !is.na(get(var1)) & !is.na(get(var2)) & !is.na(adm_1)])
  if (nrow(dat) == 0) {
    warning(sprintf("No paired %s/%s observations for survey %s", label1, label2, survid))
    return(NULL)
  }

  dat[, tp := as.integer(get(var1) == 1L & get(var2) == 1L)]
  dat[, fp := as.integer(get(var1) == 1L & get(var2) == 0L)]
  dat[, fn := as.integer(get(var1) == 0L & get(var2) == 1L)]
  dat[, tn := as.integer(get(var1) == 0L & get(var2) == 0L)]

  pp <- paste0("pct_", label1, "_pos_", label2, "_pos")
  pn <- paste0("pct_", label1, "_pos_", label2, "_neg")
  np <- paste0("pct_", label1, "_neg_", label2, "_pos")
  nn <- paste0("pct_", label1, "_neg_", label2, "_neg")

  des <- .make_design(dat)

  counts <- dat[, .(
    n_both_tested     = .N,
    n_tp = sum(tp), n_fp = sum(fp), n_fn = sum(fn), n_tn = sum(tn),
    n_eff_both_tested = sum(sample_weight)^2 / sum(sample_weight^2)
  ), by = adm_1]

  sub_est <- as.data.table(svyby(~tp + fp + fn + tn, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, old = c("tp","fp","fn","tn"), new = c(pp, pn, np, nn))
  sub_est[, pct_agree := get(pp) + get(nn)]
  sub_est[, kappa := {
    p_o <- get(pp) + get(nn)
    p_e <- (get(pp) + get(pn)) * (get(pp) + get(np)) +
           (get(np) + get(nn)) * (get(pn) + get(nn))
    (p_o - p_e) / (1 - p_e)
  }]
  keep_cols <- c("adm_1", pp, pn, np, nn, "pct_agree", "kappa")

  nat_est <- coef(svymean(~tp + fp + fn + tn, des, na.rm = TRUE))
  nat_row <- dat[, .(
    adm_1 = "national", n_both_tested = .N,
    n_tp = sum(tp), n_fp = sum(fp), n_fn = sum(fn), n_tn = sum(tn),
    n_eff_both_tested = sum(sample_weight)^2 / sum(sample_weight^2)
  )]
  nat_row[, (pp) := nat_est[["tp"]]]
  nat_row[, (pn) := nat_est[["fp"]]]
  nat_row[, (np) := nat_est[["fn"]]]
  nat_row[, (nn) := nat_est[["tn"]]]
  nat_row[, pct_agree := get(pp) + get(nn)]
  nat_row[, kappa := {
    p_o <- get(pp) + get(nn)
    p_e <- (get(pp) + get(pn)) * (get(pp) + get(np)) +
           (get(np) + get(nn)) * (get(pn) + get(nn))
    (p_o - p_e) / (1 - p_e)
  }]

  out <- rbindlist(list(merge(counts, sub_est[, ..keep_cols], by = "adm_1"), nat_row), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}

get_concordance           <- function(survid) .get_concordance_pair(survid, "rdt_result",   "micro_result", "rdt",   "micro")

# =============================================================================
# INDICATOR FUNCTIONS — KR-BASED
# get_fever(), get_fever_pr(), get_careseeking_formal(), get_act_coverage()
# =============================================================================
get_fever <- function(survid) {
  dat <- copy(extract_kr[SurveyId == survid & !is.na(fever) & !is.na(adm_1)])
  if (nrow(dat) == 0) { warning(sprintf("No fever observations for survey %s", survid)); return(NULL) }

  des    <- .make_design(dat)
  counts <- dat[, .(n_children = .N, n_fever = sum(fever == 1L)), by = adm_1]

  sub_est <- as.data.table(svyby(~fever, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, c("fever", "se"), c("prev_fever", "se_prev_fever"))

  nat_est <- svymean(~fever, des, na.rm = TRUE)
  nat_row <- dat[, .(adm_1 = "national", n_children = .N, n_fever = sum(fever == 1L))]
  nat_row[, prev_fever    := coef(nat_est)[["fever"]]]
  nat_row[, se_prev_fever := SE(nat_est)]

  out <- rbindlist(list(merge(counts, sub_est[, .(adm_1, prev_fever, se_prev_fever)], by = "adm_1"),
                        nat_row), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}


get_fever_pr <- function(survid) {
  # Fallback for surveys without a KR file: extract h22 from the PR file,
  # restricting to de-facto residents (hv103 == 1) under 5 years (hv105 < 5).
  if (!"fever_pr" %in% names(extract_pr)) return(NULL)
  dat <- copy(extract_pr[
    SurveyId == survid &
    !is.na(fever_pr) &
    !is.na(adm_1) &
    as.integer(zap_labels(hv103)) == 1L &
    as.integer(zap_labels(hv105)) <  5L
  ])
  if (nrow(dat) == 0) {
    warning(sprintf("No PR-based fever observations for survey %s", survid))
    return(NULL)
  }
  des    <- .make_design(dat)
  counts <- dat[, .(n_children = .N, n_fever = sum(fever_pr == 1L)), by = adm_1]

  sub_est <- as.data.table(svyby(~fever_pr, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, c("fever_pr", "se"), c("prev_fever", "se_prev_fever"))

  nat_est <- svymean(~fever_pr, des, na.rm = TRUE)
  nat_row <- dat[, .(adm_1 = "national", n_children = .N, n_fever = sum(fever_pr == 1L))]
  nat_row[, prev_fever    := coef(nat_est)[["fever_pr"]]]
  nat_row[, se_prev_fever := SE(nat_est)]

  out <- rbindlist(list(
    merge(counts, sub_est[, .(adm_1, prev_fever, se_prev_fever)], by = "adm_1"),
    nat_row
  ), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}

get_careseeking_formal <- function(survid) {
  dat <- copy(extract_kr[SurveyId == survid & fever == 1L & !is.na(adm_1)])
  if (nrow(dat) == 0) { warning(sprintf("No febrile children for survey %s", survid)); return(NULL) }

  des    <- .make_design(dat)
  counts <- dat[, .(n_febrile = .N, n_formal = sum(formal_care)), by = adm_1]

  sub_est <- as.data.table(svyby(~formal_care, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, c("formal_care", "se"), c("prop_formal_care", "se_prop_formal_care"))

  nat_est <- svymean(~formal_care, des, na.rm = TRUE)
  nat_row <- dat[, .(adm_1 = "national", n_febrile = .N, n_formal = sum(formal_care))]
  nat_row[, prop_formal_care    := coef(nat_est)[["formal_care"]]]
  nat_row[, se_prop_formal_care := SE(nat_est)]

  out <- rbindlist(list(merge(counts, sub_est[, .(adm_1, prop_formal_care, se_prop_formal_care)], by = "adm_1"),
                        nat_row), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}

# Proportion of febrile children who received ACT in the last two weeks.
# Denominator: febrile children with act_treatment non-missing (h37j/h37g).
get_act_coverage <- function(survid) {
  dat <- copy(extract_kr[SurveyId == survid & fever == 1L &
                           !is.na(act_treatment) & !is.na(adm_1)])
  if (nrow(dat) == 0) { warning(sprintf("No febrile children with ACT data for survey %s", survid)); return(NULL) }

  des    <- .make_design(dat)
  counts <- dat[, .(n_febrile = .N, n_act = sum(act_treatment == 1L)), by = adm_1]

  sub_est <- as.data.table(svyby(~act_treatment, ~adm_1, des, svymean, na.rm = TRUE))
  sub_est[, adm_1 := as.character(adm_1)]
  setnames(sub_est, c("act_treatment", "se"), c("prop_act", "se_prop_act"))

  nat_est <- svymean(~act_treatment, des, na.rm = TRUE)
  nat_row <- dat[, .(adm_1 = "national", n_febrile = .N, n_act = sum(act_treatment == 1L))]
  nat_row[, prop_act    := coef(nat_est)[["act_treatment"]]]
  nat_row[, se_prop_act := SE(nat_est)]

  out <- rbindlist(list(merge(counts, sub_est[, .(adm_1, prop_act, se_prop_act)], by = "adm_1"),
                        nat_row), fill = TRUE)
  out[, country     := surveys[SurveyId == survid, CountryName]]
  out[, survey_year := surveys[SurveyId == survid, SurveyYear]]
  out
}

# =============================================================================
# RUN ACROSS ALL SURVEYS
# =============================================================================
cat("Pulling RDT data...\n")
df_rdt <- rbindlist(Filter(Negate(is.null), lapply(surveys$SurveyId, get_rdt)), fill = TRUE)

cat("Pulling microscopy data...\n")
df_micro <- rbindlist(Filter(Negate(is.null), lapply(surveys$SurveyId, get_microscopy)), fill = TRUE)

cat("Pulling RDT vs microscopy concordance data...\n")
df_conc <- rbindlist(Filter(Negate(is.null), lapply(surveys$SurveyId, get_concordance)), fill = TRUE)

cat("Pulling fever data...\n")
df_fever <- rbindlist(Filter(Negate(is.null), lapply(surveys$SurveyId, get_fever)), fill = TRUE)


# For surveys without a KR file, fall back to PR-based h22 extraction
kr_survey_ids <- unique(extract_kr$SurveyId)
pr_only_ids   <- setdiff(surveys$SurveyId, kr_survey_ids)
if (length(pr_only_ids) > 0L) {
  cat(sprintf("No KR file for %d survey(s) — using PR-based fever fallback: %s\n",
              length(pr_only_ids), paste(pr_only_ids, collapse = ", ")))
  df_fever_pr <- rbindlist(Filter(Negate(is.null),
                                   lapply(pr_only_ids, get_fever_pr)), fill = TRUE)
  if (nrow(df_fever_pr) > 0L)
    df_fever <- rbindlist(list(df_fever, df_fever_pr), fill = TRUE)
}

cat("Pulling formal care-seeking data...\n")
df_careseeking <- rbindlist(Filter(Negate(is.null), lapply(surveys$SurveyId, get_careseeking_formal)), fill = TRUE)

cat("Pulling ACT coverage among febrile children...\n")
df_act_coverage <- rbindlist(Filter(Negate(is.null), lapply(surveys$SurveyId, get_act_coverage)), fill = TRUE)


# =============================================================================
# TIDY
# =============================================================================
clean <- function(dt) {
  dt[, sub_region  := stringr::str_to_title(adm_1)]
  dt[, survey_year := as.integer(survey_year)]
  dt[, adm_1       := NULL]
  dt[!grepl("refugee", sub_region, ignore.case = TRUE)]
}

df_rdt            <- clean(df_rdt)[order(sub_region, survey_year)]
df_micro          <- clean(df_micro)[order(sub_region, survey_year)]
df_conc           <- clean(df_conc)[order(sub_region, survey_year)]
df_fever         <- clean(df_fever)[order(sub_region, survey_year)]
df_careseeking   <- clean(df_careseeking)[order(sub_region, survey_year)]
df_act_coverage  <- clean(df_act_coverage)[order(sub_region, survey_year)]

stopifnot(
  all(df_rdt$n_positive   + df_rdt$n_negative   == df_rdt$n_tested),
  all(df_micro$n_positive + df_micro$n_negative == df_micro$n_tested),
  all(df_conc$n_tp + df_conc$n_fp + df_conc$n_fn + df_conc$n_tn == df_conc$n_both_tested),
  all(df_fever$n_fever        <= df_fever$n_children),
  all(df_careseeking$n_formal <= df_careseeking$n_febrile)
)

# =============================================================================
# COMBINED LONG TABLE
# One row per sub_region × survey_year × country × indicator.
# n_denominator / n_numerator are unweighted counts;
# estimate and se are design-weighted.
# =============================================================================
df_combined <- rbindlist(list(
  df_rdt[, .(sub_region, survey_year, country,
             indicator = "rdt", n_denominator = n_tested, n_numerator = n_positive,
             estimate = prev_rdt, se = se_prev_rdt)],
  df_micro[, .(sub_region, survey_year, country,
               indicator = "microscopy", n_denominator = n_tested, n_numerator = n_positive,
               estimate = prev_micro, se = se_prev_micro)],
  df_fever[, .(sub_region, survey_year, country,
               indicator = "fever", n_denominator = n_children, n_numerator = n_fever,
               estimate = prev_fever, se = se_prev_fever)],
  df_careseeking[, .(sub_region, survey_year, country,
                     indicator = "formal_care", n_denominator = n_febrile, n_numerator = n_formal,
                     estimate = prop_formal_care, se = se_prop_formal_care)],
  df_conc[, .(sub_region, survey_year, country,
              indicator = "conc_rdt_pos_micro_pos", n_denominator = n_both_tested,
              n_numerator = n_tp, estimate = pct_rdt_pos_micro_pos, se = NA_real_)],
  df_conc[, .(sub_region, survey_year, country,
              indicator = "conc_rdt_pos_micro_neg", n_denominator = n_both_tested,
              n_numerator = n_fp, estimate = pct_rdt_pos_micro_neg, se = NA_real_)],
  df_conc[, .(sub_region, survey_year, country,
              indicator = "conc_rdt_neg_micro_pos", n_denominator = n_both_tested,
              n_numerator = n_fn, estimate = pct_rdt_neg_micro_pos, se = NA_real_)],
  df_conc[, .(sub_region, survey_year, country,
              indicator = "conc_rdt_neg_micro_neg", n_denominator = n_both_tested,
              n_numerator = n_tn, estimate = pct_rdt_neg_micro_neg, se = NA_real_)],
  df_conc[, .(sub_region, survey_year, country,
              indicator = "conc_agree", n_denominator = n_both_tested,
              n_numerator = n_tp + n_tn, estimate = pct_agree, se = NA_real_)],
  df_conc[, .(sub_region, survey_year, country,
              indicator = "conc_kappa", n_denominator = n_both_tested,
              n_numerator = NA_integer_, estimate = kappa, se = NA_real_)]
))[order(sub_region, survey_year, indicator)]

} else {

  message("Loading cached outputs from DIR_PULL (set REPULL_DATA <- TRUE to re-download from DHS).")

  df_rdt          <- fread(file.path(DIR_PULL, "dhs_mis_rdt.csv"))
  df_micro        <- fread(file.path(DIR_PULL, "dhs_mis_microscopy.csv"))
  df_conc         <- fread(file.path(DIR_PULL, "dhs_mis_concordance.csv"))
  df_fever        <- fread(file.path(DIR_PULL, "dhs_mis_fever.csv"))
  df_careseeking  <- fread(file.path(DIR_PULL, "dhs_mis_careseeking_formal.csv"))
  df_combined     <- fread(file.path(DIR_PULL, "dhs_mis_combined.csv"))
  df_act_coverage <- fread(file.path(DIR_PULL, "dhs_mis_act_coverage.csv"))

}  # end REPULL_DATA / cache block

# =============================================================================
# OUTPUT
# =============================================================================
# Tables read by steps 1b and 1c (and the data stream figure). Every other
# indicator the original pull computed is in archive/src/01_pull_mis_data.R.
write_csv(df_rdt,          file.path(DIR_PULL, "dhs_mis_rdt.csv"))
write_csv(df_micro,        file.path(DIR_PULL, "dhs_mis_microscopy.csv"))
write_csv(df_conc,         file.path(DIR_PULL, "dhs_mis_concordance.csv"))
write_csv(df_fever,        file.path(DIR_PULL, "dhs_mis_fever.csv"))
write_csv(df_careseeking,  file.path(DIR_PULL, "dhs_mis_careseeking_formal.csv"))
write_csv(df_combined,     file.path(DIR_PULL, "dhs_mis_combined.csv"))
write_csv(df_act_coverage, file.path(DIR_PULL, "dhs_mis_act_coverage.csv"))

# Individual-level microdata — KR children with merged PR test results.
# Used by step 1c (individual file with deletion risk) and the step 2 regression.
# Skipped when loading from cache (file already exists on disk).
if (PULL_FRESH) {
ind_slim <- extract_kr[
  !is.na(rdt_result) | !is.na(micro_result),
  .(SurveyId,
    cluster         = as.integer(zap_labels(v001)),
    stratum         = as.integer(zap_labels(v022)),
    sample_weight,
    adm_1,
    age_months,
    interview_month, interview_year,
    rdt_result, micro_result, fever, formal_care, act_treatment)
]
if ("rdt_type" %in% names(extract_kr)) {
  ind_slim[, rdt_type := extract_kr[!is.na(rdt_result) | !is.na(micro_result), rdt_type]]
} else {
  ind_slim[, rdt_type := NA_character_]
}
ind_slim <- merge(
  ind_slim,
  surveys[, .(SurveyId, country = CountryName, survey_year = as.integer(SurveyYear))],
  by = "SurveyId"
)
ind_slim[, sub_region := stringr::str_to_title(adm_1)]
ind_slim[, adm_1      := NULL]
ind_slim <- ind_slim[!grepl("refugee", sub_region, ignore.case = TRUE)]
ind_slim <- ind_slim[, .(country, survey_year, sub_region, cluster, stratum,
                          sample_weight, age_months,
                          interview_month, interview_year, rdt_type,
                          rdt_result, micro_result, fever, formal_care, act_treatment)]
check_individual_completeness(ind_slim)
# Fever (KR, or the PR fallback) is a warning only: a survey can lack the
# fever question; the individual file must have every tested survey
check_survey_sets(df_conc, list(individual = ind_slim, fever = df_fever),
                  strict = c(individual = TRUE, fever = FALSE))
compare_survey_summary(ind_slim, DIR_PULL)
fwrite(ind_slim, file.path(DIR_PULL, "dhs_mis_individual.csv"))
write_run_manifest(DIR_PULL, RUN_START, "src/1_data_prep/01_format_dhs_mis_data.R")
}  # end ind_slim block

cat(sprintf("\nRDT written:              %d rows\n", nrow(df_rdt)))
cat(sprintf("Microscopy written:       %d rows\n", nrow(df_micro)))
cat(sprintf("Concordance written:      %d rows\n", nrow(df_conc)))
cat(sprintf("Fever written:            %d rows\n", nrow(df_fever)))
cat(sprintf("Care-seeking written:     %d rows\n", nrow(df_careseeking)))
cat(sprintf("Combined written:         %d rows\n", nrow(df_combined)))
cat(sprintf("ACT coverage written:     %d rows\n", nrow(df_act_coverage)))
