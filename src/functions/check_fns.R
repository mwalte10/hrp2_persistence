# =============================================================================
# src/functions/check_fns.R
#
# ANALYSIS COMPONENT: DATA CLEANING AND TRANSFORMATION (function library: data checks for steps 1-2)
#
# Data checks for steps 1a, 1c and 2. They exist because of the November 2026
# KR-PR bug: b19 (age in months) is absent from pre-DHS-7 KR files, so the
# under-5 filter silently dropped every child in most pre-2016 surveys, the
# individual file kept their test results with no region, weight, age or
# fever, and step 2 fitted 25 of ~55 surveys without any error. A stale
# individual file then survived a re-pull because its write was skipped.
#
# Each check works per survey, because a pooled share hides a survey that
# drops out entirely. Structural failures stop(); judgement calls warn().
# =============================================================================
library(data.table)

.check_fail <- function(msg, strict = TRUE) {
  if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
}

.fmt_surveys <- function(x, max_n = 15) {
  x <- unique(as.character(x))
  if (length(x) > max_n) x <- c(x[seq_len(max_n)], sprintf("... (%d more)", length(x) - max_n))
  paste(x, collapse = ", ")
}

# --- 1. Variable coverage by survey ------------------------------------------
# Share non-missing of each variable within each survey. A required variable
# that is entirely missing in a survey (absent from that recode, or renamed)
# stops the run; optional variables are only reported. A column absent from
# the extract altogether counts as 0% everywhere.
check_var_coverage <- function(dt, required, optional = character(), label,
                               by = "SurveyId", out_csv = NULL, strict = TRUE) {
  vars <- c(required, optional)
  cov <- dt[, lapply(setNames(vars, vars), function(v)
    if (v %in% names(.SD)) mean(!is.na(.SD[[v]])) else 0), by = by,
    .SDcols = intersect(vars, names(dt))]
  if (!is.null(out_csv)) fwrite(cov, out_csv)
  bad <- melt(cov[, c(by, required), with = FALSE], id.vars = by,
              variable.name = "variable", value.name = "share")[share == 0]
  if (nrow(bad) > 0)
    .check_fail(sprintf(
      "[%s] required variable(s) entirely missing in %d survey(s): %s. Check the recode for a renamed or absent variable before filtering on it.",
      label, uniqueN(bad[[by]]),
      .fmt_surveys(bad[, sprintf("%s (%s)", get(by), variable)])), strict)
  cat(sprintf("[%s] variable coverage OK for %d survey(s) (%d required variable(s))\n",
              label, nrow(cov), length(required)))
  invisible(cov)
}

# --- 2. Rows lost to a filter, by survey -------------------------------------
# n_before / n_after: data.tables of counts by survey (dt[, .N, by = SurveyId]).
# Expected losses for the KR under-5 alive filter are deceased children
# (~5-10%); a survey losing more than max_loss, or all of its rows, means the
# filter variable is missing or coded differently there.
check_filter_loss <- function(n_before, n_after, label, by = "SurveyId",
                              max_loss = 0.3, strict = TRUE) {
  m <- merge(n_before, n_after, by = by, all.x = TRUE, suffixes = c("_before", "_after"))
  m[is.na(N_after), N_after := 0L]
  m[, loss := 1 - N_after / N_before]
  bad <- m[loss > max_loss]
  if (nrow(bad) > 0)
    .check_fail(sprintf(
      "[%s] filter removed more than %.0f%% of rows in %d survey(s): %s",
      label, 100 * max_loss, nrow(bad),
      .fmt_surveys(bad[, sprintf("%s (%.0f%%)", get(by), 100 * loss)])), strict)
  cat(sprintf("[%s] filter loss OK: max %.0f%% across %d survey(s)\n",
              label, 100 * max(m$loss), nrow(m)))
  invisible(m)
}

# --- 3. KR-PR match rate, by survey ------------------------------------------
# Among tested under-5 children, the share with a KR record (b5 non-missing
# after the all.y merge). Children whose mother was not interviewed have no KR
# record, so 75-95% is normal. Only surveys that have a KR file are checked:
# PR-only surveys are 0% by construction.
check_kr_pr_match <- function(extract_kr, kr_ids, warn_below = 0.6, out_csv = NULL) {
  tested <- extract_kr[SurveyId %in% kr_ids &
                       (!is.na(rdt_result) | !is.na(micro_result)) &
                       (is.na(age_months) | age_months < 60L)]
  rate <- tested[, .(n_tested = .N, n_matched = sum(!is.na(b5)),
                     pct_matched = round(100 * mean(!is.na(b5)), 1)), by = SurveyId]
  if (!is.null(out_csv)) fwrite(rate[order(pct_matched)], out_csv)
  low <- rate[pct_matched < 100 * warn_below]
  if (nrow(low) > 0)
    warning(sprintf(
      "[KR-PR merge] fewer than %.0f%% of tested children matched to a KR record in %d survey(s): %s. Their fever and care seeking will be missing.",
      100 * warn_below, nrow(low),
      .fmt_surveys(low[order(pct_matched), sprintf("%s (%.0f%%)", SurveyId, pct_matched)])),
      call. = FALSE)
  cat(sprintf("[KR-PR merge] median match %.0f%% (range %.0f-%.0f%%) across %d survey(s)\n",
              median(rate$pct_matched), min(rate$pct_matched), max(rate$pct_matched), nrow(rate)))
  invisible(rate)
}

# --- 4. Same surveys across outputs ------------------------------------------
# Every survey in the reference table (by default the PR-based concordance,
# which only needs the test results) should appear in each other table.
# tables: named list of data.tables with the key columns.
check_survey_sets <- function(reference, tables, key = c("country", "survey_year"),
                              ref_label = "concordance", strict = TRUE) {
  ref_keys <- unique(reference[, ..key])
  for (nm in names(tables)) {
    missing <- fsetdiff(ref_keys, unique(tables[[nm]][, ..key]))
    if (nrow(missing) > 0)
      .check_fail(sprintf("[survey sets] %d survey(s) in %s are missing from %s: %s",
                          nrow(missing), ref_label, nm,
                          .fmt_surveys(do.call(paste, missing))),
                  strict = if (length(strict) > 1) strict[[nm]] else strict)
  }
  cat(sprintf("[survey sets] checked %s against %s (%d survey(s))\n",
              paste(names(tables), collapse = ", "), ref_label, nrow(ref_keys)))
  invisible(TRUE)
}

# --- 5. Individual-file completeness, by survey -------------------------------
# Every tested child should have sub-region, weight and age (from KR, or PR
# for children with no KR record); fever comes from KR only, so it is
# expected for most but not all children.
check_individual_completeness <- function(ind, min_core = 0.99, min_fever = 0.5,
                                          strict = TRUE) {
  comp <- ind[, .(n = .N,
                  sub_region    = mean(!is.na(sub_region) & sub_region != ""),
                  sample_weight = mean(!is.na(sample_weight)),
                  age_months    = mean(!is.na(age_months)),
                  fever         = mean(!is.na(fever))),
              by = .(country, survey_year)]
  core <- melt(comp[, .(country, survey_year, sub_region, sample_weight, age_months)],
               id.vars = c("country", "survey_year"))[value < min_core]
  if (nrow(core) > 0)
    .check_fail(sprintf(
      "[individual file] sub-region, weight or age below %.0f%% in %d survey(s): %s",
      100 * min_core, uniqueN(core[, .(country, survey_year)]),
      .fmt_surveys(core[, sprintf("%s %d %s (%.0f%%)", country, survey_year, variable, 100 * value)])),
      strict)
  low_fever <- comp[fever < min_fever]
  if (nrow(low_fever) > 0)
    warning(sprintf("[individual file] fever recorded for fewer than %.0f%% of children in %d survey(s): %s",
                    100 * min_fever, nrow(low_fever),
                    .fmt_surveys(low_fever[, sprintf("%s %d (%.0f%%)", country, survey_year, 100 * fever)])),
            call. = FALSE)
  cat(sprintf("[individual file] completeness OK for %d survey(s)\n", nrow(comp)))
  invisible(comp)
}

# --- 6. Run manifest and freshness --------------------------------------------
# Step 1a writes run_manifest.csv listing every file it wrote in this run
# (modified since run_start) with its modification time. Downstream steps call
# check_pull_fresh(): each input must be in the manifest with an unchanged
# modification time, so a file whose write was skipped, or that was changed
# by hand afterwards, is caught. Without a manifest (outputs from before this
# check existed) it falls back to requiring every input to be no more than
# tol_hours older than the reference file.
write_run_manifest <- function(dir, run_start, script) {
  f <- list.files(dir, full.names = TRUE)
  f <- f[!dir.exists(f) & basename(f) != "run_manifest.csv"]
  info <- file.info(f)
  info <- data.table(file = basename(f), size = info$size, mtime = info$mtime)[mtime >= run_start]
  info[, `:=`(run_start = format(run_start, "%Y-%m-%d %H:%M:%S"), script = script,
              mtime = format(mtime, "%Y-%m-%d %H:%M:%OS3"))]
  fwrite(info, file.path(dir, "run_manifest.csv"))
  cat(sprintf("[manifest] %d file(s) written this run recorded in %s\n",
              nrow(info), file.path(dir, "run_manifest.csv")))
  invisible(info)
}

check_pull_fresh <- function(dir, files, reference = "dhs_mis_concordance.csv",
                             tol_hours = 1, strict = TRUE) {
  paths <- file.path(dir, files)
  absent <- files[!file.exists(paths)]
  if (length(absent) > 0)
    .check_fail(sprintf("[freshness] input file(s) missing from %s: %s", dir,
                        paste(absent, collapse = ", ")), strict)
  man_path <- file.path(dir, "run_manifest.csv")
  if (file.exists(man_path)) {
    man <- fread(man_path, colClasses = list(character = "mtime"))
    now <- data.table(file = files,
                      mtime_now = format(file.info(paths)$mtime, "%Y-%m-%d %H:%M:%OS3"))
    m <- merge(now, man[, .(file, mtime)], by = "file", all.x = TRUE)
    not_in_run <- m[is.na(mtime), file]
    changed    <- m[!is.na(mtime) & mtime != mtime_now, file]
    if (length(not_in_run) > 0)
      .check_fail(sprintf(
        "[freshness] not written by the last step 1a run (%s): %s. Its write was skipped, so it predates the other outputs.",
        man$run_start[1], paste(not_in_run, collapse = ", ")), strict)
    if (length(changed) > 0)
      .check_fail(sprintf("[freshness] modified after the last step 1a run (%s): %s",
                          man$run_start[1], paste(changed, collapse = ", ")), strict)
    cat(sprintf("[freshness] %d input(s) all from the step 1a run of %s\n",
                length(files), man$run_start[1]))
  } else {
    ref_time <- file.info(file.path(dir, reference))$mtime
    age_h <- as.numeric(difftime(ref_time, file.info(paths)$mtime, units = "hours"))
    stale <- files[age_h > tol_hours]
    if (length(stale) > 0)
      .check_fail(sprintf(
        "[freshness] older than %s by more than %g h (no run_manifest.csv to check against): %s",
        reference, tol_hours, paste(sprintf("%s (%.0f h)", stale, age_h[age_h > tol_hours]),
                                    collapse = ", ")), strict)
    warning("[freshness] no run_manifest.csv in ", dir,
            "; checked file times only. Re-run step 1a to create one.", call. = FALSE)
  }
  invisible(TRUE)
}

# --- 7. Surveys reaching the regression ----------------------------------------
# brm() drops any row with a missing covariate without saying so, so the check
# is on complete cases of every variable in the formula, not on the prepared
# data. expected: survids that have paired RDT/LM results in the included
# countries. Stops if fewer than min_share of them reach the fit.
check_model_surveys <- function(model_data, formula, expected, min_share = 0.9,
                                strict = TRUE) {
  vars <- intersect(all.vars(formula), names(model_data))
  cc <- complete.cases(model_data[, ..vars])
  by_s <- model_data[, .(n_prepared = .N), by = survid]
  by_s <- merge(by_s, model_data[cc, .(n_fit = .N), by = survid], by = "survid", all.x = TRUE)
  by_s[is.na(n_fit), n_fit := 0L]
  by_s <- merge(data.table(survid = unique(expected)), by_s, by = "survid", all = TRUE)
  by_s[is.na(n_prepared), n_prepared := 0L]
  by_s[is.na(n_fit), n_fit := 0L]
  lost <- by_s[n_fit == 0]
  partial <- by_s[n_fit > 0 & n_fit < 0.9 * n_prepared]
  print(by_s[order(n_fit)], nrows = 200)
  if (nrow(partial) > 0) {
    na_by_var <- model_data[survid %in% partial$survid,
                            lapply(.SD, function(x) round(100 * mean(is.na(x)))),
                            by = survid, .SDcols = vars]
    warning(sprintf("[model data] more than 10%% of rows lost to missing covariates in %d survey(s): %s",
                    nrow(partial), .fmt_surveys(partial$survid)), call. = FALSE)
    print(na_by_var)
  }
  share <- uniqueN(by_s[n_fit > 0 & survid %in% expected, survid]) / uniqueN(expected)
  if (share < min_share) {
    # Two different causes: no rows left after data prep (missing fever,
    # treatment or tests), or rows present but every one has an NA covariate
    lost_prep <- lost[n_prepared == 0, survid]
    lost_na   <- lost[n_prepared > 0, survid]
    if (length(lost_na) > 0)
      print(model_data[survid %in% lost_na,
                       lapply(.SD, function(x) round(100 * mean(is.na(x)))),
                       by = survid, .SDcols = vars])
    .check_fail(sprintf(
      "[model data] only %d of %d surveys with paired tests reach the fit (%.0f%%, minimum %.0f%%). No rows after data prep (fever, care/ACT or tests missing): %s. Rows dropped by brm() for NA covariates (%% missing above): %s.",
      uniqueN(by_s[n_fit > 0 & survid %in% expected, survid]), uniqueN(expected),
      100 * share, 100 * min_share,
      if (length(lost_prep)) .fmt_surveys(lost_prep) else "none",
      if (length(lost_na)) .fmt_surveys(lost_na) else "none"), strict)
  }
  cat(sprintf("[model data] %d of %d surveys reach the fit; %s rows of %s prepared\n",
              uniqueN(by_s[n_fit > 0, survid]), uniqueN(expected),
              format(sum(by_s$n_fit), big.mark = ","),
              format(sum(by_s$n_prepared), big.mark = ",")))
  invisible(by_s)
}

# --- 8. Comparison with the previous run ---------------------------------------
# Per-survey summary of the individual file, compared with the copy saved by
# the previous run: flags surveys that appear, disappear, or change n tested
# or n with fever by more than rel_change. Warns only; the diff is saved for
# review after every re-pull.
compare_survey_summary <- function(ind, dir, rel_change = 0.1) {
  cur <- ind[, .(n_tested = .N,
                 n_fever  = sum(!is.na(fever)),
                 prev_rdt = round(weighted.mean(rdt_result, sample_weight, na.rm = TRUE), 4),
                 prev_lm  = round(weighted.mean(micro_result, sample_weight, na.rm = TRUE), 4)),
             by = .(country, survey_year)]
  path <- file.path(dir, "survey_summary.csv")
  if (file.exists(path)) {
    old <- fread(path)
    d <- merge(old, cur, by = c("country", "survey_year"), all = TRUE,
               suffixes = c("_previous", "_current"))
    d[, status := fcase(
      is.na(n_tested_previous), "new",
      is.na(n_tested_current),  "dropped",
      abs(n_tested_current - n_tested_previous) > rel_change * n_tested_previous |
        abs(n_fever_current - n_fever_previous) > rel_change * pmax(n_fever_previous, 1), "changed",
      default = "same")]
    fwrite(d[order(status, country, survey_year)], file.path(dir, "survey_summary_diff.csv"))
    flagged <- d[status != "same"]
    if (nrow(flagged) > 0)
      warning(sprintf(
        "[run comparison] %d survey(s) new, dropped or changed by more than %.0f%% since the previous run (see survey_summary_diff.csv): %s",
        nrow(flagged), 100 * rel_change,
        .fmt_surveys(flagged[, sprintf("%s %d (%s)", country, survey_year, status)])),
        call. = FALSE)
    else cat("[run comparison] no survey changed since the previous run\n")
  } else {
    cat("[run comparison] no previous survey_summary.csv; saving this run's as the baseline\n")
  }
  fwrite(cur, path)
  invisible(cur)
}
