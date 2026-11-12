suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sf)
  library(site)
  library(stringdist)
})

source("./src/functions/paths.R")
source("./src/functions/check_fns.R")
source("./src/functions/data_stream_fns.R")   # survey_test_years()

# Every step 1a output read below must come from the same step 1a run (see
# check_fns.R: a skipped write once left a months-old individual file here)
check_pull_fresh(DIR_PULL, tol_hours = 3, files = c(
  "dhs_mis_individual.csv", "dhs_mis_concordance.csv", "dhs_mis_rdt.csv",
  "dhs_mis_microscopy.csv", "dhs_mis_fever.csv", "dhs_mis_careseeking_formal.csv",
  "dhs_mis_combined.csv", "dhs_mis_act_coverage.csv",
  "dhs_subnational_boundaries.RDS"))

# =============================================================================
# src/1_data_prep/03_structure_mis_data.R          (was src/02_structure_mis_data.R)
#
# ANALYSIS COMPONENT: DATA CLEANING AND TRANSFORMATION (site matching, deletion and seasonality inputs)
#   Reads the restricted site files (SITE_CACHE_DIR) and DHS/MIS-derived
#   files; every output in DIR_STRUCTURE is git-ignored.
#
# Step 1c. Run after 02_survey_design_effects.R (Section 1 reads its
# survey_design_effects.RDS).
#
# Structures raw MIS/DHS data (from 01_pull_mis_data.R) for two downstream
# uses: (1) regression analysis -- cluster-level prevalence linked to
# site-file covariates + HRP2 deletion risk; (2) EIR/HRP2 calibration --
# population prevalence at the lowest representative admin unit and
# nationally. Merged from three originally-separate scripts (still run as
# three sections below, in this order, since Section 3 needs Section 2's
# own output):
#   Section 1 -- Map site prevalence            (was 02a_map_site_prev.R)
#   Section 2 -- Aggregate data + HRP2 deletion  (was 02b_aggregate_data.R)
#   Section 3 -- Extract site seasonality        (was 02c_site_seasonality.R)
# Section 1 is independent of Sections 2/3 (reads 01's outputs and
# 02_survey_design_effects.R's survey_design_effects.RDS).
#
# INPUT (paths in src/functions/paths.R):
#   DIR_PULL/*.csv, dhs_subnational_boundaries.RDS
#   HRP2_RISK_CSV, UGA_NAME1_MAP_CSV (data/)
#   UGA_SNT_SITE_FILE                    -- Uganda SNT (district-level, name_1 = DHS region)
#   DIR_STRUCTURE/survey_design_effects.RDS
#   site::fetch_site() / available_sites() -- cached under DIR_STRUCTURE/site_cache*/
#
# OUTPUT (all under DIR_STRUCTURE):
#   countries_included.RDS               -- countries with more than one paired survey (Section 1)
#   matched_site_mis_prevalence.RDS      -- site-modelled + DHS RDT/LM prevalence by name_1 x year (Section 1)
#   dhs_pfhrp2_merged.csv                -- aggregate sub-regional dataset + HRP2 deletion risk (Section 2, updated in Section 3)
#   dhs_mis_individual_hrp2.csv          -- individual-level DHS/MIS + HRP2 deletion risk (Section 2, updated in Section 3)
#
# New aggregate-level columns added in Section 3 (on top of Section 2's):
#   survey_month_mean, rdt_brand_modal, site_seasonality_type, site_peak_month_1/2,
#   site_match_level, site_months_to_peak_mean, site_smc_cov
# New individual-level columns added in Section 3:
#   site_seasonality_type, site_peak_month_1/2, site_match_level, site_smc_cov,
#   site_months_to_peak, prev_micro, prev_rdt, prop_care_or_act
# =============================================================================

# Real-world admin boundary shapefiles (DHS's here) commonly have degenerate
# vertices (zero-length edges) that GEOS tolerates but sf's default S2
# spherical engine rejects outright ("Loop N is not valid: Edge M is
# degenerate"). Disabling S2 falls back to GEOS's more tolerant planar-on-WGS84
# handling; areas are only used below for relative overlap fractions, not
# exact geodesic area, so the precision loss doesn't matter here. Set once,
# globally, since it's a session option Section 1 needs and Section 3 doesn't
# touch (harmless either way).
sf::sf_use_s2(FALSE)

OUT_DIR <- DIR_STRUCTURE
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ISO3 lookup for all DHS countries (static) -- shared by Sections 1 and 3
# (identical list in both original scripts; defined once here).
DHS_ISO3 <- c(
  "Angola"                     = "AGO",
  "Benin"                      = "BEN",
  "Burkina Faso"               = "BFA",
  "Burundi"                    = "BDI",
  "Cameroon"                   = "CMR",
  "Chad"                       = "TCD",
  "Congo Democratic Republic"  = "COD",
  "Cote d'Ivoire"              = "CIV",
  "Ethiopia"                   = "ETH",
  "Gabon"                      = "GAB",
  "Gambia"                     = "GMB",
  "Ghana"                      = "GHA",
  "Guinea"                     = "GIN",
  "Kenya"                      = "KEN",
  "Lesotho"                    = "LSO",
  "Liberia"                    = "LBR",
  "Madagascar"                 = "MDG",
  "Malawi"                     = "MWI",
  "Mali"                       = "MLI",
  "Mauritania"                 = "MRT",
  "Mozambique"                 = "MOZ",
  "Niger"                      = "NER",
  "Nigeria"                    = "NGA",
  "Rwanda"                     = "RWA",
  "Senegal"                    = "SEN",
  "Sierra Leone"               = "SLE",
  "South Africa"               = "ZAF",
  "Tanzania"                   = "TZA",
  "Togo"                       = "TGO",
  "Uganda"                     = "UGA",
  "Zambia"                     = "ZMB",
  "Zimbabwe"                   = "ZWE"
)

# Reverse of DHS_ISO3 (iso3c -> canonical/plain-ASCII country name) -- the
# CANONICAL country spelling this whole file (and every downstream
# consumer keying on `country` as a string, e.g. 04_grid_simulation.R
# Section F's own countries_included.RDS filter) should use. Needed because
# the site package's own s$prevalence$country field is NOT guaranteed to
# match: confirmed for real it returns "Côte d'Ivoire" (accented) where
# DHS_ISO3/dhs_mis_combined.csv's own spelling is "Cote d'Ivoire" (plain
# ASCII) -- same root cause behind the Section 1.5 merge-key bug and the
# region_of() accent-lookup fix already in this file, now also confirmed to
# silently drop every CIV row from 04_grid_simulation.R Section F's own
# country filter (mis[country %in% countries_included.RDS$country]) even
# after both of those other fixes. Overriding matched_site_mis_prevalence's
# own `country` column to this canonical spelling (see Section 1.5 below)
# fixes the root cause once, for every downstream consumer, instead of
# patching each string comparison site individually as they're discovered.
ISO3_TO_COUNTRY <- setNames(names(DHS_ISO3), DHS_ISO3)

# UNICEF-style ESARO (Eastern and Southern Africa) / WCARO (West and Central
# Africa) regional office grouping -- shared here (rather than defined once
# down in Section 2.7, where it used to live) so Section 1.5's own
# region-level exclusion reporting and Section 2.7's africa_region column
# both read off the SAME classification. Covers every country in DHS_ISO3
# above; extended from the original 17-country map (which was missing
# Benin, Cameroon, Chad, Congo Democratic Republic, Cote d'Ivoire, Ethiopia,
# Gabon, Gambia, Lesotho, Mauritania, Niger, Senegal, South Africa,
# Tanzania, Zimbabwe -- confirmed gap: Benin/Cote d'Ivoire/Senegal are all
# present in matched_site_mis_prevalence, so the original map would have
# left them as NA in any region-based report).
REGION_MAP <- c(
  Angola = "ESA", Burundi = "ESA", Ethiopia = "ESA", Kenya = "ESA",
  Lesotho = "ESA", Madagascar = "ESA", Malawi = "ESA", Mozambique = "ESA",
  Rwanda = "ESA", "South Africa" = "ESA", Tanzania = "ESA", Uganda = "ESA",
  Zambia = "ESA", Zimbabwe = "ESA",
  Benin = "WCA", "Burkina Faso" = "WCA", Cameroon = "WCA", Chad = "WCA",
  "Congo Democratic Republic" = "WCA", "Cote d'Ivoire" = "WCA", Gabon = "WCA",
  Gambia = "WCA", Ghana = "WCA", Guinea = "WCA", Liberia = "WCA",
  Mali = "WCA", Mauritania = "WCA", Niger = "WCA", Nigeria = "WCA",
  Senegal = "WCA", "Sierra Leone" = "WCA", Togo = "WCA"
)

# Looks up REGION_MAP through an ASCII transliteration on BOTH sides (its
# own names are already plain ASCII, so this only transliterates the
# incoming country string) -- a country column's actual spelling isn't
# guaranteed to match REGION_MAP's own accent-free keys exactly. Confirmed
# needed for real: the site package's own s$prevalence$country field for
# Cote d'Ivoire is "Côte d'Ivoire" (accented), which a raw REGION_MAP[country]
# lookup would silently return NA for (same failure class as the
# country-string merge-key bug this file's own Section 1.5 header
# documents) even though REGION_MAP already has a "Cote d'Ivoire" entry.
region_of <- function(country) REGION_MAP[iconv(country, to = "ASCII//TRANSLIT")]

# =============================================================================
# SECTION 1: Map site prevalence
# (formerly src/02a_map_site_prev.R)
#
# Spatially matches each DHS survey's sub-national (adm_1) boundary polygons
# (outputs/01_pull_mis_data/dhs_subnational_boundaries.RDS) against the site
# package's GADM admin-1 polygons ($shape$level_1), for every country that has
# more than one year of paired RDT + microscopy (LM) prevalence estimates.
#
# This is a spatial match rather than the fuzzy name-string crosswalk used in
# Section 3 because DHS sub-regions are not always simple GADM admin-1 units
# -- some surveys group several admin-1 areas into custom epidemiological
# zones (e.g. Angola's MIS "stable mesoendemic" / "hyperendemic" strata),
# which have no sensible admin-1 name to fuzzy-match against but do have real
# geometry that can be overlaid on the GADM boundaries.
#
# OUTPUT: dhs_site_name1_overlap.csv/.RDS -- one row per (DHS region x
# overlapping site name_1) pair, area-weighted in both directions, so
# downstream code can pick whichever matching rule it needs (best single
# match, or the full many-to-many overlap).
# =============================================================================

CACHE_DIR <- file.path(OUT_DIR, "site_cache")
dir.create(CACHE_DIR, showWarnings = FALSE)

# Overlaps smaller than this fraction of the DHS region's area are dropped as
# edge-touch / sliver artefacts rather than real matches.
MIN_OVERLAP_FRAC <- 0.01

# -----------------------------------------------------------------------------
# 1.1 COUNTRIES WITH >1 YEAR OF PAIRED RDT + MICROSCOPY (LM) ESTIMATES
# A "paired" year is a country x survey_year where both indicators were
# collected (in any sub_region) -- i.e. the underlying DHS survey ran both
# tests, not just that the country has both indicators in different years.
# -----------------------------------------------------------------------------
# survey_test_years() (data_stream_fns.R) drops single-child rows and the 2008
# Senegal survey; the appendix data stream figure uses the same function.
paired_years <- survey_test_years(fread(file.path(DIR_PULL, "dhs_mis_combined.csv")))[has_rdt & has_micro]

n_paired_years <- paired_years[, .N, by = country]
countries_keep <- sort(n_paired_years[N > 1, country])

cat(sprintf("%d countries have >1 year of paired RDT + microscopy estimates:\n  %s\n\n",
            length(countries_keep), paste(countries_keep, collapse = ", ")))

iso3_keep    <- DHS_ISO3[countries_keep]
missing_iso3 <- countries_keep[is.na(iso3_keep)]
if (length(missing_iso3) > 0)
  warning("No iso3c mapping for: ", paste(missing_iso3, collapse = ", "))
countries_keep <- countries_keep[!is.na(iso3_keep)]
iso3_keep      <- iso3_keep[!is.na(iso3_keep)]
saveRDS(list(iso = iso3_keep, country = countries_keep), 
        file.path(OUT_DIR, "countries_included.RDS"))

# -----------------------------------------------------------------------------
# 1.2 LOAD DHS SUBNATIONAL BOUNDARIES
# -----------------------------------------------------------------------------
dhs_boundaries <- readRDS(file.path(DIR_PULL, "dhs_subnational_boundaries.RDS"))
dhs_boundaries$survey_year <- as.integer(dhs_boundaries$survey_year)   # stored as character in the source RDS
dhs_boundaries <- st_make_valid(dhs_boundaries)

# Reproject to an equal-area CRS before any area/intersection/union work. With
# S2 off (see above), st_area() on unprojected lon/lat geometry requires the
# lwgeom package (not installed); projecting sidesteps that and gives genuine
# metric areas instead of degree^2 approximations -- all countries here are in
# Africa, so a single continental equal-area projection covers every one.
AFRICA_EQAREA_CRS <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
dhs_boundaries <- st_transform(dhs_boundaries, AFRICA_EQAREA_CRS)

# The raw SDR export has some (SurveyId, REGCODE) regions split across more
# than one feature row -- for some surveys these are exact duplicate polygons
# (e.g. RW2013MIS), for others they're legitimate disjoint multi-part regions
# stored as separate rows instead of one MULTIPOLYGON (e.g. ET2016DHS regcode
# 8). Dissolving by (SurveyId, REGCODE) fixes both: true duplicates collapse
# with no area change, and genuine multi-part regions get correctly unioned.
# Left un-dissolved, either case double-counts overlap area in the matching
# below (observed as pct_of_dhs_region > 1 for RW2013MIS before this fix).
n_before <- nrow(dhs_boundaries)
dhs_boundaries <- dhs_boundaries %>%
  group_by(SurveyId, survey_year, country, REGCODE, REGNAME) %>%
  summarise(.groups = "drop")
cat(sprintf("Dissolved DHS boundaries: %d rows -> %d rows (multi-part/duplicate regions merged)\n",
            n_before, nrow(dhs_boundaries)))

# SYNTHETIC RW2017 boundary, reusing RW2015DHS's own geometry -- the DHS
# Spatial Data Repository has no boundary release at all for Rwanda's 2017
# MIS round (confirmed: RW's real releases are 2005/2008/2010/2013/2015/
# 2019/2023/2025, no 2017), even though the raw indicator data
# (dhs_mis_rdt.csv/dhs_mis_microscopy.csv) has real RW 2017 RDT+LM
# prevalence. With no boundary row at all, Section 1.3's spatial overlay
# below has nothing to match RW/2017 site name_1s against, so those survey
# rows would otherwise never reach matched_site_mis_prevalence no matter
# how good the region-name matching logic gets downstream -- this is a
# genuine gap in the pulled geometry, not a name-matching problem.
# Substitution justified directly: RW2017's own sub_region labels (East,
# Kigali City, National, North, South, West) are IDENTICAL to RW2015DHS's
# boundary REGNAME values -- Rwanda's 5-province admin1 structure was
# unchanged between those rounds (RW2019DHS renamed "Kigali City" to
# "Kigali", confirming a real -- but LATER -- boundary change, not one that
# already applied by 2017). SurveyId is a clearly-marked synthetic string
# (not a real DHS SurveyId), so this is never mistaken for genuine 2017
# boundary data if inspected downstream.
rw2017_synthetic <- dhs_boundaries[dhs_boundaries$SurveyId == "RW2015DHS", ]
if (nrow(rw2017_synthetic) == 0) {
  warning("RW2015DHS boundary not found -- cannot synthesize RW2017 boundary geometry.")
} else {
  rw2017_synthetic$SurveyId    <- "RW2017MIS_SYNTH_FROM_RW2015DHS"
  rw2017_synthetic$survey_year <- 2017L
  dhs_boundaries <- rbind(dhs_boundaries, rw2017_synthetic)
  cat(sprintf("Synthesized %d RW2017 boundary rows from RW2015DHS geometry (no real RW2017 boundary release exists).\n",
              nrow(rw2017_synthetic)))
}

# -----------------------------------------------------------------------------
# 1.3 FETCH SITE ADMIN-1 SHAPES + SPATIALLY MATCH TO DHS BOUNDARIES
# -----------------------------------------------------------------------------
fetch_site_cached <- function(iso) {
  path <- file.path(CACHE_DIR, paste0(iso, "_site.RDS"))
  if (file.exists(path)) return(readRDS(path))
  s <- tryCatch(
    site::fetch_site(iso3c = iso, admin_level = 1, urban_rural = TRUE),
    error = function(e) {
      warning(sprintf("fetch_site('%s') failed: %s", iso, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(s)) saveRDS(s, path)
  s
}

# One country's worth of DHS-region x site-name_1 area overlaps.
match_country <- function(country_name, iso) {
  cat(sprintf("Matching %s (%s)...\n", country_name, iso))

  s <- fetch_site_cached(iso)
  if (is.null(s) || is.null(s$shape$level_1)) {
    warning(sprintf("  No site admin-1 shape for %s (%s) -- skipping", country_name, iso))
    return(NULL)
  }
  site_shp <- st_make_valid(st_transform(s$shape$level_1, st_crs(dhs_boundaries)))

  dhs_shp <- dhs_boundaries[dhs_boundaries$country == country_name, ]
  if (nrow(dhs_shp) == 0) {
    warning(sprintf("  No DHS boundaries found for %s -- skipping", country_name))
    return(NULL)
  }

  dhs_area_dt <- as.data.table(st_drop_geometry(dhs_shp))[
    , .(SurveyId, survey_year, dhs_regname = REGNAME, dhs_regcode = REGCODE,
        dhs_area = as.numeric(st_area(dhs_shp)))
  ]
  site_area_dt <- data.table(site_name_1 = site_shp$name_1,
                              site_area   = as.numeric(st_area(site_shp)))

  # Bulk overlay (not a per-pair loop): st_intersection.sf on two full layers
  # is the supported n x m path and tracks attribute duplication correctly
  # when a single pair's intersection splits into multiple geometries -- the
  # per-pair 1x1 form used to hit that case as an sf internal bug (mismatched
  # geometry/attribute row counts from a shared-boundary sliver).
  dhs_sub  <- dhs_shp[, c("SurveyId", "survey_year", "REGNAME", "REGCODE")]
  names(dhs_sub)[names(dhs_sub) %in% c("REGNAME", "REGCODE")] <- c("dhs_regname", "dhs_regcode")
  site_sub <- setNames(site_shp[, "name_1"], c("site_name_1", attr(site_shp, "sf_column")))

  inter <- tryCatch(
    suppressWarnings(st_intersection(dhs_sub, site_sub)),
    error = function(e) {
      warning(sprintf("  st_intersection failed for %s: %s", country_name, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(inter) || nrow(inter) == 0) {
    warning(sprintf("  No overlapping site/DHS polygons found for %s", country_name))
    return(NULL)
  }

  # Drop stray point/line slivers from shared-boundary topology noise, keep
  # only the polygonal overlap.
  inter <- suppressWarnings(st_collection_extract(inter, "POLYGON"))
  inter$overlap_area <- as.numeric(st_area(inter))

  rows <- as.data.table(st_drop_geometry(inter))[
    , .(overlap_area = sum(overlap_area)),
    by = .(SurveyId, survey_year, dhs_regname, dhs_regcode, site_name_1)
  ]
  rows <- merge(rows, dhs_area_dt,  by = c("SurveyId", "survey_year", "dhs_regname", "dhs_regcode"))
  rows <- merge(rows, site_area_dt, by = "site_name_1")
  rows[, `:=`(pct_of_dhs_region  = overlap_area / dhs_area,
              pct_of_site_region = overlap_area / site_area,
              dhs_area = NULL, site_area = NULL)]

  rows <- rows[pct_of_dhs_region >= MIN_OVERLAP_FRAC | pct_of_site_region >= MIN_OVERLAP_FRAC]
  if (nrow(rows) == 0) {
    warning(sprintf("  No overlaps above %.0f%% found for %s", MIN_OVERLAP_FRAC * 100, country_name))
    return(NULL)
  }
  rows[, `:=`(country = country_name, iso3c = iso)]
  rows
}

site_dhs_match_list <- Map(match_country, countries_keep, iso3_keep)
site_dhs_matches     <- rbindlist(Filter(Negate(is.null), site_dhs_match_list), fill = TRUE)

# Best single site name_1 for each DHS region (by area covered), and best
# single DHS region for each site name_1 -- kept as two flags so either
# direction of the crosswalk can be read straight off one table.
setorder(site_dhs_matches, SurveyId, dhs_regcode, -pct_of_dhs_region)
site_dhs_matches[, best_match_for_dhs_region  := seq_len(.N) == 1L, by = .(SurveyId, dhs_regcode)]
setorder(site_dhs_matches, SurveyId, site_name_1, -pct_of_site_region)
site_dhs_matches[, best_match_for_site_region := seq_len(.N) == 1L, by = .(SurveyId, site_name_1)]
setorder(site_dhs_matches, country, survey_year, dhs_regcode, -pct_of_dhs_region)

# -----------------------------------------------------------------------------
# 1.4 SAVE
# -----------------------------------------------------------------------------

n_countries_matched <- length(unique(site_dhs_matches$country))
cat(sprintf(
  "\nWrote %d DHS-region x site-name_1 overlap rows (%d best-for-DHS-region matches) for %d/%d countries to %s\n",
  nrow(site_dhs_matches), sum(site_dhs_matches$best_match_for_dhs_region), n_countries_matched,
  length(countries_keep), OUT_DIR
))

# -----------------------------------------------------------------------------
# 1.5 matched_site_mis_prevalence -- site-package modelled prevalence + DHS
# RDT/LM, by site name_1 x year
#
# pfpr / pvpr are the site package's own modelled Pf/Pv prevalence surface
# (s$prevalence, already at name_1 x year -- independent of the DHS
# crosswalk). pfpr_rdt_num/denom and pfpr_lm_num/denom are the raw DHS survey
# counts (n_positive / n_tested) for whichever DHS region best overlaps that
# site name_1 (best_match_for_site_region), mapped in via country + survey
# year + normalised region name. Only site name_1 x year combos with an
# actual DHS survey that year get non-NA RDT/LM columns; pfpr/pvpr are
# populated for every year the site package models.
# -----------------------------------------------------------------------------
# adm_1 (sub_region, from the PR/KR microdata's hv024/v024 factor label) and
# REGNAME (from the SDR boundary shapefile) name the same DHS-defined regions
# but are exported inconsistently in three distinct ways:
#  - accents/hyphenation ("Boké" vs "boke", "Centre Est" vs "centre-est")
#    -> fixed by transliterating + stripping punctuation below.
#  - bilingual/adjectival forms ("Nord"/"Sud" vs "north"/"south", "Northern"
#    vs "north", "Fct" vs "Federal Capital Territory") -> fixed by the alias
#    table, mirroring Section 3's francophone SITE_ALIASES.
#  - the boundary shapefile grouping site name_1 units into a COARSER region
#    than the survey actually reports at (Nigeria's boundary has 6
#    geopolitical zones, e.g. "North West", but RDT/microscopy prevalence is
#    reported per state, e.g. "Kaduna" -- the same granularity as site
#    name_1). Fixed below by trying site_name_1 itself as the join key first.
# A residual few surveys assign REGNAME to a genuinely different
# classification with no lexical relationship to either name (Kenya's
# malaria epidemiological zones vs its administrative provinces, Madagascar's
# climate zones vs its 22 official regions) -- no name-based match, fuzzy or
# aliased, can resolve those; they stay NA and are counted below.
REGION_ALIASES <- c(
  "nord" = "north", "sud" = "south", "est" = "east", "ouest" = "west",
  "nord ouest" = "north west", "sud ouest" = "south west",
  "centre est" = "centre east", "centre ouest" = "centre west",
  "centre nord" = "centre north", "centre sud" = "centre south",
  "northern" = "north", "southern" = "south", "eastern" = "east", "western" = "west",
  "fct" = "federal capital territory"
)
norm_region <- function(x) {
  x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- trimws(gsub("\\s+", " ", x))
  # Alias lookup applied PER TOKEN, not to the whole string -- REGION_ALIASES
  # was previously keyed on the full normalised string, so e.g. "northern" ->
  # "north" only fired when a region's ENTIRE name was "Northern"; it silently
  # did nothing for "Mid Northern" (stays "mid northern", never reducing to
  # "mid north" the way a same-round "Mid North" survey label needed to
  # match). Per-token application still reproduces every existing multi-word
  # alias correctly (e.g. "nord ouest" -> tokens "nord"/"ouest" -> "north"/
  # "west" -> "north west", identical to the old whole-phrase entry), so no
  # existing match changes -- it only additionally fixes compound names built
  # from an aliased word plus something else, e.g. Uganda's 2014 MIS
  # "Mid Northern" boundary region vs. that round's own "Mid North" RDT/
  # microscopy source label.
  vapply(strsplit(x, " "), function(toks) {
    hit <- REGION_ALIASES[toks]
    paste(ifelse(is.na(hit), toks, hit), collapse = " ")
  }, character(1))
}

# Word-order-invariant key: same tokens, sorted -- catches compound region
# names that differ only in word order (e.g. Uganda 2024's boundary shapefile
# "Buganda South" vs. that round's own RDT/microscopy source label "South
# Buganda"). A token swap moves every character after it, so plain edit
# distance is roughly the shorter token's own length -- nowhere near
# fuzzy_key()'s tight character-distance threshold -- so this needs its own
# tier, not a looser fuzzy threshold.
sort_tokens <- function(x) {
  vapply(strsplit(x, " "), function(toks) paste(sort(toks), collapse = " "), character(1))
}

site_match <- site_dhs_matches[best_match_for_site_region == TRUE,
                       .(country, iso3c, survey_year, site_name_1,
                         site_key       = norm_region(site_name_1),
                         dhs_region_key = norm_region(dhs_regname))]

rdt   <- fread(file.path(DIR_PULL, "dhs_mis_rdt.csv"))
micro <- fread(file.path(DIR_PULL, "dhs_mis_microscopy.csv"))

rdt_j <- unique(rdt[, .(country, survey_year, dhs_region_key = norm_region(sub_region),
                        pfpr_rdt_num = n_positive, pfpr_rdt_denom = n_tested)])
micro_j <- unique(micro[, .(country, survey_year, dhs_region_key = norm_region(sub_region),
                            pfpr_lm_num = n_positive, pfpr_lm_denom = n_tested)])

dup_key <- function(dt) dt[, .N, by = .(country, survey_year, dhs_region_key)][N > 1]
if (nrow(dup_key(rdt_j)) > 0 || nrow(dup_key(micro_j)) > 0)
  warning("Normalised (country, survey_year, region) key is not unique in RDT/microscopy data -- check for near-duplicate region names.")

# Known false-positive fuzzy pairs: two REAL but DIFFERENT DHS regions that
# happen to land within fuzzy_key()'s own edit-distance threshold of each
# other, confirmed by checking that country/year's actual source region list
# directly (not just the string distance) -- an audit of every fuzzy match
# in the full dataset found exactly these two:
#   - Madagascar: "antsiranana" (the site's own name) vs "atsinanana" -- a
#     real but geographically UNRELATED region (near Toamasina, nowhere
#     close to historical Antsiranana/Diana in the far north). Antsiranana's
#     crosswalk key for 2011/2013/2016 (survey rounds whose boundary
#     shapefile encodes eco-climatic zones, not admin regions -- see this
#     section's own header) never has a genuine name-based match, so the
#     fuzzy fallback was silently pulling in Atsinanana's real prevalence
#     data instead.
#   - Liberia: "south central" vs "north central" -- confirmed Liberia's own
#     2013 RDT/microscopy source tables only ever contain "North Central",
#     never "South Central" at all that year; the site's own crosswalk key
#     ("south central", from a different survey round's own regional
#     scheme) was matching the wrong region purely on string distance.
# Every OTHER fuzzy match found in that same audit (Burkina Faso "Haut-
# Bassins"/"Hauts Bassins", Mali "Timbuktu"/"Tombouctou", Senegal "Thies"/
# "Thics" -- the latter a mangled-encoding artifact in the raw CSV, not a
# different place -- Uganda "Bugishu"/"Bugisu") is a genuine spelling/
# transliteration variant of the SAME region, not a cross-region false
# positive, and is left alone. Kept as an explicit pair list (checked in
# EITHER direction, not tied to a specific country) rather than a
# per-country special case, so it generalises if the same failure mode ever
# shows up for a different pair.
FUZZY_MATCH_BLOCKLIST <- list(
  c("antsiranana", "atsinanana"),
  c("south central", "north central")
)
is_blocked_fuzzy_pair <- function(a, b) {
  any(vapply(FUZZY_MATCH_BLOCKLIST, function(pair) setequal(c(a, b), pair), logical(1)))
}

# Fuzzy fallback restricted to candidates from the SAME country + survey year
# (never cross-matches a different region), within a tight edit-distance
# threshold. Walks candidates in order of increasing distance so a blocked
# pair (see FUZZY_MATCH_BLOCKLIST above) can still fall through to the next-
# best real candidate rather than just refusing to match at all. Used as the
# last step within each tier below.
fuzzy_key <- function(target, cand_keys) {
  if (is.na(target) || length(cand_keys) == 0) return(NA_character_)
  d       <- adist(target, cand_keys)[1, ]
  thresh  <- max(2L, floor(0.2 * nchar(target)))
  ord     <- order(d)
  for (idx in ord) {
    if (d[idx] > thresh) break
    if (!is_blocked_fuzzy_pair(target, cand_keys[idx])) return(cand_keys[idx])
  }
  NA_character_
}

# Tier 1: join on site_name_1 itself (handles surveys that report at a finer
# granularity than the boundary shapefile's polygon, e.g. Nigeria). Tier 2:
# fall back to the boundary's own (possibly coarser) dhs_regname -- correct
# for surveys where dhs_regname genuinely IS the reporting level (e.g.
# Angola's epidemiological-zone MIS rounds).
# value_cols: every source column carried onto the matched site row (counts
# plus any per-source extras such as a source-row id or design effect);
# value_cols[1L] doubles as the "matched yet?" indicator column.
match_prev_source <- function(source_dt, value_cols, label) {
  probe <- value_cols[1L]
  keyed <- copy(site_match)
  keyed[, resolved_key := site_key]
  keyed <- merge(keyed, source_dt, by.x = c("country", "survey_year", "resolved_key"),
                 by.y = c("country", "survey_year", "dhs_region_key"), all.x = TRUE)
  n_tier1 <- sum(!is.na(keyed[[probe]]))

  need <- is.na(keyed[[probe]])
  if (any(need)) {
    keyed[need, resolved_key := dhs_region_key]
    keyed[, (value_cols) := NULL]
    keyed <- merge(keyed, source_dt, by.x = c("country", "survey_year", "resolved_key"),
                   by.y = c("country", "survey_year", "dhs_region_key"), all.x = TRUE)
  }
  n_tier2 <- sum(!is.na(keyed[[probe]])) - n_tier1

  # Tier 3: word-order-invariant match (see sort_tokens()'s own header) --
  # row-wise via set(), same pattern as the fuzzy tier below, so it only
  # ever touches rows still unmatched after tiers 1-2 and can't disturb an
  # already-resolved row's values.
  still_need <- which(is.na(keyed[[probe]]))
  n_tier3 <- 0L
  if (length(still_need) > 0) {
    source_sorted <- copy(source_dt)
    source_sorted[, dhs_region_key := sort_tokens(dhs_region_key)]
    for (i in still_need) {
      for (candidate_key in unique(sort_tokens(c(keyed$site_key[i], keyed$dhs_region_key[i])))) {
        cand <- source_sorted[country == keyed$country[i] & survey_year == keyed$survey_year[i]
                               & dhs_region_key == candidate_key]
        if (nrow(cand) > 0) {
          for (vc in value_cols) set(keyed, i, vc, cand[[vc]][1])
          n_tier3 <- n_tier3 + 1L
          break
        }
      }
    }
  }

  still_need <- which(is.na(keyed[[probe]]))
  n_fuzzy <- 0L
  for (i in still_need) {
    for (candidate_key in unique(c(keyed$site_key[i], keyed$dhs_region_key[i]))) {
      cand <- source_dt[country == keyed$country[i] & survey_year == keyed$survey_year[i]]
      hit  <- fuzzy_key(candidate_key, cand$dhs_region_key)
      if (!is.na(hit)) {
        row <- cand[dhs_region_key == hit][1]
        for (vc in value_cols) set(keyed, i, vc, row[[vc]])
        n_fuzzy <- n_fuzzy + 1L
        break
      }
    }
  }
  cat(sprintf("  %s: %d matched via site name_1 directly, %d more via the DHS boundary region, %d more via word-order-invariant match, %d more via fuzzy match\n",
              label, n_tier1, n_tier2, n_tier3, n_fuzzy))
  keyed[, c("country", "iso3c", "survey_year", "site_name_1", value_cols), with = FALSE]
}

# Source-row ids (src_*) identify WHICH (survey, DHS region) row a site-year's
# counts came from -- several sites can match the same row (one DHS region
# overlapping several admin1s), and 04_grid_simulation.R Section F uses these
# ids to down-weight shared rows so each surveyed child counts once in total.
# deff_* (02b_survey_design_effects.R, per survey) rides along with its source
# row. Joint RDT x microscopy counts (children with both tests) come from
# 01's concordance table: tp = RDT+/LM+, fp = RDT+/LM-, fn = RDT-/LM+, tn = RDT-/LM-.
survey_deff <- readRDS(file.path(OUT_DIR, "survey_design_effects.RDS"))[
  , .(country, survey_year, deff_rdt, deff_lm, deff_joint)]
conc <- fread(file.path(DIR_PULL, "dhs_mis_concordance.csv"))
conc_j <- unique(conc[, .(country, survey_year, dhs_region_key = norm_region(sub_region),
                          joint_n = n_both_tested, joint_tp = n_tp, joint_fp = n_fp,
                          joint_fn = n_fn, joint_tn = n_tn)])
add_source_extras <- function(dt, src_col, deff_col) {
  dt <- merge(dt, survey_deff[, c("country", "survey_year", deff_col), with = FALSE],
              by = c("country", "survey_year"), all.x = TRUE)
  dt[, (src_col) := sprintf("%s|%d|%s", country, survey_year, dhs_region_key)]
  dt
}
rdt_j   <- add_source_extras(rdt_j,   "src_rdt",   "deff_rdt")
micro_j <- add_source_extras(micro_j, "src_lm",    "deff_lm")
conc_j  <- add_source_extras(conc_j,  "src_joint", "deff_joint")
if (nrow(dup_key(conc_j)) > 0)
  warning("Normalised (country, survey_year, region) key is not unique in concordance data -- check for near-duplicate region names.")

rdt_matched   <- match_prev_source(rdt_j,   c("pfpr_rdt_num", "pfpr_rdt_denom", "src_rdt", "deff_rdt"), "RDT")
micro_matched <- match_prev_source(micro_j, c("pfpr_lm_num",  "pfpr_lm_denom",  "src_lm",  "deff_lm"),  "microscopy")
joint_matched <- match_prev_source(conc_j,  c("joint_n", "joint_tp", "joint_fp", "joint_fn", "joint_tn",
                                              "src_joint", "deff_joint"), "RDT x microscopy (joint)")

dhs_prev_by_site1 <- Reduce(function(a, b) merge(a, b, by = c("country", "iso3c", "survey_year", "site_name_1")),
                            list(rdt_matched, micro_matched, joint_matched))

n_unmatched <- dhs_prev_by_site1[is.na(pfpr_rdt_denom) & is.na(pfpr_lm_denom), .N]
cat(sprintf("  %d/%d site name_1 x survey-year crosswalk rows had no matching DHS RDT/microscopy row (incl. survey years with no RDT/microscopy test at all, and a residual few surveys where REGNAME is a non-administrative classification such as Kenya's epidemiological zones)\n",
            n_unmatched, nrow(dhs_prev_by_site1)))

# site package's modelled prevalence surface, reloaded from the per-country
# cache written above (no re-fetch). s$prevalence is stratified by
# urban_rural -- confirmed directly (e.g. MDG/Antananarivo/2013: rural pfpr
# = 0.0316, urban = 0.0243) -- dropping that column without aggregating used
# to leave TWO rows per (country, iso3c, name_1, year), one per stratum, not
# one. Collapsed here to a single row via population-weighted mean (same
# convention cluster_fns.R::get_adm1_prevalence_weighted_from_site() already
# uses elsewhere for this identical urban/rural duplication problem),
# weighted by each species' own at-risk population
# (s$population$population_total's par_pf/par_pv, same urban_rural keying).
# weighted_mean_safe() falls back to a plain mean when a group's weights sum
# to zero (confirmed this happens for real -- several countries have
# zero/NA par_pf|par_pv for some strata, e.g. KEN/UGA/BDI) rather than
# propagating a divide-by-zero NaN into pfpr/pvpr.
weighted_mean_safe <- function(x, w) {
  w[is.na(w)] <- 0
  if (sum(w) == 0) mean(x, na.rm = TRUE) else stats::weighted.mean(x, w, na.rm = TRUE)
}
site_prev_list <- lapply(unique(site_dhs_matches$iso3c), function(iso) {
  path <- file.path(CACHE_DIR, paste0(iso, "_site.RDS"))
  if (!file.exists(path)) return(NULL)
  s <- readRDS(path)
  if (is.null(s$prevalence)) return(NULL)

  prev <- as.data.table(s$prevalence)[, .(country, iso3c, name_1, urban_rural, year, pfpr, pvpr)]
  pop  <- as.data.table(s$population$population_total)[, .(country, iso3c, name_1, urban_rural, year, par_pf, par_pv)]
  prev <- merge(prev, pop, by = c("country", "iso3c", "name_1", "urban_rural", "year"), all.x = TRUE)

  prev[, .(
    pfpr = weighted_mean_safe(pfpr, par_pf),
    pvpr = weighted_mean_safe(pvpr, par_pv)
  ), by = .(country, iso3c, name_1, year)]
})
site_prev <- rbindlist(Filter(Negate(is.null), site_prev_list), fill = TRUE)

# Joined on (iso3c, name_1, year) only -- NOT country, even though both
# sides have a country column. Confirmed by direct byte comparison: the
# site package's own s$prevalence$country for Cote d'Ivoire is "Côte
# d'Ivoire" (accented, UTF-8 0xc3 0xb4), while the raw DHS RDT/microscopy
# CSV's country column (ultimately from this file's own DHS_ISO3 list) is
# "Cote d'Ivoire" (plain ASCII) -- identical iso3c/name_1/year on both
# sides, but requiring an exact country-STRING match on top of that
# silently dropped every CIV row from this merge (all.x = TRUE just filled
# them with NA RDT/LM columns, which the "must have both diagnostics"
# exclusion below then removed with no visible error at all). iso3c already
# uniquely identifies the country on both sides -- country is redundant as
# a join key and, as this incident shows, fragile to accent/encoding
# variants a free-text field is never guaranteed to agree on.
matched_site_mis_prevalence <- merge(
  site_prev,
  dhs_prev_by_site1[, .(iso3c, site_name_1, survey_year,
                         pfpr_rdt_num, pfpr_rdt_denom, pfpr_lm_num, pfpr_lm_denom,
                         joint_n, joint_tp, joint_fp, joint_fn, joint_tn,
                         src_rdt, src_lm, src_joint, deff_rdt, deff_lm, deff_joint)],
  by.x = c("iso3c", "name_1", "year"),
  by.y = c("iso3c", "site_name_1", "survey_year"),
  all.x = TRUE
)
# Override with the canonical DHS_ISO3 spelling (see ISO3_TO_COUNTRY's own
# header at the top of this file) rather than trust site_prev's own
# (site-package-sourced) country string, which isn't guaranteed to agree
# with it -- fixes the root cause of the CIV accent-mismatch bug for every
# downstream consumer of matched_site_mis_prevalence.RDS's country column,
# not just this file's own merges.
matched_site_mis_prevalence[, country := ISO3_TO_COUNTRY[iso3c]]
setcolorder(matched_site_mis_prevalence, c("country", "iso3c", "name_1", "year", "pfpr", "pvpr",
                              "pfpr_rdt_num", "pfpr_rdt_denom", "pfpr_lm_num", "pfpr_lm_denom"))
setorder(matched_site_mis_prevalence, country, name_1, year)

# Every exclusion below is logged (site-years removed AND remaining, plus
# the underlying survey INDIVIDUALS removed AND remaining -- summed
# separately for RDT and LM, since a row's two denominators aren't always
# identical or fully overlapping) to matched_site_mis_prevalence_exclusions.
# csv, so the effect of each filter on the data actually feeding 04/05/06 is
# auditable without re-deriving it by hand. keep_mask is evaluated against
# the CURRENT state of matched_site_mis_prevalence at each step, so
# "remaining" is a running total and rows already dropped by an earlier step
# never appear in a later step's "removed" count.
# ALL_ISO3C: every country present BEFORE any exclusion -- fixed for the
# whole log so a country excluded entirely at an early step (e.g. Madagascar)
# still gets its own excluded_sites_<ISO> column (reading 0 thereafter,
# nothing left to exclude), rather than the column silently disappearing
# once that country has no rows left to compute a count from.
ALL_ISO3C <- sort(unique(matched_site_mis_prevalence$iso3c))
excluded_site_cols <- function(dt, kept) {
  # Per-country count of sites FULLY dropped at this step -- present
  # somewhere in dt (before this step's filter) but with NO row left in
  # kept (after it). Deliberately NOT "sites that lost >=1 row" -- a step
  # like the RDT-or-LM filter routinely removes SOME years from a site
  # while other years survive, which isn't "excluded" in any reasonable
  # sense (the site is still in matched_site_mis_prevalence). This
  # definition is the one that reconciles exactly with the
  # n_sites_remaining column's own step-to-step drop (verified: their sum
  # equals -diff(n_sites_remaining) at every step).
  before_sites <- unique(dt[,   .(iso3c, name_1)])
  after_sites  <- unique(kept[, .(iso3c, name_1)])
  dropped <- before_sites[!after_sites, on = .(iso3c, name_1)]
  counts <- setNames(rep(0L, length(ALL_ISO3C)), ALL_ISO3C)
  if (nrow(dropped) > 0) {
    by_iso <- dropped[, .(n = .N), by = iso3c]
    counts[by_iso$iso3c] <- by_iso$n
  }
  as.data.table(as.list(counts))[, setNames(.SD, paste0("excluded_sites_", ALL_ISO3C))]
}

exclusion_log <- list()
log_exclusion <- function(dt, keep_mask, label) {
  removed <- dt[!keep_mask]
  kept    <- dt[keep_mask]
  exclusion_log[[length(exclusion_log) + 1L]] <<- cbind(
    data.table::data.table(
      step                        = label,
      n_site_years_removed        = nrow(removed),
      n_individuals_rdt_removed   = sum(removed$pfpr_rdt_denom, na.rm = TRUE),
      n_individuals_lm_removed    = sum(removed$pfpr_lm_denom,  na.rm = TRUE),
      n_site_years_remaining      = nrow(kept),
      n_individuals_rdt_remaining = sum(kept$pfpr_rdt_denom, na.rm = TRUE),
      n_individuals_lm_remaining  = sum(kept$pfpr_lm_denom,  na.rm = TRUE),
      n_countries_remaining       = uniqueN(kept$iso3c),
      n_years_remaining           = uniqueN(kept$year),
      n_sites_remaining           = uniqueN(kept[, .(iso3c, name_1)]),
      n_country_years_remaining   = uniqueN(kept[, .(iso3c, year)])
    ),
    excluded_site_cols(dt, kept)
  )
  kept
}
exclusion_log[[1L]] <- cbind(
  data.table::data.table(
    step = "Before any exclusions", n_site_years_removed = 0L,
    n_individuals_rdt_removed = 0L, n_individuals_lm_removed = 0L,
    n_site_years_remaining = nrow(matched_site_mis_prevalence),
    n_individuals_rdt_remaining = sum(matched_site_mis_prevalence$pfpr_rdt_denom, na.rm = TRUE),
    n_individuals_lm_remaining  = sum(matched_site_mis_prevalence$pfpr_lm_denom,  na.rm = TRUE),
    n_countries_remaining       = uniqueN(matched_site_mis_prevalence$iso3c),
    n_years_remaining           = uniqueN(matched_site_mis_prevalence$year),
    n_sites_remaining           = uniqueN(matched_site_mis_prevalence[, .(iso3c, name_1)]),
    n_country_years_remaining   = uniqueN(matched_site_mis_prevalence[, .(iso3c, year)])
  ),
  excluded_site_cols(matched_site_mis_prevalence, matched_site_mis_prevalence)
)

##Removing instances with denominator == 1 -- run FIRST (order requested),
##before Nairobi/Madagascar and the RDT/LM pairing filter below, so this
##needs to be NA-safe on its own: pfpr_rdt_denom is still NA for plenty of
##rows at this point (any LM-only row, not yet excluded by the pairing
##filter further down) -- `NA != 1` is NA, not TRUE, which data.table's row
##subsetting would otherwise silently drop/mishandle. is.na(...) | ... != 1
##keeps an NA-RDT row here (nothing to flag as "denom == 1" yet) and lets
##the pairing filter below remove it on its own, real grounds instead.
matched_site_mis_prevalence <- log_exclusion(
  matched_site_mis_prevalence,
  is.na(matched_site_mis_prevalence$pfpr_rdt_denom) | matched_site_mis_prevalence$pfpr_rdt_denom != 1,
  "RDT denominator == 1 (unreliable single-respondent estimate)"
)
##Removing the one kenya location that has geographic bounds
matched_site_mis_prevalence <- log_exclusion(
  matched_site_mis_prevalence, matched_site_mis_prevalence$name_1 != 'Nairobi',
  "Nairobi (Kenya) -- geographic bounds issue"
)
##Removing Madagascar -- several survey rounds report by eco-climatic zone
##rather than admin region (no lexical/geometric match to site name_1 at
##all), and the fuzzy-match fallback that used to paper over this produced
##at least one confirmed wrong cross-region match (Antsiranana <- Atsinanana,
##a real but unrelated region) -- excluding the whole country rather than
##trusting the remaining per-year matches site by site.
matched_site_mis_prevalence <- log_exclusion(
  matched_site_mis_prevalence, matched_site_mis_prevalence$iso3c != 'MDG',
  "Madagascar -- eco-climatic zone survey rounds, no reliable admin1 match"
)
##Removing every (site, year) that only has ONE of RDT/microscopy -- 04/05/06
##need a paired RDT+LM year to fit against, not a single-diagnostic one.
##Supersedes (and is a strict superset of) the old Liberia-2016-only
##exclusion this file used to have here -- that round is RDT-only and gets
##caught here too, along with Senegal 2008, Uganda 2016, and Mali/Malawi's
##own recent RDT-only/LM-only rounds. Filtered here (the table 04/05/06
##actually read), not upstream in the raw RDT/microscopy pull, so
##03_model_rdt_individual.R's own inputs (dhs_mis_individual_hrp2.csv,
##countries_included.RDS) are untouched.
matched_site_mis_prevalence <- log_exclusion(
  matched_site_mis_prevalence,
  !is.na(matched_site_mis_prevalence$pfpr_rdt_denom) & !is.na(matched_site_mis_prevalence$pfpr_lm_denom),
  "Site-years missing RDT or LM (both diagnostics required)"
)
##Removing sites left with FEWER than 2 paired RDT+LM survey years -- a
##single time point carries no information about DECAY (04/05/06 fit HRP2
##deletion's effect on RDT-vs-LM divergence over time within a site), so a
##site with only one qualifying year contributes nothing to that model even
##though the row itself is individually valid. Site-level (not site-year),
##by construction: a site failing this drops ALL its remaining rows, not
##just the "extra" one.
site_paired_year_counts <- matched_site_mis_prevalence[
  , .(n_paired_years = uniqueN(year)), by = .(iso3c, name_1)
]
sites_with_2plus <- site_paired_year_counts[n_paired_years >= 2, .(iso3c, name_1)]
matched_site_mis_prevalence <- log_exclusion(
  matched_site_mis_prevalence,
  matched_site_mis_prevalence[, paste(iso3c, name_1)] %in% sites_with_2plus[, paste(iso3c, name_1)],
  "Sites with fewer than 2 years of paired RDT+LM measurements"
)

exclusion_summary <- data.table::rbindlist(exclusion_log)
cat("\nExclusion summary:\n")
print(exclusion_summary)

# --- report by geographical region (ESA/WCA, REGION_MAP -- see its own
# header at the top of this file) ---
# 1. Region rollup of each step's own per-country excluded-site columns
#    (sum the ISO columns belonging to each region) -- same steps/rows as
#    exclusion_summary above, just regrouped from 16 country columns down
#    to 2 region columns for a coarser read.
excl_iso_cols <- grep("^excluded_sites_", names(exclusion_summary), value = TRUE)
excl_iso      <- sub("^excluded_sites_", "", excl_iso_cols)
excl_region   <- REGION_MAP[names(DHS_ISO3)[match(excl_iso, DHS_ISO3)]]
region_exclusion_summary <- exclusion_summary[, c("step", "n_sites_remaining")]
for (r in unique(excl_region)) {
  cols <- excl_iso_cols[excl_region == r]
  region_exclusion_summary[[paste0("excluded_sites_", r)]] <- rowSums(exclusion_summary[, ..cols])
}
cat("\nExclusion summary by region (ESA/WCA):\n")
print(region_exclusion_summary)

# 2. Final-state summary: the dataset actually feeding 04/05/06, split by
#    region -- sites/site-years/individuals remaining, not excluded counts.
final_by_region <- matched_site_mis_prevalence[, .(
  n_countries     = uniqueN(iso3c),
  n_sites         = uniqueN(paste(iso3c, name_1)),
  n_site_years    = .N,
  n_individuals_rdt = sum(pfpr_rdt_denom, na.rm = TRUE),
  n_individuals_lm  = sum(pfpr_lm_denom,  na.rm = TRUE)
), by = .(africa_region = region_of(country))]
cat("\nFinal matched_site_mis_prevalence, by region:\n")
print(final_by_region)

saveRDS(matched_site_mis_prevalence, file.path(OUT_DIR, "matched_site_mis_prevalence.RDS"))

cat(sprintf(
  "\nWrote matched_site_mis_prevalence: %d rows (%d name_1 x year rows with matched DHS RDT and/or LM data) to %s\n",
  nrow(matched_site_mis_prevalence),
  matched_site_mis_prevalence[!is.na(pfpr_rdt_denom) | !is.na(pfpr_lm_denom), .N],
  OUT_DIR
))

# =============================================================================
# SECTION 2: Aggregate data and merge HRP2 deletion
# (formerly src/02b_aggregate_data.R)
#
# Reads every sub-regional aggregate output, merges into one wide table, then
# joins pfhrp2 deletion-risk scores by admin-1 name matching.
#
# Matching strategy (applied in order):
#   1. Exact (case-normalised, stripped punctuation)
#   2. Fuzzy  (Jaro-Winkler distance <= 0.12, closest match wins)
#   3. Unmatched rows receive NA for all pfhrp2 columns + match_flag = "no_match"
#
# Uganda caveat: pfhrp2 is at district level (~130 districts) while DHS
# sub_regions are aggregated regions ("Acholi", "Central 1", etc.). Fuzzy
# matching would find the closest-sounding district, which is often wrong.
# All Uganda joins therefore use an explicit district->region crosswalk
# instead (built below), not string distance.
#
# OUTPUT: dhs_mis_individual_hrp2.csv (dhs_pfhrp2_merged.csv after Section 3)
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 READ AND MERGE ALL DHS OUTPUTS
# -----------------------------------------------------------------------------
rdt2      <- fread(file.path(DIR_PULL, "dhs_mis_rdt.csv"))
micro2    <- fread(file.path(DIR_PULL, "dhs_mis_microscopy.csv"))
fever     <- fread(file.path(DIR_PULL, "dhs_mis_fever.csv"))
conc      <- fread(file.path(DIR_PULL, "dhs_mis_concordance.csv"))
careseek  <- fread(file.path(DIR_PULL, "dhs_mis_careseeking_formal.csv"))
act_cov   <- fread(file.path(DIR_PULL, "dhs_mis_act_coverage.csv"))

keys <- c("country", "survey_year", "sub_region")

# RDT -- base table
base <- copy(rdt2)
setnames(base, c("n_tested", "n_positive", "n_negative", "prev_rdt", "se_prev_rdt"),
               c("n_tested_rdt", "n_positive_rdt", "n_negative_rdt", "prev_rdt", "se_prev_rdt"))

# Microscopy
micro_j <- micro2[, .(country, survey_year, sub_region,
                      n_tested_micro  = n_tested,
                      n_positive_micro = n_positive,
                      prev_micro, se_prev_micro)]
base <- merge(base, micro_j, by = keys, all = TRUE)

# Fever prevalence
fever_j <- fever[, .(country, survey_year, sub_region,
                      n_children, n_fever, prev_fever, se_prev_fever)]
base <- merge(base, fever_j, by = keys, all = TRUE)

# RDT x Microscopy concordance (all febrile + non-febrile)
conc_j <- conc[, .(country, survey_year, sub_region,
                    n_both_tested,
                    pct_rdt_pos_micro_pos, pct_rdt_pos_micro_neg,
                    pct_rdt_neg_micro_pos, pct_rdt_neg_micro_neg)]
base <- merge(base, conc_j, by = keys, all = TRUE)

# Formal care-seeking
care_j <- careseek[, .(country, survey_year, sub_region,
                        n_febrile, n_formal, prop_formal_care, se_prop_formal_care)]
base <- merge(base, care_j, by = keys, all = TRUE)

# ACT coverage among febrile children
act_cov_j <- act_cov[, .(country, survey_year, sub_region,
                           n_act, prop_act, se_prop_act)]
base <- merge(base, act_cov_j, by = keys, all = TRUE)

cat(sprintf("Merged DHS table: %d rows x %d columns\n", nrow(base), ncol(base)))

# -----------------------------------------------------------------------------
# 2.2 LOAD pfhrp2 DELETION RISK AND FILTER TO STUDY COUNTRIES
# -----------------------------------------------------------------------------
iso_to_country <- c(
  KEN = "Kenya",     UGA = "Uganda",    RWA = "Rwanda",
  BDI = "Burundi",   MDG = "Madagascar",MWI = "Malawi",
  AGO = "Angola",    MOZ = "Mozambique",ZMB = "Zambia",
  NGA = "Nigeria",   GHA = "Ghana",     GIN = "Guinea",
  LBR = "Liberia",   SLE = "Sierra Leone", TGO = "Togo",
  MLI = "Mali",      BFA = "Burkina Faso", TZA = 'Tanzania',
  SEN = 'Senegal', NER = "Niger", CMR = "Cameroon",
  BEN = 'Benin', BWA = 'Botswana', TCD = 'Chad', COD = 'Congo Democratic Republic',
  CIV = "Cote d'Ivoire", ETH = 'Ethiopia', GAB = 'Gabon', GMB = 'Gambia',
  LSO = 'Lesotho', MRT = 'Mauritania', ZAF = 'South Africa', ZWE = 'Zimbabwe'
)

hrp2_raw <- fread(HRP2_RISK_CSV)
hrp2 <- hrp2_raw[iso %in% names(iso_to_country),
                  .(iso, name_1, hrp2_risk, hrp2_prospective_risk)]
hrp2[, country := iso_to_country[iso]]
cat(sprintf("pfhrp2 rows for study countries: %d across %d countries\n",
            nrow(hrp2), uniqueN(hrp2$country)))

# -----------------------------------------------------------------------------
# 2.3 NAME NORMALISATION
# -----------------------------------------------------------------------------
normalise <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9 ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

hrp2[,  name_1_norm   := normalise(name_1)]
base[,  sub_reg_norm  := normalise(sub_region)]

# Uganda explicit crosswalk: 2009/2014 DHS sub_region names -> name_1_2025 region names
uga_map_raw <- fread(UGA_NAME1_MAP_CSV)
uga_xwalk   <- unique(uga_map_raw[, .(sub_reg_norm = normalise(name_1),
                                       region_2025  = name_1_2025)])

# Uganda district -> region lookup (pfhrp2 is district-level; this maps each
# district to its name_1_2025 region so pfhrp2 risk can be aggregated per DHS region)
uga_district_region <- rbindlist(list(
  data.table(district = c("Agago","Amuru","Gulu","Kitgum","Lamwo","Nwoya","Omoro","Pader"),              region_2025 = "Acholi"),
  data.table(district = c("Buhweju","Bushenyi","Ibanda","Isingiro","Mbarara","Mitooma","Ntungamo","Rubirizi","Sheema"), region_2025 = "Ankole"),
  data.table(district = c("Bududa","Bulambuli","Manafwa","Mbale","Sironko"),                             region_2025 = "Bugisu"),
  data.table(district = c("Busia","Bukedea","Butaleja","Kibuku","Pallisa","Tororo"),                     region_2025 = "Bukedi"),
  data.table(district = c("Buliisa","Hoima","Kagadi","Kakumiro","Kibaale","Kiryandongo","Masindi"),      region_2025 = "Bunyoro"),
  data.table(district = c("Bugiri","Buyende","Iganga","Jinja","Kaliro","Kamuli","Luuka","Mayuge","Namayingo","Namutumba"), region_2025 = "Busoga"),
  data.table(district = "Kampala",                                                                       region_2025 = "Kampala"),
  data.table(district = c("Abim","Amudat","Kaabong","Kotido","Moroto","Nakapiripirit","Napak"),          region_2025 = "Karamoja"),
  data.table(district = c("Kabale","Kanungu","Kisoro","Rubanda","Rukungiri"),                            region_2025 = "Kigezi"),
  data.table(district = c("Alebtong","Amolatar","Apac","Dokolo","Kole","Lira","Otuke","Oyam"),           region_2025 = "Lango"),
  data.table(district = c("Kayunga","Kyankwanzi","Luwero","Mityana","Mubende","Nakaseke","Nakasongola","Wakiso"), region_2025 = "North Buganda"),
  data.table(district = c("Butambala","Buvuma","Gomba","Kalangala","Kalungu","Lwengo","Lyantonde","Masaka","Mpigi","Mukono","Rakai","Ssembabule"), region_2025 = "South Buganda"),
  data.table(district = c("Amuria","Kaberamaido","Katakwi","Kumi","Ngora","Serere","Soroti"),            region_2025 = "Teso"),
  data.table(district = c("Bundibugyo","Kabarole","Kamwenge","Kyegegwa","Kyenjojo","Ntoroko"),           region_2025 = "Tooro"),
  data.table(district = c("Adjumani","Arua","Koboko","Maracha","Moyo","Nebbi","Yumbe","Zombo"),          region_2025 = "West Nile"),
  # Elgon: 2024 DHS sub-region carved from eastern Bugisu
  data.table(district = c("Bukwo","Kapchorwa","Kween"),                                                  region_2025 = "Elgon")
))
uga_district_region[, district_norm := normalise(district)]
uga_district_region[, region_norm   := normalise(region_2025)]

# Kenya province -> county lookup (DHS 2015/2020 uses pre-2013 provinces;
# pfhrp2 is at county level -- aggregate risk across constituent counties)
ken_province_county <- rbindlist(list(
  data.table(county = c("Kiambu","Kirinyaga","Murang'a","Nyandarua","Nyeri"),
             province_norm = "central"),
  data.table(county = c("Kilifi","Kwale","Lamu","Mombasa","Taita Taveta","Tana River"),
             province_norm = "coast"),
  data.table(county = c("Embu","Isiolo","Kitui","Machakos","Makueni","Marsabit","Meru","Tharaka-Nithi"),
             province_norm = "eastern"),
  data.table(county = "Nairobi",
             province_norm = "nairobi"),
  data.table(county = c("Garissa","Mandera","Wajir"),
             province_norm = "north eastern"),
  data.table(county = c("Homa Bay","Kisii","Kisumu","Migori","Nyamira","Siaya"),
             province_norm = "nyanza"),
  data.table(county = c("Baringo","Bomet","Elgeyo-Marakwet","Kajiado","Kericho","Laikipia",
                         "Nakuru","Nandi","Narok","Samburu","Trans Nzoia","Turkana",
                         "Uasin Gishu","West Pokot"),
             province_norm = "rift valley"),
  data.table(county = c("Bungoma","Busia","Kakamega","Vihiga"),
             province_norm = "western")
))
ken_province_county[, county_norm := normalise(county)]

# Explicit overrides: DHS sub_region name variants that fail normalised exact
# matching (accents stripped to spaces, parenthetical suffixes, hyphens, or
# DHS city-level splits that map back to their parent region in pfhrp2).
# Keyed on (country, sub_reg_norm) -- checked before all other matching logic.
sub_region_overrides <- rbindlist(list(
  # Angola -- province names; "Bié" accent strips to "bi e" under normalise()
  data.table(country = "Angola",
             sub_reg_norm  = c("huambo","cuando cubango","lunda norte","bi e","bie","moxico"),
             pfhrp2_name_1 = c("Huambo","Cuando Cubango","Lunda Norte","Bie","Bie","Moxico")),
  # Burundi -- province names with encoding differences
  data.table(country = "Burundi",
             sub_reg_norm  = c("cibitoke","ngozi","rutana"),
             pfhrp2_name_1 = c("Cibitoke","Ngozi","Rutana")),
  # Cameroon -- DHS adds "(Excludes City)" suffixes, uses hyphens, or French/English name variants
  data.table(country = "Cameroon",
             sub_reg_norm  = c("adamaoua","littoral","littoral excludes douala",
                                "north west","north west","south west",
                                "centre excludes yaounde","yaounde","douala",
                                "east","north","south","west","far north"),
             pfhrp2_name_1 = c("Adamaoua","Littoral","Littoral",
                                "North West","North West","South West",
                                "Centre","Centre","Littoral",
                                "Est","Nord","Sud","Ouest","Extreme Nord")),
  # Nigeria -- state names with encoding differences
  data.table(country = "Nigeria",
             sub_reg_norm  = c("anambra","rivers","sokoto","delta",
                                "nassarawa","nasarawa","fct","abuja"),
             pfhrp2_name_1 = c("Anambra","Rivers","Sokoto","Delta",
                                "Nassarawa","Nassarawa","Abuja","Abuja")),
  # Togo -- region name with encoding difference
  data.table(country = "Togo",
             sub_reg_norm  = "kara",
             pfhrp2_name_1 = "Kara")
))

# Risk ordering for aggregation when multiple pfhrp2 districts map to one DHS region
RISK_ORDER <- c(Slight = 1L, Marginal = 2L, Moderate = 3L, High = 4L)

aggregate_risk <- function(risks) {
  risks <- risks[!is.na(risks)]
  if (length(risks) == 0) return(NA_character_)
  # Modal risk; break ties toward higher risk
  tab  <- sort(table(risks), decreasing = TRUE)
  top  <- names(tab)[tab == max(tab)]
  top[which.max(RISK_ORDER[top])]
}

# -----------------------------------------------------------------------------
# 2.4 MATCH pfhrp2 -> DHS sub_regions
# -----------------------------------------------------------------------------
# For each unique country x sub_region in base, find the best pfhrp2 match.

unique_regions <- unique(base[sub_region != "National",
                               .(country, sub_region, sub_reg_norm)])

# Angola: abbreviated province names used inside parenthetical zone labels
angola_abbrev <- c(
  "K. Norte"   = "Cuanza Norte",  "K. Sul"    = "Cuanza Sul",
  "K. Kubango" = "Cuando Cubango","L. Norte"  = "Lunda Norte",
  "L. Sul"     = "Lunda Sul",     "Malange"   = "Malanje",
  "Bié"        = "Bie",           "Bie"       = "Bie",
  "Cunene"     = "Cunene"
)

# Burundi: DHS zone -> province list (pfhrp2 is province-level)
# "rutana" is in Sud; user notes it sometimes appears in Ouest -- kept in Sud only.
# "bujumbura mairie" is a standalone DHS sub_region that exact-matches pfhrp2 directly,
# but also included in Ouest for when the DHS groups it there.
bdi_zone_provinces <- list(
  "centre est" = c("Gitega","Karuzi","Mwaro","Muramvya","Ruyigi","Cankuzo"),
  "nord"       = c("Kayanza","Kirundo","Muyinga","Ngozi"),
  "ouest"      = c("Bubanza","Bujumbura Mairie","Bujumbura Rural","Cibitoke","Rumonge"),
  "sud"        = c("Bururi","Makamba","Rutana")
)

# Nigeria: geopolitical zone -> state list (2010/2015 DHS uses zones; pfhrp2 is state-level)
# State names use pfhrp2 spellings (e.g. "Nassarawa", "Abuja")
nga_zone_states <- list(
  "north central" = c("Abuja","Benue","Kogi","Kwara","Nassarawa","Niger","Plateau"),
  "north east"    = c("Adamawa","Bauchi","Borno","Gombe","Taraba","Yobe"),
  "north west"    = c("Jigawa","Kaduna","Kano","Katsina","Kebbi","Sokoto","Zamfara"),
  "south east"    = c("Abia","Anambra","Ebonyi","Enugu","Imo"),
  "south south"   = c("Akwa Ibom","Bayelsa","Cross River","Delta","Edo","Rivers"),
  "south west"    = c("Ekiti","Lagos","Ogun","Ondo","Osun","Oyo")
)

# Portuguese zone names (2011 survey) -> constituent pfhrp2 province names
angola_zones <- list(
  "hiperendemica"         = c("Cabinda","Uige","Cuanza Norte","Malanje","Lunda Norte","Lunda Sul"),
  "mesoendemica estavel"  = c("Zaire","Luanda","Bengo","Benguela","Cuanza Sul","Huambo","Bie"),
  "mesoendemica instavel" = c("Moxico","Cuando Cubango","Cunene","Huila","Namibe")
)

match_pfhrp2 <- function(ctry, sreg_orig, sreg_norm) {
  pool <- hrp2[country == ctry]
  if (nrow(pool) == 0)
    return(data.table(hrp2_name_1 = NA_character_,
                      hrp2_risk = NA_character_,
                      hrp2_prospective_risk = NA_character_,
                      match_dist = NA_real_,
                      match_flag = "country_not_in_pfhrp2"))

  # Angola: parse province names from parenthetical zone labels
  if (ctry == "Angola") {
    paren <- regmatches(sreg_orig, regexpr("(?<=\\()([^)]+)(?=\\))", sreg_orig, perl = TRUE))
    if (length(paren) > 0) {
      # Split comma-separated province tokens, expand abbreviations
      parts    <- trimws(strsplit(paren, ",")[[1]])
      resolved <- sapply(parts, function(p) {
        if (p %in% names(angola_abbrev)) angola_abbrev[[p]]
        else {
          m <- pool[name_1_norm == normalise(p)]
          if (nrow(m) > 0) m$name_1[1] else p
        }
      }, USE.NAMES = FALSE)
      matched <- pool[name_1 %in% resolved]
      if (nrow(matched) > 0)
        return(data.table(
          hrp2_name_1           = paste(matched$name_1, collapse = ";"),
          hrp2_risk             = aggregate_risk(matched$hrp2_risk),
          hrp2_prospective_risk = aggregate_risk(matched$hrp2_prospective_risk),
          match_dist            = 0,
          match_flag            = "angola_zone_aggregated"
        ))
    }
    # Portuguese zone name without parentheses (2011 survey)
    if (sreg_norm %in% names(angola_zones)) {
      matched <- pool[name_1 %in% angola_zones[[sreg_norm]]]
      if (nrow(matched) > 0)
        return(data.table(
          hrp2_name_1           = paste(matched$name_1, collapse = ";"),
          hrp2_risk             = aggregate_risk(matched$hrp2_risk),
          hrp2_prospective_risk = aggregate_risk(matched$hrp2_prospective_risk),
          match_dist            = 0,
          match_flag            = "angola_zone_aggregated"
        ))
    }
    # Fall through to explicit override / exact / fuzzy for plain province names
  }

  # Explicit override -- checked before all other logic
  ov <- sub_region_overrides[country == ctry & sub_reg_norm == sreg_norm]
  if (nrow(ov) > 0) {
    matched <- pool[name_1 == ov$pfhrp2_name_1[1]]
    if (nrow(matched) > 0)
      return(data.table(hrp2_name_1           = matched$name_1[1],
                        hrp2_risk             = matched$hrp2_risk[1],
                        hrp2_prospective_risk = matched$hrp2_prospective_risk[1],
                        match_dist            = 0,
                        match_flag            = "explicit_override"))
  }

  # Uganda: crosswalk -> district lookup -> aggregate pfhrp2 risk by region
  if (ctry == "Uganda") {
    # Step 1: resolve DHS sub_region to name_1_2025 region name(s)
    xw <- uga_xwalk[sub_reg_norm == sreg_norm]
    if (nrow(xw) > 0) {
      # 2009/2014 DHS name (e.g. "Central 1") -- crosswalk gives region_2025
      target_regions <- normalise(xw$region_2025)
    } else if (sreg_norm %in% uga_district_region$region_norm) {
      # 2018/2024 DHS name already IS a name_1_2025 region (e.g. "Acholi")
      target_regions <- sreg_norm
    } else {
      # Unknown Uganda sub_region
      dists <- stringdist(sreg_norm, pool$name_1_norm, method = "jw")
      best  <- which.min(dists)
      return(data.table(hrp2_name_1 = pool$name_1[best],
                        hrp2_risk = pool$hrp2_risk[best],
                        hrp2_prospective_risk = pool$hrp2_prospective_risk[best],
                        match_dist = round(dists[best], 4),
                        match_flag = "uncertain_aggregation"))
    }

    # Step 2: find all pfhrp2 districts in those region(s) and aggregate risk
    districts <- uga_district_region[region_norm %in% target_regions, district_norm]
    matched   <- pool[name_1_norm %in% districts]

    if (nrow(matched) == 0) {
      dists <- stringdist(sreg_norm, pool$name_1_norm, method = "jw")
      best  <- which.min(dists)
      return(data.table(hrp2_name_1 = pool$name_1[best],
                        hrp2_risk = pool$hrp2_risk[best],
                        hrp2_prospective_risk = pool$hrp2_prospective_risk[best],
                        match_dist = round(dists[best], 4),
                        match_flag = "uncertain_aggregation"))
    }

    flag <- if (nrow(matched) > 1) "crosswalk_aggregated" else "crosswalk"
    return(data.table(
      hrp2_name_1           = paste(matched$name_1, collapse = ";"),
      hrp2_risk             = aggregate_risk(matched$hrp2_risk),
      hrp2_prospective_risk = aggregate_risk(matched$hrp2_prospective_risk),
      match_dist            = 0,
      match_flag            = flag
    ))
  }

  # Kenya: DHS uses pre-2013 provinces; pfhrp2 is at county level
  if (ctry == "Kenya") {
    kpc <- ken_province_county[province_norm == sreg_norm]
    if (nrow(kpc) > 0) {
      matched <- pool[name_1_norm %in% kpc$county_norm]
      if (nrow(matched) > 0) {
        return(data.table(
          hrp2_name_1           = paste(matched$name_1, collapse = ";"),
          hrp2_risk             = aggregate_risk(matched$hrp2_risk),
          hrp2_prospective_risk = aggregate_risk(matched$hrp2_prospective_risk),
          match_dist            = 0,
          match_flag            = "crosswalk_aggregated"
        ))
      }
    }
    # Province not in lookup -- fall back to fuzzy, flag uncertain
    dists <- stringdist(sreg_norm, pool$name_1_norm, method = "jw")
    best  <- which.min(dists)
    return(data.table(hrp2_name_1 = pool$name_1[best],
                      hrp2_risk = pool$hrp2_risk[best],
                      hrp2_prospective_risk = pool$hrp2_prospective_risk[best],
                      match_dist = round(dists[best], 4),
                      match_flag = "fuzzy_uncertain"))
  }

  # All other countries: exact first
  exact <- pool[name_1_norm == sreg_norm]
  if (nrow(exact) >= 1) {
    return(data.table(hrp2_name_1 = exact$name_1[1],
                      hrp2_risk = exact$hrp2_risk[1],
                      hrp2_prospective_risk = exact$hrp2_prospective_risk[1],
                      match_dist = 0,
                      match_flag = "exact"))
  }

  # Nigeria: if exact match failed, check geopolitical zone -> state aggregation
  if (ctry == "Nigeria" && sreg_norm %in% names(nga_zone_states)) {
    matched <- pool[name_1 %in% nga_zone_states[[sreg_norm]]]
    if (nrow(matched) > 0)
      return(data.table(
        hrp2_name_1           = paste(matched$name_1, collapse = ";"),
        hrp2_risk             = aggregate_risk(matched$hrp2_risk),
        hrp2_prospective_risk = aggregate_risk(matched$hrp2_prospective_risk),
        match_dist            = 0,
        match_flag            = "nga_zone_aggregated"
      ))
  }

  # Burundi: if exact match failed, check DHS zone -> province aggregation
  if (ctry == "Burundi" && sreg_norm %in% names(bdi_zone_provinces)) {
    matched <- pool[name_1 %in% bdi_zone_provinces[[sreg_norm]]]
    if (nrow(matched) > 0)
      return(data.table(
        hrp2_name_1           = paste(matched$name_1, collapse = ";"),
        hrp2_risk             = aggregate_risk(matched$hrp2_risk),
        hrp2_prospective_risk = aggregate_risk(matched$hrp2_prospective_risk),
        match_dist            = 0,
        match_flag            = "bdi_zone_aggregated"
      ))
  }

  # Fuzzy match (Jaro-Winkler)
  dists <- stringdist(sreg_norm, pool$name_1_norm, method = "jw")
  best  <- which.min(dists)
  flag  <- if (dists[best] < 0.30) "fuzzy_confident" else "fuzzy_uncertain"

  data.table(hrp2_name_1 = pool$name_1[best],
             hrp2_risk = pool$hrp2_risk[best],
             hrp2_prospective_risk = pool$hrp2_prospective_risk[best],
             match_dist = round(dists[best], 4),
             match_flag = flag)
}

hrp2_matches <- rbindlist(
  mapply(match_pfhrp2,
         unique_regions$country,
         unique_regions$sub_region,
         unique_regions$sub_reg_norm,
         SIMPLIFY = FALSE)
)
unique_regions <- cbind(unique_regions, hrp2_matches)

# Lesotho: entirely landlocked within South Africa -- override with ZAF modal risk
zaf_risks <- hrp2[country == "South Africa" & !is.na(hrp2_risk)]
if (nrow(zaf_risks) > 0) {
  zaf_modal       <- aggregate_risk(zaf_risks$hrp2_risk)
  zaf_modal_prosp <- aggregate_risk(zaf_risks$hrp2_prospective_risk)
  unique_regions[country == "Lesotho",
                 `:=`(hrp2_risk             = zaf_modal,
                      hrp2_prospective_risk  = zaf_modal_prosp,
                      hrp2_name_1            = "South Africa (surrounds Lesotho)",
                      match_dist             = NA_real_,
                      match_flag             = "lesotho_zaf_surrounds")]
}

# For non-exact matches: if every pfhrp2 region in that country shares the same
# risk category, the sub-region assignment is irrelevant -- use the uniform
# country value and promote the flag to confirm it.
country_uniform <- hrp2[!is.na(hrp2_risk),
                         .(n_distinct = uniqueN(hrp2_risk), uniform_risk = hrp2_risk[1]),
                         by = country][n_distinct == 1L]

NON_EXACT <- c("fuzzy_confident", "fuzzy_uncertain",
               "crosswalk_aggregated", "uncertain_aggregation")

unique_regions <- merge(unique_regions,
                        country_uniform[, .(country, uniform_risk)],
                        by = "country", all.x = TRUE)

unique_regions[match_flag %in% NON_EXACT & !is.na(uniform_risk),
               `:=`(hrp2_risk  = uniform_risk,
                    match_flag = paste0(match_flag, "+country_uniform"))]
unique_regions[, uniform_risk := NULL]

cat("\n--- Match flag summary ---\n")
print(unique_regions[, .N, by = match_flag][order(-N)])

cat("\n--- Fuzzy uncertain matches (review recommended) ---\n")
print(unique_regions[match_flag == "fuzzy_uncertain",
                     .(country, sub_region, hrp2_name_1, match_dist)])

# -----------------------------------------------------------------------------
# 2.5 JOIN MATCH TABLE ONTO BASE
# -----------------------------------------------------------------------------
join_cols <- c("hrp2_name_1", "hrp2_risk", "hrp2_prospective_risk",
               "match_dist", "match_flag")

# National rows get country-level summary (modal risk category)
country_modal_risk <- hrp2[, {
  tab <- table(hrp2_risk[!is.na(hrp2_risk)])
  if (length(tab) == 0) .(modal_risk = NA_character_)
  else .(modal_risk = names(which.max(tab)))
}, by = country]

# Lesotho national row: use South Africa modal (mirrors subnational override above)
if (exists("zaf_modal") && nrow(country_modal_risk[country == "Lesotho"]) > 0)
  country_modal_risk[country == "Lesotho", modal_risk := zaf_modal]

base <- merge(base,
              unique_regions[, c(keys[c(1,3)], join_cols), with = FALSE],
              by = c("country", "sub_region"), all.x = TRUE)

# Fill National rows with country modal risk
base <- merge(base, country_modal_risk, by = "country", all.x = TRUE)
base[sub_region == "National" & is.na(hrp2_risk),
     `:=`(hrp2_risk             = modal_risk,
          hrp2_prospective_risk  = NA_character_,
          hrp2_name_1            = "national_modal",
          match_dist             = NA_real_,
          match_flag             = "national_modal_risk")]
base[, modal_risk := NULL]
base[, sub_reg_norm := NULL]

# -----------------------------------------------------------------------------
# 2.6 NATIONAL / SUBNATIONAL FLAG
# -----------------------------------------------------------------------------
base[, level := fifelse(sub_region == "National", "national", "subnational")]

# -----------------------------------------------------------------------------
# 2.7 REGION LABEL
# -----------------------------------------------------------------------------
# REGION_MAP is a shared top-of-file constant now (see its own header there
# for the full ESARO/WCARO coverage and why it moved) -- was previously
# redefined here as a narrower, 17-country local copy.
base[, africa_region := region_of(country)]

# -----------------------------------------------------------------------------
# 2.8 CARE-OR-ACT AGGREGATE + INDIVIDUAL-LEVEL HRP2 JOIN
# Proportion of febrile children who sought formal care OR received ACT in the
# past two weeks. Uses survey sample weights but NOT full svydesign (no
# stratification variable available in individual CSV) -- weighted mean is an
# approximation relative to the design-based estimates for other indicators.
# Denominator: febrile children with BOTH formal_care and act_treatment non-NA.
#
# Read once and reused for both the care-or-act aggregate AND the individual-
# level hrp2_risk join below (previously two separate freads of the same
# 500k-row file).
# -----------------------------------------------------------------------------
ind_raw <- fread(file.path(DIR_PULL, "dhs_mis_individual.csv"))
ind_raw[, care_or_act := as.integer(formal_care == 1L | act_treatment == 1L)]

care_or_act_sr <- ind_raw[
  fever == 1L & !is.na(formal_care) & !is.na(act_treatment),
  .(
    n_febrile_with_tx = .N,
    n_care_or_act     = sum(care_or_act),
    prop_care_or_act  = weighted.mean(care_or_act, w = sample_weight)
  ),
  by = .(country, sub_region, survey_year)
]

base <- merge(base, care_or_act_sr[, .(country, sub_region, survey_year,
                                        n_febrile_with_tx, n_care_or_act,
                                        prop_care_or_act)],
              by = keys, all.x = TRUE)

n_coa_matched <- base[sub_region != "National" & !is.na(prop_care_or_act), .N]
cat(sprintf("prop_care_or_act: %d / %d subnational rows have data\n",
            n_coa_matched, base[sub_region != "National", .N]))

# dhs_pfhrp2_merged.csv is written once, at the end of Section 3

# Re-uses the unique_regions lookup already built above -- same matching
# logic, no extra name matching needed.
hrp2_lookup <- unique(unique_regions[,
  .(country, sub_region, hrp2_risk, hrp2_prospective_risk, match_flag)])

ind_hrp2 <- merge(ind_raw, hrp2_lookup,
                  by = c("country", "sub_region"), all.x = TRUE)

n_matched <- sum(!is.na(ind_hrp2$hrp2_risk))
cat(sprintf(
  "\nIndividual-level + HRP2: %d children  |  %d (%.1f%%) matched to hrp2_risk\n",
  nrow(ind_hrp2), n_matched, 100 * n_matched / nrow(ind_hrp2)
))
if (anyNA(ind_hrp2$hrp2_risk)) {
  unmatched <- ind_hrp2[is.na(hrp2_risk), unique(paste(country, sub_region))]
  cat(sprintf("  Unmatched sub-regions (%d): %s\n",
              length(unmatched), paste(unmatched, collapse = "; ")))
}

fwrite(ind_hrp2, file.path(OUT_DIR, "dhs_mis_individual_hrp2.csv"))
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "dhs_mis_individual_hrp2.csv")))

# -----------------------------------------------------------------------------
# 2.9 WRITE UGANDA DISTRICT-LEVEL HRP2 RISK FILE
# For use by 09_malariasimulation_integration.R (path_hrp2_risk).
# The raw file has four levels: Slight, Marginal, Moderate, High.
# Three output levels: "Slight" | "Marginal/Moderate" | "High".
# Marginal and Moderate are collapsed; Slight is kept separate.
#
# 32 districts in the SNT site file are absent from the pfhrp2 source.
# Fallback rules applied in order:
#   1. City districts -> "Slight/Marginal" (urban areas have lower deletion burden)
#   2. Carved-out districts -> risk of the parent district they split from
#   3. Kiruhura, Kazo -> Ankole regional modal (all 9 known Ankole districts = High)
# -----------------------------------------------------------------------------
uga_hrp2 <- hrp2_raw[iso == "UGA", .(district = name_1, hrp2_risk)]
uga_hrp2[, hrp2_risk := fcase(
  hrp2_risk == "Slight",                   "Slight",
  hrp2_risk %in% c("Marginal", "Moderate"), "Marginal/Moderate",
  hrp2_risk == "High",                      "High"
)]

# Districts missing from pfhrp2 source -- fallback assignments
uga_fallback <- rbindlist(list(
  # --- Cities -> Slight (urban areas have lower deletion burden) ---
  data.table(district  = c("Arua City", "Fort Portal City", "Gulu City",
                            "Hoima City", "Jinja City", "Lira City",
                            "Masaka City", "Mbale City", "Mbarara City",
                            "Soroti City"),
             hrp2_risk = "Slight",
             fallback  = "city"),

  # --- Carved-out districts: risk from parent district (pfhrp2 lookup) ---
  # Parent -> risk (raw) -> mapped:
  #   Iganga=Marginal->M/M, Kabarole=High, Kaberamaido=Marginal->M/M, Amuria=Marginal->M/M,
  #   Kotido=Marginal->M/M, Mubende=High, Kamwenge=High, Apac=Marginal->M/M, Rakai=High,
  #   Moyo=High, Nakapiripirit=Marginal->M/M, Kabale=High, Mbarara=High,
  #   Nebbi=Marginal->M/M
  data.table(district  = c("Bugweri",     # <- Iganga      (Marginal)
                            "Bunyangabu",  # <- Kabarole    (High)
                            "Kalaki",      # <- Kaberamaido (Marginal)
                            "Kapelebyong", # <- Amuria      (Marginal)
                            "Karenga",     # <- Kotido      (Marginal)
                            "Kassanda",    # <- Mubende     (High)
                            "Kitagwenda",  # <- Kamwenge    (High)
                            "Kwania",      # <- Apac        (Marginal)
                            "Kyotera",     # <- Rakai       (High)
                            "Nabilatuk",   # <- Nakapiripirit (Marginal)
                            "Obongi",      # <- Moyo        (High)
                            "Pakwach",     # <- Nebbi       (Marginal)
                            "Rukiga",      # <- Kabale      (High)
                            "Rwampara"),   # <- Mbarara     (High)
             hrp2_risk = c("Marginal/Moderate", "High", "Marginal/Moderate", "Marginal/Moderate",
                            "Marginal/Moderate", "High", "High",              "Marginal/Moderate",
                            "High", "Marginal/Moderate", "High",              "Marginal/Moderate",
                            "High", "High"),
             fallback  = "parent_district"),

  # Parent district maps to Slight or Moderate
  data.table(district  = c("Butebo",      # <- Pallisa     (Slight)
                            "Kikuube",     # <- Hoima + Kibaale (Moderate)
                            "Madi-Okollo", # <- Arua        (Moderate)
                            "Namisindwa",  # <- Manafwa     (Moderate)
                            "Terego"),     # <- Arua        (Moderate)
             hrp2_risk = c("Slight", "Marginal/Moderate", "Marginal/Moderate",
                            "Marginal/Moderate", "Marginal/Moderate"),
             fallback  = "parent_district"),

  # --- Ankole regional modal (Kiruhura pre-dates pfhrp2 era; all 9 known
  #     Ankole districts are High; Kazo was carved from Kiruhura) ---
  data.table(district  = c("Kiruhura", "Kazo"),
             hrp2_risk = "High",
             fallback  = "region_modal_Ankole")
))

uga_hrp2 <- rbind(uga_hrp2[, fallback := "pfhrp2_source"],
                  uga_fallback,
                  fill = TRUE)

cat(sprintf(
  "Uganda HRP2 risk: %d districts (%d from pfhrp2 source, %d fallback)\n",
  nrow(uga_hrp2),
  uga_hrp2[fallback == "pfhrp2_source", .N],
  uga_hrp2[fallback != "pfhrp2_source", .N]
))
cat("Fallback breakdown:\n")
print(uga_hrp2[fallback != "pfhrp2_source", .N, by = fallback])

cat(sprintf(
  "Saved: %d Slight  %d Marginal/Moderate  %d High\n",
  uga_hrp2[hrp2_risk == "Slight",            .N],
  uga_hrp2[hrp2_risk == "Marginal/Moderate", .N],
  uga_hrp2[hrp2_risk == "High",              .N]
))

# -----------------------------------------------------------------------------
# Bridge into Section 3 -- reuse Section 2's in-memory objects rather than
# re-reading the CSVs just written back off disk.
# -----------------------------------------------------------------------------
base_merged <- base
ind         <- ind_hrp2

# =============================================================================
# SECTION 3: Extract site seasonality
# (formerly src/02c_site_seasonality.R)
#
# Adds survey timing, RDT brand, site-package seasonality, and site-package
# SMC coverage to the datasets Section 2 just built.
#
# New individual-level columns:
#   site_seasonality_type, site_peak_month_1/2, site_match_level,
#   site_months_to_peak, site_smc_cov (area-level, NOT individual receipt --
#   every child in the same sub_region-survey shares the same value),
#   prev_micro, prev_rdt, prop_care_or_act (sub-regional, for
#   03_model_rdt_individual.R)
#
# New aggregate-level columns (in addition to the above site_* columns):
#   survey_month_mean, rdt_brand_modal, site_months_to_peak_mean
# =============================================================================

# -----------------------------------------------------------------------------
# 3.1 HELPERS
# -----------------------------------------------------------------------------

# Circular mean of month integers 1-12.
# Encoding: theta = 2*pi*(m-1)/12, so January(1) -> 0, July(7) -> pi.
# Inverse:  m = round(angle/(2*pi)*12) mod 12 + 1
# Returns the nearest integer month (1-12).
circ_mean_month <- function(months) {
  m <- as.integer(months[!is.na(months)])
  if (length(m) == 0L) return(NA_integer_)
  theta <- 2 * pi * (m - 1L) / 12
  angle <- atan2(mean(sin(theta)), mean(cos(theta)))   # -pi to pi
  raw   <- angle / (2 * pi) * 12                       # -6 to 6
  (as.integer(round(raw)) %% 12L) + 1L                 # 1-12
}

# Circular difference: peak_month - survey_month, wrapped to [-6, 6].
# Positive = survey is before the peak; negative = survey is after.
circ_diff_month <- function(peak, survey) {
  d <- as.integer(peak) - as.integer(survey)
  ifelse(d >  6L, d - 12L,
  ifelse(d < -6L, d + 12L, d))
}

# Modal value (first of ties).
modal <- function(x) {
  x <- x[!is.na(x) & x != "" & x != "na"]
  if (length(x) == 0L) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1L]
}

# -----------------------------------------------------------------------------
# 3.2 Classify seasonality from monthly_rainfall
# monthly_avg: named numeric of length 12 (months 1-12), mm.
# CV threshold for "no seasonal trend": if rainfall is too uniform month-to-month.
# Peak significance threshold: a local max is kept only if >= 65% of global max.
# -----------------------------------------------------------------------------

classify_seasonality <- function(monthly_avg) {
  stopifnot(length(monthly_avg) == 12L)

  # Coefficient of variation
  mn <- mean(monthly_avg, na.rm = TRUE)
  if (is.na(mn) || mn == 0) {
    return(list(type = "no_seasonal_trend", peak1 = NA_integer_, peak2 = NA_integer_))
  }
  cv <- sd(monthly_avg, na.rm = TRUE) / mn
  if (cv < 0.25) {
    return(list(type = "no_seasonal_trend", peak1 = NA_integer_, peak2 = NA_integer_))
  }

  # Local maxima (circular neighbours)
  n <- 12L
  is_local_max <- sapply(seq_len(n), function(i) {
    prev_i <- ((i - 2L) %% n) + 1L
    next_i <- (i        %% n) + 1L
    monthly_avg[i] > monthly_avg[prev_i] & monthly_avg[i] > monthly_avg[next_i]
  })

  peak_idx <- which(is_local_max)
  if (length(peak_idx) == 0L) peak_idx <- which.max(monthly_avg)

  # Retain only peaks >= 65% of global maximum
  max_rain  <- max(monthly_avg)
  sig_peaks <- peak_idx[monthly_avg[peak_idx] >= 0.65 * max_rain]
  if (length(sig_peaks) == 0L) sig_peaks <- which.max(monthly_avg)

  if (length(sig_peaks) >= 2L) {
    # Keep two strongest peaks
    top2 <- sig_peaks[order(-monthly_avg[sig_peaks])][1:2]
    top2 <- sort(top2)
    return(list(type = "bimodal", peak1 = top2[1L], peak2 = top2[2L]))
  }

  return(list(type = "unimodal", peak1 = sig_peaks[1L], peak2 = NA_integer_))
}

# WHO SMC-style seasonality index: the max 3-consecutive-month (~91-day) share
# of a site's annual rainfall total, computed from the site package's own
# Fourier-reconstructed DAILY rainfall profile (s$seasonality$fourier_prediction
# -- distinct from classify_seasonality()'s 12-point monthly average above).
# >=0.60 is the standard WHO SMC "highly seasonal transmission" threshold.
# Confirmed sensible against every currently-cached admin-1 unit (494 units,
# 30 countries): the top of the real range is exactly the classic Sahelian SMC
# belt (Chad/Mauritania/Niger/Senegal/Mali, index 0.82-0.86) and the bottom is
# Congo-basin/Lake-Victoria equatorial rainforest (DRC/Uganda/Kenya, index
# 0.32-0.34) -- the WHO threshold cleanly separates the two.
# A raw peak/trough RATIO was tried first and rejected: many admin-1 units'
# fitted daily profile touches ~0 rainfall at the annual minimum, which blows
# the ratio up to 5-6 orders of magnitude and makes it useless as a
# discriminator -- a share-of-annual-total index doesn't have that
# divide-by-near-zero failure mode.
seasonality_index <- function(t, profile) {
  ord <- order(t)
  t <- t[ord]; profile <- pmax(profile[ord], 0)
  n <- length(t)
  if (n < 12L) return(NA_real_)
  win   <- max(1L, round(n * 91 / 365))
  total <- sum(profile)
  if (total <= 0) return(NA_real_)
  ext  <- c(profile, profile[seq_len(win - 1L)])   # circular wrap
  roll <- vapply(seq_len(n), function(i) sum(ext[i:(i + win - 1L)]), numeric(1))
  max(roll) / total
}

# -----------------------------------------------------------------------------
# 3.3 Identify, fetch, and load site data for all DHS countries
# Strategy:
#   1. DHS_ISO3 (defined at top of file) already covers all DHS countries.
#   2. Call available_sites() to find which ISO3 codes have site package data.
#   3. For each matched country: load from cache if present, else fetch_site()
#      and save to cache.
#   4. Uganda: use UGA_SNT_SITE_FILE (SNT file) which has name_1 = DHS region names.
#
# Uses its own cache directory (site_cache_seasonality), separate from
# Section 1's site_cache: Section 1 fetches with urban_rural = TRUE (needed
# for its admin-1 shape matching); this section doesn't care about the
# urban/rural split, but keeping the caches separate avoids one section's
# fetch_site() call silently overwriting the other's differently-shaped
# cached object under the same iso3c filename.
# -----------------------------------------------------------------------------
CACHE_DIR_SEASONALITY <- file.path(OUT_DIR, "site_cache_seasonality")
dir.create(CACHE_DIR_SEASONALITY, showWarnings = FALSE, recursive = TRUE)

# Drop any columns written by a previous run so all merges start from a clean slate
SITE_COLS_ALL <- c("survey_month_mean", "rdt_brand_modal",
                   "site_seasonality_type", "site_peak_month_1", "site_peak_month_2",
                   "site_match_level", "site_months_to_peak_mean", "site_smc_cov",
                   "site_prev_pfpr", "site_seasonality_index")
stale <- intersect(names(base_merged), SITE_COLS_ALL)
if (length(stale) > 0L) base_merged[, (stale) := NULL]
dhs_countries <- sort(unique(base_merged$country))
dhs_iso3      <- DHS_ISO3[dhs_countries]
dhs_iso3      <- dhs_iso3[!is.na(dhs_iso3)]   # drop any not in lookup

# Query the site package registry
cat("Querying available_sites()...\n")
avail <- tryCatch(
  available_sites(),
  error = function(e) {
    warning("available_sites() failed (", conditionMessage(e),
            ") -- will use cached files only.")
    NULL
  }
)

if (!is.null(avail)) {
  avail_iso3 <- unique(as.character(avail$iso3c))
  to_fetch   <- intersect(dhs_iso3, avail_iso3)
  no_data    <- setdiff(dhs_countries, names(DHS_ISO3[DHS_ISO3 %in% to_fetch]))
  cat(sprintf("  Site data available for %d / %d DHS countries\n",
              length(to_fetch), length(dhs_iso3)))
  if (length(no_data) > 0)
    cat("  No site data for:", paste(sort(no_data), collapse = ", "), "\n")
} else {
  # Fallback: use whatever is already cached
  cached_iso3 <- sub("_site\\.RDS$", "", list.files(CACHE_DIR_SEASONALITY, pattern="\\.RDS$"))
  to_fetch    <- intersect(dhs_iso3, cached_iso3)
  cat(sprintf("  Falling back to %d cached files\n", length(to_fetch)))
}

# Load each country: from cache if present, fetch_site() + save otherwise
snt_path     <- UGA_SNT_SITE_FILE
uga_from_snt <- file.exists(snt_path)

# SMC (seasonal malaria chemoprevention) coverage extraction -- the site
# package returns interventions data in two different shapes depending on
# the cache:
#   (a) flat annual data frame with a smc_cov column directly (name_1, year,
#       smc_cov one row each) -- e.g. KEN, BFA, GHA
#   (b) nested list ($smc$implementation) with one row per round (and
#       sometimes per urban/rural stratum) -- e.g. NER, SEN, TCD, CMR, and
#       Uganda's SNT file. Collapsed to one annual mean per name_1 x year.
# Not Uganda-specific: about half the cached RDS files use the nested format.
extract_smc <- function(s, iso, ctry_name) {
  if (!is.null(s$interventions) && is.data.frame(s$interventions) &&
      "smc_cov" %in% names(s$interventions)) {
    smc_dt <- as.data.table(s$interventions)
    return(smc_dt[, .(country = ctry_name, iso3c = iso, name_1, year, smc_cov)])
  }
  impl <- s$interventions$smc$implementation
  if (!is.null(impl)) {
    smc_dt <- as.data.table(impl)[, .(smc_cov = mean(smc_cov, na.rm = TRUE)),
                                  by = .(name_1, year)]
    smc_dt[, `:=`(country = ctry_name, iso3c = iso)]
    return(smc_dt[, .(country, iso3c, name_1, year, smc_cov)])
  }
  NULL
}

# Site-package modelled prevalence (PfPR/PvPR) -- a flat data frame in every
# cached file seen so far (unlike interventions' two shapes above), keyed by
# name_1 x year, with an extra name_2 (district) level in Uganda's SNT file
# only. Aggregating by name_1 x year handles both: it's a no-op where name_2
# doesn't exist, and a mean-across-districts collapse where it does.
extract_prev <- function(s, iso, ctry_name) {
  df <- s$prevalence
  if (is.null(df) || !is.data.frame(df) || !"pfpr" %in% names(df)) return(NULL)
  dt <- as.data.table(df)[, .(pfpr = mean(pfpr, na.rm = TRUE)), by = .(name_1, year)]
  dt[, `:=`(country = ctry_name, iso3c = iso)]
  dt[, .(country, iso3c, name_1, year, pfpr)]
}

# Seasonality amplitude (site_seasonality_index, see seasonality_index() above)
# -- climatological, no year dimension, same as rain/classify_seasonality's own
# site_seasonality_type. Uses fourier_prediction (t, profile), a different
# sub-table of s$seasonality than monthly_rainfall.
extract_seas_idx <- function(s, iso, ctry_name) {
  df <- s$seasonality$fourier_prediction
  if (is.null(df) || !is.data.frame(df) || !all(c("t", "profile") %in% names(df))) return(NULL)
  dt <- as.data.table(df)[, .(site_seasonality_index = seasonality_index(t, profile)), by = .(name_1)]
  dt[, `:=`(country = ctry_name, iso3c = iso)]
  dt[, .(country, iso3c, name_1, site_seasonality_index)]
}

cat("Loading site data...\n")
# Extracts rainfall (for seasonality), SMC coverage (for site_smc_cov), and
# modelled prevalence (for site_prev_pfpr) from the same loaded/cached site
# object, so each RDS is only fetched/read once.
site_list <- lapply(to_fetch, function(iso) {
  # Uganda: SNT file preferred (name_1 aligns with DHS regions)
  if (iso == "UGA" && uga_from_snt) return(list(rain = NULL, smc = NULL, prev = NULL))

  cache_path <- file.path(CACHE_DIR_SEASONALITY, paste0(iso, "_site.RDS"))
  if (file.exists(cache_path)) {
    s <- readRDS(cache_path)
  } else {
    ctry <- names(DHS_ISO3)[DHS_ISO3 == iso]
    cat(sprintf("  Fetching %s (%s) from site package...\n", ctry, iso))
    s <- tryCatch(fetch_site(iso), error = function(e) {
      warning(sprintf("fetch_site('%s') failed: %s", iso, conditionMessage(e)))
      NULL
    })
    if (is.null(s)) return(list(rain = NULL, smc = NULL, prev = NULL))
    saveRDS(s, cache_path)
  }

  rain <- if (!is.null(s$seasonality$monthly_rainfall)) {
    rain_dt <- as.data.table(s$seasonality$monthly_rainfall)
    rain_dt[, .(country = names(DHS_ISO3)[DHS_ISO3 == iso], iso3c = iso,
                name_1, year, month, rainfall)]
  } else NULL

  smc     <- extract_smc(s, iso, names(DHS_ISO3)[DHS_ISO3 == iso])
  prev    <- extract_prev(s, iso, names(DHS_ISO3)[DHS_ISO3 == iso])
  seas_idx <- extract_seas_idx(s, iso, names(DHS_ISO3)[DHS_ISO3 == iso])

  list(rain = rain, smc = smc, prev = prev, seas_idx = seas_idx)
})
rain_list     <- lapply(site_list, `[[`, "rain")
smc_list      <- lapply(site_list, `[[`, "smc")
prev_list     <- lapply(site_list, `[[`, "prev")
seas_idx_list <- lapply(site_list, `[[`, "seas_idx")

# Uganda SNT
if (uga_from_snt && "UGA" %in% to_fetch) {
  snt     <- readRDS(snt_path)
  rain_uga <- as.data.table(snt$seasonality$monthly_rainfall)
  rain_uga <- rain_uga[, .(country = "Uganda", iso3c = "UGA",
                            rainfall = mean(rainfall, na.rm = TRUE)),
                        by = .(name_1, year, month)]
  rain_list <- c(rain_list, list(rain_uga))
  cat("  Uganda: using SNT file (name_1 = DHS region level)\n")

  smc_uga <- extract_smc(snt, "UGA", "Uganda")
  if (!is.null(smc_uga)) smc_list <- c(smc_list, list(smc_uga))

  prev_uga <- extract_prev(snt, "UGA", "Uganda")
  if (!is.null(prev_uga)) prev_list <- c(prev_list, list(prev_uga))

  seas_idx_uga <- extract_seas_idx(snt, "UGA", "Uganda")
  if (!is.null(seas_idx_uga)) seas_idx_list <- c(seas_idx_list, list(seas_idx_uga))
}

rain_all <- rbindlist(rain_list[!sapply(rain_list, is.null)], fill = TRUE)
cat(sprintf("  Loaded rainfall data for %d countries\n", uniqueN(rain_all$iso3c)))

smc_all <- rbindlist(smc_list[!sapply(smc_list, is.null)], fill = TRUE)
cat(sprintf("  Loaded SMC coverage data for %d countries\n", uniqueN(smc_all$iso3c)))

prev_all <- rbindlist(prev_list[!sapply(prev_list, is.null)], fill = TRUE)
cat(sprintf("  Loaded modelled prevalence data for %d countries\n", uniqueN(prev_all$iso3c)))

seas_idx_all <- rbindlist(seas_idx_list[!sapply(seas_idx_list, is.null)], fill = TRUE)
cat(sprintf("  Loaded seasonality index data for %d countries\n", uniqueN(seas_idx_all$iso3c)))

# Average across years (2000-2024) -> mean monthly profile per name_1
mean_rain <- rain_all[, .(rainfall = mean(rainfall, na.rm = TRUE)), by = .(country, iso3c, name_1, month)]
setorder(mean_rain, country, name_1, month)

# Classify each name_1
cat("Classifying seasonality patterns...\n")
seasonality_admin1 <- mean_rain[, {
  avg <- rainfall[order(month)]   # ensure month 1..12 order
  res <- classify_seasonality(avg)
  .(site_seasonality_type = res$type,
    site_peak_month_1     = as.integer(res$peak1),
    site_peak_month_2     = as.integer(res$peak2))
}, by = .(country, iso3c, name_1)]

# Merge in the seasonality index (site_seasonality_index) computed above --
# climatological, no year dimension, same key as site_seasonality_type.
seasonality_admin1 <- merge(seasonality_admin1, seas_idx_all,
                            by = c("country", "iso3c", "name_1"), all.x = TRUE)

# Country-level fallback: aggregate rainfall across all name_1s, reclassify
mean_rain_country <- mean_rain[, .(rainfall = mean(rainfall, na.rm = TRUE)), by = .(country, iso3c, month)]
seasonality_country <- mean_rain_country[, {
  avg <- rainfall[order(month)]
  res <- classify_seasonality(avg)
  .(site_seasonality_type = res$type,
    site_peak_month_1     = as.integer(res$peak1),
    site_peak_month_2     = as.integer(res$peak2))
}, by = .(country, iso3c)]
seasonality_country[, site_match_level := "country"]

# Country-level seasonality index fallback: mean across that country's own
# admin-1 indices (no country-level Fourier fit exists to recompute directly
# from, unlike rainfall which can just be re-averaged and reclassified above).
seas_idx_country <- seas_idx_all[, .(site_seasonality_index = mean(site_seasonality_index, na.rm = TRUE)),
                                 by = .(country, iso3c)]
seasonality_country <- merge(seasonality_country, seas_idx_country,
                             by = c("country", "iso3c"), all.x = TRUE)

cat(sprintf("  Classified %d admin-1 units across %d countries\n",
            nrow(seasonality_admin1), uniqueN(seasonality_admin1$country)))
cat(sprintf("  Seasonality types: %s\n",
            paste(sort(table(seasonality_admin1$site_seasonality_type)), names(sort(table(seasonality_admin1$site_seasonality_type))), sep = " ", collapse = " | ")))

# -----------------------------------------------------------------------------
# 3.4 Build DHS sub_region -> site name_1 crosswalk
# Matching hierarchy (applied in order, stops at first success):
#   1. Exact clean match           -> "admin1"
#   2. Contains / substring match  -> "admin1_partial"
#   3. Levenshtein fuzzy match     -> "admin1_fuzzy"
#   4. Country fallback, tagged by within-country variation:
#        "country_uniform"  -- peak months consistent across admin-1 units (SD < 1 mo)
#        "country_variable" -- meaningful spatial variation; fallback loses information
# Language alias tables translate French (Cameroon) and Kinyarwanda (Rwanda) site
# names to English before matching.
# -----------------------------------------------------------------------------

clean_str <- function(x) {
  x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))   # transliterate accents (e.g. é->e, â->a, ô->o)
  x <- gsub("[^a-z0-9 ]", " ", x)                  # strip hyphens, parentheses, etc.
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

# Language alias tables: site name_1 -> English equivalent
# Applied before any matching step; covers French (CMR) and Kinyarwanda (RWA).
SITE_ALIASES <- list(
  # Cameroon -- French -> English
  "est"           = "east",
  "nord"          = "north",
  "ouest"         = "west",
  "sud"           = "south",
  "nord ouest"    = "north west",
  "sud ouest"     = "south west",
  "extreme nord"  = "far north",
  # Rwanda -- Kinyarwanda -> English
  "amajyaruguru"    = "north",
  "amajyepfo"       = "south",
  "iburasirazuba"   = "east",
  "iburengerazuba"  = "west",
  "umujyi wa kigali" = "kigali city"
)

translate_site <- function(x) {
  cl    <- clean_str(x)
  alias <- SITE_ALIASES[[cl]]
  if (is.null(alias)) cl else alias
}

# Circular SD of peak months (in months; range 0-6).
# Uses mean resultant length formula: SD_circ = sqrt(-2 ln R_bar) * 12/(2*pi).
circ_sd_month <- function(months) {
  m <- as.integer(months[!is.na(months)])
  if (length(m) < 2L) return(0)
  theta <- 2 * pi * (m - 1L) / 12
  R_bar <- sqrt(mean(cos(theta))^2 + mean(sin(theta))^2)
  sqrt(-2 * log(max(R_bar, 1e-10))) * 12 / (2 * pi)
}

# Within-country variation: circular SD across all admin-1 peak months
country_var <- seasonality_admin1[, .(
  peak_circ_sd = circ_sd_month(site_peak_month_1),
  n_admin1     = .N
), by = country]
country_var[, country_tag := ifelse(peak_circ_sd < 1.0, "country_uniform", "country_variable")]

# Per-country matching function
match_dhs_to_site <- function(sub, site_n1, site_clean_translated) {
  sub_cl  <- clean_str(sub)
  n       <- length(site_n1)
  if (n == 0L) return(list(site_name_1 = NA_character_, level = "country"))

  # 1. Exact match (on translated site names)
  idx <- which(site_clean_translated == sub_cl)
  if (length(idx) > 0L) return(list(site_name_1 = site_n1[idx[1L]], level = "admin1"))

  # 2. Contains match: site name is substring of DHS name, or DHS name is substring of site name
  contains_idx <- which(sapply(site_clean_translated, function(s)
    nchar(s) > 2L && (grepl(s, sub_cl, fixed = TRUE) || grepl(sub_cl, s, fixed = TRUE))))
  if (length(contains_idx) > 0L) {
    # Pick the one whose length is closest to the query
    best <- contains_idx[which.min(abs(nchar(site_clean_translated[contains_idx]) - nchar(sub_cl)))]
    return(list(site_name_1 = site_n1[best], level = "admin1_partial"))
  }

  # 3. Fuzzy: Levenshtein distance <= max(2, 25% of query length)
  dists     <- as.integer(adist(sub_cl, site_clean_translated))
  min_dist  <- min(dists)
  threshold <- max(2L, floor(0.25 * nchar(sub_cl)))
  if (min_dist <= threshold) {
    best <- which(dists == min_dist)[1L]
    return(list(site_name_1 = site_n1[best], level = "admin1_fuzzy"))
  }

  list(site_name_1 = NA_character_, level = "country")
}

# All DHS sub_region x country combinations
dhs_regions <- unique(base_merged[, .(country, sub_region)])
dhs_regions <- dhs_regions[sub_region != "National" & sub_region != ""]

cat("Building DHS sub_region -> site name_1 crosswalk...\n")

crosswalk <- rbindlist(lapply(unique(dhs_regions$country), function(ctry) {
  dhs_subs <- unique(dhs_regions[country == ctry, sub_region])
  site_n1  <- unique(seasonality_admin1[country == ctry, name_1])
  site_tr  <- sapply(site_n1, translate_site, USE.NAMES = FALSE)   # translated

  rbindlist(lapply(dhs_subs, function(sub) {
    res <- match_dhs_to_site(sub, site_n1, site_tr)
    data.table(country = ctry, sub_region = sub,
               site_name_1 = res$site_name_1, site_match_level = res$level)
  }))
}))

# For country-fallback rows, tag as uniform or variable
crosswalk <- merge(crosswalk, country_var[, .(country, country_tag)], by = "country", all.x = TRUE)
crosswalk[site_match_level == "country" & !is.na(country_tag),
          site_match_level := country_tag]
crosswalk[, country_tag := NULL]

# Report match quality
match_summary <- crosswalk[, .N, by = site_match_level][order(-N)]
cat("  Match level breakdown:\n")
print(match_summary)

cat("\n  Admin-1 matches by country (all tiers):\n")
print(crosswalk[grepl("^admin1", site_match_level), .N, by = .(country, site_match_level)][order(country, site_match_level)])

cat("\n  Country-variable fallbacks (high within-country variation, SD >= 1 month):\n")
var_countries <- country_var[country_tag == "country_variable", .(country, n_admin1, peak_circ_sd = round(peak_circ_sd, 2))]
unmatched_var <- merge(
  crosswalk[site_match_level == "country_variable", .(country, sub_region)],
  var_countries, by = "country")
if (nrow(unmatched_var) > 0L) print(unique(unmatched_var[, .(country, n_admin1, peak_circ_sd)])[order(-peak_circ_sd)]) else cat("  None\n")

# -----------------------------------------------------------------------------
# 3.5 Merge seasonality onto crosswalk
# -----------------------------------------------------------------------------

# Admin-1 matched rows (exact, partial, or fuzzy -- all have a site_name_1)
xw_admin1 <- crosswalk[grepl("^admin1", site_match_level)]
xw_admin1 <- merge(xw_admin1,
  seasonality_admin1[, .(country, name_1, site_seasonality_type, site_peak_month_1, site_peak_month_2,
                          site_seasonality_index)],
  by.x = c("country", "site_name_1"), by.y = c("country", "name_1"), all.x = TRUE)

# Country fallback rows (uniform or variable -- no admin-1 match)
xw_country <- crosswalk[grepl("^country", site_match_level)]
match_levels_backup <- xw_country$site_match_level   # preserve uniform/variable tag
xw_country <- merge(xw_country,
  seasonality_country[, .(country, site_seasonality_type, site_peak_month_1, site_peak_month_2,
                           site_seasonality_index)],
  by = "country", all.x = TRUE)
xw_country[, site_match_level := match_levels_backup]   # restore uniform/variable tag

site_cols <- c("site_seasonality_type", "site_peak_month_1", "site_peak_month_2", "site_match_level",
               "site_seasonality_index")
xw_full <- rbindlist(
  list(xw_admin1[, c("country", "sub_region", site_cols), with = FALSE],
       xw_country[, c("country", "sub_region", site_cols), with = FALSE]),
  fill = TRUE)

# Include a row for National sub_region (use country-level seasonality)
national_rows <- merge(
  unique(base_merged[sub_region == "National", .(country)]),
  seasonality_country[, .(country, site_seasonality_type, site_peak_month_1, site_peak_month_2,
                           site_seasonality_index)],
  by = "country", all.x = TRUE)
national_rows[, sub_region       := "National"]
national_rows[, site_match_level := "country"]
xw_full <- rbindlist(list(xw_full, national_rows[, c("country","sub_region", site_cols), with = FALSE]),
                     fill = TRUE)

# -----------------------------------------------------------------------------
# 3.6 SMC coverage by sub_region x year
#
# Unlike seasonality (a single classification per admin-1, averaged across all
# years), SMC coverage varies by year, so this keeps a year dimension and joins
# it onto each survey's own survey_year below, rather than collapsing to one
# static value per sub_region as 3.5 does.
# -----------------------------------------------------------------------------

# One row per admin-1 x year (mean handles any accidental duplicates; standard
# site-package format already has exactly one row per name_1 x year)
smc_admin1_year_raw <- smc_all[, .(smc_cov = mean(smc_cov, na.rm = TRUE)),
                               by = .(country, iso3c, name_1, year)]

# Zero-fill admin-1 x year combinations that have rainfall/seasonality data
# (i.e. are known real admin-1 units) but no smc_cov row at all. Some sources
# (e.g. Uganda's SNT file) list ONLY the admin-1 units actually enrolled in an
# SMC programme (e.g. just "Karamoja") rather than including every district
# with an explicit 0 -- unlike the standard site-package format, which does
# list non-implementing districts at 0 (see e.g. Kenya, Niger). Left as-is,
# the missing districts would be silently excluded from any country-level
# average, which would then be biased toward whichever few districts DO have
# an SMC programme. Only zero-fills years the SMC source actually covers for
# that country -- it does not assume 0 coverage for years the source has no
# data for at all (e.g. years before a country's SMC programme existed in the
# source data, which should stay NA, not 0).
admin1_years_seen <- unique(rain_all[, .(country, iso3c, name_1, year)])
smc_admin1_year <- merge(admin1_years_seen, smc_admin1_year_raw,
                          by = c("country", "iso3c", "name_1", "year"), all.x = TRUE)

smc_year_range <- smc_admin1_year_raw[, .(min_year = min(year), max_year = max(year)), by = country]
smc_admin1_year <- merge(smc_admin1_year, smc_year_range, by = "country", all.x = TRUE)
smc_admin1_year[is.na(smc_cov) & !is.na(min_year) & year >= min_year & year <= max_year,
                smc_cov := 0]
smc_admin1_year[, c("min_year", "max_year") := NULL]

# Country x year fallback (mean across that country's admin-1 units) --
# mirrors seasonality_country's role in 3.5, but keeping the year dimension
# since SMC coverage (unlike seasonality) genuinely varies by year.
smc_country_year <- smc_admin1_year[, .(smc_cov = mean(smc_cov, na.rm = TRUE)),
                                    by = .(country, iso3c, year)]

# Admin-1 matched crosswalk rows: fan out to every year that admin-1 has an
# SMC value for (allow.cartesian: a genuine one-to-many fan-out on year, not
# a runaway duplicate-key bug -- each (country, sub_region) crosswalk row has
# exactly one site_name_1, matched here against every year that name_1 has).
xw_admin1_smc <- crosswalk[grepl("^admin1", site_match_level)]
xw_admin1_smc <- merge(xw_admin1_smc,
  smc_admin1_year[, .(country, name_1, year, smc_cov)],
  by.x = c("country", "site_name_1"), by.y = c("country", "name_1"),
  all.x = TRUE, allow.cartesian = TRUE)

# Country-fallback crosswalk rows (no admin-1 match): same fan-out, using the
# country x year mean instead.
xw_country_smc <- crosswalk[grepl("^country", site_match_level)]
xw_country_smc <- merge(xw_country_smc,
  smc_country_year[, .(country, year, smc_cov)],
  by = "country", all.x = TRUE, allow.cartesian = TRUE)

smc_xw <- rbindlist(list(
  xw_admin1_smc[, .(country, sub_region, year, smc_cov)],
  xw_country_smc[, .(country, sub_region, year, smc_cov)]
), fill = TRUE)

# National sub_region rows: country x year fallback, same convention as
# 3.5's national_rows for seasonality.
national_smc <- merge(
  unique(base_merged[sub_region == "National", .(country)]),
  smc_country_year[, .(country, year, smc_cov)],
  by = "country", all.x = TRUE, allow.cartesian = TRUE)
national_smc[, sub_region := "National"]

smc_xw_full <- rbindlist(
  list(smc_xw, national_smc[, .(country, sub_region, year, smc_cov)]),
  fill = TRUE)
setnames(smc_xw_full, c("year", "smc_cov"), c("survey_year", "site_smc_cov"))

cat(sprintf("  SMC coverage: %d (country, sub_region, survey_year) rows built (%d non-NA)\n",
            nrow(smc_xw_full), smc_xw_full[!is.na(site_smc_cov), .N]))

# -----------------------------------------------------------------------------
# 3.6b Site-package modelled PfPR2-10 (site_prev_pfpr) by sub_region x year
#
# Same fan-out-by-year pattern as 3.6's SMC block (prev_all already has a year
# dimension, unlike seasonality). No zero-fill step: unlike SMC coverage,
# "no prevalence row for this admin-1/year" has no natural zero -- it just
# means the site package has no modelled estimate there, so it stays NA.
# -----------------------------------------------------------------------------

prev_admin1_year <- prev_all[, .(pfpr = mean(pfpr, na.rm = TRUE)), by = .(country, iso3c, name_1, year)]
prev_country_year <- prev_admin1_year[, .(pfpr = mean(pfpr, na.rm = TRUE)), by = .(country, iso3c, year)]

xw_admin1_prev <- crosswalk[grepl("^admin1", site_match_level)]
xw_admin1_prev <- merge(xw_admin1_prev,
  prev_admin1_year[, .(country, name_1, year, pfpr)],
  by.x = c("country", "site_name_1"), by.y = c("country", "name_1"),
  all.x = TRUE, allow.cartesian = TRUE)

xw_country_prev <- crosswalk[grepl("^country", site_match_level)]
xw_country_prev <- merge(xw_country_prev,
  prev_country_year[, .(country, year, pfpr)],
  by = "country", all.x = TRUE, allow.cartesian = TRUE)

prev_xw <- rbindlist(list(
  xw_admin1_prev[, .(country, sub_region, year, pfpr)],
  xw_country_prev[, .(country, sub_region, year, pfpr)]
), fill = TRUE)

national_prev <- merge(
  unique(base_merged[sub_region == "National", .(country)]),
  prev_country_year[, .(country, year, pfpr)],
  by = "country", all.x = TRUE, allow.cartesian = TRUE)
national_prev[, sub_region := "National"]

prev_xw_full <- rbindlist(
  list(prev_xw, national_prev[, .(country, sub_region, year, pfpr)]),
  fill = TRUE)
setnames(prev_xw_full, c("year", "pfpr"), c("survey_year", "site_prev_pfpr"))

cat(sprintf("  Modelled PfPR2-10: %d (country, sub_region, survey_year) rows built (%d non-NA)\n",
            nrow(prev_xw_full), prev_xw_full[!is.na(site_prev_pfpr), .N]))

# -----------------------------------------------------------------------------
# 3.7 Survey timing aggregation from individual data
# -----------------------------------------------------------------------------

cat("Computing survey timing aggregates...\n")
stale_ind <- intersect(names(ind), c(SITE_COLS_ALL, "site_months_to_peak"))
if (length(stale_ind) > 0L) ind[, (stale_ind) := NULL]

# Check that interview_month is present (needs 01 re-pull)
if (!"interview_month" %in% names(ind)) {
  warning(paste(
    "interview_month not in the individual dataset.",
    "Re-run 01_pull_mis_data.R (with REPULL_DATA = TRUE) to populate it.",
    "site_months_to_peak will be NA for all rows."))
  ind[, interview_month := NA_integer_]
}
if (!"rdt_type" %in% names(ind)) ind[, rdt_type := NA_character_]

# Admin-1: circular mean month and modal RDT brand
admin1_timing <- ind[, .(
  survey_month_mean = circ_mean_month(interview_month),
  rdt_brand_modal   = modal(rdt_type)
), by = .(country, survey_year, sub_region)]

# Admin-0 (national): circular mean across all sub-regions within country x survey
admin0_timing <- ind[, .(
  sub_region        = "National",
  survey_month_mean = circ_mean_month(interview_month),
  rdt_brand_modal   = modal(rdt_type)
), by = .(country, survey_year)]

timing_all <- rbindlist(list(admin1_timing, admin0_timing), fill = TRUE)

# -----------------------------------------------------------------------------
# 3.8 Proxy column: nearest_peak_month - survey_month
# -----------------------------------------------------------------------------

compute_proxy <- function(interview_month_vec, peak1_vec, peak2_vec, type_vec) {
  mapply(function(im, p1, p2, tp) {
    if (is.na(im) || is.na(p1) || tp == "no_seasonal_trend") return(NA_integer_)
    d1 <- circ_diff_month(p1, im)
    if (!is.na(p2)) {
      d2 <- circ_diff_month(p2, im)
      return(if (abs(d1) <= abs(d2)) d1 else d2)
    }
    d1
  }, interview_month_vec, peak1_vec, peak2_vec, type_vec, SIMPLIFY = TRUE)
}

# -----------------------------------------------------------------------------
# 3.9 Update aggregate CSV (dhs_pfhrp2_merged.csv)
# -----------------------------------------------------------------------------

cat(sprintf("Updating %s...\n", file.path(OUT_DIR, "dhs_pfhrp2_merged.csv")))

# Merge timing
base_merged <- merge(base_merged, timing_all, by = c("country", "survey_year", "sub_region"), all.x = TRUE)

# Fallback: where survey_month_mean is still NA, use national (country x survey_year) circular mean
nat_fallback <- admin0_timing[, .(country, survey_year, nat_month_mean = survey_month_mean)]
base_merged  <- merge(base_merged, nat_fallback, by = c("country", "survey_year"), all.x = TRUE)
base_merged[is.na(survey_month_mean), survey_month_mean := nat_month_mean]
base_merged[, nat_month_mean := NULL]
n_filled <- base_merged[sub_region != "National", sum(!is.na(survey_month_mean))]
cat(sprintf("  survey_month_mean: %d/%d non-NA sub_region rows after national fallback\n",
            n_filled, base_merged[sub_region != "National", .N]))

# Merge seasonality crosswalk
base_merged <- merge(base_merged, xw_full, by = c("country", "sub_region"), all.x = TRUE)

# Merge SMC coverage (country x sub_region x survey_year -- see 3.6)
base_merged <- merge(base_merged, smc_xw_full, by = c("country", "sub_region", "survey_year"), all.x = TRUE)

# Merge modelled PfPR2-10 (country x sub_region x survey_year -- see 3.6b)
base_merged <- merge(base_merged, prev_xw_full, by = c("country", "sub_region", "survey_year"), all.x = TRUE)

# Proxy using circular-mean survey month
base_merged[, site_months_to_peak_mean := as.integer(compute_proxy(
  survey_month_mean, site_peak_month_1, site_peak_month_2, site_seasonality_type))]

fwrite(base_merged, file.path(OUT_DIR, "dhs_pfhrp2_merged.csv"))
cat(sprintf("  Saved: %d rows, %d columns\n", nrow(base_merged), ncol(base_merged)))

# -----------------------------------------------------------------------------
# 3.10 Update individual CSV (dhs_mis_individual_hrp2.csv)
# -----------------------------------------------------------------------------

IND_PATH <- file.path(OUT_DIR, "dhs_mis_individual_hrp2.csv")
cat(sprintf("Updating %s...\n", IND_PATH))

# Merge seasonality onto individuals via country x sub_region
ind <- merge(ind, xw_full, by = c("country", "sub_region"), all.x = TRUE)

# Merge SMC coverage via country x sub_region x survey_year (see 3.6) --
# every child in the same sub_region-survey shares the same area-level value,
# not individual receipt (see this section's header).
ind <- merge(ind, smc_xw_full, by = c("country", "sub_region", "survey_year"), all.x = TRUE)

# Merge modelled PfPR2-10 via country x sub_region x survey_year (see 3.6b) --
# area-level MAP/site-package baseline, not the DHS-survey microscopy
# prevalence sr_micro_10pp already derives from below.
ind <- merge(ind, prev_xw_full, by = c("country", "sub_region", "survey_year"), all.x = TRUE)

# Merge sub-regional microscopy prevalence, RDT prevalence, and care-or-ACT
# prevalence from the aggregate table -- 03_model_rdt_individual.R's own
# sr_micro_10pp/sr_rdt_10pp/sr_care_or_act_10pp derive from these and expect
# them already present. prev_rdt is pulled through analogously to prev_micro
# (both are sub-regional aggregates from Section 2.1's base table, carried
# unchanged into base_merged).
# hrp2_risk is NOT merged here -- unlike prev_micro/prev_rdt/prop_care_or_act,
# it's already a native column of ind (merged in Section 2.8).
micro_care_lookup <- unique(base_merged[, .(country, sub_region, survey_year, prev_micro, prev_rdt, prop_care_or_act)])
ind <- merge(ind, micro_care_lookup, by = c("country", "sub_region", "survey_year"), all.x = TRUE)

# Individual-level proxy (uses per-row interview_month)
ind[, site_months_to_peak := as.integer(compute_proxy(
  interview_month, site_peak_month_1, site_peak_month_2, site_seasonality_type))]

fwrite(ind, IND_PATH)
cat(sprintf("  Saved: %d rows, %d columns\n", nrow(ind), ncol(ind)))

# -----------------------------------------------------------------------------
# 3.11 Summary diagnostics
# -----------------------------------------------------------------------------

cat("\n--- Seasonality classification summary (admin-1 units) ---\n")
print(seasonality_admin1[, .N, by = site_seasonality_type])

cat("\n--- Countries with site data ---\n")
print(sort(unique(seasonality_admin1$country)))

cat("\n--- Countries WITHOUT site data (will have NA site_* columns) ---\n")
all_dhs_countries <- unique(base_merged$country)
print(sort(setdiff(all_dhs_countries, unique(seasonality_admin1$country))))

cat("\n--- Admin-1 match rate by country (all tiers) ---\n")
print(crosswalk[, .(
  total         = .N,
  exact_match   = sum(site_match_level == "admin1"),
  partial_match = sum(site_match_level == "admin1_partial"),
  fuzzy_match   = sum(site_match_level == "admin1_fuzzy"),
  country_fb    = sum(grepl("^country", site_match_level)),
  pct_admin1    = round(100 * mean(grepl("^admin1", site_match_level)))
), by = country][order(-pct_admin1)])

cat("\n--- site_months_to_peak distribution (aggregate) ---\n")
print(summary(base_merged$site_months_to_peak_mean))

cat("\n--- site_smc_cov coverage (aggregate) ---\n")
cat(sprintf("  %d/%d sub_region rows non-NA\n",
            base_merged[sub_region != "National", sum(!is.na(site_smc_cov))],
            base_merged[sub_region != "National", .N]))
print(summary(base_merged$site_smc_cov))

cat("\n--- prev_micro/prev_rdt/prop_care_or_act coverage on the individual-level merge ---\n")
cat(sprintf("  prev_micro: %d/%d individuals non-NA after join\n",
            sum(!is.na(ind$prev_micro)), nrow(ind)))
cat(sprintf("  prev_rdt: %d/%d individuals non-NA after join\n",
            sum(!is.na(ind$prev_rdt)), nrow(ind)))
cat(sprintf("  prop_care_or_act: %d/%d individuals non-NA after join\n",
            sum(!is.na(ind$prop_care_or_act)), nrow(ind)))

cat("\n--- site_prev_pfpr (modelled PfPR2-10) coverage on the individual-level merge ---\n")
cat(sprintf("  %d/%d individuals non-NA after join\n",
            sum(!is.na(ind$site_prev_pfpr)), nrow(ind)))
print(summary(ind$site_prev_pfpr))

cat("\n--- site_seasonality_index coverage on the individual-level merge ---\n")
cat(sprintf("  %d/%d individuals non-NA after join\n",
            sum(!is.na(ind$site_seasonality_index)), nrow(ind)))
print(summary(ind$site_seasonality_index))

cat("\nDone.\n")
