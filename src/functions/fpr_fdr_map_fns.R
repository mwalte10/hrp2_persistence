# =============================================================================
# src/functions/fpr_fdr_map_fns.R
#
# ANALYSIS COMPONENT: REGRESSION ANALYSIS (function library: FPR/FDR maps)
#
# Sub-regional FPR/FDR maps from the step 2 regression (src/plots.R
# PLOT_FPR_FDR_MAPS). Ported from mal_proj src/03k_brms_country_fpr_fdr_maps.R
# and src/functions/map_helpers.R.
#
# For each country's most recent survey, one page of seven maps:
#   A  SR RDT prevalence x PfHRP2 deletion risk (all 9 levels in the legend)
#   B  SMC (has_smc, the binary the model sees)
#   C  Fever prevalence (all children)
#   D  Care seeking (formal care or ACT) among febrile children -- the DHS/MIS
#      treatment-seeking questions are only asked about children with fever
#      in the past 2 weeks, so care_or_act is ~0 among non-febrile children
#      (0-2.5% in every survey); an all-children denominator would just
#      reproduce panel C.
#   E  FPR = P(RDT+ | LM-)
#   F  FDR = P(LM- | RDT+)
#   G  P(RDT+ and LM-) -- absolute false positives (= FPR x P(LM-) = FDR x P(RDT+))
# plus one page with every country on a map of sub-Saharan Africa.
#
# E-G are posterior predictions at a fixed age, during/after peak season,
# INCLUDING the random effect of that country's most recent survey
# (re_formula = NULL). This departs from the FPR/FDR lollipop (population-
# level) on purpose: these maps describe each country, and survey
# intercepts (SD ~0.55 on each logit) move sub-regional FDR by a median
# 6 pp, and up to ~40 pp (Madagascar 2016), relative to the population-level
# prediction.
#
# fever and care_or_act are individual-level terms, so each sub-region's
# rates marginalise over its own observed (survey-weighted) joint
# distribution of covariates -- joint rather than the product of the
# marginals, because care seeking is almost entirely conditional on fever.
# sr_rdt_cat/hrp2_risk/has_smc are in that joint too rather than assumed
# constant within sub-region (panels A/B show the modal level). The mixture
# is taken on the four outcome-cell probabilities within each posterior draw
# BEFORE forming the ratios, so the rates are population rates for the
# sub-region's children of that age, not an average of per-stratum ratios.
#
# Boundaries: Natural Earth admin-1 via the survey-specific lookups below
# (SURVEY_CONFIG), except Burundi, Cote d'Ivoire, Ghana and Togo, which use
# geoBoundaries gbOpen (GEOBOUNDARIES_DIR, see paths.R).
# =============================================================================

# ---- Helpers ----
# norm_name() rather than mal_proj's norm(), which would mask base::norm.
norm_name <- function(x) iconv(tolower(trimws(x)), to = "ASCII//TRANSLIT")

apply_fixes <- function(x, fix_vec) {
  idx <- match(x, names(fix_vec))
  ifelse(!is.na(idx), fix_vec[idx], x)
}

# Join sub-region data onto boundary polygons by normalised name, optionally
# dissolving polygons into survey regions first (dissolve_lookup: normalised
# polygon name -> normalised survey sub_region).
make_map_sf <- function(geo_sf, dt_year, dissolve_lookup = NULL) {
  geo_sf$ne_norm <- norm_name(geo_sf$name)

  if (!is.null(dissolve_lookup)) {
    geo_sf$dhs_region <- dissolve_lookup[geo_sf$ne_norm]
    geo_sf <- geo_sf[!is.na(geo_sf$dhs_region), ] |>
      dplyr::group_by(dhs_region) |>
      dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
    join_key <- "dhs_region"
  } else {
    join_key <- "ne_norm"
  }

  dt_norm <- copy(dt_year)
  dt_norm[, join_key := norm_name(sub_region)]

  merged <- merge(
    as.data.frame(geo_sf), dt_norm,
    by.x = join_key, by.y = "join_key",
    all.x = TRUE
  )
  sf::st_as_sf(merged)
}

# ---- Country-specific name fixes ----
fixes_sierra_leone <- c(
  "eastern province" = "eastern", "northern province" = "northern",
  "southern province" = "southern", "western province" = "western",
  "eastern area" = "eastern", "northern area" = "northern"
)
fixes_rwanda <- c(
  "east" = "eastern", "north" = "northern",
  "south" = "southern", "west" = "western",
  "kigali" = "kigali city"
)
fixes_kenya        <- c("north eastern" = "north-eastern")
fixes_nigeria_states <- c("fct" = "federal capital territory", "nasarawa" = "nassarawa")
fixes_togo         <- c("centrale" = "centre")
fixes_mali         <- c("tombouctou" = "timbuktu")
fixes_madagascar   <- c(
  "amoron i mania"      = "amoron'i mania",
  "atsimo andrefana"    = "atsimo-andrefana",
  "atsimo atsinanana"   = "atsimo-atsinanana",
  "vatovavy fitovinany" = "vatovavy-fitovinany",
  "alaotra mangoro"     = "alaotra-mangoro"
)

# ---- Dissolve lookups: NE-normalised polygon name → DHS-normalised region ----

dissolve_bf <- c(
  "bale"="boucle du mouhoun","banwa"="boucle du mouhoun","kossi"="boucle du mouhoun",
  "mou houn"="boucle du mouhoun","nayala"="boucle du mouhoun","sourou"="boucle du mouhoun",
  "komoe"="cascades","leraba"="cascades",
  "kadiogo"="centre",
  "boulgou"="centre-est","koulpelogo"="centre-est","kouritenga"="centre-est",
  "bam"="centre-nord","namentenga"="centre-nord","sanmatenga"="centre-nord",
  "boulkiemde"="centre-ouest","sanguie"="centre-ouest","sissili"="centre-ouest","ziro"="centre-ouest",
  "bazega"="centre-sud","nahouri"="centre-sud","zoundweogo"="centre-sud",
  "gnagna"="est","gourma"="est","kompienga"="est","komondjari"="est","tapoa"="est",
  "houet"="hauts-bassins","kenedougou"="hauts-bassins","tuy"="hauts-bassins",
  "loroum"="nord","passore"="nord","yatenga"="nord","zondoma"="nord",
  "ganzourgou"="plateau central","kourweogo"="plateau central","oubritenga"="plateau central",
  "oudalan"="sahel","seno"="sahel","soum"="sahel","yagha"="sahel",
  "bougouriba"="sud-ouest","ioba"="sud-ouest","noumbiel"="sud-ouest","poni"="sud-ouest"
)
fixes_bf <- c(
  "boucle de mouhoun" = "boucle du mouhoun",
  "centre est"        = "centre-est",
  "centre nord"       = "centre-nord",
  "centre ouest"      = "centre-ouest",
  "centre sud"        = "centre-sud",
  "hauts basins"      = "hauts-bassins",
  "sud ouest"         = "sud-ouest"
)

dissolve_mw <- c(
  "chitipa"="northern","rumphi"="northern","mzimba"="northern",
  "nkhata bay"="northern","likoma"="northern",
  "kasungu"="central","nkhotakota"="central","ntchisi"="central",
  "dowa"="central","salima"="central","lilongwe"="central",
  "mchinji"="central","dedza"="central","ntcheu"="central",
  "mangochi"="southern","machinga"="southern","balaka"="southern",
  "zomba"="southern","chiradzulu"="southern","blantyre"="southern",
  "mwanza"="southern","neno"="southern","thyolo"="southern",
  "mulanje"="southern","phalombe"="southern","chikwawa"="southern",
  "nsanje"="southern"
)
fixes_mw <- c("north" = "northern", "south" = "southern")

dissolve_gn <- c(
  "boke"="boke","boffa"="boke","fria"="boke","gaoual"="boke","koundara"="boke",
  "conakry"="conakry",
  "faranah"="faranah","dabola"="faranah","dinguiraye"="faranah","kissidougou"="faranah",
  "kankan"="kankan","kerouane"="kankan","kouroussa"="kankan","mandiana"="kankan","siguiri"="kankan",
  "kindia"="kindia","coyah"="kindia","dubreka"="kindia","forecariah"="kindia","telimele"="kindia",
  "labe"="labe","koubia"="labe","lelouma"="labe","mali"="labe","pita"="labe","tougue"="labe",
  "mamou"="mamou","dalaba"="mamou",
  "nzerekore"="nzerekore","beyla"="nzerekore","gueckedou"="nzerekore","lola"="nzerekore",
  "macenta"="nzerekore","yomou"="nzerekore"
)

dissolve_lr <- c(
  "montserrado"="monrovia",
  "bong"="north central","nimba"="north central","lofa"="north central",
  "bomi"="north western","grand cape mount"="north western","gbapolu"="north western",
  "margibi"="south central","grand bassa"="south central","river cess"="south central",
  "sinoe"="south eastern a","grand gedeh"="south eastern a","river gee"="south eastern a",
  "maryland"="south eastern b","grand kru"="south eastern b"
)
fixes_lr <- c("greater monrovia" = "monrovia")

dissolve_ng_zones <- c(
  "benue"="north central","federal capital territory"="north central",
  "kogi"="north central","kwara"="north central","nassarawa"="north central",
  "niger"="north central","plateau"="north central",
  "adamawa"="north east","bauchi"="north east","borno"="north east",
  "gombe"="north east","taraba"="north east","yobe"="north east",
  "jigawa"="north west","kaduna"="north west","kano"="north west",
  "katsina"="north west","kebbi"="north west","sokoto"="north west","zamfara"="north west",
  "abia"="south east","anambra"="south east","ebonyi"="south east",
  "enugu"="south east","imo"="south east",
  "akwa ibom"="south south","bayelsa"="south south","cross river"="south south",
  "delta"="south south","edo"="south south","rivers"="south south",
  "ekiti"="south west","lagos"="south west","ogun"="south west",
  "ondo"="south west","osun"="south west","oyo"="south west"
)

dissolve_ug_2009 <- c(
  "kiboga"="central 1","nakaseke"="central 1","luweero"="central 1",
  "mubende"="central 1","mityana"="central 1","gomba"="central 1","mpigi"="central 1",
  "nakasongola"="central 1","kyankwanzi"="central 1","butambala"="central 1",
  "bukomansimbi"="central 1",
  "masaka"="central 2","rakai"="central 2","kalangala"="central 2",
  "mukono"="central 2","buikwe"="central 2","kalungu"="central 2",
  "lwengo"="central 2","lyantonde"="central 2","sembabule"="central 2",
  "kayunga"="central 2","buvuma"="central 2",
  "jinja"="east central","iganga"="east central","kamuli"="east central",
  "buyende"="east central","luuka"="east central","kaliro"="east central",
  "mayuge"="east central","namutumba"="east central","bugiri"="east central",
  "namayingo"="east central",
  "kampala"="kampala","wakiso"="kampala",
  "mbale"="mid eastern","sironko"="mid eastern","bududa"="mid eastern",
  "manafwa"="mid eastern","bulambuli"="mid eastern","kapchorwa"="mid eastern",
  "bukwa"="mid eastern","kween"="mid eastern",
  "tororo"="mid eastern","busia"="mid eastern","butaleja"="mid eastern",
  "pallisa"="mid eastern","kibuku"="mid eastern","budaka"="mid eastern",
  "soroti"="mid eastern","kumi"="mid eastern","ngora"="mid eastern",
  "katakwi"="mid eastern","amuria"="mid eastern","serere"="mid eastern",
  "kaberamaido"="mid eastern","bukedea"="mid eastern",
  "gulu"="mid northern","kitgum"="mid northern","lamwo"="mid northern",
  "amuru"="mid northern","nwoya"="mid northern","agago"="mid northern",
  "pader"="mid northern","lira"="mid northern","apac"="mid northern",
  "dokolo"="mid northern","amolatar"="mid northern","kole"="mid northern",
  "otuke"="mid northern","oyam"="mid northern","alebtong"="mid northern",
  "hoima"="mid western","masindi"="mid western","kibale"="mid western",
  "kiryandongo"="mid western","buliisa"="mid western",
  "kabarole"="mid western","kamwenge"="mid western","kyenjojo"="mid western",
  "ntoroko"="mid western","kyegegwa"="mid western",
  "kasese"="mid western","bundibugyo"="mid western",
  "moroto"="north east","kotido"="north east","kaabong"="north east",
  "abim"="north east","nakapiripirit"="north east","napak"="north east",
  "amudat"="north east",
  "kabale"="south","kisoro"="south","kanungu"="south","rukungiri"="south",
  "ibanda"="south","isingiro"="south","mbarara"="south","bushenyi"="south",
  "ntungamo"="south","mitooma"="south","buhweju"="south","sheema"="south",
  "rubirizi"="south",
  "arua"="west nile","adjumani"="west nile","moyo"="west nile",
  "nebbi"="west nile","zombo"="west nile","koboko"="west nile",
  "yumbe"="west nile","maracha"="west nile"
)

dissolve_ug_2014 <- dissolve_ug_2009
dissolve_ug_2014[dissolve_ug_2014 == "mid northern"] <- "mid north"
dissolve_ug_2014[dissolve_ug_2014 == "south"]        <- "south western"

dissolve_ug_2018 <- c(
  "gulu"="acholi","kitgum"="acholi","lamwo"="acholi","amuru"="acholi",
  "nwoya"="acholi","agago"="acholi","pader"="acholi","omoro"="acholi",
  "mbarara"="ankole","bushenyi"="ankole","ntungamo"="ankole","isingiro"="ankole",
  "ibanda"="ankole","mitooma"="ankole","buhweju"="ankole","sheema"="ankole",
  "rubirizi"="ankole",
  "mbale"="bugisu","sironko"="bugisu","bududa"="bugisu","manafwa"="bugisu",
  "bulambuli"="bugisu",
  "tororo"="bukedi","busia"="bukedi","butaleja"="bukedi",
  "pallisa"="bukedi","kibuku"="bukedi","budaka"="bukedi",
  "hoima"="bunyoro","masindi"="bunyoro","kibale"="bunyoro",
  "kiryandongo"="bunyoro","buliisa"="bunyoro",
  "jinja"="busoga","iganga"="busoga","kamuli"="busoga","buyende"="busoga",
  "luuka"="busoga","kaliro"="busoga","mayuge"="busoga","namutumba"="busoga",
  "bugiri"="busoga","namayingo"="busoga",
  "kampala"="kampala","wakiso"="kampala",
  "moroto"="karamoja","kotido"="karamoja","kaabong"="karamoja",
  "abim"="karamoja","nakapiripirit"="karamoja","napak"="karamoja","amudat"="karamoja",
  "kabale"="kigezi","kisoro"="kigezi","kanungu"="kigezi","rukungiri"="kigezi",
  "lira"="lango","apac"="lango","dokolo"="lango","amolatar"="lango",
  "kole"="lango","otuke"="lango","oyam"="lango","alebtong"="lango",
  "luweero"="north buganda","nakaseke"="north buganda","nakasongola"="north buganda",
  "kiboga"="north buganda","kyankwanzi"="north buganda","kayunga"="north buganda",
  "mubende"="north buganda","mityana"="north buganda","gomba"="north buganda",
  "masaka"="south buganda","rakai"="south buganda","kalangala"="south buganda",
  "mukono"="south buganda","buikwe"="south buganda","kalungu"="south buganda",
  "lwengo"="south buganda","lyantonde"="south buganda","sembabule"="south buganda",
  "buvuma"="south buganda","butambala"="south buganda","bukomansimbi"="south buganda",
  "mpigi"="south buganda",
  "soroti"="teso","kumi"="teso","ngora"="teso","katakwi"="teso",
  "amuria"="teso","serere"="teso","kaberamaido"="teso","bukedea"="teso",
  "kabarole"="tooro","kamwenge"="tooro","kyenjojo"="tooro",
  "ntoroko"="tooro","kyegegwa"="tooro","kasese"="tooro","bundibugyo"="tooro",
  "arua"="west nile","adjumani"="west nile","moyo"="west nile",
  "nebbi"="west nile","zombo"="west nile","koboko"="west nile",
  "yumbe"="west nile","maracha"="west nile"
)

dissolve_ug_2024 <- dissolve_ug_2018
dissolve_ug_2024[dissolve_ug_2024 == "tooro"] <- "toro"
dissolve_ug_2024[c("mbale","sironko","bududa","manafwa","bulambuli",
                    "kapchorwa","bukwa","kween")] <- "elgon"

# ---- Survey-level config ----
SURVEY_CONFIG <- list(
  "Angola_2011"        = list(skip = TRUE),
  "Burundi_2012"       = list(skip = TRUE),
  "Ghana_2016"         = list(ne = "Ghana"),
  "Ghana_2019"         = list(ne = "Ghana"),
  "Guinea_2021"        = list(ne = "Guinea",       dissolve = dissolve_gn),
  "Kenya_2015"         = list(ne = "Kenya",         fixes   = fixes_kenya),
  "Kenya_2020"         = list(ne = "Kenya",         fixes   = fixes_kenya),
  "Liberia_2009"       = list(ne = "Liberia",       dissolve = dissolve_lr, fixes = fixes_lr),
  "Liberia_2011"       = list(ne = "Liberia",       dissolve = dissolve_lr, fixes = fixes_lr),
  "Liberia_2022"       = list(ne = "Liberia",       dissolve = dissolve_lr, fixes = fixes_lr),
  "Madagascar_2011"    = list(ne = "Madagascar",    fixes   = fixes_madagascar),
  "Madagascar_2013"    = list(skip = TRUE),
  "Madagascar_2016"    = list(ne = "Madagascar",    fixes   = fixes_madagascar),
  "Malawi_2012"        = list(ne = "Malawi",        dissolve = dissolve_mw, fixes = fixes_mw),
  "Malawi_2014"        = list(ne = "Malawi",        dissolve = dissolve_mw, fixes = fixes_mw),
  "Malawi_2017"        = list(ne = "Malawi",        dissolve = dissolve_mw, fixes = fixes_mw),
  "Mali_2015"          = list(ne = "Mali",          fixes   = fixes_mali),
  "Nigeria_2010"       = list(ne = "Nigeria",       dissolve = dissolve_ng_zones),
  "Nigeria_2015"       = list(ne = "Nigeria",       dissolve = dissolve_ng_zones),
  "Nigeria_2021"       = list(ne = "Nigeria",       fixes   = fixes_nigeria_states),
  "Rwanda_2017"        = list(ne = "Rwanda",        fixes   = fixes_rwanda),
  "Rwanda_2023"        = list(ne = "Rwanda",        fixes   = fixes_rwanda),
  "Sierra Leone_2016"  = list(ne = "Sierra Leone",  fixes   = fixes_sierra_leone),
  "Togo_2017"          = list(ne = "Togo",          fixes   = fixes_togo),
  "Uganda_2009"        = list(ne = "Uganda",        dissolve = dissolve_ug_2009),
  "Uganda_2014"        = list(ne = "Uganda",        dissolve = dissolve_ug_2014),
  "Uganda_2018"        = list(ne = "Uganda",        dissolve = dissolve_ug_2018),
  "Uganda_2024"        = list(ne = "Uganda",        dissolve = dissolve_ug_2024),
  "Burkina Faso_2014"  = list(ne = "Burkina Faso",  dissolve = dissolve_bf, fixes = fixes_bf),
  "Burkina Faso_2017"  = list(ne = "Burkina Faso",  dissolve = dissolve_bf, fixes = fixes_bf)
)


# ---- Boundary source per country ----
NE_NAME <- c("Cote d'Ivoire" = "Ivory Coast")
# Spelling-only differences between survey sub_region and boundary names
# (applied after any configured fixes).
EXTRA_FIXES <- list(
  "Benin"   = c("atacora" = "atakora", "atlantic" = "atlantique", "couffo" = "kouffo"),
  "Burundi" = c("karusi" = "karuzi"),
  "Guinea"  = c("n'zerekore" = "nzerekore")
)

# Natural Earth's admin-1 doesn't match these surveys' sub-regions, so they
# use geoBoundaries gbOpen, cached in GEOBOUNDARIES_DIR (downloaded on first
# use if missing). name_fn cleans shapeName to the survey's spelling.
GEOBOUNDARIES <- list(
  # 14 districts (post-2011 scheme); NE has the older 19 regions
  "Cote d'Ivoire" = list(iso = "CIV", licence = "CC BY 4.0", name_fn = function(x) {
    x <- sub("^District Autonome De? ?D?'?", "", x)
    x <- sub("^Valle Du", "Vallee Du", x)
    sub("^Bas-Sassandra$", "Bas Sassandra", x)
  }),
  # 16 regions (2018 reform); NE has the older 10
  "Ghana" = list(iso = "GHA", licence = "ODbL/CC BY-SA (OpenStreetMap)",
                 name_fn = function(x) sub(" Region$", "", x)),
  # NE predates Rumonge province (split from Bururi, 2015)
  "Burundi" = list(iso = "BDI", licence = "CC0", name_fn = identity),
  # The survey reports Lome separately from Maritime but no ADM1 source does,
  # so Togo is built from ADM2 prefectures: Lome Commune + Golfe form the
  # "Agglomeration de Lome" (DHS/MIS Grand Lome definition), every other
  # prefecture is dissolved into its ADM1 region (point-on-surface join).
  "Togo" = list(iso = "TGO", licence = "ODbL/CC BY-SA (OpenStreetMap)", build_fn = function() {
    adm1 <- gb_read("TGO", "ADM1")
    adm2 <- gb_read("TGO", "ADM2")
    # Planar operations on lon/lat are fine at prefecture scale; silence sf's notes
    s2_was <- sf::sf_use_s2()
    on.exit(suppressMessages(sf::sf_use_s2(s2_was)))
    suppressMessages(suppressWarnings({
      sf::sf_use_s2(FALSE)
      region <- adm1$shapeName[unlist(sf::st_intersects(sf::st_point_on_surface(adm2), adm1))]
      adm2$name <- ifelse(adm2$shapeName %in% c("Lome Commune", "Golfe"),
                          "Agglomeration De Lome", sub(" Region$", "", region))
      adm2 |> dplyr::group_by(name) |>
        dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
    }))
  })
)
gb_read <- function(iso, adm = "ADM1") {
  path <- file.path(GEOBOUNDARIES_DIR, sprintf("geoBoundaries-%s-%s_simplified.geojson", iso, adm))
  if (!file.exists(path)) {
    dir.create(GEOBOUNDARIES_DIR, showWarnings = FALSE, recursive = TRUE)
    meta <- jsonlite::fromJSON(sprintf("https://www.geoboundaries.org/api/current/gbOpen/%s/%s/", iso, adm))
    utils::download.file(meta$simplifiedGeometryGeoJSON, path, mode = "wb", quiet = TRUE)
  }
  sf::st_transform(sf::st_read(path, quiet = TRUE), 4326)
}
read_geoboundaries <- function(gb) {
  if (!is.null(gb$build_fn)) return(gb$build_fn())
  geo <- gb_read(gb$iso)
  geo$name <- gb$name_fn(geo$shapeName)
  geo
}

# SURVEY_CONFIG was hand-built for specific survey waves. If a country's most
# recent survey has no entry, reuse its latest configured wave's lookup
# (admin-1 schemes are usually stable between waves), else match names
# directly against Natural Earth.
resolve_map_config <- function(ctry, yr) {
  if (ctry %in% names(GEOBOUNDARIES)) return(list(ne = NA_character_, source = "geoBoundaries"))
  cfg <- SURVEY_CONFIG[[paste0(ctry, "_", yr)]]
  if (!is.null(cfg) && !isTRUE(cfg$skip)) return(c(cfg, source = "exact"))
  keys <- grep(paste0("^", ctry, "_"), names(SURVEY_CONFIG), value = TRUE)
  keys <- keys[!vapply(SURVEY_CONFIG[keys], function(x) isTRUE(x$skip), logical(1))]
  if (length(keys) > 0) {
    k <- keys[which.max(as.integer(sub(".*_", "", keys)))]
    return(c(SURVEY_CONFIG[[k]], source = k))
  }
  list(ne = if (ctry %in% names(NE_NAME)) NE_NAME[[ctry]] else ctry, source = "direct")
}

# One country's sub-region data joined onto its boundaries. Returns the sf,
# the unmatched sub-regions and a caption note on the boundary source.
country_map_sf <- function(ctry, yr, map_dt) {
  cfg <- resolve_map_config(ctry, yr)
  geo <- if (cfg$source == "geoBoundaries") {
    read_geoboundaries(GEOBOUNDARIES[[ctry]])
  } else {
    rnaturalearth::ne_states(country = cfg$ne, returnclass = "sf")
  }

  dt_c <- copy(map_dt[country == ctry])
  dt_c[, sub_region_label := sub_region]
  dt_c[, sub_region := norm_name(sub_region)]
  if (!is.null(cfg$fixes)) dt_c[, sub_region := apply_fixes(sub_region, cfg$fixes)]
  if (!is.null(EXTRA_FIXES[[ctry]])) dt_c[, sub_region := apply_fixes(sub_region, EXTRA_FIXES[[ctry]])]

  geo_c <- sf::st_transform(make_map_sf(geo, dt_c, dissolve_lookup = cfg$dissolve), 4326)
  geo_c$country <- ctry
  unmatched <- dt_c[!sub_region %in% geo_c$sub_region, sub_region_label]

  note <- ""
  if (cfg$source != "exact") {
    gb <- GEOBOUNDARIES[[ctry]]
    note <- paste0(note, "\nBoundaries: ",
                   switch(cfg$source,
                          direct        = "direct name match to Natural Earth admin-1",
                          geoBoundaries = paste0("geoBoundaries gbOpen (", gb$licence, ")",
                                                 if (!is.null(gb$build_fn))
                                                   "; ADM2 prefectures dissolved to survey regions, Lome = Lome Commune + Golfe"),
                          paste0("lookup built for ", sub("_", " ", cfg$source))), ".")
  }
  if (length(unmatched) > 0) {
    note <- paste0(note, "\nNot mapped (no matching boundary): ",
                   paste(strwrap(paste(unmatched, collapse = ", "), 120), collapse = "\n"), ".")
  }
  list(sf = geo_c, unmatched = unmatched, note = note, source = cfg$source)
}

# Sub-Saharan Africa outline (World Bank grouping), cropped to mainland +
# Madagascar: far-flung islands (e.g. South Africa's Prince Edward Islands
# at ~47S, Mauritius/Seychelles) otherwise stretch the frame.
ssa_background <- function() {
  bg <- rnaturalearth::ne_countries(continent = "Africa", scale = "medium", returnclass = "sf")
  s2_was <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(s2_was)))
  suppressMessages(suppressWarnings({
    sf::sf_use_s2(FALSE)
    sf::st_crop(bg[bg$region_wb == "Sub-Saharan Africa", ],
                c(xmin = -18, ymin = -35.5, xmax = 52, ymax = 28))
  }))
}

# ---- Predictions ----
# Children in each country's most recent survey, prepared as for the fitted
# model so sr_rdt_cat/hrp2_risk/has_smc are coded identically. age_mean_months
# (attribute) is the mean prepare_brms_full_data() centres on, taken before
# its treatment-completeness filter.
latest_survey_children <- function(mis_dir = DIR_STRUCTURE) {
  ind <- prepare_brms_full_data(mis_dir)
  ind_raw <- fread(file.path(mis_dir, "dhs_mis_individual_hrp2.csv"))
  include <- readRDS(file.path(mis_dir, "countries_included.RDS"))
  ind_raw <- ind_raw[country %in% include$country]
  ind_raw <- ind_raw[!is.na(rdt_result) & !is.na(micro_result) & !is.na(fever)]

  latest <- ind[, .(survey_year = max(survey_year)), by = country]
  dat <- merge(ind, latest, by = c("country", "survey_year"))
  dat[is.na(sample_weight), sample_weight := 1]
  setattr(dat, "age_mean_months", mean(ind_raw$age_months))
  dat[]
}

# Sub-region covariates and posterior FPR/FDR/P(RDT+ and LM-) for every
# country's most recent survey, at age `age_months`, during/after peak
# season, with that survey's random effect. Returns one row per sub-region
# (mean and 95% CrI of each rate, in %).
fpr_fdr_subregion_rates <- function(fit, age_months, mis_dir = DIR_STRUCTURE) {
  dat <- latest_survey_children(mis_dir)
  age_mean_months <- attr(dat, "age_mean_months")
  dat[, sr_key := paste(country, sub_region, sep = "__")]

  mode_of <- function(x, w) {
    tab <- tapply(w, as.character(x), sum)
    names(tab)[which.max(tab)]
  }
  sr_dt <- dat[, .(
    country     = country[1],
    survey_year = survey_year[1],
    sub_region  = sub_region[1],
    sr_rdt_cat  = mode_of(sr_rdt_cat, sample_weight),
    hrp2_risk   = mode_of(hrp2_risk,  sample_weight),
    has_smc     = as.integer(mode_of(has_smc, sample_weight)),
    mixed_cov   = uniqueN(sr_rdt_cat) > 1 | uniqueN(hrp2_risk) > 1 | uniqueN(has_smc) > 1,
    prev_fever  = weighted.mean(fever, sample_weight),
    prop_care   = weighted.mean(care_or_act[fever == 1], sample_weight[fever == 1])
  ), by = sr_key]

  cov_cols <- c("sr_rdt_cat", "hrp2_risk", "has_smc", "fever", "care_or_act")
  grid_dt <- dat[, .(w = sum(sample_weight)), by = c("sr_key", "survid", cov_cols)]
  grid_dt[, w := w / sum(w), by = sr_key]
  grid_dt[, sr_rdt_cat := factor(sr_rdt_cat, levels = levels(fit$data$sr_rdt_cat))]
  grid_dt[, hrp2_risk  := factor(hrp2_risk,  levels = levels(fit$data$hrp2_risk))]
  grid_dt[, during_or_after_peak := 1L]
  grid_dt[, age_c := age_months - age_mean_months]
  stopifnot(all(grid_dt$survid %in% fit$data$survid))

  grid_epred <- brms::posterior_epred(fit, newdata = grid_dt, re_formula = NULL)

  # Per draw: mix the four cell probabilities over each sub-region's joint
  # covariate weights (W: grid rows x sub-regions), then form the rates.
  W <- matrix(0, nrow(grid_dt), nrow(sr_dt), dimnames = list(NULL, sr_dt$sr_key))
  W[cbind(seq_len(nrow(grid_dt)), match(grid_dt$sr_key, sr_dt$sr_key))] <- grid_dt$w
  cell_mix <- lapply(setNames(nm = c("p_p", "n_p", "p_n", "n_n")),
                     function(k) grid_epred[, , k] %*% W)
  rate_draws <- list(
    fpr    = cell_mix$n_p / (cell_mix$n_p + cell_mix$n_n) * 100,
    fdr    = cell_mix$n_p / (cell_mix$n_p + cell_mix$p_p) * 100,
    abs_fp = cell_mix$n_p * 100
  )
  rates_dt <- Reduce(function(a, b) merge(a, b, by = "sr_key"), lapply(names(rate_draws), function(v) {
    m <- rate_draws[[v]]
    out <- data.table(sr_key = colnames(m), colMeans(m),
                      apply(m, 2, quantile, 0.025), apply(m, 2, quantile, 0.975))
    setnames(out, c("sr_key", v, paste0(v, "_lo"), paste0(v, "_hi")))
  }))

  map_dt <- merge(sr_dt, rates_dt, by = "sr_key")
  sr_levels   <- c("<10%", "10-20%", "≥20%")
  hrp2_levels <- c("Slight", "Medium", "High")
  map_dt[, prev_hrp2 := factor(
    paste0(gsub(">=", "≥", sr_rdt_cat), ", ", gsub("Marginal/Moderate", "Medium", hrp2_risk)),
    levels = as.vector(t(outer(sr_levels, hrp2_levels, paste, sep = ", ")))
  )]
  map_dt[, smc_lab := factor(ifelse(has_smc == 1, "Present", "Absent"), levels = c("Absent", "Present"))]
  map_dt[]
}

# ---- Contribution of each fixed effect to a country's FPR/FDR ----
# Shapley decomposition of each country's rate (most recent survey, same
# children, weights, age and season as the maps) into the reference rate and
# one contribution per factor. The rates are non-linear in the coefficients,
# so a factor's contribution depends on which others are already set; the
# Shapley value averages it over every order of adding the factors, which
# makes the contributions sum exactly to rate - reference rate and splits
# the prevalence x deletion-risk interaction between those two factors.
# Reference: <10% prevalence, slight deletion risk, no SMC, fever or care,
# before peak season, mean age, no survey effect. Contributions can be
# negative. Posterior means only: v(S) is averaged over draws before the
# decomposition (the Shapley value is linear in v, so this is the posterior
# mean contribution).
SHAPLEY_FACTORS <- c(sr_rdt_cat = "Sub-regional prevalence", hrp2_risk = "Deletion risk",
                     has_smc = "SMC", fever = "Fever", care_or_act = "Formal care or ACT",
                     during_or_after_peak = "During/after peak season", age_c = "Age (4 years vs mean)",
                     survey = "Survey effect")

fpr_fdr_shapley_country <- function(fit, age_months, mis_dir = DIR_STRUCTURE) {
  dat <- latest_survey_children(mis_dir)
  cov_cols <- c("sr_rdt_cat", "hrp2_risk", "has_smc", "fever", "care_or_act")
  grid_dt <- dat[, .(w = sum(sample_weight)), by = c("country", "survey_year", "survid", cov_cols)]
  grid_dt[, w := w / sum(w), by = country]
  grid_dt[, sr_rdt_cat := factor(sr_rdt_cat, levels = levels(fit$data$sr_rdt_cat))]
  grid_dt[, hrp2_risk  := factor(hrp2_risk,  levels = levels(fit$data$hrp2_risk))]
  grid_dt[, during_or_after_peak := 1L]
  grid_dt[, age_c := age_months - attr(dat, "age_mean_months")]
  stopifnot(all(grid_dt$survid %in% fit$data$survid))

  ref <- list(sr_rdt_cat = factor(levels(fit$data$sr_rdt_cat)[1], levels = levels(fit$data$sr_rdt_cat)),
              hrp2_risk  = factor(levels(fit$data$hrp2_risk)[1],  levels = levels(fit$data$hrp2_risk)),
              has_smc = 0L, fever = 0L, care_or_act = 0L, during_or_after_peak = 0L, age_c = 0)
  fac <- names(SHAPLEY_FACTORS)
  n_fac <- length(fac)
  surveys <- unique(grid_dt[, .(country, survey_year)])
  W <- matrix(0, nrow(grid_dt), nrow(surveys))
  W[cbind(seq_len(nrow(grid_dt)), match(grid_dt$country, surveys$country))] <- grid_dt$w

  # v[S, country]: posterior mean rate with the factors in S at their observed
  # values and the rest at the reference
  n_sets <- 2^n_fac
  v <- list(fpr = matrix(NA_real_, n_sets, nrow(surveys)), fdr = matrix(NA_real_, n_sets, nrow(surveys)))
  for (s in seq_len(n_sets) - 1) {
    in_s <- bitwAnd(s, 2^(seq_len(n_fac) - 1)) > 0
    nd <- copy(grid_dt)
    for (f in setdiff(fac[!in_s], "survey")) set(nd, j = f, value = rep(ref[[f]], nrow(nd)))
    ep <- brms::posterior_epred(fit, newdata = nd, re_formula = if (in_s[n_fac]) NULL else NA)
    mix <- lapply(setNames(nm = c("p_p", "n_p", "n_n")), function(k) ep[, , k] %*% W)
    v$fpr[s + 1, ] <- colMeans(mix$n_p / (mix$n_p + mix$n_n)) * 100
    v$fdr[s + 1, ] <- colMeans(mix$n_p / (mix$n_p + mix$p_p)) * 100
  }

  set_size <- vapply(seq_len(n_sets) - 1, function(s) sum(bitwAnd(s, 2^(seq_len(n_fac) - 1)) > 0), numeric(1))
  out <- rbindlist(lapply(names(v), function(o) {
    phi <- vapply(seq_len(n_fac), function(i) {
      bit <- 2^(i - 1)
      without <- which(bitwAnd(seq_len(n_sets) - 1, bit) == 0)
      wt <- factorial(set_size[without]) * factorial(n_fac - set_size[without] - 1) / factorial(n_fac)
      colSums(wt * (v[[o]][without + bit, , drop = FALSE] - v[[o]][without, , drop = FALSE]))
    }, numeric(nrow(surveys)))
    phi <- matrix(phi, nrow(surveys))
    stopifnot(max(abs(rowSums(phi) - (v[[o]][n_sets, ] - v[[o]][1, ]))) < 1e-8)
    rbind(
      data.table(surveys, outcome = toupper(o), component = "Reference", value = v[[o]][1, ]),
      data.table(surveys[rep(seq_len(nrow(surveys)), n_fac)], outcome = toupper(o),
                 component = rep(unname(SHAPLEY_FACTORS), each = nrow(surveys)), value = as.vector(phi)),
      data.table(surveys, outcome = toupper(o), component = "Total", value = v[[o]][n_sets, ])
    )
  }))
  out[]
}

# ---- Observed RDT:LM ratio ----
# Sub-region ratio of survey-weighted RDT to LM prevalence among children
# tested by both (step 1a's PR-based concordance table), for each included
# country's most recent survey with both tests, joined onto the same
# boundaries as the FPR/FDR maps. Ratio is NA where LM prevalence is 0.
rdt_lm_ratio_sf <- function(pull_dir = DIR_PULL, mis_dir = DIR_STRUCTURE) {
  conc <- fread(file.path(pull_dir, "dhs_mis_concordance.csv"))
  include <- readRDS(file.path(mis_dir, "countries_included.RDS"))
  conc <- conc[country %in% include$country & tolower(sub_region) != "national"]
  conc[, `:=`(prev_rdt = pct_rdt_pos_micro_pos + pct_rdt_pos_micro_neg,
              prev_lm  = pct_rdt_pos_micro_pos + pct_rdt_neg_micro_pos)]
  latest <- conc[, .(survey_year = max(survey_year)), by = country][order(country)]
  dt <- merge(conc, latest, by = c("country", "survey_year"))
  dt[, ratio := fifelse(prev_lm > 0, prev_rdt / prev_lm, NA_real_)]
  geo <- lapply(seq_len(nrow(latest)), function(i) {
    m <- country_map_sf(latest$country[i], latest$survey_year[i],
                        dt[country == latest$country[i],
                           .(country, survey_year, sub_region, prev_rdt, prev_lm, ratio, n_both_tested)])
    if (length(m$unmatched) > 0)
      message(sprintf("%s %d: no boundary for %s", latest$country[i], latest$survey_year[i],
                      paste(m$unmatched, collapse = ", ")))
    m$sf[, c("country", "survey_year", "sub_region", "prev_rdt", "prev_lm", "ratio", "n_both_tested")]
  })
  list(sf = do.call(rbind, geo), surveys = latest, n_lm_zero = dt[prev_lm == 0, .N])
}

# ---- Page ----
# Prevalence x deletion: green = Slight, yellow = Medium, orange/red = High
# risk; darker = higher SR prevalence.
COLS_PREV_HRP2 <- c(
  "<10%, Slight"       = "#c7e9c0", "<10%, Medium"       = "#fff3a0", "<10%, High"       = "#fdae6b",
  "10-20%, Slight"     = "#74c476", "10-20%, Medium"     = "#f0d43a", "10-20%, High"     = "#e6550d",
  "≥20%, Slight"  = "#238b45", "≥20%, Medium"  = "#b59a00", "≥20%, High"  = "#8c2d04"
)

MAP_THEME <- theme_void(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(size = 10, face = "bold", margin = margin(b = 3)),
    plot.margin      = margin(6, 6, 6, 6),
    legend.position  = "bottom",
    legend.title.position = "top",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 7),
    legend.key.width  = unit(0.4, "cm"),
    legend.key.height = unit(0.3, "cm")
  )

# Pad a map's extent to a fixed display aspect (height:width, after
# coord_sf's 1/cos(latitude) correction) so every country gets the same
# panel shape -- otherwise patchwork shrinks the whole page around a tall,
# narrow country like Togo.
fixed_frame <- function(geo_sf, aspect = 1) {
  bb    <- sf::st_bbox(geo_sf)
  cos_y <- cos(mean(bb[c("ymin", "ymax")]) * pi / 180)
  dx <- bb[["xmax"]] - bb[["xmin"]]
  dy <- bb[["ymax"]] - bb[["ymin"]]
  if (dy / (dx * cos_y) > aspect) dx <- dy / (aspect * cos_y) else dy <- aspect * dx * cos_y
  cx <- mean(bb[c("xmin", "xmax")]); cy <- mean(bb[c("ymin", "ymax")])
  coord_sf(xlim = cx + c(-1, 1) * dx * 0.52, ylim = cy + c(-1, 1) * dy * 0.52, expand = FALSE)
}

# show.legend = TRUE on panels A/B forces keys for levels absent from that
# country (ggplot2 >= 3.5 otherwise only draws keys present in the layer).
# bg_sf: optional background (country outlines) drawn under and over the data.
# Every polygon is first drawn in MAP_NO_DATA: the discrete scales (A/B) drop
# NA polygons, so without this layer unsurveyed regions (e.g. northern Mali)
# vanished from A/B but showed as NA fill in C-G.
MAP_NO_DATA <- "#E6DFCF"
map_panel <- function(show_all_keys = NA, bg_sf = NULL, border_lwd = 0.3) {
  list(
    if (!is.null(bg_sf)) geom_sf(data = bg_sf, inherit.aes = FALSE,
                                 fill = MAP_NO_DATA, colour = "grey75", linewidth = 0.2),
    geom_sf(data = function(d) d, inherit.aes = FALSE, fill = MAP_NO_DATA,
            colour = if (border_lwd > 0) "white" else NA, linewidth = border_lwd),
    geom_sf(colour = if (border_lwd > 0) "white" else NA, linewidth = border_lwd,
            show.legend = show_all_keys),
    if (!is.null(bg_sf)) geom_sf(data = bg_sf, inherit.aes = FALSE,
                                 fill = NA, colour = "grey40", linewidth = 0.2),
    facet_wrap(~ strip),
    MAP_THEME
  )
}
pct_0_100 <- list(limits = c(0, 100), breaks = seq(0, 100, 25),
                  guide = guide_colourbar(theme = theme(legend.key.width = unit(5, "cm"))),
                  labels = scales::label_percent(accuracy = 1, scale = 1),
                  na.value = MAP_NO_DATA)
blues_0_100 <- function(name) {
  do.call(scale_fill_distiller, c(list(palette = "Blues", direction = 1, name = name), pct_0_100))
}
# E-G use a pink-magenta ramp rather than magma: magma's yellow-orange-red range repeated
# the panel A deletion-risk colours
COLS_RATE <- c("#FDE0EF", "#F1B6DA", "#DE77AE", "#C51B7D", "#7A0177")
rate_0_100 <- function(name) {
  do.call(scale_fill_gradientn, c(list(colours = COLS_RATE, name = name), pct_0_100))
}
# SMC No/Yes take the two ends of the C/D Blues bar
COLS_SMC <- c("Absent" = "#F7FBFF", "Present" = "#08306B")

# One 2 x 4 page (A-G + empty slot). age_label goes into the caption.
# publication = TRUE (SSA figure): no title, subtitle or caption on the image
# (the caption is returned as attr(p, "caption") for the manuscript) and
# short markdown panel titles.
build_fpr_fdr_page <- function(geo_sf, title, subtitle, note, age_label,
                               bg_sf = NULL, border_lwd = 0.3, publication = FALSE) {
  with_strip <- function(lab) { geo_sf$strip <- lab; geo_sf }
  frame <- fixed_frame(if (is.null(bg_sf)) geo_sf else bg_sf)
  panel <- function(show_all_keys = NA) {
    c(map_panel(show_all_keys, bg_sf = bg_sf, border_lwd = border_lwd), list(frame))
  }
  # Strip labels are markdown on the publication page (italic PfHRP2)
  lab <- if (publication) {
    c(a = "Sub-regional RDT prevalence x<br>*PfHRP2* deletion risk",
      b = "Seasonal malaria<br>chemoprevention",
      c = "Reported fever in<br>past two weeks",
      d = "Reported seeking care<br>in past two weeks",
      e = "False positivity rate,<br>P(RDT+ | LM-)",
      f = "False discovery rate,<br>P(LM- | RDT+)",
      g = "RDT positive, LM negative")
  } else {
    c(a = "SR RDT prevalence x\nPfHRP2 deletion risk",
      b = "Seasonal malaria chemoprevention\n(coverage > 1%)",
      c = "Fever in past 2 weeks\n(all children)",
      d = "Care seeking, febrile children only\n(only asked after recent fever)",
      e = "False positivity rate,\nP(RDT+ | LM-)",
      f = "False discovery rate,\nP(LM- | RDT+)",
      g = "Absolute false positives,\nP(RDT+ and LM-)")
  }
  # Each panel keeps its own legend underneath: Blues for the survey covariates
  # (C-D, and B's No/Yes ends), the pink-magenta ramp for the model outputs (E-G)
  p_a <- ggplot(with_strip(lab[["a"]]), aes(fill = prev_hrp2)) +
    panel(show_all_keys = TRUE) +
    scale_fill_manual(values = COLS_PREV_HRP2,
                      name = if (publication) "Sub-regional RDT prevalence,\ndeletion risk" else "SR prevalence, deletion risk",
                      drop = FALSE, na.value = MAP_NO_DATA, na.translate = FALSE) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))
  p_b <- ggplot(with_strip(lab[["b"]]), aes(fill = smc_lab)) +
    panel(show_all_keys = TRUE) +
    scale_fill_manual(values = COLS_SMC,
                      name = if (publication) "Seasonal malaria\nchemoprevention" else "SMC",
                      drop = FALSE, na.value = MAP_NO_DATA, na.translate = FALSE) +
    guides(fill = guide_legend(override.aes = list(colour = "grey50", linewidth = 0.2)))
  p_c <- ggplot(with_strip(lab[["c"]]), aes(fill = prev_fever * 100)) +
    panel() + blues_0_100("Fever (%)")
  p_d <- ggplot(with_strip(lab[["d"]]), aes(fill = prop_care * 100)) +
    panel() + blues_0_100(if (publication) "Formal care or ACT (%)" else "Formal care or ACT (%)")
  p_e <- ggplot(with_strip(lab[["e"]]), aes(fill = fpr)) +
    panel() + rate_0_100("FPR (%)")
  p_f <- ggplot(with_strip(lab[["f"]]), aes(fill = fdr)) +
    panel() + rate_0_100("FDR (%)")
  p_g <- ggplot(with_strip(lab[["g"]]), aes(fill = abs_fp)) +
    panel() + rate_0_100("RDT+ and LM- (%)")

  caption <- paste0("E-G: children ", age_label, ", surveyed during/after peak season; ",
                    "predictions include the survey random effect. Beige: not surveyed or no data.", note)
  tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16))

  if (!publication) {
    return(wrap_plots(p_a, p_b, p_c, p_d, p_e, p_f, p_g, plot_spacer(), nrow = 2, ncol = 4) +
      plot_annotation(
        title = title, subtitle = subtitle, caption = caption, tag_levels = "A",
        theme = theme(plot.title    = element_text(face = "bold", size = 22, hjust = 0.5),
                      plot.subtitle = element_text(size = 13, colour = "grey30", hjust = 0.5),
                      plot.caption  = element_text(size = 9, colour = "grey40", hjust = 0))
      ) & tag_theme)
  }

  # One grid for the whole page so rows line up across the vertical rule:
  #   left (survey covariates): A B / C D, A and B with their own legends,
  #   C and D sharing one Blues bar; right (model outputs): E over F, sharing
  #   one pink-magenta bar (P(RDT+ and LM-), panel G on the country pages, is left
  #   off the publication figure). The bars are pulled out of C and F and
  #   placed in their own bottom row. Tags are set per plot because
  #   tag_levels would also letter the headers, rule and bars.
  # Text, keys and tags are ~1.5x MAP_THEME's (larger again for print, 2026-11-10).
  pub_theme <- theme(strip.text   = ggtext::element_markdown(size = 16, face = "bold", margin = margin(b = 4)),
                     plot.tag     = element_text(face = "bold", size = 22),
                     legend.title = element_text(size = 15),
                     legend.text  = element_text(size = 13),
                     legend.key.width = unit(0.7, "cm"), legend.key.height = unit(0.5, "cm"))
  no_leg   <- theme(legend.position = "none")
  tagged   <- function(p, tag, leg = TRUE) p + labs(tag = tag) + pub_theme + if (!leg) no_leg
  pct_name <- "Percentage (%)"
  header   <- function(txt) {
    wrap_elements(full = grid::textGrob(txt, gp = grid::gpar(fontsize = 21, fontface = "bold")))
  }
  bar_leg <- function(p, width_cm) {
    g <- ggplotGrob(p + pub_theme +
                      guides(fill = guide_colourbar(theme = theme(legend.key.width  = unit(width_cm, "cm"),
                                                                   legend.key.height = unit(0.54, "cm")))))
    wrap_elements(full = g$grobs[[which(g$layout$name == "guide-box-bottom")]])
  }
  rule <- wrap_elements(full = grid::segmentsGrob(0.5, 0.01, 0.5, 0.99, gp = grid::gpar(col = "grey40", lwd = 1.5)))
  design <- "
    HHZI
    ABZE
    CDZF
    ccZf
  "
  p <- wrap_plots(
    H = header("Fixed-effect covariates,\nmost recent survey"),
    I = header("Estimates,\nmost recent survey"),
    Z = rule,
    A = tagged(p_a, "A"), B = tagged(p_b, "B"),
    C = tagged(p_c, "C", leg = FALSE), D = tagged(p_d, "D", leg = FALSE),
    E = tagged(p_e, "E", leg = FALSE), F = tagged(p_f, "F", leg = FALSE),
    c = bar_leg(p_c + blues_0_100(pct_name), 10.8),
    f = bar_leg(p_f + rate_0_100("Rate (%)"), 8),
    design = design, widths = c(1, 1, 0.04, 1), heights = c(0.16, 1, 1, 0.14)
  )
  caption <- sub("^E-G:", "E-F:", caption)
  attr(p, "caption") <- paste(c(subtitle, caption), collapse = "\n")
  p
}
