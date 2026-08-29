# Set working directory to this script's location.
# Works in: RStudio (interactive/Source), VSCode source(), Rscript --file=
.script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),          # RStudio
  error = function(e) tryCatch(
    dirname(normalizePath(sys.frames()[[1]]$ofile)),             # source() in VSCode/R terminal
    error = function(e) {
      f <- sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])
      if (length(f) && nzchar(f)) dirname(normalizePath(f))     # Rscript --file=
      else getwd()
    }
  )
)
setwd(.script_dir)

# =============================================================================
# taxotox_install.R
# -----------------------------------------------------------------------------
# Purpose : Build the toxicity reference dataset used by the TaxoTox Shiny app.
#           Queries the local ECOTOX SQLite cache, filters to aquatic organisms
#           (fish / algae / crustacean) and acute lethal/effective-concentration
#           endpoints (LC50 / EC50), normalises all concentrations to ng/L, and
#           writes the result to Data/final_ecotox_data.fst.
#
# Output  : ../Data/final_ecotox_data.fst
#           One row per test result; key columns used by app.R:
#             cas_number   – CAS Registry Number (chemical identity)
#             ecotox_group – taxonomic group: "fish" | "algae" | "crustacean"
#             conc_ng_L    – concentration converted to ng/L
#           app.R takes median(conc_ng_L) per (cas_number, ecotox_group) as the
#           species-sensitivity denominator when computing Toxic Units (TU).
#
# Prerequisites:
#   1. ECOTOX data must be downloaded at least once -- Step 0a below offers to
#      do this interactively, or does it automatically (only when a newer
#      release exists) on non-interactive runs.
#   2. The SQLite path used below (search for dbConnect(RSQLite::SQLite())
#      must match the ECOTOXr cache location.
#      Check with: ECOTOXr::get_ecotox_sqlite_file()
#   3. Step Z at the end of this file cleans up superseded ECOTOX releases
#      left behind by Step 0a (prompted interactively, auto with a sanity
#      guard when non-interactive) -- see that section for details.
#
# Frequency : Re-run whenever the ECOTOX database is updated (≈ quarterly).
#
# Authors  : Yair Suari & Noam Gridish, 2025
# =============================================================================

# =============================================================================
# Section checkpoint system
# -----------------------------------------------------------------------------
# .ask_run(section_id, label) — called before each slow section.
#   • Reads ../Data/install_timestamps.csv to find when this section last ran.
#   • Prints a prompt: "Last run: <date>  Run again? [y/n]"
#   • Returns TRUE (run) or FALSE (skip).
#   • Always returns TRUE when running non-interactively (Rscript).
# .mark_done(section_id) — call at the END of each section to record completion.
# =============================================================================

.ts_file <- "../Data/install_timestamps.csv"

.load_timestamps <- function() {
  if (file.exists(.ts_file))
    tryCatch(read.csv(.ts_file, stringsAsFactors = FALSE),
             error = function(e) data.frame(section = character(), timestamp = character()))
  else
    data.frame(section = character(), timestamp = character())
}

.ask_run <- function(section_id, label) {
  # Non-interactive mode (Rscript): always run
  if (!interactive()) return(TRUE)

  ts <- .load_timestamps()
  last <- ts$timestamp[ts$section == section_id]
  last_str <- if (length(last) == 1 && nzchar(last)) last else "never"

  cat(sprintf(
    "\n\033[1m[%s]\033[0m %s\n  Last completed: %s\n  Run this section? [y/n]: ",
    section_id, label, last_str
  ))
  ans <- tolower(trimws(readline()))
  ans == "y" || ans == "yes"
}

.mark_done <- function(section_id) {
  ts  <- .load_timestamps()
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (section_id %in% ts$section)
    ts$timestamp[ts$section == section_id] <- now
  else
    ts <- rbind(ts, data.frame(section = section_id, timestamp = now,
                               stringsAsFactors = FALSE))
  write.csv(ts, .ts_file, row.names = FALSE)
  message(sprintf("[%s] completed at %s", section_id, now))
}

#######################################################################

# Libraries

library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(tidyverse)
library(ggpattern)
library(ggpubr)
library(ggpmisc)
library(hrbrthemes)
library(remotes)
library(purrr)  
library(readxl)
library(stringr)
library(RSQLite)
library(data.table)
library(fst)
library(tcltk)
library(ECOTOXr)
library(ssdtools)

#######################################################################
# Baseline row counts (pre-run), used by Step Z's cleanup safety check
# -----------------------------------------------------------------------
# Cheap metadata-only reads of whatever final_ecotox_data.fst /
# taxotox_reference.fst already exist on disk, captured before Step 0a or
# Step 1 can overwrite them. NA if this is the first-ever run.
.prev_final_ecotox_n <- tryCatch(nrow(fst::fst("../Data/final_ecotox_data.fst")),
                                  error = function(e) NA_integer_)
.prev_reference_n    <- tryCatch(nrow(fst::fst("../Data/taxotox_reference.fst")),
                                  error = function(e) NA_integer_)

#######################################################################
# Step 0a: ECOTOX database update
# -----------------------------------------------------------------------
# Merged in from the former standalone update_ecotox.R (deleted -- this was
# the only thing it did beyond re-running this file, and having two entry
# points was more confusing than one script with an extra optional step).
# Downloads the latest EPA ECOTOX Knowledgebase release (hundreds of MB,
# 10-30 min) and rebuilds the local SQLite cache Step 1 below reads from.
#
# Interactive runs: manual y/n prompt, same as always -- a human is present
# to watch a multi-hundred-MB download happen.
#
# Non-interactive runs (Rscript): unlike every other section in this file,
# this does NOT use .ask_run() -- that helper auto-answers "yes" for every
# section when run non-interactively, which is right for the comparatively
# cheap, idempotent steps elsewhere in this file but would be wrong here if
# applied blindly (unconditionally re-downloading hundreds of MB on every
# unattended run). Instead, it checks the release version first via
# check_ecotox_version() and only downloads when a newer release actually
# exists -- so scheduled/automated installs stay current without ever
# re-downloading data that's already up to date. A failed/offline version
# check fails open (treated as up to date) so an unattended run never
# hard-blocks on a network hiccup.
if (interactive()) {
  ans <- toupper(trimws(readline(
    "Check for a new ECOTOX release before building the reference table? [Data download, 10-30 min] [y/n]: "
  )))
  if (ans == "Y") {
    message("Downloading latest ECOTOX release...")
    ECOTOXr::download_ecotox_data()
    message("ECOTOX database update complete.")
  }
} else {
  message("Step 0a: non-interactive run -- checking ECOTOX release version...")
  up_to_date <- tryCatch(ECOTOXr::check_ecotox_version(), error = function(e) {
    message("Step 0a: version check failed (", conditionMessage(e),
            ") -- using cached database.")
    TRUE
  })
  if (isFALSE(up_to_date)) {
    message("Step 0a: newer ECOTOX release available -- downloading automatically...")
    ECOTOXr::download_ecotox_data(ask = FALSE)
    message("Step 0a: ECOTOX database update complete.")
  } else {
    message("Step 0a: already on the latest ECOTOX release -- skipping download.")
  }
}

#######################################################################
# Step 0b: Pending curation check
# -----------------------------------------------------------------------
# Compounds encountered during app sessions but not yet in Known_CAS.fst
# are logged to a Google Sheet (see taxotox_curate.R). Curation itself
# requires human judgement per compound (Add/Skip/Reject) and can't be
# automated, so -- same reasoning as Step 0a -- this only ever prompts when
# interactive; a non-interactive run just prints how many items are
# waiting and continues, so nothing blocks or hangs waiting for input that
# will never arrive on an unattended Rscript run.
{
  .curation_sheet_id <- Sys.getenv("TAXOTOX_SHEET_ID",
                                   unset = "1pftfWQNfIStasPqH1CvpDSAH3yHxYmjhggvAwlBaxDA")
  .n_pending_curation <- tryCatch({
    library(googlesheets4)
    key_path <- "taxotox-service.json"
    if (file.exists(key_path)) gs4_auth(path = key_path) else gs4_auth()
    nrow(read_sheet(.curation_sheet_id, col_types = "ccc"))
  }, error = function(e) {
    message("Could not check the curation Google Sheet: ", conditionMessage(e))
    NA_integer_
  })

  if (!is.na(.n_pending_curation) && .n_pending_curation > 0) {
    if (interactive()) {
      ans <- toupper(trimws(readline(sprintf(
        "%d compound(s) pending curation in the Google Sheet. Run taxotox_curate.R now? [y/n]: ",
        .n_pending_curation
      ))))
      if (ans == "Y") {
        # taxotox_curate.R stop()s intentionally when there's nothing left to
        # review (e.g. every pending row turns out to already be in
        # Known_CAS) -- fine when it's run standalone, but a bare source()
        # here would let that stop() abort this entire install run. Catch it
        # and continue instead.
        tryCatch(
          source("taxotox_curate.R"),
          error = function(e) message("Step 0b: taxotox_curate.R stopped early (",
                                       conditionMessage(e), ") -- continuing install.")
        )
      }
    } else {
      message(sprintf(
        "Step 0b: %d compound(s) pending curation in the Google Sheet -- run taxotox_curate.R interactively to review them.",
        .n_pending_curation
      ))
    }
  }
}

#######################################################################
# Step 1–6: Query ECOTOX, normalise units, write final_ecotox_data.fst
# -----------------------------------------------------------------------
# Slow step (~5–10 min): runs a large SQL query against the local ECOTOX
# SQLite database and writes final_ecotox_data.fst (~1–2 M rows).
# Re-run only when the ECOTOX database has been updated (≈ quarterly).

if (.ask_run("S1_ecotox", "ECOTOX query + unit conversion + write final_ecotox_data.fst")) {

conn <- dbConnect(RSQLite::SQLite(), ECOTOXr::get_ecotox_sqlite_file())

#combining all relevant columns from result, test and chemical df

mysearch <- paste0("
-- ============================================================
-- TaxoTox ECOTOX extraction query
-- Purpose: pull all aquatic toxicity results (fish, algae,
--   crustacean) with concentration values converted to a
--   single numeric field (min_concentration) ready for unit
--   conversion in R.
-- Source tables (ECOTOXr SQLite cache):
--   results   - one row per measured endpoint result
--   tests     - experimental conditions, links to species/chemical
--   species   - taxonomic group (ecotox_group)
--   chemicals - chemical name and CAS number
-- ============================================================

-- CTE: enrich each result row with species, chemical, and
-- test metadata, and compute calc_conc1_mean.
-- calc_conc1_mean: use conc1_mean when valid (>0);
--   fall back to (conc1_min + conc1_max)/2 when mean is
--   missing or zero but range bounds are available.
WITH enriched_results AS (
    SELECT
        r.result_id,
        t.test_id,
        t.reference_number,
        t.exposure_type,
        t.exposure_duration_mean,
        t.exposure_duration_unit,
        t.media_type,
        r.endpoint,
        r.trend,
        r.effect,
        r.conc1_unit,
        r.measurement,
        r.measurement_comments,
        r.conc1_mean,
        r.conc2_mean,
        r.conc3_mean,
        r.conc1_min,
        r.conc1_max,
        s.species,
        s.genus,
        s.family,
        s.ecotox_group,
        c.chemical_name,
        c.cas_number,
        -- Step 1: compute adjusted conc1_mean
        CASE
            WHEN (r.conc1_mean IS NULL OR r.conc1_mean <= 0)
                 AND r.conc1_min > 0 AND r.conc1_max > 0
            THEN (r.conc1_min + r.conc1_max) / 2.0
            ELSE r.conc1_mean
        END AS calc_conc1_mean
    FROM results r
    JOIN tests t ON r.test_id = t.test_id
    JOIN species s ON t.species_number = s.species_number
    JOIN chemicals c ON t.test_cas = c.cas_number
)

SELECT
    er.result_id,

    -- Step 2: min_concentration = lowest available concentration
    -- across the up to three concentration columns (conc1/2/3).
    -- Using the minimum follows the conservative approach for
    -- mixture toxicity: we prefer the most sensitive value when
    -- multiple concentrations are reported for a single result.
    -- Only positive (non-zero) values are considered.
    CASE
        WHEN er.calc_conc1_mean > 0 AND er.conc2_mean > 0 AND er.conc3_mean > 0 THEN
            CASE
                WHEN er.calc_conc1_mean <= er.conc2_mean AND er.calc_conc1_mean <= er.conc3_mean THEN er.calc_conc1_mean
                WHEN er.conc2_mean <= er.conc3_mean THEN er.conc2_mean
                ELSE er.conc3_mean
            END
        WHEN er.calc_conc1_mean > 0 AND er.conc2_mean > 0 THEN
            CASE WHEN er.calc_conc1_mean <= er.conc2_mean THEN er.calc_conc1_mean ELSE er.conc2_mean END
        WHEN er.calc_conc1_mean > 0 AND er.conc3_mean > 0 THEN
            CASE WHEN er.calc_conc1_mean <= er.conc3_mean THEN er.calc_conc1_mean ELSE er.conc3_mean END
        WHEN er.conc2_mean > 0 AND er.conc3_mean > 0 THEN
            CASE WHEN er.conc2_mean <= er.conc3_mean THEN er.conc2_mean ELSE er.conc3_mean END
        WHEN er.calc_conc1_mean > 0 THEN er.calc_conc1_mean
        WHEN er.conc2_mean > 0 THEN er.conc2_mean
        WHEN er.conc3_mean > 0 THEN er.conc3_mean
        ELSE NULL
    END AS min_concentration,

    er.endpoint,
    er.cas_number,
    er.trend,
    er.effect,
    er.exposure_type,
    er.exposure_duration_mean,
    er.exposure_duration_unit,
    er.media_type,
    er.measurement,
    er.measurement_comments,
    er.test_id,
    er.reference_number,
    er.conc1_unit,
    t.test_cas,
    er.species,
    er.genus,
    er.family,

    -- Normalise ecotox_group labels: ECOTOX uses varied strings
    -- (e.g. 'Fish', 'Saltwater Fish'); map all to lowercase tokens
    -- so downstream R code can filter with simple equality checks.
    CASE
        WHEN LOWER(er.ecotox_group) LIKE '%fish%' THEN 'fish'
        WHEN LOWER(er.ecotox_group) LIKE '%algae%' THEN 'algae'
        WHEN LOWER(er.ecotox_group) LIKE '%crustacean%' THEN 'crustacean'
        ELSE er.ecotox_group
    END AS ecotox_group,

    er.chemical_name

FROM enriched_results er
JOIN tests t ON er.test_id = t.test_id

-- Keep only the three aquatic organism groups relevant to TaxoTox.
-- All other ecotox_groups (birds, mammals, etc.) are excluded.
WHERE LOWER(er.ecotox_group) LIKE '%fish%'
   OR LOWER(er.ecotox_group) LIKE '%algae%'
   OR LOWER(er.ecotox_group) LIKE '%crustacean%';
")


filterd_ecotox_data <- dbGetQuery(conn, mysearch)

# Step 2: Normalise unit strings to the four supported SI forms.
# "AI" prefix (Active Ingredient) is stripped — we treat lab concentrations as
# active ingredient concentrations regardless of formulation.
# Unit aliases used in ECOTOX:
#   ppt  → ng/L  (parts-per-trillion)
#   ppb  → µg/L  (parts-per-billion)
#   ppm  → mg/L  (parts-per-million)
#   0/00 → g/L   (per-mille, i.e. g/kg ≈ g/L in dilute aqueous solutions)
#   dm³  → L     (cubic decimetre = litre)
# Rows with any other unit (e.g. nmol/L, % body weight) are dropped.

relevent_units <- c("ng/L", "ug/L", "mg/L", "g/L")

filterd_ecotox_data_conc_unit <- filterd_ecotox_data %>%
  mutate(
    conc1_unit = str_remove(conc1_unit, "^AI\\s*"),
    conc1_unit = str_replace_all(conc1_unit, "\\bdm3\\b", "L"),
    conc1_unit = str_replace_all(conc1_unit, "ppt",  "ng/L"),
    conc1_unit = str_replace_all(conc1_unit, "ppm",  "mg/L"),
    conc1_unit = str_replace_all(conc1_unit, "ppb",  "ug/L"),
    conc1_unit = str_replace_all(conc1_unit, "0/00", "g/L")
  ) %>%
  filter(conc1_unit %in% relevent_units)

# Step 3: Build a unit-to-ng/L conversion look-up table.
# factor_to_ng_L is multiplied by min_concentration (in the original unit)
# to produce conc_ng_L — the common scale used throughout app.R.
conversion_df <- tibble(
  conc1_unit     = relevent_units,
  factor_to_ng_L = c(1, 1e3, 1e6, 1e9)
)

# Diagnostic: count rows per endpoint label (informational only, not used downstream).
endpoint_count <- filterd_ecotox_data_conc_unit %>%
  group_by(endpoint) %>%
  count(endpoint)

# Step 4: Apply conversion, back-transform log-reported values, filter endpoints.
final_ecotox_data <- filterd_ecotox_data_conc_unit %>%
  left_join(conversion_df, by = "conc1_unit") %>%
  mutate(min_concentration = as.numeric(min_concentration)) %>%
  # Some ECOTOX entries report log10(concentration) under labels "(log)LC50" /
  # "(log)EC50". Back-transform to the linear scale before unit conversion.
  mutate(min_concentration = if_else(
    str_detect(endpoint, "^\\(log\\)"),
    10 ^ min_concentration,
    min_concentration
  )) %>%
  # Compute final concentration in ng/L.
  mutate(conc_ng_L = min_concentration * factor_to_ng_L) %>%
  # Strip label decorators so endpoint values become canonical "EC50" / "LC50" / "IC50".
  #   "(log)EC50" → "EC50",  "LC50*" → "LC50",  "LC50/" → "LC50",  "IC50/" → "IC50"
  mutate(endpoint = str_remove(endpoint, "^\\(log\\)"),
         endpoint = str_replace_all(endpoint, "[*/]", "")) %>%
  # Keep only acute lethal / effective-/inhibitory-concentration endpoints.
  # These represent the most widely reported, standardised benchmarks in ECOTOX
  # and are the denominators in the Toxic Unit (TU) framework. IC50 is included
  # alongside EC50 for algae in particular: ECOTOX's IC50 rows share the exact
  # same effect/measurement code profile as its EC50 rows (population growth
  # rate, biomass, chlorophyll, photosynthesis) -- different source studies
  # label the identical growth-inhibition assay "IC50" vs "EC50" depending on
  # publication convention, not a biological difference. Confirmed empirically:
  # 234 compounds have BOTH an EC50 and an IC50 algae result in ECOTOX. See
  # Docs/TaxoTox_Technical_Methods.md Section 5.4.
  filter(endpoint %in% c("LC50", "EC50", "IC50")) %>%
  mutate(effect = str_replace_all(effect, "[~/]", ""))

# Step 5: Summary table (median per chemical × taxonomic group).
# Not used by app.R directly — kept for diagnostics.
taxotox_data <- final_ecotox_data %>%
  group_by(cas_number, ecotox_group) %>%
  summarize(median_min_conc  = median(min_concentration, na.rm = TRUE),
            observation_count = n())

# Step 6: Write to FST for fast loading in app.R.
# app.R reads this file at startup and computes median(conc_ng_L) per
# (cas_number, ecotox_group) as the species-sensitivity denominator.
write_fst(final_ecotox_data, "../Data/final_ecotox_data.fst", compress = 50)
.mark_done("S1_ecotox")

} else {
  # Load existing file so downstream steps have final_ecotox_data in memory
  message("S1_ecotox skipped — loading existing final_ecotox_data.fst...")
  final_ecotox_data <- read.fst("../Data/final_ecotox_data.fst", as.data.table = FALSE)
}

# =============================================================================
# Step 7 (I-1a): Build the skeleton of taxotox_reference.fst
# -----------------------------------------------------------------------------
# taxotox_reference.fst is the pre-aggregated denominator table used by app.R
# instead of the runtime median computation.  One row per (cas_number,
# ecotox_group).  This step populates:
#   median_lc50_ng_L  – ECOTOX median (current TU denominator)
#   n_ecotox          – number of ECOTOX test results behind the median
#   hc5_ssd_ng_L      – HC5 from Species Sensitivity Distribution (n >= 5 only)
#   hc5_model_ng_L    – HC5 from external model (populated in step I-1b; NA here)
#   hc5_method        – "SSD" | "Model" | "Both" | NA
#
# Subsequent steps add: predicted_lc50_ng_L (I-2), moa_group (I-3),
# benchmark column (I-6), stc_nowell_ng_L (I-10, separate script
# taxotox_install_nowell.R), then write_fst() is called at step I-4.
# =============================================================================

# --- 7a: Per-(cas_number, ecotox_group) aggregates --------------------------
reference <- final_ecotox_data %>%
  group_by(cas_number, ecotox_group) %>%
  summarise(
    chemical_name    = first(chemical_name),   # one name per CAS (multiple synonyms → first)
    n_ecotox         = n(),
    median_lc50_ng_L = median(conc_ng_L, na.rm = TRUE),
    .groups = "drop"
  )

# --- 7a-ii: Nowell et al. (2014) "Taxon-Sensitive PTI" denominator ----------
# stc_nowell_ng_L / n_stc_nowell are NOT computed here. They used to be
# re-derived from this script's own ECOTOX pull (5th-percentile/minimum STC,
# duration/measurement filtering, cladoceran-genus restriction) but that
# diverged substantially from Nowell et al.'s actual published values (e.g.
# Diazinon fish STC: ~45 ug/L recomputed here vs 85 ug/L published) because
# Nowell's STC also draws on OPP aquatic-life benchmarks and the Pesticide
# Properties Database, not ECOTOX alone, and because ECOTOX has been updated
# since 2014. Run taxotox_install_nowell.R (step I-10) after this script to
# join Nowell's own published Appendix B values directly instead.

# --- 7b: SSD-based HC5 (hc5_ssd_ng_L) --------------------------------------
# SLOW step (~15–30 min): fits a log-normal SSD for every compound×group pair
# with n ≥ 5 ECOTOX results. Re-run only when final_ecotox_data.fst changes.
# If skipped, hc5_ssd_ng_L is loaded from the existing taxotox_reference.fst.

if (.ask_run("S7b_ssd", "SSD fitting for all compound\u00d7group pairs with n \u2265 5 (slow ~15\u201330 min)")) {

# Fits a log-normal Species Sensitivity Distribution to all LC50/EC50 values
# for each compound × group pair with at least 5 data points.
# HC5 = 5th percentile of the fitted distribution = the concentration at which
# 5% of species would be affected (protecting 95% of species).
#
# Log-normal is the standard regulatory default (ANZG 2018, EU TGD, REACH).
# Fitting failures (non-convergence, degenerate data) return NA silently.
#
# NOTE: hc5_ssd_ng_L covers ~15–20% of compound×group rows in ECOTOX.
# Rows with fewer data points receive hc5_model_ng_L in step I-1b.
# Both columns are kept separate so the Excel output can report which method
# was used for each compound (see plan items A-13, A-14, A-16).

fit_hc5_ssd_full <- function(values) {
  # values: numeric vector of conc_ng_L for one compound x group pair
  # Returns a 1-row tibble: hc5_ssd_ng_L (ng/L), sigma_ssd_lnorm (audit only),
  # and sigma_ssd_log10 (the one actually consumed downstream). All NA_real_
  # if fitting fails or n < 5.
  #
  # UNITS WARNING: ssdtools' ssd_fit_dists(dists="lnorm") reports `sdlog` in
  # the same convention as base R's plnorm/qlnorm -- a NATURAL-log SD (i.e.
  # log(X) ~ N(meanlog, sdlog), confirmed empirically against a known-sigma
  # simulation: fitted sdlog matches log10_sigma * log(10), not log10_sigma
  # itself). sigma_group_fallback (step 7c, below) is a log10 SD by
  # construction (it inverts a log10 relationship). msPAF_EC50's formula
  # (app.R's .calc_hc5(), Docs/TaxoTox_Technical_Methods.md Section 5.5) is
  # pnorm(log10(S)/sigma_mix) -- a log10-based formula -- so sigma_ec50 (the
  # coalesce of these two in step 7c) MUST be log10 throughout. Converting
  # here, once, at the source, rather than trusting every downstream
  # consumer to remember the conversion.
  #
  # sigma_ssd_log10 is sigma_i in the multi-substance Potentially Affected
  # Fraction (msPAF_EC50) formula (de Zwart & Posthuma, 2005, Environ.
  # Toxicol. Chem. 24(11):2665-2679, Eq. 4-5) -- extracted from this exact
  # SSD fit (not re-fitted separately) so it stays consistent with
  # hc5_ssd_ng_L above. See Docs/TaxoTox_Technical_Methods.md Section 5.5.
  # Warnings (e.g. optimiser non-convergence) are suppressed -- they are
  # expected for extreme data distributions and handled by returning NA.
  .na_result <- tibble(hc5_ssd_ng_L = NA_real_, sigma_ssd_lnorm = NA_real_,
                       sigma_ssd_log10 = NA_real_)
  withCallingHandlers(
    tryCatch({
      df <- data.frame(Conc = values[values > 0 & is.finite(values)])
      if (nrow(df) < 5) {
        return(.na_result)
      }
      fit <- ssd_fit_dists(df, dists = "lnorm")
      hc  <- ssd_hc(fit, proportion = 0.05, ci = FALSE)

      # coef(fit) returns a tibble (dist, term, est, se); for a single
      # "lnorm" fit this has two rows, term %in% c("meanlog", "sdlog").
      co        <- coef(fit)
      sdlog_val <- co$est[co$dist == "lnorm" & co$term == "sdlog"]
      sdlog_val <- if (length(sdlog_val) == 1) sdlog_val else NA_real_

      tibble(
        hc5_ssd_ng_L    = as.numeric(hc$est[[1]]),
        sigma_ssd_lnorm = sdlog_val,                    # natural-log, audit trail only
        sigma_ssd_log10 = sdlog_val / log(10)            # log10 -- what sigma_ec50 needs
      )
    }, error = function(e) .na_result),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

eligible <- final_ecotox_data %>%
  filter(!is.na(conc_ng_L), conc_ng_L > 0) %>%
  group_by(cas_number, ecotox_group) %>%
  filter(n() >= 5) %>%
  ungroup()

n_eligible <- n_distinct(eligible[c("cas_number", "ecotox_group")])
message(sprintf("Step 7b: fitting SSD for %d compound\u00d7group pairs (n \u2265 5)...",
                n_eligible))

hc5_ssd <- eligible %>%
  group_by(cas_number, ecotox_group) %>%
  group_modify(~ fit_hc5_ssd_full(.x$conc_ng_L)) %>%
  ungroup()

# Join HC5 onto reference; initialise model column and method flag
reference <- reference %>%
  left_join(hc5_ssd, by = c("cas_number", "ecotox_group")) %>%
  mutate(
    hc5_model_ng_L = NA_real_,          # populated in step I-1b
    hc5_method     = case_when(
      !is.na(hc5_ssd_ng_L) ~ "SSD",
      TRUE                 ~ NA_character_
    )
  )

message(sprintf(
  "Step 7b complete: hc5_ssd_ng_L computed for %d rows; %d rows remain NA (will use scaling, step I-1b).",
  sum(!is.na(reference$hc5_ssd_ng_L)),
  sum(is.na(reference$hc5_ssd_ng_L))
))
.mark_done("S7b_ssd")

} else {
  message("S7b_ssd skipped -- loading hc5_ssd_ng_L / sigma_ssd_log10 from existing taxotox_reference.fst...")
  if (file.exists("../Data/taxotox_reference.fst")) {
    prev <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
      select(cas_number, ecotox_group, any_of(c("hc5_ssd_ng_L", "sigma_ssd_lnorm", "sigma_ssd_log10"))) %>%
      mutate(cas_number = as.character(cas_number))
    reference <- reference %>%
      mutate(cas_number = as.character(cas_number)) %>%
      left_join(prev, by = c("cas_number", "ecotox_group")) %>%
      mutate(
        hc5_model_ng_L = NA_real_,
        hc5_method     = case_when(!is.na(hc5_ssd_ng_L) ~ "SSD", TRUE ~ NA_character_)
      )
    # sigma_ssd_lnorm (natural-log sdlog, audit only) was added to
    # taxotox_reference.fst's schema alongside this change; sigma_ssd_log10
    # (what step 7c's coalesce() actually needs -- msPAF_EC50 is a log10
    # formula) was added later still. any_of() above silently omits whichever
    # column a pre-existing file predates, so backfill both: derive
    # sigma_ssd_log10 from sigma_ssd_lnorm / log(10) if a file has the old
    # natural-log-only column but not yet the log10 one.
    if (!"sigma_ssd_lnorm" %in% names(reference)) reference$sigma_ssd_lnorm <- NA_real_
    if (!"sigma_ssd_log10" %in% names(reference)) {
      reference$sigma_ssd_log10 <- reference$sigma_ssd_lnorm / log(10)
    }
    message(sprintf("  Loaded hc5_ssd_ng_L for %d rows (sigma_ssd_log10 for %d).",
                    sum(!is.na(reference$hc5_ssd_ng_L)),
                    sum(!is.na(reference$sigma_ssd_log10))))
  } else {
    warning("No existing taxotox_reference.fst found — hc5_ssd_ng_L will be NA. Run S7b_ssd at least once.")
    reference <- reference %>%
      mutate(hc5_ssd_ng_L = NA_real_, sigma_ssd_lnorm = NA_real_, sigma_ssd_log10 = NA_real_,
             hc5_model_ng_L = NA_real_, hc5_method = NA_character_)
  }
}

# --- 7c (I-1b): Empirical scaling factor HC5 (hc5_model_ng_L) --------------
# Standard practice (EU TGD, REACH) uses an Assessment Factor (AF) when fewer
# than 5 species are available — typically AF = 10 for a 3-taxon dataset.
# Here we replace the fixed AF with an empirically calibrated factor derived
# from data-rich compounds in ECOTOX (those with n_ecotox >= 5 where a proper
# SSD was fitted in step 7b).
#
# k = median(hc5_ssd_ng_L / median_lc50_ng_L)  per ecotox_group
#
# k is computed separately per group because the HC5/median ratio differs
# between fish, algae, and crustacean due to differences in inter-species
# variability within each group.
#
# hc5_model_ng_L is populated for ALL rows with a median_lc50_ng_L (i.e.
# any row with ECOTOX data), including data-rich rows — this allows direct
# comparison between SSD-derived and scaled estimates.
#
# For synthetic rows with n_ecotox = 0 (added at step I-2 onwards),
# hc5_model_ng_L = predicted_lc50_ng_L × k; that assignment is deferred
# until predicted_lc50_ng_L is available.
#
# The scaling factors are saved to Data/hc5_scaling_factors.csv for
# audit trails and reuse without re-running the full install.

scaling_factors <- reference %>%
  filter(!is.na(hc5_ssd_ng_L), !is.na(median_lc50_ng_L), median_lc50_ng_L > 0) %>%
  mutate(ratio = hc5_ssd_ng_L / median_lc50_ng_L) %>%
  group_by(ecotox_group) %>%
  summarise(
    scaling_factor = median(ratio, na.rm = TRUE),
    n_calibration  = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Group-level fallback sigma for msPAF_EC50 (de Zwart & Posthuma, 2005,
    # Environ. Toxicol. Chem. 24(11):2665-2679, Eq. 4), used for the ~80% of
    # compound x group rows with no direct SSD fit (sigma_ssd_lnorm NA).
    # Inverts the log-normal SSD relationship already used above for
    # hc5_model_ng_L:
    #   HC5 = median_LC50 * 10^(qnorm(0.05) * sigma)
    #       = median_LC50 * scaling_factor   [scaling_factor == k, this group]
    # so   sigma_group_fallback = log10(scaling_factor) / qnorm(0.05)
    # qnorm(0.05) is used directly (not a hardcoded -1.6449) so this stays
    # tied to the proportion = 0.05 passed to ssd_hc() in step 7b above --
    # if that changes, this does too. See Posthuma & de Zwart (2012),
    # Environ. Toxicol. Chem. 31(10):2175-2181, for the fitted-over-fallback
    # priority used when this is joined back onto `reference` below.
    sigma_group_fallback = log10(scaling_factor) / qnorm(0.05)
  )

message("Step 7c: HC5/median scaling factors per taxonomic group:")
for (i in seq_len(nrow(scaling_factors))) {
  message(sprintf("  %-12s  k = %.4f  sigma_fallback = %.4f  (calibrated from %d compounds)",
                  scaling_factors$ecotox_group[i],
                  scaling_factors$scaling_factor[i],
                  scaling_factors$sigma_group_fallback[i],
                  scaling_factors$n_calibration[i]))
}

write.csv(scaling_factors, "../Data/hc5_scaling_factors.csv", row.names = FALSE)

reference <- reference %>%
  left_join(scaling_factors, by = "ecotox_group") %>%
  mutate(
    hc5_model_ng_L = if_else(
      !is.na(median_lc50_ng_L) & median_lc50_ng_L > 0,
      median_lc50_ng_L * scaling_factor,
      NA_real_
    ),
    hc5_method = case_when(
      !is.na(hc5_ssd_ng_L) & !is.na(hc5_model_ng_L) ~ "Both",
      !is.na(hc5_ssd_ng_L)                            ~ "SSD",
      !is.na(hc5_model_ng_L)                          ~ "Scaled",
      TRUE                                             ~ NA_character_
    ),
    # sigma_ec50: the sigma actually consumed by app.R's msPAF_EC50
    # calculation (de Zwart & Posthuma, 2005, Eq. 4-6). Prefers this
    # compound's own fitted SSD sdlog when available, otherwise falls back
    # to this taxonomic group's calibrated sigma -- the same priority order
    # already used for hc5_denom = coalesce(hc5_ssd_ng_L, hc5_model_ng_L),
    # so any row with a usable HC5 denominator also has a usable sigma.
    # MUST use sigma_ssd_log10, not sigma_ssd_lnorm: sigma_group_fallback is
    # a log10 SD (log10(k)/qnorm(0.05)) and app.R's msPAF_EC50 formula
    # (pnorm(log10(S)/sigma_mix)) is log10-based throughout, but
    # sigma_ssd_lnorm is ssdtools' natural-log sdlog -- coalescing that
    # straight in here mixed two log bases (a factor of log(10) ~ 2.303) and
    # badly overestimated msPAF_EC50 for every SSD-fitted compound. See
    # fit_hc5_ssd_full() above for the empirical confirmation.
    sigma_ec50 = dplyr::coalesce(sigma_ssd_log10, sigma_group_fallback)
  ) %>%
  select(-scaling_factor, -n_calibration, -sigma_group_fallback)

message(sprintf(
  "Step 7c complete: hc5_model_ng_L populated for %d rows (%d SSD+Scaled, %d Scaled-only, %d SSD-only).",
  sum(!is.na(reference$hc5_model_ng_L)),
  sum(reference$hc5_method == "Both",  na.rm = TRUE),
  sum(reference$hc5_method == "Scaled", na.rm = TRUE),
  sum(reference$hc5_method == "SSD",   na.rm = TRUE)
))

# Step I-4 (partial): Write taxotox_reference.fst after I-1a/I-1b.
# Overwritten again after I-2 if the CompTox API key is available.
write_fst(reference, "../Data/taxotox_reference.fst", compress = 50)
message(sprintf("taxotox_reference.fst written (pre-API): %d rows, %d columns.",
                nrow(reference), ncol(reference)))

# =============================================================================
# Step I-2: CompTox/OPERA gap-fill
# -----------------------------------------------------------------------------
# Queries the EPA CompTox Dashboard (CCTE API) for OPERA-predicted aquatic
# LC50/EC50 values for all compounds in Known_CAS.fst.
#
# What this step does:
#   a) Resolves CAS → DTXSID via search endpoint (one call per CAS)
#   b2) Fetches molecular weight via batch detail endpoint (averageMass field)
#   c) Fetches OPERA aquatic toxicity predictions per DTXSID (batch)
#   d) Converts predictions from mol/L to ng/L via MW
#   d) Adds predicted_lc50_ng_L to existing reference rows (ECOTOX-covered)
#   e) Creates new rows (n_ecotox = 0) for Known_CAS compounds missing from ECOTOX
#   f) Applies hc5_model_ng_L = predicted_lc50_ng_L × k for n_ecotox = 0 rows
#   g) Sets lc50_source: "ECOTOX_median" | "CompTox_predicted" | "Both"
#
# API key setup (free, register at https://www.epa.gov/comptox-tools):
#   Option A — environment variable (add to ~/.Renviron):
#              COMPTOX_API_KEY=your-key-here
#   Option B — key file (excluded from git):
#              save key as a single line in Code/comptox_api_key.txt
#
# CTX API base URL : https://comptox.epa.gov/ctx-api
# OPERA properties used (units: mol/L):
#   Fish      : "96 Hour Fathead Minnow LC50"
#   Crustacean: "48 Hour Daphnia Magna LC50"
#   Algae     : not available (no suitable OPERA endpoint)
#
# If no API key is present, this step is skipped gracefully and the reference
# file written above (ECOTOX-only) is kept as-is.
# =============================================================================

library(httr)
library(jsonlite)

# Null-coalescing helper (avoids rlang dependency)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# --- API key -----------------------------------------------------------------

.get_api_key <- function() {
  key <- Sys.getenv("COMPTOX_API_KEY", unset = "")
  if (nzchar(key)) return(key)
  key_file <- "comptox_api_key.txt"
  if (file.exists(key_file)) return(trimws(readLines(key_file, n = 1L)))
  ""
}

COMPTOX_KEY <- .get_api_key()

if (!nzchar(COMPTOX_KEY)) {
  message(paste(
    "Step I-2: SKIPPED — CompTox API key not found.",
    "  To enable gap-fill, register free at https://www.epa.gov/comptox-tools then either:",
    "    a) add  COMPTOX_API_KEY=your-key  to ~/.Renviron, or",
    "    b) save the key as Code/comptox_api_key.txt  (git-ignored).",
    sep = "\n"
  ))
} else if (!.ask_run("S_i2_comptox",
    "CompTox/OPERA LC50 gap-fill — CAS\u2192DTXSID lookup + batch property fetch (~5\u201315 min)")) {
  message("S_i2_comptox skipped — existing predicted_lc50_ng_L values (if any) preserved in taxotox_reference.fst.")
} else {

  COMPTOX_BASE <- "https://comptox.epa.gov/ctx-api"
  BATCH_SIZE   <- 200L    # search endpoint: one call per CAS (sequential)
                          # property endpoint: max 1000 DTXSIDs per POST

  # --- Safe HTTP helpers ---------------------------------------------------

  # A stalled TCP connection (no response, no RST) hangs httr's underlying
  # curl call indefinitely if no timeout is set -- observed directly: a
  # rebuild sat blocked for 5.7+ hours on a single POST with ~0 CPU usage.
  # 30s is generous for a JSON API response; a genuinely stalled connection
  # will not resolve by waiting longer, so this timeout is what lets
  # max_tries actually retry instead of the whole pipeline hanging forever.
  .ctox_get <- function(path, ..., max_tries = 3L) {
    url <- paste0(COMPTOX_BASE, path)
    for (i in seq_len(max_tries)) {
      resp <- tryCatch(
        GET(url, add_headers("x-api-key" = COMPTOX_KEY,
                             "accept" = "application/json"),
            httr::timeout(30), ...),
        error = function(e) NULL
      )
      if (!is.null(resp) && status_code(resp) == 200L) return(resp)
      if (i < max_tries) Sys.sleep(i)
    }
    NULL
  }

  # 120s (not 30s): the /chemical/property/predicted/search/by-dtxsid/
  # endpoint returns EVERY predicted property per compound (hundreds of QSAR
  # model rows), not just the 2 we filter for client-side. Measured directly:
  # 200 DTXSIDs -> ~19MB / ~94s; 1000 DTXSIDs didn't finish in 120s (10.9MB
  # received and still going). 120s covers .fetch_props's smaller batch size
  # (see below) with margin, and is still generous for the much lighter
  # /chemical/detail/search/by-dtxsid/ (MW) endpoint's 1000-item batches.
  .ctox_post <- function(path, body_vec, ..., max_tries = 3L) {
    url  <- paste0(COMPTOX_BASE, path)
    body <- toJSON(body_vec, auto_unbox = FALSE)
    for (i in seq_len(max_tries)) {
      resp <- tryCatch(
        POST(url,
             add_headers("x-api-key"    = COMPTOX_KEY,
                         "accept"       = "application/json",
                         "content-type" = "application/json"),
             body = body, encode = "raw", httr::timeout(120), ...),
        error = function(e) NULL
      )
      if (!is.null(resp) && status_code(resp) == 200L) return(resp)
      if (i < max_tries) Sys.sleep(i)
    }
    NULL
  }

  # --- I-2a: Load Known_CAS ------------------------------------------------

  message("Step I-2a: Loading Known_CAS...")
  known_cas <- read.fst("../Data/Known_CAS.fst", as.data.table = FALSE)

  # Normalise the CAS column name (handles cas_number / casrn / CAS)
  cas_col <- intersect(c("cas_number", "casrn", "CAS", "CASRN"), names(known_cas))[1]
  if (is.na(cas_col))
    stop("Known_CAS.fst has no recognisable CAS column (expected: cas_number or casrn).")
  target_cas <- unique(as.character(known_cas[[cas_col]]))
  message(sprintf("  %d unique CAS numbers in Known_CAS", length(target_cas)))

  # --- I-2b: CAS → DTXSID via /chemical/search/equal/{casrn} (one call per CAS)
  # Note: the search endpoint does NOT return molWeight — MW is fetched separately
  # in I-2b2 via the batch detail endpoint.

  message("Step I-2b: CAS \u2192 DTXSID lookup (", length(target_cas), " compounds)...")

  # CAS->DTXSID is a static/definitional mapping (not time-varying), so a
  # prior run's cache is safe to reuse indefinitely -- only CAS numbers not
  # already in the cache need a live API call. This turns a ~2,600-call,
  # multi-hour sequential lookup into a few-second cache load on every rerun
  # except when Known_CAS has genuinely new CAS numbers.
  .dtxsid_cache_path <- "../Data/dtxsid_map.csv"
  dtxsid_cache <- if (file.exists(.dtxsid_cache_path)) {
    tryCatch(read.csv(.dtxsid_cache_path, stringsAsFactors = FALSE, colClasses = "character"),
             error = function(e) NULL)
  } else NULL

  cached_cas <- if (!is.null(dtxsid_cache)) intersect(target_cas, dtxsid_cache$cas_number) else character(0)
  new_cas    <- setdiff(target_cas, cached_cas)

  message(sprintf("  %d / %d CAS reused from cache; %d need a live lookup",
                  length(cached_cas), length(target_cas), length(new_cas)))

  dtxsid_new <- purrr::map_dfr(seq_along(new_cas), function(i) {
    cas <- new_cas[i]
    if (i %% 100 == 0)
      message(sprintf("  %d / %d", i, length(new_cas)))

    resp <- .ctox_get(paste0("/chemical/search/equal/", utils::URLencode(cas, reserved = TRUE)))
    if (is.null(resp))
      return(tibble(cas_number = cas, dtxsid = NA_character_, preferred_name = NA_character_))

    parsed <- tryCatch(
      fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE),
      error = function(e) NULL
    )
    if (is.null(parsed) || nrow(parsed) == 0)
      return(tibble(cas_number = cas, dtxsid = NA_character_, preferred_name = NA_character_))

    # Prefer exact CASRN match; fall back to rank-1 result
    hit <- parsed[!is.na(parsed$casrn) & parsed$casrn == cas, ]
    if (nrow(hit) == 0) hit <- parsed[1, ]

    tibble(
      cas_number     = cas,
      dtxsid         = as.character(hit$dtxsid[1]        %||% NA),
      preferred_name = as.character(hit$preferredName[1] %||% NA)
    )
  })

  dtxsid_map <- bind_rows(
    if (length(cached_cas) > 0)
      dtxsid_cache %>% filter(cas_number %in% cached_cas) %>%
        mutate(mw_g_mol = suppressWarnings(as.numeric(mw_g_mol)))
    else tibble(),
    dtxsid_new
  )

  message(sprintf("  DTXSID resolved : %d / %d",
                  sum(!is.na(dtxsid_map$dtxsid)), length(target_cas)))

  # --- I-2b2: Batch MW from /chemical/detail/search/by-dtxsid/ ---------------
  # The search endpoint does not return molWeight; the detail endpoint does
  # (field: averageMass). One batch POST per 1000 DTXSIDs.

  message("Step I-2b2: Fetching molecular weights (batch detail)...")

  .fetch_detail_mw <- function(dtxsids, batch_size = 1000L) {
    batches <- split(dtxsids, ceiling(seq_along(dtxsids) / batch_size))
    purrr::map_dfr(seq_along(batches), function(i) {
      batch <- batches[[i]]
      message(sprintf("  MW batch %d / %d (%d DTXSIDs)", i, length(batches), length(batch)))
      resp <- .ctox_post("/chemical/detail/search/by-dtxsid/", batch)
      if (is.null(resp)) { warning(sprintf("  MW batch %d failed.", i)); return(tibble()) }
      parsed <- tryCatch(
        fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE),
        error = function(e) NULL
      )
      if (is.null(parsed) || nrow(parsed) == 0) return(tibble())
      parsed %>%
        select(dtxsid, mw_g_mol = averageMass) %>%
        mutate(mw_g_mol = as.numeric(mw_g_mol))
    })
  }

  # Cached rows already carry mw_g_mol from a prior run's I-2b2 -- only fetch
  # MW for dtxsids that don't have it yet (new compounds, or cached ones that
  # never resolved a value).
  if (!"mw_g_mol" %in% names(dtxsid_map)) dtxsid_map$mw_g_mol <- NA_real_
  needs_mw    <- dtxsid_map$dtxsid[!is.na(dtxsid_map$dtxsid) & is.na(dtxsid_map$mw_g_mol)]
  message(sprintf("  %d / %d DTXSIDs already have MW from cache; %d need a fetch",
                  sum(!is.na(dtxsid_map$dtxsid) & !is.na(dtxsid_map$mw_g_mol)),
                  sum(!is.na(dtxsid_map$dtxsid)), length(needs_mw)))
  mw_table    <- .fetch_detail_mw(needs_mw)

  dtxsid_map  <- dtxsid_map %>%
    left_join(mw_table, by = "dtxsid", suffix = c("", "_new")) %>%
    mutate(mw_g_mol = dplyr::coalesce(mw_g_mol, mw_g_mol_new)) %>%
    select(-any_of("mw_g_mol_new"))

  message(sprintf("  MW available    : %d / %d DTXSIDs",
                  sum(!is.na(dtxsid_map$mw_g_mol)), sum(!is.na(dtxsid_map$dtxsid))))

  all_dtxsids <- dtxsid_map$dtxsid[!is.na(dtxsid_map$dtxsid)]

  # --- I-2c: OPERA aquatic toxicity predictions (batch POST) ---------------
  #
  # Endpoint  : POST /chemical/property/predicted/search/by-dtxsid/
  # Body      : JSON array of DTXSIDs (max 1000 per request)
  # Units     : mol/L  (NOT -log10 — confirmed empirically)
  # Properties used:
  #   Fish      : "96 Hour Fathead Minnow LC50"
  #   Crustacean: "48 Hour Daphnia Magna LC50"
  #   Algae     : not available — predicted_lc50_ng_L left NA for algae
  #               (Tetrahymena IGC50 is the only ecotoxicity proxy but is
  #                not an algal endpoint and is scientifically inappropriate)

  PROP_FISH  <- "96 Hour Fathead Minnow LC50"
  PROP_CRUST <- "48 Hour Daphnia Magna LC50"

  # Helper: batch-fetch a property endpoint and pivot to wide format
  # batch_size = 150, not 1000: this endpoint returns every predicted
  # property per compound (hundreds of rows each), so a 1000-item batch's
  # response is enormous and unreliable -- measured directly: 200 items took
  # ~94s/~19MB, 1000 items didn't complete in 120s. 150 keeps each request
  # comfortably inside .ctox_post's 120s timeout.
  .fetch_props <- function(dtxsids, endpoint, prop_names, batch_size = 150L) {
    batches <- split(dtxsids, ceiling(seq_along(dtxsids) / batch_size))
    purrr::map_dfr(seq_along(batches), function(i) {
      batch <- batches[[i]]
      message(sprintf("  Batch %d / %d (%d DTXSIDs)", i, length(batches), length(batch)))
      resp <- .ctox_post(endpoint, batch)
      if (is.null(resp)) { warning(sprintf("  Batch %d failed.", i)); return(tibble()) }
      parsed <- tryCatch(
        fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE),
        error = function(e) NULL
      )
      if (is.null(parsed) || nrow(parsed) == 0) return(tibble())
      parsed %>%
        filter(propName %in% prop_names) %>%
        select(dtxsid, propName, propValue) %>%
        mutate(propValue = as.numeric(propValue)) %>%
        tidyr::pivot_wider(names_from = propName, values_from = propValue,
                           values_fn = first)
    })
  }

  # --- I-2c: OPERA ecotoxicity predictions ---------------------
  message("Step I-2c: Fetching OPERA ecotoxicity predictions...")
  raw_tox <- .fetch_props(all_dtxsids,
                          "/chemical/property/predicted/search/by-dtxsid/",
                          c(PROP_FISH, PROP_CRUST))

  fish_col  <- PROP_FISH
  crust_col <- PROP_CRUST

  message(sprintf("  OPERA fish predictions       : %d / %d DTXSIDs",
                  if (fish_col  %in% names(raw_tox)) sum(!is.na(raw_tox[[fish_col]]))  else 0L,
                  nrow(raw_tox)))
  message(sprintf("  OPERA crustacean predictions : %d / %d DTXSIDs",
                  if (crust_col %in% names(raw_tox)) sum(!is.na(raw_tox[[crust_col]])) else 0L,
                  nrow(raw_tox)))

  if (nrow(raw_tox) > 0 && !fish_col %in% names(raw_tox)) {
    warning(paste(
      "Fish LC50 property not found — property name may have changed.",
      "Check available names with:",
      "  r <- httr::POST('https://comptox.epa.gov/ctx-api/chemical/property/predicted/search/by-dtxsid/',",
      "    httr::add_headers('x-api-key'=key,'content-type'='application/json'),",
      "    body=jsonlite::toJSON(c('DTXSID8021482')), encode='raw')",
      "  unique(jsonlite::fromJSON(httr::content(r,'text'))$propName)",
      sep = "\n"
    ))
  }

  # --- I-2d: Unit conversion — mol/L → ng/L --------------------------------
  # Units confirmed as mol/L (not -log10).
  # ng/L = mol/L × MW (g/mol) × 1e9
  # MW comes from step I-2b2 (averageMass field in detail endpoint).

  opera_joined <- dtxsid_map %>%       # has mw_g_mol from detail endpoint (step I-2b2)
    left_join(raw_tox, by = "dtxsid")

  # Extract tox columns safely (NA vector if column absent after join)
  fish_vals  <- if (fish_col  %in% names(opera_joined)) opera_joined[[fish_col]]
                else rep(NA_real_, nrow(opera_joined))
  crust_vals <- if (crust_col %in% names(opera_joined)) opera_joined[[crust_col]]
                else rep(NA_real_, nrow(opera_joined))

  opera_wide <- opera_joined %>%
    mutate(
      pred_fish_ng_L  = if_else(
        !is.na(fish_vals)  & !is.na(mw_g_mol) & mw_g_mol > 0,
        fish_vals  * mw_g_mol * 1e9, NA_real_),
      pred_crust_ng_L = if_else(
        !is.na(crust_vals) & !is.na(mw_g_mol) & mw_g_mol > 0,
        crust_vals * mw_g_mol * 1e9, NA_real_)
      # algae: no suitable OPERA endpoint — left NA by design
    )

  # Pivot to long format — fish and crustacean only
  opera_long <- bind_rows(
    opera_wide %>% transmute(cas_number, chemical_name = preferred_name,
                             ecotox_group = "fish",       predicted_lc50_ng_L = pred_fish_ng_L),
    opera_wide %>% transmute(cas_number, chemical_name = preferred_name,
                             ecotox_group = "crustacean", predicted_lc50_ng_L = pred_crust_ng_L)
  ) %>%
    filter(!is.na(predicted_lc50_ng_L))

  message(sprintf("Step I-2d: %d compound×group pairs with predicted LC50 (ng/L)",
                  nrow(opera_long)))

  # --- I-2e: Merge into reference ------------------------------------------

  # Normalise cas_number to character in both tables before joining
  reference  <- reference  %>% mutate(cas_number = as.character(cas_number))
  opera_long <- opera_long %>% mutate(cas_number = as.character(cas_number))

  # Add predicted_lc50_ng_L to existing ECOTOX rows
  reference <- reference %>%
    left_join(
      opera_long %>% select(cas_number, ecotox_group, predicted_lc50_ng_L),
      by = c("cas_number", "ecotox_group")
    )

  # Create n_ecotox = 0 rows for Known_CAS compounds absent from ECOTOX
  existing_keys <- reference %>% select(cas_number, ecotox_group)

  new_rows <- opera_long %>%
    anti_join(existing_keys, by = c("cas_number", "ecotox_group")) %>%
    mutate(
      n_ecotox         = 0L,
      median_lc50_ng_L = NA_real_,
      hc5_ssd_ng_L     = NA_real_,
      sigma_ssd_lnorm  = NA_real_,   # no SSD fit possible with 0 ECOTOX points
      sigma_ssd_log10  = NA_real_,   # no SSD fit possible with 0 ECOTOX points
      hc5_model_ng_L   = NA_real_,
      hc5_method       = NA_character_
    ) %>%
    left_join(scaling_factors, by = "ecotox_group") %>%
    mutate(
      hc5_model_ng_L = if_else(
        !is.na(predicted_lc50_ng_L) & predicted_lc50_ng_L > 0,
        predicted_lc50_ng_L * scaling_factor,
        NA_real_
      ),
      hc5_method = if_else(!is.na(hc5_model_ng_L), "Scaled", NA_character_),
      # sigma_ssd_log10 is NA for all n_ecotox = 0 rows, so this coalesces
      # straight to the group fallback (de Zwart & Posthuma, 2005) -- see
      # the sigma_ec50 comment at step 7c above.
      sigma_ec50 = dplyr::coalesce(sigma_ssd_log10, sigma_group_fallback)
    ) %>%
    select(-any_of(c("scaling_factor", "n_calibration", "sigma_group_fallback")))

  reference <- bind_rows(reference, new_rows)

  # --- I-2f: lc50_source column --------------------------------------------

  reference <- reference %>%
    mutate(
      lc50_source = case_when(
        !is.na(median_lc50_ng_L) & !is.na(predicted_lc50_ng_L) ~ "Both",
        !is.na(median_lc50_ng_L)                                ~ "ECOTOX_median",
        !is.na(predicted_lc50_ng_L)                             ~ "CompTox_predicted",
        TRUE                                                     ~ NA_character_
      )
    )

  message(sprintf(paste(
    "Step I-2 complete: %d total reference rows",
    "  ECOTOX only     : %d",
    "  CompTox only    : %d",
    "  Both            : %d",
    "  No denominator  : %d", sep = "\n"),
    nrow(reference),
    sum(reference$lc50_source == "ECOTOX_median",    na.rm = TRUE),
    sum(reference$lc50_source == "CompTox_predicted", na.rm = TRUE),
    sum(reference$lc50_source == "Both",              na.rm = TRUE),
    sum(is.na(reference$lc50_source))
  ))

  # Overwrite reference with API-enriched version
  write_fst(reference, "../Data/taxotox_reference.fst", compress = 50)
  message(sprintf("taxotox_reference.fst updated (post-API): %d rows, %d columns.",
                  nrow(reference), ncol(reference)))

  # Save dtxsid_map for reuse by I-3 (avoids re-querying the search endpoint)
  write.csv(dtxsid_map, "../Data/dtxsid_map.csv", row.names = FALSE)
  message("dtxsid_map.csv written for reuse by I-3.")

  .mark_done("S_i2_comptox")

}  # end if (nzchar(COMPTOX_KEY))

# =============================================================================
# Step I-3: Mode of Action (MoA) classification
# -----------------------------------------------------------------------------
# moa_group / moa_source are NOT computed here. They used to be a ~100-name
# regex heuristic (organophosphate/carbamate/triazine/phenylurea/pyrethroid
# name patterns) plus a LogKow/Verhaar check that contributed nothing in
# practice, defaulting everything else to "narcosis" -- meaning ~97% of
# compounds got a conservative placeholder, not a real classification. Run
# taxotox_install_moa.R after this script to join a real, externally curated
# Mode-of-Action classification instead: Kramer et al. (2024, Sci. Data,
# https://doi.org/10.1038/s41597-023-02904-7), freely available at
# https://doi.org/10.5281/zenodo.10071824. Unclassified compounds get an
# explicit moa_group = "unknown" rather than being folded into "narcosis".
# =============================================================================

# =============================================================================
# Step I-4: Final write of taxotox_reference.fst
# -----------------------------------------------------------------------------
# This step consolidates everything into the definitive reference table.
# It loads the current taxotox_reference.fst (written by whichever previous
# sections ran), enforces canonical column order, fills any missing columns
# with NA, and writes the final file.
#
# Running this step is safe at any point — it never overwrites data, only
# reorders and canonicalises columns. Re-running after adding the benchmark
# column (step I-6) will produce a file with that column too.
#
# Canonical column order:
#   Identity    : cas_number, ecotox_group, chemical_name
#   ECOTOX data : n_ecotox, median_lc50_ng_L
#   HC5         : hc5_ssd_ng_L, hc5_model_ng_L, hc5_method
#   CompTox     : predicted_lc50_ng_L, lc50_source
#   Mode of action: moa_group, moa_source
#   Benchmarks  : benchmark_usepa_fish_acute_ng_L, _crust_acute_ng_L,
#                 _algae_acute_ng_L (EU EQS / AU ANZG / CA CCME were removed
#                 due to low compound coverage -- see
#                 taxotox_install_benchmarks.R)
# =============================================================================

if (.ask_run("S_i4_finalise",
    "Finalise taxotox_reference.fst — enforce column order and print coverage report")) {

  if (!file.exists("../Data/taxotox_reference.fst")) {
    stop("taxotox_reference.fst not found — run steps I-1a through I-3 first.")
  }

  ref_final <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
    mutate(cas_number = as.character(cas_number))

  # Canonical column list — columns not yet produced are added as NA
  canonical_cols <- c(
    # Identity
    "cas_number", "ecotox_group", "chemical_name",
    # ECOTOX aggregates
    "n_ecotox", "median_lc50_ng_L",
    # HC5
    "hc5_ssd_ng_L", "hc5_model_ng_L", "hc5_method",
    # msPAF_EC50 (de Zwart & Posthuma, 2005, Environ. Toxicol. Chem.
    # 24(11):2665-2679). sigma_ssd_lnorm: SSD-fit sdlog in ssdtools' native
    # NATURAL-log units, retained for audit alongside hc5_ssd_ng_L only --
    # NOT consumed downstream. sigma_ssd_log10: the same value converted to
    # log10 (/ log(10)), which is what sigma_ec50 actually needs, since
    # sigma_group_fallback and app.R's msPAF_EC50 formula are both log10-based.
    # sigma_ec50: coalesce(sigma_ssd_log10, group-level fallback) -- the value
    # app.R actually consumes for the msPAF_EC50 calculation. NOT the same
    # metric as msPAF-NOEC (chronic SSDs; Posthuma et al., 2019, ETC
    # 38(4):905-917) -- see Docs/TaxoTox_Technical_Methods.md Section 5.5.
    "sigma_ssd_lnorm", "sigma_ssd_log10", "sigma_ec50",
    # CompTox gap-fill
    "predicted_lc50_ng_L", "lc50_source",
    # Mode of Action
    "moa_group", "moa_source",
    # Benchmarks (I-6 — NA until that step runs)
    "benchmark_usepa_fish_acute_ng_L",
    "benchmark_usepa_crust_acute_ng_L",
    "benchmark_usepa_algae_acute_ng_L"
  )

  # Add any canonical columns not yet present as NA
  for (col in canonical_cols) {
    if (!col %in% names(ref_final))
      ref_final[[col]] <- NA_real_
  }

  # Reorder: canonical columns first, then any extras in existing order
  extra_cols <- setdiff(names(ref_final), canonical_cols)
  ref_final  <- ref_final[, c(canonical_cols, extra_cols)]

  write_fst(ref_final, "../Data/taxotox_reference.fst", compress = 50)

  # Coverage report
  message(sprintf("\ntaxotox_reference.fst finalised: %d rows, %d columns",
                  nrow(ref_final), ncol(ref_final)))
  message("\nColumn coverage:")
  for (col in canonical_cols) {
    n_present <- sum(!is.na(ref_final[[col]]))
    pct       <- round(100 * n_present / nrow(ref_final), 1)
    status    <- if (pct == 0) "  [NOT RUN]" else if (pct == 100) "  [complete]" else ""
    message(sprintf("  %-42s  %5d / %d  (%5.1f%%)%s",
                    col, n_present, nrow(ref_final), pct, status))
  }

  .mark_done("S_i4_finalise")

}  # end I-4

# =============================================================================
# Step I-12: Join satellite reference tables
# -----------------------------------------------------------------------------
# Runs the three standalone scripts that join external published datasets
# into taxotox_reference.fst, so a single taxotox_install.R run produces a
# complete reference table without separately remembering to run each one.
# Each is also fully runnable on its own afterward (see its own header) --
# e.g. to refresh just one source without re-running everything above.
# Order between the three doesn't matter: each only adds/overwrites its own
# columns, and none of them touches the others' columns.
#   taxotox_install_benchmarks.R -- US EPA aquatic benchmarks (I-6)
#   taxotox_install_nowell.R     -- Nowell et al. (2014) Taxon-Sensitive PTI
#   taxotox_install_moa.R        -- Kramer et al. (2024) Mode of Action
# =============================================================================

if (.ask_run("S_i12_benchmarks",
    "Join US EPA aquatic benchmarks (taxotox_install_benchmarks.R)")) {
  source("taxotox_install_benchmarks.R")
  .mark_done("S_i12_benchmarks")
} else {
  message("S_i12_benchmarks skipped.")
}

if (.ask_run("S_i12_nowell",
    "Join Nowell et al. (2014) Taxon-Sensitive PTI values (taxotox_install_nowell.R)")) {
  source("taxotox_install_nowell.R")
  .mark_done("S_i12_nowell")
} else {
  message("S_i12_nowell skipped.")
}

if (.ask_run("S_i12_moa",
    "Join Kramer et al. (2024) Mode of Action classification (taxotox_install_moa.R)")) {
  source("taxotox_install_moa.R")
  .mark_done("S_i12_moa")
} else {
  message("S_i12_moa skipped.")
}

#######################################################################
# Step Z: ECOTOX SQLite cache housekeeping
# -----------------------------------------------------------------------
# download_ecotox_data() (Step 0a) never deletes superseded releases --
# every dated ecotox_ascii_MM_DD_YYYY.sqlite (+ .log, _cit.txt, extracted
# ASCII folder) accumulates in the ECOTOXr cache dir
# (ECOTOXr::get_ecotox_path()) forever unless removed by hand.
# get_ecotox_sqlite_file() (Step 1 above, and validate_benchmarks.R) always
# auto-selects the newest one, so a stale copy costs disk space (~1-1.5 GB
# each) rather than correctness -- but it's easy to forget about.
#
# Deliberately the LAST step in this script: only clean up a superseded
# release once the full pipeline above has successfully rebuilt
# final_ecotox_data.fst / taxotox_reference.fst and execution reached this
# point without error.
#
# Interactive: list superseded release(s) and prompt before deleting, same
# as every other manual step in this file -- a human is present to review.
#
# Non-interactive: auto-delete, but only if this run's rebuild looks sane --
# compare the row counts just produced against the pre-run baseline
# captured near the top of this file (.prev_final_ecotox_n /
# .prev_reference_n). Skip (and warn) if either dropped below 90% of its
# previous value, or if there was no previous baseline to compare against
# (first-ever run) -- in both cases, leaving the old release in place as a
# fallback is safer than guessing.
.ecotox_cache_dir <- ECOTOXr::get_ecotox_path()
.ecotox_local <- attr(ECOTOXr::check_ecotox_availability(.ecotox_cache_dir), "files")

if (!is.null(.ecotox_local) && nrow(.ecotox_local) > 1) {
  .ecotox_keep  <- .ecotox_local$database[which.max(.ecotox_local$date)]
  .ecotox_stale <- setdiff(.ecotox_local$database, .ecotox_keep)
  .ecotox_stale_stub <- sub("\\.sqlite$", "", .ecotox_stale)

  .delete_stale_ecotox_releases <- function(stubs, cache_dir) {
    for (s in stubs) {
      f_sqlite <- file.path(cache_dir, paste0(s, ".sqlite"))
      f_log    <- file.path(cache_dir, paste0(s, ".log"))
      f_cit    <- file.path(cache_dir, paste0(s, "_cit.txt"))
      f_dir    <- file.path(cache_dir, s)
      if (file.exists(f_sqlite)) file.remove(f_sqlite)
      if (file.exists(f_log))    file.remove(f_log)
      if (file.exists(f_cit))    file.remove(f_cit)
      if (dir.exists(f_dir))     unlink(f_dir, recursive = TRUE)
      message("  Deleted: ", s)
    }
  }

  message(sprintf(
    "\nStep Z: %d superseded ECOTOX release(s) found (keeping %s):",
    length(.ecotox_stale_stub), .ecotox_keep
  ))
  for (s in .ecotox_stale_stub) message("  - ", s, " (.sqlite, .log, _cit.txt, extracted folder)")

  if (interactive()) {
    ans <- toupper(trimws(readline(
      "Delete superseded release(s) now to free disk space? [y/n]: "
    )))
    if (ans == "Y") {
      .delete_stale_ecotox_releases(.ecotox_stale_stub, .ecotox_cache_dir)
    } else {
      message("Step Z: skipped -- superseded release(s) left in place.")
    }
  } else {
    .new_final_ecotox_n <- tryCatch(nrow(fst::fst("../Data/final_ecotox_data.fst")),
                                     error = function(e) NA_integer_)
    .new_reference_n    <- tryCatch(nrow(fst::fst("../Data/taxotox_reference.fst")),
                                     error = function(e) NA_integer_)
    .rebuild_looks_sane <-
      !is.na(.prev_final_ecotox_n) && !is.na(.prev_reference_n) &&
      !is.na(.new_final_ecotox_n) && !is.na(.new_reference_n) &&
      .new_final_ecotox_n >= 0.9 * .prev_final_ecotox_n &&
      .new_reference_n    >= 0.9 * .prev_reference_n

    if (.rebuild_looks_sane) {
      message("Step Z: non-interactive run, rebuild looks sane -- auto-deleting superseded release(s)...")
      .delete_stale_ecotox_releases(.ecotox_stale_stub, .ecotox_cache_dir)
    } else {
      message("Step Z: SKIPPED -- non-interactive run and rebuild size could not be confirmed sane ",
              "(no prior baseline, or row counts dropped >10%). Superseded release(s) left in place; ",
              "review manually before deleting.")
    }
  }
}

#######################################################################


