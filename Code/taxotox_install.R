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
#   1. Run update_ecotox.R once to download ECOTOX data from the EPA website.
#   2. The SQLite path on line 42 must match the ECOTOXr cache location.
#      Check with: ECOTOXr::get_ecotox_sqlite_file()
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
    er.measurement,
    er.measurement_comments,
    er.test_id,
    er.reference_number,
    er.conc1_unit,
    t.test_cas,
    er.species,

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
  # Strip label decorators so endpoint values become canonical "EC50" / "LC50".
  #   "(log)EC50" → "EC50",  "LC50*" → "LC50",  "LC50/" → "LC50"
  mutate(endpoint = str_remove(endpoint, "^\\(log\\)"),
         endpoint = str_replace_all(endpoint, "[*/]", "")) %>%
  # Keep only acute lethal / effective-concentration endpoints.
  # These represent the most widely reported, standardised benchmarks in ECOTOX
  # and are the denominators in the Toxic Unit (TU) framework.
  filter(endpoint %in% c("LC50", "EC50")) %>%
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
# benchmark columns (I-6 to I-9), then write_fst() is called at step I-4.
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

fit_hc5_ssd <- function(values) {
  # values: numeric vector of conc_ng_L for one compound × group pair
  # Returns HC5 in ng/L, or NA_real_ if fitting fails or n < 5.
  # Warnings (e.g. optimiser non-convergence) are suppressed — they are
  # expected for extreme data distributions and handled by returning NA.
  withCallingHandlers(
    tryCatch({
      df <- data.frame(Conc = values[values > 0 & is.finite(values)])
      if (nrow(df) < 5) return(NA_real_)
      fit <- ssd_fit_dists(df, dists = "lnorm")
      hc  <- ssd_hc(fit, proportion = 0.05, ci = FALSE)
      as.numeric(hc$est[[1]])
    }, error = function(e) NA_real_),
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
  summarise(
    hc5_ssd_ng_L = fit_hc5_ssd(conc_ng_L),
    .groups = "drop"
  )

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
  message("S7b_ssd skipped — loading hc5_ssd_ng_L from existing taxotox_reference.fst...")
  if (file.exists("../Data/taxotox_reference.fst")) {
    prev <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
      select(cas_number, ecotox_group, any_of("hc5_ssd_ng_L")) %>%
      mutate(cas_number = as.character(cas_number))
    reference <- reference %>%
      mutate(cas_number = as.character(cas_number)) %>%
      left_join(prev, by = c("cas_number", "ecotox_group")) %>%
      mutate(
        hc5_model_ng_L = NA_real_,
        hc5_method     = case_when(!is.na(hc5_ssd_ng_L) ~ "SSD", TRUE ~ NA_character_)
      )
    message(sprintf("  Loaded hc5_ssd_ng_L for %d rows.", sum(!is.na(reference$hc5_ssd_ng_L))))
  } else {
    warning("No existing taxotox_reference.fst found — hc5_ssd_ng_L will be NA. Run S7b_ssd at least once.")
    reference <- reference %>%
      mutate(hc5_ssd_ng_L = NA_real_, hc5_model_ng_L = NA_real_, hc5_method = NA_character_)
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
  )

message("Step 7c: HC5/median scaling factors per taxonomic group:")
for (i in seq_len(nrow(scaling_factors))) {
  message(sprintf("  %-12s  k = %.4f  (calibrated from %d compounds)",
                  scaling_factors$ecotox_group[i],
                  scaling_factors$scaling_factor[i],
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
    )
  ) %>%
  select(-scaling_factor, -n_calibration)

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

  .ctox_get <- function(path, ..., max_tries = 3L) {
    url <- paste0(COMPTOX_BASE, path)
    for (i in seq_len(max_tries)) {
      resp <- tryCatch(
        GET(url, add_headers("x-api-key" = COMPTOX_KEY,
                             "accept" = "application/json"), ...),
        error = function(e) NULL
      )
      if (!is.null(resp) && status_code(resp) == 200L) return(resp)
      if (i < max_tries) Sys.sleep(i)
    }
    NULL
  }

  .ctox_post <- function(path, body_vec, ..., max_tries = 3L) {
    url  <- paste0(COMPTOX_BASE, path)
    body <- toJSON(body_vec, auto_unbox = FALSE)
    for (i in seq_len(max_tries)) {
      resp <- tryCatch(
        POST(url,
             add_headers("x-api-key"    = COMPTOX_KEY,
                         "accept"       = "application/json",
                         "content-type" = "application/json"),
             body = body, encode = "raw", ...),
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

  dtxsid_map <- purrr::map_dfr(seq_along(target_cas), function(i) {
    cas <- target_cas[i]
    if (i %% 100 == 0)
      message(sprintf("  %d / %d", i, length(target_cas)))

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

  all_dtxsids <- dtxsid_map$dtxsid[!is.na(dtxsid_map$dtxsid)]
  mw_table    <- .fetch_detail_mw(all_dtxsids)

  dtxsid_map  <- dtxsid_map %>%
    left_join(mw_table, by = "dtxsid")

  message(sprintf("  MW available    : %d / %d DTXSIDs",
                  sum(!is.na(dtxsid_map$mw_g_mol)), sum(!is.na(dtxsid_map$dtxsid))))

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
  .fetch_props <- function(dtxsids, endpoint, prop_names, batch_size = 1000L) {
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
      hc5_method = if_else(!is.na(hc5_model_ng_L), "Scaled", NA_character_)
    ) %>%
    select(-any_of(c("scaling_factor", "n_calibration")))

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
# Populates two columns in taxotox_reference.fst:
#   moa_group  — broad mechanistic group used to partition compounds in CAMA:
#                  "narcosis"          baseline toxicity (non-polar or polar narcosis)
#                  "AChE_inhibition"   acetylcholinesterase inhibitors
#                                      (organophosphates, carbamates)
#                  "PSII_inhibition"   Photosystem II inhibitors
#                                      (triazines, phenylureas, diazinones)
#                  "pyrethroid"        sodium channel modulators
#                  "reactive"          unspecific reactive chemicals
#   moa_source — how the group was assigned:
#                  "name_heuristic"    matched a chemical name pattern
#                  "Verhaar_LogKow"    Verhaar classification via LogKow
#                  "default_narcosis"  no information — conservative default
#
# Classification strategy (applied in priority order):
#   1. Chemical name pattern matching against preferred_name from CompTox
#      and chemical_name from ECOTOX — catches major pesticide classes reliably
#   2. LogKow-based Verhaar classification (Verhaar et al. 1992) for residuals
#      Polar narcotics (LogKow < 2) and non-polar narcotics (LogKow ≥ 2)
#      both map to "narcosis"; reactive class not classified without SMARTS
#   3. Default: "narcosis" (most common aquatic MoA, most conservative for CAMA)
#
# LogKow is fetched from the CompTox predicted property endpoint
# (propName = "LogKow: Octanol-Water"). The dtxsid_map is reloaded from
# Data/dtxsid_map.csv written at the end of I-2.
#
# Reference: Verhaar, H.J.M., van Leeuwen, C.J., Hermens, J.L.M. (1992).
#   Classifying environmental pollutants. Chemosphere, 25(4), 471-491.
# =============================================================================

if (!nzchar(tryCatch(Sys.getenv("COMPTOX_API_KEY", unset = ""),
                     error = function(e) ""))) {
  message("Step I-3: SKIPPED — CompTox API key not found (needed for LogKow fetch).")
} else if (!file.exists("../Data/taxotox_reference.fst")) {
  message("Step I-3: SKIPPED — taxotox_reference.fst not found. Run I-2 first.")
} else if (!.ask_run("S_i3_moa",
    "Mode of Action classification — LogKow fetch + name heuristics (~2\u20135 min)")) {
  message("S_i3_moa skipped — existing moa_group values (if any) preserved.")
} else {

  library(httr)
  library(jsonlite)

  # Reuse helpers and key from I-2 scope if available; otherwise rebuild
  if (!exists("COMPTOX_KEY") || !nzchar(COMPTOX_KEY)) {
    COMPTOX_KEY <- Sys.getenv("COMPTOX_API_KEY", unset = "")
    if (!nzchar(COMPTOX_KEY) && file.exists("comptox_api_key.txt"))
      COMPTOX_KEY <- trimws(readLines("comptox_api_key.txt", n = 1L))
  }
  if (!exists("COMPTOX_BASE")) COMPTOX_BASE <- "https://comptox.epa.gov/ctx-api"
  if (!exists(".ctox_get")) {
    .ctox_get <- function(path, ..., max_tries = 3L) {
      url <- paste0(COMPTOX_BASE, path)
      for (i in seq_len(max_tries)) {
        resp <- tryCatch(
          GET(url, add_headers("x-api-key" = COMPTOX_KEY, "accept" = "application/json"), ...),
          error = function(e) NULL
        )
        if (!is.null(resp) && status_code(resp) == 200L) return(resp)
        if (i < max_tries) Sys.sleep(i)
      }
      NULL
    }
  }
  if (!exists(".ctox_post")) {
    .ctox_post <- function(path, body_vec, ..., max_tries = 3L) {
      url  <- paste0(COMPTOX_BASE, path)
      body <- toJSON(body_vec, auto_unbox = FALSE)
      for (i in seq_len(max_tries)) {
        resp <- tryCatch(
          POST(url,
               add_headers("x-api-key"    = COMPTOX_KEY,
                           "accept"       = "application/json",
                           "content-type" = "application/json"),
               body = body, encode = "raw", ...),
          error = function(e) NULL
        )
        if (!is.null(resp) && status_code(resp) == 200L) return(resp)
        if (i < max_tries) Sys.sleep(i)
      }
      NULL
    }
  }
  if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

  # --- I-3a: Load reference and dtxsid_map -----------------------------------

  reference <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
    mutate(cas_number = as.character(cas_number))

  if (file.exists("../Data/dtxsid_map.csv")) {
    dtxsid_map_i3 <- read.csv("../Data/dtxsid_map.csv", stringsAsFactors = FALSE)
    message(sprintf("  Loaded dtxsid_map: %d entries (%d with DTXSID)",
                    nrow(dtxsid_map_i3), sum(!is.na(dtxsid_map_i3$dtxsid))))
  } else {
    # dtxsid_map.csv is written at the end of I-2 but may be absent if I-2 ran
    # before that save was added.  Re-derive DTXSID from Known_CAS inline.
    message("dtxsid_map.csv not found — running CAS\u2192DTXSID lookup inline...")
    known_cas_i3 <- read.fst("../Data/Known_CAS.fst", as.data.table = FALSE)
    cas_col_i3   <- intersect(c("cas_number", "casrn", "CAS", "CASRN"), names(known_cas_i3))[1]
    target_cas_i3 <- unique(as.character(known_cas_i3[[cas_col_i3]]))
    message(sprintf("  %d unique CAS numbers to resolve", length(target_cas_i3)))

    dtxsid_map_i3 <- purrr::map_dfr(seq_along(target_cas_i3), function(i) {
      cas <- target_cas_i3[i]
      if (i %% 100 == 0) message(sprintf("  %d / %d", i, length(target_cas_i3)))
      resp <- .ctox_get(paste0("/chemical/search/equal/",
                               utils::URLencode(cas, reserved = TRUE)))
      if (is.null(resp))
        return(tibble(cas_number = cas, dtxsid = NA_character_, preferred_name = NA_character_))
      parsed <- tryCatch(
        fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE),
        error = function(e) NULL
      )
      if (is.null(parsed) || nrow(parsed) == 0)
        return(tibble(cas_number = cas, dtxsid = NA_character_, preferred_name = NA_character_))
      hit <- parsed[!is.na(parsed$casrn) & parsed$casrn == cas, ]
      if (nrow(hit) == 0) hit <- parsed[1, ]
      tibble(cas_number     = cas,
             dtxsid         = as.character(hit$dtxsid[1]        %||% NA),
             preferred_name = as.character(hit$preferredName[1] %||% NA))
    })
    write.csv(dtxsid_map_i3, "../Data/dtxsid_map.csv", row.names = FALSE)
    message(sprintf("  DTXSID resolved: %d / %d — saved to dtxsid_map.csv",
                    sum(!is.na(dtxsid_map_i3$dtxsid)), nrow(dtxsid_map_i3)))
  }

  all_dtxsids_i3 <- dtxsid_map_i3$dtxsid[!is.na(dtxsid_map_i3$dtxsid)]

  # --- I-3b: Fetch LogKow from CompTox predicted property endpoint -----------

  PROP_LOGKOW <- "LogKow: Octanol-Water"

  message("Step I-3b: Fetching LogKow from CompTox predicted properties...")

  .fetch_logkow <- function(dtxsids, batch_size = 1000L) {
    batches <- split(dtxsids, ceiling(seq_along(dtxsids) / batch_size))
    purrr::map_dfr(seq_along(batches), function(i) {
      batch <- batches[[i]]
      message(sprintf("  LogKow batch %d / %d (%d DTXSIDs)", i, length(batches), length(batch)))
      resp <- .ctox_post("/chemical/property/predicted/search/by-dtxsid/", batch)
      if (is.null(resp)) { warning(sprintf("  LogKow batch %d failed.", i)); return(tibble()) }
      parsed <- tryCatch(
        fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE),
        error = function(e) NULL
      )
      if (is.null(parsed) || nrow(parsed) == 0) return(tibble())
      rows <- parsed[parsed$propName == PROP_LOGKOW, c("dtxsid", "propValue")]
      if (nrow(rows) == 0) return(tibble())
      # Take first value per DTXSID (multiple model sources may be returned)
      tibble(dtxsid = rows$dtxsid, logkow = as.numeric(rows$propValue)) %>%
        group_by(dtxsid) %>% slice(1) %>% ungroup()
    })
  }

  if (length(all_dtxsids_i3) == 0) {
    logkow_table <- tibble(dtxsid = character(), logkow = numeric())
    message("  No DTXSIDs available — LogKow table empty.")
  } else {
    logkow_table <- .fetch_logkow(all_dtxsids_i3)
    message(sprintf("  LogKow retrieved for %d / %d unique DTXSIDs",
                    n_distinct(logkow_table$dtxsid), length(all_dtxsids_i3)))
  }

  # Merge LogKow back to dtxsid_map, then to cas_number.
  # Deduplicate by cas_number (keep first): dtxsid_map_i3 may have multiple
  # CAS numbers mapping to the same DTXSID, creating many-to-many with logkow_table.
  cas_logkow <- dtxsid_map_i3 %>%
    left_join(logkow_table, by = "dtxsid") %>%
    select(cas_number, preferred_name, any_of("logkow")) %>%
    { if (!"logkow" %in% names(.)) mutate(., logkow = NA_real_) else . } %>%
    group_by(cas_number) %>% slice(1) %>% ungroup()

  # --- I-3c: Chemical name pattern matching ----------------------------------
  # Matches preferred_name (CompTox) and chemical_name (ECOTOX) against
  # known chemical class patterns. Applied before LogKow-based classification.
  # Patterns are case-insensitive and target the most common aquatic toxicants.

  .classify_by_name <- function(name) {
    if (is.na(name) || !nzchar(name)) return(NA_character_)
    n <- tolower(name)

    # Acetylcholinesterase (AChE) inhibitors
    # -- Organophosphates: phosphate/phosphorothioate/phosphonothioate esters
    if (grepl("chlorpyrifos|malathion|diazinon|parathion|azinphos|dimethoate|
               phosmet|fenthion|methidathion|profenofos|triazophos|phorate|
               disulfoton|terbufos|acephate|omethoate|methamidophos|
               dichlorvos|trichlorfon|naled|ethoprop|cadusafos|
               pirimiphos|fenitrothion|temephos|isazofos|quinalphos|
               monocrotophos|phosphamidon|vamidothion|demeton",
              n, perl = TRUE)) return("AChE_inhibition")
    # -- Carbamates
    if (grepl("carbofuran|carbaryl|methomyl|aldicarb|oxamyl|pirimicarb|
               thiodicarb|fenobucarb|propoxur|bendiocarb|aminocarb|
               formetanate|alanycarb|furathiocarb",
              n, perl = TRUE)) return("AChE_inhibition")

    # Photosystem II (PSII) inhibitors
    # -- Triazines
    if (grepl("atrazine|simazine|cyanazine|terbuthylazine|prometryn|
               ametryn|terbutryn|propazine|desmetryn|atraton|
               hexazinone|metribuzin",
              n, perl = TRUE)) return("PSII_inhibition")
    # -- Phenylureas
    if (grepl("diuron|isoproturon|linuron|chlortoluron|metobromuron|
               monolinuron|metoxuron|fluometuron|tebuthiuron|
               chlorbromuron|buturon",
              n, perl = TRUE)) return("PSII_inhibition")
    # -- Other PSII inhibitors
    if (grepl("bromacil|terbacil|lenacil|bentazon|ioxynil|bromoxynil|
               metribuzin|amicarbazone|tembotrione",
              n, perl = TRUE)) return("PSII_inhibition")

    # Pyrethroids (sodium channel modulators)
    if (grepl("permethrin|cypermethrin|deltamethrin|cyhalothrin|
               bifenthrin|fenvalerate|esfenvalerate|cyfluthrin|
               flucythrinate|fenpropathrin|tralomethrin|
               tefluthrin|acrinathrin|etofenprox|tau-fluvalinate|
               resmethrin|allethrin|tetramethrin|prallethrin|
               pyrethrin|pyrethroid|pyrethr",
              n, perl = TRUE)) return("pyrethroid")

    # Reactive chemicals (unspecific electrophilic reactivity)
    # Epoxides, Michael acceptors, alpha-beta unsaturated carbonyls
    if (grepl("epoxide|acrylate|acrolein|formaldehyde|
               chloroacetamide|methyl vinyl|isocyanate|
               2-chloroacet|alpha-chloro",
              n, perl = TRUE)) return("reactive")

    NA_character_  # no pattern matched
  }

  # --- I-3d: Combine name patterns and LogKow into moa assignments ----------

  # Build a per-cas_number MoA table (one entry per compound, not per group)
  moa_by_cas <- reference %>%
    select(cas_number, chemical_name) %>%
    distinct() %>%
    left_join(cas_logkow %>% select(cas_number, preferred_name, logkow),
              by = "cas_number") %>%
    rowwise() %>%
    mutate(
      # Try name patterns on both ECOTOX name and CompTox preferred name
      moa_from_ecotox_name  = .classify_by_name(chemical_name),
      moa_from_comptox_name = .classify_by_name(preferred_name),

      # Priority: CompTox preferred name > ECOTOX chemical name
      moa_from_name = dplyr::coalesce(moa_from_comptox_name, moa_from_ecotox_name),

      # LogKow-based Verhaar classification for residuals
      # Verhaar class 1 (LogKow ≥ 2): non-polar narcosis
      # Verhaar class 2 (LogKow < 2):  polar narcosis
      # Both map to "narcosis" in CAMA (same mechanism, different potency)
      moa_from_logkow = case_when(
        !is.na(logkow) ~ "narcosis",
        TRUE            ~ NA_character_
      ),

      # Final assignment (priority: name > LogKow > default)
      moa_group = dplyr::coalesce(moa_from_name, moa_from_logkow, "narcosis"),

      moa_source = case_when(
        !is.na(moa_from_name)   ~ "name_heuristic",
        !is.na(moa_from_logkow) ~ "Verhaar_LogKow",
        TRUE                    ~ "default_narcosis"
      )
    ) %>%
    ungroup() %>%
    select(cas_number, moa_group, moa_source) %>%
    # Deduplicate: one row per cas_number (guards against duplicate cas_logkow rows)
    group_by(cas_number) %>% slice(1) %>% ungroup()

  # Summarise coverage
  message("Step I-3d: Mode of Action classification summary:")
  moa_by_cas %>%
    count(moa_group, moa_source) %>%
    arrange(desc(n)) %>%
    { message(paste(capture.output(print(as.data.frame(.))), collapse = "\n")) }

  # --- I-3e: Join MoA into reference and write fst --------------------------

  # moa_group and moa_source apply per compound (same across all ecotox_groups).
  # Drop any existing moa columns before joining (handles re-runs cleanly).
  reference <- reference %>%
    select(-any_of(c("moa_group", "moa_source"))) %>%
    left_join(moa_by_cas, by = "cas_number")

  message(sprintf(
    "Step I-3 complete: moa_group populated for %d / %d reference rows.",
    sum(!is.na(reference$moa_group)), nrow(reference)
  ))

  write_fst(reference, "../Data/taxotox_reference.fst", compress = 50)
  message(sprintf("taxotox_reference.fst updated (post-MoA): %d rows, %d columns.",
                  nrow(reference), ncol(reference)))

  .mark_done("S_i3_moa")

}  # end I-3

# =============================================================================
# Step I-4: Final write of taxotox_reference.fst
# -----------------------------------------------------------------------------
# This step consolidates everything into the definitive reference table.
# It loads the current taxotox_reference.fst (written by whichever previous
# sections ran), enforces canonical column order, fills any missing columns
# with NA, and writes the final file.
#
# Running this step is safe at any point — it never overwrites data, only
# reorders and canonicalises columns. Re-running after adding benchmark
# columns (steps I-6 to I-11) will produce a file with those columns too.
#
# Canonical column order:
#   Identity    : cas_number, ecotox_group, chemical_name
#   ECOTOX data : n_ecotox, median_lc50_ng_L
#   HC5         : hc5_ssd_ng_L, hc5_model_ng_L, hc5_method
#   CompTox     : predicted_lc50_ng_L, lc50_source
#   Mode of action: moa_group, moa_source
#   Benchmarks  : benchmark_usepa_fish_acute_ng_L, _crust_acute_ng_L,
#                 _algae_acute_ng_L, benchmark_eu_eqs_aa_marine_ng_L,
#                 benchmark_au_anzg_fresh_ng_L, _marine_ng_L,
#                 benchmark_ca_ccme_fresh_lt_ng_L
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
    # CompTox gap-fill
    "predicted_lc50_ng_L", "lc50_source",
    # Mode of Action
    "moa_group", "moa_source",
    # Benchmarks (I-6 to I-9 — NA until those steps run)
    "benchmark_usepa_fish_acute_ng_L",
    "benchmark_usepa_crust_acute_ng_L",
    "benchmark_usepa_algae_acute_ng_L",
    "benchmark_eu_eqs_aa_marine_ng_L",
    "benchmark_au_anzg_fresh_ng_L",
    "benchmark_au_anzg_marine_ng_L",
    "benchmark_ca_ccme_fresh_lt_ng_L"
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

#######################################################################


