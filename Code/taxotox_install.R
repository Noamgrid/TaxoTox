setwd("c:/Users/Yair/Documents/yair/research/estuary_rehabilitation/yairs_stuff/Students/Noam/TaxoTox/Code")

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

#######################################################################
# Step 1: Connect to the local ECOTOX SQLite database and extract data.
# The database is built by update_ecotox.R (ECOTOXr::download_ecotox_data()).
# To find the correct path on your machine:  ECOTOXr::get_ecotox_sqlite_file()

conn <- dbConnect(RSQLite::SQLite(), "C:/Users/Yair/AppData/Local/ECOTOXr/ECOTOXr/Cache/ecotox_ascii_03_13_2025.sqlite")

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

#######################################################################


