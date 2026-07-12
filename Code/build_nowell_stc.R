# One-off: compute Nowell et al. (2014) "Taxon-Sensitive PTI" STC columns and
# merge them into the existing taxotox_reference.fst, without re-running the
# full taxotox_install.R pipeline (which would also re-trigger the slow SSD
# fit and CompTox API steps under non-interactive Rscript execution, since
# .ask_run() always returns TRUE when !interactive()).
#
# Both fish and cladoceran STC are queried fresh here (rather than reusing
# final_ecotox_data.fst) because that file predates the genus/family and
# exposure-duration fields this script needs. See taxotox_install.R for the
# full pipeline; this script keeps `.nowell_stc()`/`.duration_ok()` and the
# conversion constants in sync with it by design, and should be retired once
# final_ecotox_data.fst is rebuilt with those fields (at which point
# taxotox_install.R alone is sufficient).

library(dplyr)
library(stringr)
library(fst)
library(RSQLite)
library(ECOTOXr)
library(bit64)  # required for correct integer64 -> character CAS-number conversion

stopifnot("Run this script from the Code/ directory (paths are relative to it)" =
  file.exists("../Data/taxotox_reference.fst"))

NOWELL_CLADOCERAN_GENERA <- c(
  "Acantholeberis", "Acroperus", "Alona", "Alonella", "Bosmina",
  "Ceriodaphnia", "Chydorus", "Daphnia", "Diaphanosoma", "Disparalona",
  "Eurycercus", "Moina", "Moinodaphnia", "Pleuroxus", "Pseudosida",
  "Scapholeberis", "Simocephalus"
)

.nowell_stc <- function(values) {
  values <- values[!is.na(values) & values > 0]
  n <- length(values)
  if (n == 0) return(NA_real_)
  if (n > 12) as.numeric(quantile(values, probs = 0.05, names = FALSE))
  else min(values)
}

# Nowell et al. (2014) main paper Table 1 ("Standard criteria for toxicity
# test duration, endpoint, and measured effect, by taxonomic group") specifies
# 96-hour LC50 for fish and 48-hour EC50/LC50 for cladocerans as the standard.
# Most ECOTOX records don't have exposure duration coded at all (75.5% "NC"
# for fish LC50, 44% for cladoceran IMBL-EC50 in this database) -- excluding
# those outright would discard the majority of usable data for a metadata gap,
# not because the test was actually non-standard. So: keep records with
# unrecorded duration ("NC"/NA) or matching the standard; exclude only records
# where a NON-standard duration is explicitly logged.
.duration_ok <- function(unit, mean_val, standard_hours) {
  is.na(unit) | unit %in% c("NC", "") | (unit == "h" & mean_val == standard_hours)
}

# NOTE: a freshwater-only filter (media_type "FW" vs "SW") was tried and
# reverted -- Nowell et al.'s title says "...to Freshwater Aquatic
# Organisms", but empirically restricting to media_type=="FW" made the fish
# fit WORSE (R² 0.813 -> 0.738) rather than better, unlike every other fix in
# this file. Likely explanation: excluding saltwater records pushed some
# compounds below the n>12 percentile threshold, forcing a noisier
# minimum-based STC. Not reapplying without further investigation.

relevent_units <- c("ng/L", "ug/L", "mg/L", "g/L")
conversion_df <- tibble::tibble(
  conc1_unit     = relevent_units,
  factor_to_ng_L = c(1, 1e3, 1e6, 1e9)
)

.convert_units <- function(raw) {
  raw %>%
    mutate(
      conc1_unit = str_remove(conc1_unit, "^AI\\s*"),
      conc1_unit = str_replace_all(conc1_unit, "\\bdm3\\b", "L"),
      conc1_unit = str_replace_all(conc1_unit, "ppt",  "ng/L"),
      conc1_unit = str_replace_all(conc1_unit, "ppm",  "mg/L"),
      conc1_unit = str_replace_all(conc1_unit, "ppb",  "ug/L"),
      conc1_unit = str_replace_all(conc1_unit, "0/00", "g/L")
    ) %>%
    filter(conc1_unit %in% relevent_units) %>%
    left_join(conversion_df, by = "conc1_unit") %>%
    mutate(min_concentration = as.numeric(min_concentration)) %>%
    mutate(min_concentration = if_else(
      str_detect(endpoint, "^\\(log\\)"),
      10 ^ min_concentration,
      min_concentration
    )) %>%
    mutate(conc_ng_L = min_concentration * factor_to_ng_L) %>%
    mutate(endpoint = str_remove(endpoint, "^\\(log\\)"),
           endpoint = str_replace_all(endpoint, "[*/]", ""))
}

.calc_conc1_mean_sql <- "
        CASE
            WHEN (r.conc1_mean IS NULL OR r.conc1_mean <= 0)
                 AND r.conc1_min > 0 AND r.conc1_max > 0
            THEN (r.conc1_min + r.conc1_max) / 2.0
            ELSE r.conc1_mean
        END AS calc_conc1_mean"

.min_concentration_sql <- "
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
    END AS min_concentration"

conn <- dbConnect(RSQLite::SQLite(), ECOTOXr::get_ecotox_sqlite_file())

# --- Fish STC: fresh query, LC50 only, 96h-or-unrecorded duration -----------
message("Querying ECOTOX for fish LC50 data...")
fish_query <- paste0("
WITH enriched_results AS (
    SELECT
        r.result_id, r.endpoint, r.conc1_unit, r.conc1_mean, r.conc2_mean,
        r.conc3_mean, r.conc1_min, r.conc1_max,
        t.exposure_duration_mean, t.exposure_duration_unit,
        c.cas_number,", .calc_conc1_mean_sql, "
    FROM results r
    JOIN tests t ON r.test_id = t.test_id
    JOIN species s ON t.species_number = s.species_number
    JOIN chemicals c ON t.test_cas = c.cas_number
    WHERE LOWER(s.ecotox_group) LIKE '%fish%'
      AND r.endpoint = 'LC50'
)
SELECT
    er.result_id, er.endpoint, er.conc1_unit, er.cas_number,
    er.exposure_duration_mean, er.exposure_duration_unit,", .min_concentration_sql, "
FROM enriched_results er;
")
fish_raw <- dbGetQuery(conn, fish_query)
message(sprintf("  %d raw fish LC50 bioassay rows fetched", nrow(fish_raw)))

fish_data <- .convert_units(fish_raw) %>%
  filter(endpoint == "LC50")
message(sprintf("  %d fish rows after unit conversion (before duration fix)", nrow(fish_data)))

fish_data <- fish_data %>%
  filter(.duration_ok(exposure_duration_unit, exposure_duration_mean, 96))
message(sprintf("  %d fish rows after 96h-or-unrecorded duration fix", nrow(fish_data)))

nowell_fish <- fish_data %>%
  mutate(cas_number = as.character(cas_number)) %>%
  group_by(cas_number) %>%
  summarise(
    ecotox_group    = "fish",
    stc_nowell_ng_L = .nowell_stc(conc_ng_L),
    n_stc_nowell    = sum(!is.na(conc_ng_L) & conc_ng_L > 0),
    .groups = "drop"
  )
message(sprintf("  Fish STC computed for %d compounds", nrow(nowell_fish)))

# --- Cladoceran STC: fresh, genus-restricted ECOTOX query -------------------
message("\nQuerying ECOTOX for cladoceran-genus-restricted crustacean data...")
genus_list_sql <- paste0("'", NOWELL_CLADOCERAN_GENERA, "'", collapse = ", ")

clad_query <- paste0("
WITH enriched_results AS (
    SELECT
        r.result_id, r.endpoint, r.conc1_unit, r.conc1_mean, r.conc2_mean,
        r.conc3_mean, r.conc1_min, r.conc1_max, r.measurement,
        t.exposure_duration_mean, t.exposure_duration_unit,
        c.cas_number,", .calc_conc1_mean_sql, "
    FROM results r
    JOIN tests t ON r.test_id = t.test_id
    JOIN species s ON t.species_number = s.species_number
    JOIN chemicals c ON t.test_cas = c.cas_number
    WHERE LOWER(s.ecotox_group) LIKE '%crustacean%'
      AND s.genus IN (", genus_list_sql, ")
)
SELECT
    er.result_id, er.endpoint, er.conc1_unit, er.cas_number, er.measurement,
    er.exposure_duration_mean, er.exposure_duration_unit,", .min_concentration_sql, "
FROM enriched_results er;
")

clad_raw <- dbGetQuery(conn, clad_query)
dbDisconnect(conn)
message(sprintf("  %d raw cladoceran bioassay rows fetched", nrow(clad_raw)))

clad_data <- .convert_units(clad_raw) %>%
  filter(endpoint %in% c("LC50", "EC50"))
message(sprintf("  %d cladoceran rows after unit conversion/endpoint filter (before measurement/duration fix)", nrow(clad_data)))

# Nowell et al. (2014) AppA.pdf sec. 2 explicitly restricts cladoceran EC50 data
# to the immobilization endpoint ("measurement" == "IMBL") -- LC50 is kept as-is
# since mortality is unambiguous. Without this, chronic/sub-lethal EC50 records
# (e.g. reproduction "PROG"/"FCND", population "GPOP") ~1000x more sensitive than
# acute immobilization leak into the STC. Confirmed empirically for Diazinon:
# stc_nowell_ng_L went from 0.323 ng/L (n=112, contaminated) to 239 ng/L (n=80,
# IMBL-only) after this fix -- a ~740x correction toward a plausible acute value.
clad_data <- clad_data %>%
  filter(endpoint == "LC50" | (endpoint == "EC50" & measurement == "IMBL"))
message(sprintf("  %d cladoceran rows after measurement=='IMBL' fix for EC50", nrow(clad_data)))

clad_data <- clad_data %>%
  filter(.duration_ok(exposure_duration_unit, exposure_duration_mean, 48))
message(sprintf("  %d cladoceran rows after 48h-or-unrecorded duration fix", nrow(clad_data)))

nowell_cladoceran <- clad_data %>%
  mutate(cas_number = as.character(cas_number)) %>%
  group_by(cas_number) %>%
  summarise(
    ecotox_group    = "crustacean",
    stc_nowell_ng_L = .nowell_stc(conc_ng_L),
    n_stc_nowell    = sum(!is.na(conc_ng_L) & conc_ng_L > 0),
    .groups = "drop"
  )
message(sprintf("  Cladoceran STC computed for %d compounds", nrow(nowell_cladoceran)))

# --- Merge into taxotox_reference.fst ---------------------------------------
nowell_stc_table <- bind_rows(nowell_fish, nowell_cladoceran)

file.copy("../Data/taxotox_reference.fst", "../Data/taxotox_reference.fst.bak_pre_freshwater_fix", overwrite = TRUE)
message("\nBacked up taxotox_reference.fst -> taxotox_reference.fst.bak_pre_freshwater_fix")

ref <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
  select(-any_of(c("stc_nowell_ng_L", "n_stc_nowell"))) %>%  # safe re-run
  left_join(nowell_stc_table, by = c("cas_number", "ecotox_group"))

write_fst(ref, "../Data/taxotox_reference.fst", compress = 50)
message(sprintf("\ntaxotox_reference.fst updated: %d rows, %d columns.", nrow(ref), ncol(ref)))
message(sprintf("  stc_nowell_ng_L populated: %d rows", sum(!is.na(ref$stc_nowell_ng_L))))
