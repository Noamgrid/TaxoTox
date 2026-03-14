setwd("c:/Users/Yair/Documents/yair/research/estuary_rehabilitation/yairs_stuff/Students/Noam/TaxoTox/Code")

# title: "TaxoToxInstall_AFB"
# author: "Yair & Noam"
# purpose: Build final_ecotox_data_AFB.fst using EPA Aquatic Life Benchmarks
#          (AFB/AIB/ANVPB) as denominators instead of ECOTOX medians.
#          Use this file for validation against Covert et al. (2020) PTI values.
#          The main TaxoTox pipeline (taxotox_install.R / final_ecotox_data.fst)
#          is NOT modified.

library(dplyr)
library(tidyr)
library(stringr)
library(fst)

# ── 1. Load the EPA benchmark table ──────────────────────────────────────────
# Source file: Testing_Taxotox/aquatic_benchmarks.csv
# Column layout (duplicate headers renamed by R with .1 suffix):
#   Col 1  "Pesticide"    — chemical name
#   Col 3  "CAS.number"   — CAS Registry Number
#   Col 4  "AcuteA"       — Freshwater Fish Acute (AFB)       µg/L
#   Col 6  "AcuteA.1"     — Freshwater Invertebrate Acute (AIB) µg/L
#   Col 8  "AcuteC"       — Freshwater Nonvascular Plant Acute (ANVPB) µg/L

benchmark_path <- "../../Testing_Taxotox/aquatic_benchmarks.csv"

if (!file.exists(benchmark_path))
  stop("Benchmark file not found: ", normalizePath(benchmark_path, mustWork = FALSE))

raw <- read.csv(benchmark_path, stringsAsFactors = FALSE, check.names = TRUE)

message("Loaded ", nrow(raw), " rows.")
message("Columns: ", paste(names(raw), collapse = ", "))

# Verify expected columns exist
needed <- c("Pesticide", "CAS.number", "AcuteA", "AcuteA.1", "AcuteC")
missing <- setdiff(needed, names(raw))
if (length(missing) > 0)
  stop("Unexpected column names. Missing: ", paste(missing, collapse = ", "),
       "\nActual columns: ", paste(names(raw), collapse = ", "))

# ── 2. Extract and clean benchmark values ─────────────────────────────────────

clean_num <- function(x) {
  # Remove footnote letters, spaces, "<" signs; keep digits and decimal point
  suppressWarnings(as.numeric(str_remove_all(as.character(x), "[^0-9.]")))
}

benchmarks <- raw %>%
  transmute(
    chemical_name = str_trim(Pesticide),
    cas_number    = str_remove_all(str_trim(CAS.number), "-"),  # strip dashes
    AFB           = clean_num(AcuteA),    # Freshwater Fish Acute   (µg/L)
    AIB           = clean_num(AcuteA.1),  # Freshwater Invert Acute (µg/L)
    ANVPB         = clean_num(AcuteC)     # Freshwater Plant Acute  (µg/L)
  ) %>%
  filter(
    cas_number != "",
    !is.na(cas_number),
    !is.na(AFB) | !is.na(AIB) | !is.na(ANVPB)
  ) %>%
  distinct(cas_number, .keep_all = TRUE)

message("\nRetained ", nrow(benchmarks), " unique chemicals with at least one benchmark.")
message("  Fish (AFB)             : ", sum(!is.na(benchmarks$AFB)))
message("  Invertebrate/AIB       : ", sum(!is.na(benchmarks$AIB)))
message("  Nonvascular Plant/ANVPB: ", sum(!is.na(benchmarks$ANVPB)))

# ── 4. Pivot to long format matching final_ecotox_data.fst schema ─────────────
# EPA benchmarks are in µg/L → multiply by 1000 to convert to ng/L.
# Each CAS gets one row per group. app.R takes median(conc_ng_L) per group,
# so a single row yields median = that value directly.

afb_long <- benchmarks %>%
  pivot_longer(
    cols      = c(AFB, AIB, ANVPB),
    names_to  = "benchmark_type",
    values_to = "value_ug_L"
  ) %>%
  filter(!is.na(value_ug_L), value_ug_L > 0) %>%
  mutate(
    ecotox_group = case_when(
      benchmark_type == "AFB"   ~ "fish",
      benchmark_type == "AIB"   ~ "crustacean",
      benchmark_type == "ANVPB" ~ "algae"
    ),
    conc_ng_L         = value_ug_L * 1000,   # µg/L → ng/L
    min_concentration = value_ug_L,           # µg/L (consistent with conc1_unit below)
    conc1_unit        = "ug/L",
    factor_to_ng_L    = 1000,
    endpoint          = "EC50",
    # Stub columns to satisfy read.fst and any downstream column checks
    result_id              = NA_integer_,
    test_id                = NA_integer_,
    reference_number       = NA_integer_,
    trend                  = NA_character_,
    effect                 = NA_character_,
    exposure_type          = NA_character_,
    measurement            = NA_character_,
    measurement_comments   = NA_character_,
    species                = NA_character_,
    test_cas               = cas_number
  ) %>%
  select(
    result_id, min_concentration, endpoint, cas_number, trend, effect,
    exposure_type, measurement, measurement_comments, test_id,
    reference_number, conc1_unit, test_cas, species, ecotox_group,
    chemical_name, conc_ng_L, factor_to_ng_L
  )

message("\nFinal AFB dataset: ", nrow(afb_long), " rows")
message(table(afb_long$ecotox_group))

# ── 5. Write output ────────────────────────────────────────────────────────────
# Back up the existing ECOTOX-based file (once), then overwrite with AFB version.
# To restore: rename final_ecotox_data_ECOTOX_backup.fst back to final_ecotox_data.fst

out_path    <- "../Data/final_ecotox_data.fst"
backup_path <- "../Data/final_ecotox_data_ECOTOX_backup.fst"

if (file.exists(out_path) && !file.exists(backup_path)) {
  file.copy(out_path, backup_path)
  message("Backed up existing file to: ", backup_path)
} else if (file.exists(backup_path)) {
  message("Backup already exists (", backup_path, ") — skipping backup.")
}

write_fst(afb_long, out_path, compress = 50)
message("\nSaved AFB benchmarks as: ", out_path)
message("App.R will now use EPA AFB benchmarks — no path change needed.")
message("\nTo restore original ECOTOX data, run:")
message("  file.copy('../Data/final_ecotox_data_ECOTOX_backup.fst', '../Data/final_ecotox_data.fst', overwrite = TRUE)")
