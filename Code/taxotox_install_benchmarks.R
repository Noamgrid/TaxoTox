# Set working directory to this script's location.
.script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) tryCatch(
    dirname(normalizePath(sys.frames()[[1]]$ofile)),
    error = function(e) {
      f <- sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])
      if (length(f) && nzchar(f)) dirname(normalizePath(f)) else getwd()
    }
  )
)
setwd(.script_dir)

# =============================================================================
# taxotox_install_benchmarks.R  (steps I-5 through I-11)
# -----------------------------------------------------------------------------
# Purpose : Parse national aquatic benchmark tables and join them into
#           taxotox_reference.fst as denominator columns for the
#           Benchmark Hazard Index method.
#
# Sources :
#   I-6  Data/USEPA_aquatic_benchmarks.csv   — US EPA Aquatic Life Benchmarks
#   I-7  Data/EUEPA_aquatic_benchmarks.csv   — EU WFD Directive 2013/39/EU
#   I-8  Data/Australia_aquatic_benchmarks.csv — ANZG 2018 / ANZECC 2000
#   I-9  Data/Canada_aquatic_benchmarks.csv  — CCME Water Quality Guidelines
#
# All source values are in µg/L; output stored as ng/L (× 1000).
#
# Output  : Updates ../Data/taxotox_reference.fst in place with 7 new columns:
#   benchmark_usepa_fish_acute_ng_L
#   benchmark_usepa_crust_acute_ng_L
#   benchmark_usepa_algae_acute_ng_L
#   benchmark_eu_eqs_aa_marine_ng_L
#   benchmark_au_anzg_fresh_ng_L
#   benchmark_au_anzg_marine_ng_L
#   benchmark_ca_ccme_fresh_lt_ng_L
#
# Denominator selection rationale (documented in TaxoTox_Technical_Methods.md):
#   US EPA  : acute columns only (fish = raw col 4, invertebrate = raw col 8,
#             algae = raw col 12 IC50 — see column-mapping note at I-6 below;
#             verified against EPA's live benchmarks page, NOT the same as
#             the positions the duplicate-header auto-rename suggests)
#             Chronic columns are retained in source but not used as TU denominators
#   EU EQS  : Annual Average EQS — the only value provided in this CSV extract
#             Reflects estuarine/coastal monitoring (primary TaxoTox use case)
#   AU ANZG : LOSP 95% (95% species protection level) — recommended default
#   CA CCME : Freshwater Long-Term (chronic) guideline
#
# Prerequisites : taxotox_reference.fst must exist (run taxotox_install.R first)
# Authors : Yair Suari & Noam Gridish, 2025
# =============================================================================

library(dplyr)
library(fst)

UG_TO_NG <- 1000  # all source values in µg/L; convert to ng/L

if (!file.exists("../Data/taxotox_reference.fst")) {
  stop("taxotox_reference.fst not found — run taxotox_install.R (through I-4) first.")
}

# =============================================================================
# Helper: parse a character column that may contain:
#   - plain numbers: "0.3", "12"
#   - scientific notation with × : "1.7 × 10–4", "8 × 10–6"
#   - sum notation: "Σ = 0.005"
#   - qualifiers: "≤ 0.08 (Class 1)"
#   - non-numeric: "No data", "NRG", "Insufficient data", "see footnote",
#                  "not applicable", blank
# Returns a numeric vector (NA for non-numeric inputs).
# =============================================================================
.parse_benchmark <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_

  sapply(x, function(v) {
    if (is.na(v)) return(NA_real_)

    # Scientific notation with Unicode × and – (e.g. "1.7 × 10–4")
    if (grepl("\u00d7", v)) {
      parts <- strsplit(v, "\u00d7")[[1]]
      if (length(parts) == 2) {
        mantissa <- as.numeric(trimws(parts[1]))
        exp_part <- trimws(parts[2])
        # Handle Unicode minus sign – (U+2013) and regular hyphen
        exp_part <- gsub("\u2013", "-", exp_part)
        exp_val  <- suppressWarnings(as.numeric(exp_part))
        if (!is.na(mantissa) && !is.na(exp_val))
          return(mantissa * 10^exp_val)
      }
    }

    # Sum notation: "Σ = 0.005"
    if (grepl("\u03a3", v)) {
      num <- suppressWarnings(as.numeric(gsub(".*=\\s*", "", v)))
      if (!is.na(num)) return(num)
    }

    # Strip leading inequality qualifiers: "≤ 0.08 (Class 1)" → "0.08"
    v2 <- gsub("[\u2264\u2265<>]=?\\s*", "", v)   # strip ≤ ≥ < >
    v2 <- gsub("\\s*\\(.*\\)", "", v2)             # strip parenthetical notes
    v2 <- trimws(v2)

    num <- suppressWarnings(as.numeric(v2))
    if (!is.na(num)) return(num)

    NA_real_  # "No data", "NRG", "see footnote", "not applicable", etc.
  }, USE.NAMES = FALSE)
}

# =============================================================================
# I-6  US EPA Aquatic Life Benchmarks
# =============================================================================
message("\n--- I-6: US EPA Aquatic Life Benchmarks ---")

usepa_path <- "../Data/USEPA_aquatic_benchmarks.csv"
if (!file.exists(usepa_path)) stop("Missing: ", usepa_path)

usepa_raw <- read.csv(usepa_path, check.names = TRUE, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Column mapping — READ BEFORE EDITING
# ---------------------------------------------------------------------------
# EPA's source table uses superscript footnote letters as column headers
# (e.g. "AcuteA", "ChronicB"), and BOTH the Fish and Invertebrate columns
# share the same footnote letters, as do BOTH the Nonvascular and Vascular
# Plant columns. The raw CSV therefore has genuinely duplicated header text,
# and read.csv(check.names = TRUE) auto-renames duplicates by appending
# ".1" in column order — that suffix reflects POSITION only, not organism
# identity. A previous version of this script assumed "AcuteA.1" meant
# "invertebrate acute", which is wrong and silently used the wrong
# denominator for the crust/algae Benchmark Hazard Index (confirmed via
# comparison with EPA's live benchmarks page on 2026-07-05):
#
#                    Fish Acute (ug/L)   Invert Acute (ug/L)   Algae IC50 (ug/L)
#   Chlorpyrifos           0.9                  0.05                  140
#   Dichlorvos             50                   0.0334              14000
#
# Matching these known values against raw column position (1-indexed,
# including the leading Pesticide/Year/CAS columns) gives the TRUE mapping:
#   raw col  4 ("AcuteA")   -> Freshwater fish, acute
#   raw col  8 ("AcuteC")   -> Freshwater invertebrate, acute   (NOT col 6)
#   raw col 12 ("IC50E")    -> Nonvascular plant (algae), IC50  (NOT col 8)
#
# If you download a fresh copy of USEPA_aquatic_benchmarks.csv and EPA has
# changed the column layout, the sanity check below (which re-verifies these
# same two compounds after parsing) will fail loudly with instructions,
# instead of silently reproducing this bug. Do not "fix" a failing check by
# loosening the tolerance — re-derive the column positions instead, by
# comparing a few well-known compounds against EPA's live table:
# https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/aquatic-life-benchmarks-and-ecological-risk
# ---------------------------------------------------------------------------
FISH_ACUTE_COL   <- 4
INVERT_ACUTE_COL <- 8
ALGAE_IC50_COL   <- 12

.usepa_names <- names(usepa_raw)
.fish_col    <- .usepa_names[FISH_ACUTE_COL]
.invert_col  <- .usepa_names[INVERT_ACUTE_COL]
.algae_col   <- .usepa_names[ALGAE_IC50_COL]
message(sprintf("  Column mapping: fish = raw col %d (%s), invertebrate = raw col %d (%s), algae = raw col %d (%s)",
                FISH_ACUTE_COL, .fish_col, INVERT_ACUTE_COL, .invert_col, ALGAE_IC50_COL, .algae_col))

usepa <- usepa_raw %>%
  rename(cas_number = CAS.number) %>%
  # Keep only rows with a plausible CAS number (skip "NR", blank, footnote rows)
  filter(grepl("^[0-9]+-[0-9]+-[0-9]+$", trimws(cas_number))) %>%
  mutate(
    cas_number = trimws(cas_number),
    benchmark_usepa_fish_acute_ng_L  = .parse_benchmark(.data[[.fish_col]])   * UG_TO_NG,
    benchmark_usepa_crust_acute_ng_L = .parse_benchmark(.data[[.invert_col]]) * UG_TO_NG,
    benchmark_usepa_algae_acute_ng_L = .parse_benchmark(.data[[.algae_col]])  * UG_TO_NG
  ) %>%
  select(cas_number,
         benchmark_usepa_fish_acute_ng_L,
         benchmark_usepa_crust_acute_ng_L,
         benchmark_usepa_algae_acute_ng_L) %>%
  # If a CAS appears more than once, keep the row with the most non-NA values
  group_by(cas_number) %>%
  slice_max(order_by = (!is.na(benchmark_usepa_fish_acute_ng_L)) +
                       (!is.na(benchmark_usepa_crust_acute_ng_L)) +
                       (!is.na(benchmark_usepa_algae_acute_ng_L)),
            n = 1, with_ties = FALSE) %>%
  ungroup()

message(sprintf("  US EPA: %d compounds parsed", nrow(usepa)))
message(sprintf("    fish acute  : %d non-NA", sum(!is.na(usepa$benchmark_usepa_fish_acute_ng_L))))
message(sprintf("    crust acute : %d non-NA", sum(!is.na(usepa$benchmark_usepa_crust_acute_ng_L))))
message(sprintf("    algae acute : %d non-NA", sum(!is.na(usepa$benchmark_usepa_algae_acute_ng_L))))

# ---------------------------------------------------------------------------
# Sanity check — guards against silently reintroducing the column-swap bug
# above if a future USEPA_aquatic_benchmarks.csv download reorders columns.
# Tolerance is a factor of 3 either way: wide enough to absorb an ordinary
# EPA revision to a benchmark value, tight enough that a column swap
# (which produced >100x errors in the original bug) still fails loudly.
# ---------------------------------------------------------------------------
.check_usepa_benchmark <- function(cas, name, expected_fish, expected_invert, expected_algae_ic50, tol = 3) {
  row <- usepa[usepa$cas_number == cas, ]
  if (nrow(row) == 0) {
    warning(sprintf("USEPA benchmark sanity check skipped: CAS %s (%s) not found in source table.", cas, name))
    return(invisible())
  }
  parsed  <- c(fish = row$benchmark_usepa_fish_acute_ng_L[1]  / UG_TO_NG,
              invert = row$benchmark_usepa_crust_acute_ng_L[1] / UG_TO_NG,
              algae  = row$benchmark_usepa_algae_acute_ng_L[1] / UG_TO_NG)
  expected <- c(fish = expected_fish, invert = expected_invert, algae = expected_algae_ic50)
  ratio <- parsed / expected
  bad <- is.na(ratio) | ratio < 1 / tol | ratio > tol
  if (any(bad)) {
    stop(sprintf(paste0(
      "USEPA benchmark sanity check FAILED for %s (CAS %s), column(s) %s: ",
      "parsed = %s ug/L, expected ~ %s ug/L (within %gx). ",
      "This usually means USEPA_aquatic_benchmarks.csv's column layout has ",
      "changed since this script was written. Re-derive FISH_ACUTE_COL / ",
      "INVERT_ACUTE_COL / ALGAE_IC50_COL above by comparing known compounds ",
      "against EPA's live benchmarks page (see comment above) — do not ",
      "loosen this check instead."),
      name, cas, paste(names(bad)[bad], collapse = ", "),
      paste(round(parsed[bad], 5), collapse = ", "),
      paste(expected[bad], collapse = ", "), tol))
  }
}

.check_usepa_benchmark("2921-88-2", "Chlorpyrifos", expected_fish = 0.9, expected_invert = 0.05,   expected_algae_ic50 = 140)
.check_usepa_benchmark("62-73-7",   "Dichlorvos",   expected_fish = 50,  expected_invert = 0.0334, expected_algae_ic50 = 14000)
message("  Sanity check passed: USEPA column mapping verified against Chlorpyrifos and Dichlorvos reference values.")

# =============================================================================
# I-7  EU WFD Environmental Quality Standards (Directive 2013/39/EU)
# =============================================================================
message("\n--- I-7: EU WFD Environmental Quality Standards ---")

eu_path <- "../Data/EUEPA_aquatic_benchmarks.csv"
if (!file.exists(eu_path)) stop("Missing: ", eu_path)

# File structure: 4 rows of legend, then row 5 = column header, data from row 6.
# Only 4 columns present: name, sub-name/blank, CAS, AA-EQS value (µg/L).
# This CSV contains only the single AA-EQS value — no separate inland/marine split.
eu_raw <- read.csv(eu_path, skip = 4, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("name", "subname", "cas_number", "aa_eqs_raw"))

eu <- eu_raw %>%
  slice(-1) %>%   # drop the column-header row now stored as data
  mutate(
    cas_number = trimws(cas_number),
    # Prefer sub-name for grouped entries (e.g. aldrin, dieldrin listed under parent)
    cas_number = ifelse(nzchar(trimws(subname)) & !nzchar(cas_number),
                        NA_character_, cas_number)
  ) %>%
  filter(grepl("^[0-9]+-[0-9]+-[0-9]+$", cas_number)) %>%
  mutate(
    benchmark_eu_eqs_aa_marine_ng_L = .parse_benchmark(aa_eqs_raw) * UG_TO_NG
  ) %>%
  select(cas_number, benchmark_eu_eqs_aa_marine_ng_L) %>%
  group_by(cas_number) %>%
  slice_max(order_by = !is.na(benchmark_eu_eqs_aa_marine_ng_L),
            n = 1, with_ties = FALSE) %>%
  ungroup()

message(sprintf("  EU EQS: %d compounds parsed, %d with AA-EQS value",
                nrow(eu), sum(!is.na(eu$benchmark_eu_eqs_aa_marine_ng_L))))

# =============================================================================
# I-8  Australia ANZG / ANZECC Default Guideline Values
# =============================================================================
message("\n--- I-8: Australia ANZG Guideline Values ---")

au_path <- "../Data/Australia_aquatic_benchmarks.csv"
if (!file.exists(au_path)) stop("Missing: ", au_path)

au_raw <- read.csv(au_path, stringsAsFactors = FALSE, check.names = TRUE)

# Each compound has two rows: Freshwater and Marine Water.
# Use LOSP 95% (95% species protection) — column Tox.LOSP.95.
# Unit column (Tox.LOSP.Unit) should be µg/L for all rows.

au_fresh <- au_raw %>%
  filter(trimws(Tox.Medium) == "Freshwater") %>%
  mutate(
    cas_number = trimws(Tox.CAS),
    benchmark_au_anzg_fresh_ng_L = .parse_benchmark(Tox.LOSP.95) * UG_TO_NG
  ) %>%
  filter(grepl("^[0-9]+-[0-9]+-[0-9]+$", cas_number)) %>%
  select(cas_number, benchmark_au_anzg_fresh_ng_L) %>%
  group_by(cas_number) %>%
  slice_max(order_by = !is.na(benchmark_au_anzg_fresh_ng_L),
            n = 1, with_ties = FALSE) %>%
  ungroup()

au_marine <- au_raw %>%
  filter(trimws(Tox.Medium) == "Marine Water") %>%
  mutate(
    cas_number = trimws(Tox.CAS),
    benchmark_au_anzg_marine_ng_L = .parse_benchmark(Tox.LOSP.95) * UG_TO_NG
  ) %>%
  filter(grepl("^[0-9]+-[0-9]+-[0-9]+$", cas_number)) %>%
  select(cas_number, benchmark_au_anzg_marine_ng_L) %>%
  group_by(cas_number) %>%
  slice_max(order_by = !is.na(benchmark_au_anzg_marine_ng_L),
            n = 1, with_ties = FALSE) %>%
  ungroup()

au <- full_join(au_fresh, au_marine, by = "cas_number")

message(sprintf("  AU ANZG: %d compounds total", nrow(au)))
message(sprintf("    freshwater LOSP95 : %d non-NA", sum(!is.na(au$benchmark_au_anzg_fresh_ng_L))))
message(sprintf("    marine LOSP95     : %d non-NA", sum(!is.na(au$benchmark_au_anzg_marine_ng_L))))

# =============================================================================
# I-9  Canada CCME Water Quality Guidelines
# =============================================================================
message("\n--- I-9: Canada CCME Water Quality Guidelines ---")

ca_path <- "../Data/Canada_aquatic_benchmarks.csv"
if (!file.exists(ca_path)) stop("Missing: ", ca_path)

ca_raw <- read.csv(ca_path, stringsAsFactors = FALSE, check.names = TRUE)

# Use freshwater long-term (chronic) guideline.
# Column name after check.names: X.Freshwater..Concentration..µg.L..Long.Term.
# Non-numeric values ("No data", "NRG", "Insufficient data") → NA via .parse_benchmark.

lt_col <- grep("Freshwater.*Long.Term", names(ca_raw), value = TRUE)[1]
if (is.na(lt_col)) stop("Cannot find Freshwater Long Term column in Canada CSV. Names: ",
                         paste(names(ca_raw), collapse = ", "))
message(sprintf("  Using column: %s", lt_col))

ca <- ca_raw %>%
  mutate(
    cas_number = trimws(CASRN),
    benchmark_ca_ccme_fresh_lt_ng_L = .parse_benchmark(.data[[lt_col]]) * UG_TO_NG
  ) %>%
  filter(grepl("^[0-9]+-[0-9]+-[0-9]+$", cas_number)) %>%
  select(cas_number, benchmark_ca_ccme_fresh_lt_ng_L) %>%
  group_by(cas_number) %>%
  slice_max(order_by = !is.na(benchmark_ca_ccme_fresh_lt_ng_L),
            n = 1, with_ties = FALSE) %>%
  ungroup()

message(sprintf("  CA CCME: %d compounds parsed, %d with long-term value",
                nrow(ca), sum(!is.na(ca$benchmark_ca_ccme_fresh_lt_ng_L))))

# =============================================================================
# I-11  Join all benchmarks into taxotox_reference.fst
# =============================================================================
message("\n--- I-11: Joining benchmarks into taxotox_reference.fst ---")

ref <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
  mutate(cas_number = as.character(cas_number))

# Drop any existing benchmark columns (safe re-run)
benchmark_cols <- c("benchmark_usepa_fish_acute_ng_L",
                    "benchmark_usepa_crust_acute_ng_L",
                    "benchmark_usepa_algae_acute_ng_L",
                    "benchmark_eu_eqs_aa_marine_ng_L",
                    "benchmark_au_anzg_fresh_ng_L",
                    "benchmark_au_anzg_marine_ng_L",
                    "benchmark_ca_ccme_fresh_lt_ng_L")
ref <- ref %>% select(-any_of(benchmark_cols))

# Join one at a time — all are one-row-per-cas_number tables
ref <- ref %>%
  left_join(usepa, by = "cas_number") %>%
  left_join(eu,    by = "cas_number") %>%
  left_join(au,    by = "cas_number") %>%
  left_join(ca,    by = "cas_number")

# Coverage report
message("\nBenchmark coverage in reference table:")
for (col in benchmark_cols) {
  n   <- sum(!is.na(ref[[col]]))
  pct <- round(100 * n / nrow(ref), 1)
  message(sprintf("  %-42s  %4d / %d  (%4.1f%%)", col, n, nrow(ref), pct))
}

write_fst(ref, "../Data/taxotox_reference.fst", compress = 50)
message(sprintf("\ntaxotox_reference.fst updated: %d rows, %d columns.",
                nrow(ref), ncol(ref)))
