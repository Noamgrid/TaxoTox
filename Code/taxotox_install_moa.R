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
# taxotox_install_moa.R
# -----------------------------------------------------------------------------
# Purpose : Join a real, externally curated Mode-of-Action (MoA) classification
#           into taxotox_reference.fst, for the CAMA method (Concentration
#           Addition within MoA groups, Independent Action between groups).
#
# Source  : Kramer, L., Schulze, T., Kluever, N., Altenburger, R.,
#           Hackermueller, J., Krauss, M., & Busch, W. (2024). Curated
#           mode-of-action data and effect concentrations for chemicals
#           relevant for the aquatic environment. Scientific Data, 11, 26.
#           https://doi.org/10.1038/s41597-023-02904-7
#           Data: https://doi.org/10.5281/zenodo.10071824 (free, no login)
#           Files used: Table_B_chemical_information.csv (chemical identity,
#           incl. cas_rn), Table_C_use_groups_and_mode_of_action_information_
#           for_chemicals.csv (MoA_broad/MoA_specific, keyed by the same `ID`
#           as Table B). Tables D-F (ecotoxicity data) and the R workflow /
#           graphics files in the Zenodo record are not needed -- TaxoTox has
#           its own ECOTOX pipeline for toxicity values.
#
# Why this replaces the old Step I-3 heuristic (removed from
# taxotox_install.R): that heuristic was a ~100-name regex pattern list
# defaulting everything else to "narcosis" (see the comment left in its
# place in taxotox_install.R for detail). This dataset is a real,
# systematically curated classification -- 3,387 chemicals, of which 2,149
# have an actual MoA assignment (the rest are honestly labelled "unknown" in
# the source data itself, not defaulted to narcosis).
#
# moa_group is taken verbatim from Kramer et al.'s own MoA_broad column,
# including multi-mechanism compounds where MoA_broad is itself a comma-joined
# string (e.g. "Cardiovascular system, Neuromuscular system") -- CAMA only
# needs a consistent group label to partition compounds for Concentration
# Addition, not a single "primary" mechanism, so no information is discarded
# by keeping the source label as-is.
#
# Output  : Updates ../Data/taxotox_reference.fst in place.
#   moa_group  -- Kramer et al.'s MoA_broad value, or "unknown" if the
#                 compound has no MoA classification in the source (either
#                 because Kramer et al. couldn't classify it, or because it
#                 isn't in their dataset at all).
#   moa_source -- "Kramer2024" (real classification), "unknown_in_Kramer2024"
#                 (compound is in the source dataset but MoA_broad ==
#                 "unknown" there), or "not_in_Kramer2024" (compound isn't in
#                 the source dataset at all).
#
# Prerequisites : taxotox_reference.fst must exist (run taxotox_install.R first)
# Authors : Yair Suari & Noam Gridish, 2025
# =============================================================================

library(dplyr)
library(readr)
library(fst)

tableB_path <- "../Data/Kramer2024_TableB_chemicals.csv"
tableC_path <- "../Data/Kramer2024_TableC_moa.csv"
if (!file.exists(tableB_path)) stop("Missing: ", tableB_path)
if (!file.exists(tableC_path)) stop("Missing: ", tableC_path)

if (!file.exists("../Data/taxotox_reference.fst")) {
  stop("taxotox_reference.fst not found -- run taxotox_install.R first.")
}

# =============================================================================
# Parse and join Table B (chemical identity, has cas_rn) with Table C (MoA),
# both keyed by `ID`. Confirmed during development: both tables have exactly
# 3,387 rows, one row per ID, no duplicate cas_rn in Table B and no duplicate
# ID in Table C -- a clean 1:1:1 join between ID / cas_rn / MoA row, no
# tie-breaking needed.
# =============================================================================
table_b <- read_csv(tableB_path, show_col_types = FALSE)
table_c <- read_csv(tableC_path, show_col_types = FALSE)

stopifnot(
  "Kramer2024 Table B layout changed: expected an ID column" = "ID" %in% names(table_b),
  "Kramer2024 Table B layout changed: expected a cas_rn column" = "cas_rn" %in% names(table_b),
  "Kramer2024 Table C layout changed: expected an ID column" = "ID" %in% names(table_c),
  "Kramer2024 Table C layout changed: expected a MoA_broad column" = "MoA_broad" %in% names(table_c)
)

kramer_moa <- table_c %>%
  select(ID, MoA_broad) %>%
  inner_join(table_b %>% select(ID, cas_rn), by = "ID") %>%
  # Drop rows with no assigned CAS number ("n/a" -- mixtures/generic entries
  # in Kramer et al.'s table; cannot join without a CAS).
  filter(!is.na(cas_rn), cas_rn != "n/a") %>%
  transmute(
    cas_number = gsub("-", "", trimws(cas_rn)),
    moa_group  = if_else(is.na(MoA_broad) | MoA_broad == "", "unknown", MoA_broad),
    moa_source = if_else(moa_group == "unknown", "unknown_in_Kramer2024", "Kramer2024")
  )

message(sprintf("Kramer 2024 MoA table: %d compounds parsed (%d with a real MoA classification)",
                nrow(kramer_moa), sum(kramer_moa$moa_source == "Kramer2024")))

# ---------------------------------------------------------------------------
# Sanity check -- guards against silently misparsing a future revision of
# the Kramer2024 CSVs (e.g. reordered/renamed columns). Values below are
# confirmed directly against the downloaded Table C rows for these compounds.
# ---------------------------------------------------------------------------
.check_kramer <- function(cas, name, expected_moa_broad) {
  row <- kramer_moa[kramer_moa$cas_number == gsub("-", "", cas), ]
  if (nrow(row) == 0) {
    stop(sprintf("Kramer2024 sanity check FAILED: %s (CAS %s) not found in parsed table.",
                 name, cas))
  }
  if (row$moa_group[1] != expected_moa_broad) {
    stop(sprintf(paste(
      "Kramer2024 sanity check FAILED for %s (CAS %s): parsed moa_group = %s,",
      "expected %s. Kramer2024_TableC_moa.csv structure may have changed."),
      name, cas, row$moa_group[1], expected_moa_broad))
  }
}
.check_kramer("2921-88-2", "Chlorpyrifos", "Neuromuscular system")
.check_kramer("1912-24-9", "Atrazine",     "Photosynthesis inhibition")
.check_kramer("330-54-1",  "Diuron",       "Photosynthesis inhibition")
message("Sanity check passed: Kramer2024 column mapping verified against Chlorpyrifos, Atrazine, and Diuron reference values.")

# =============================================================================
# Join into taxotox_reference.fst
# =============================================================================
ref <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
  mutate(cas_number = as.character(cas_number)) %>%
  select(-any_of(c("moa_group", "moa_source")))  # safe re-run

# moa_group/moa_source apply per compound (same across all ecotox_groups), so
# join by cas_number alone -- no new rows added (a compound with MoA data but
# no concentration/LC50 denominator is inert for CAMA either way).
ref <- ref %>%
  left_join(kramer_moa, by = "cas_number") %>%
  mutate(
    moa_group  = if_else(is.na(moa_group), "unknown", moa_group),
    moa_source = if_else(is.na(moa_source), "not_in_Kramer2024", moa_source)
  )

message("\nmoa_group coverage in reference table:")
ref %>%
  distinct(cas_number, moa_group, moa_source) %>%
  count(moa_source) %>%
  arrange(desc(n)) %>%
  { message(paste(capture.output(print(as.data.frame(.))), collapse = "\n")) }

write_fst(ref, "../Data/taxotox_reference.fst", compress = 50)
message(sprintf("\ntaxotox_reference.fst updated: %d rows, %d columns.", nrow(ref), ncol(ref)))
