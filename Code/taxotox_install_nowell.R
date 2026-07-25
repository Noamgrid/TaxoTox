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
# taxotox_install_nowell.R  (step I-10)
# -----------------------------------------------------------------------------
# Purpose : Join Nowell et al. (2014)'s own PUBLISHED Sensitive/Median Toxicity
#           Concentration values (STC/MTC) into taxotox_reference.fst, for the
#           Taxon-Sensitive PTI method.
#
# Source  : Data/Nowell2014_AppB.xlsx -- official supplementary Appendix B from
#           Nowell, L.H., Norman, J.E., Moran, P.W., Martin, J.D., & Stone, W.W.
#           (2014), Sci. Total Environ. 476-477:144-157,
#           https://doi.org/10.1016/j.scitotenv.2013.12.088
#           Downloaded from:
#           https://water.usgs.gov/nawqa/pnsp/pubs/Nowell2014_STOTEN_PTI/Nowell2014_SuppInfo_PTI.zip
#           Sheets used: "Table B.1 - Fish", "Table B.2 - Cladocerans".
#           (Table B.3 -- Benthic invertebrates -- is out of scope for TaxoTox,
#           which does not model that taxon.)
#
# Why published values instead of recomputing from ECOTOX (previous approach,
# taxotox_install.R's old "7a-ii" section):
#   - TaxoTox's own re-derivation from a local ECOTOX snapshot diverged
#     substantially from Nowell's published numbers -- e.g. Diazinon fish STC:
#     ~45 ug/L recomputed vs 85 ug/L published; Atrazine fish STC: ~2099 ug/L
#     recomputed vs 4500 ug/L published. Nowell's STC also draws on USEPA OPP
#     aquatic-life benchmarks and the Pesticide Properties Database (PPDB),
#     not ECOTOX alone, and the local ECOTOX snapshot has drifted from the
#     2012-2013 pull Nowell et al. used.
#   - Using the table Nowell published makes TaxoTox's Taxon-Sensitive PTI
#     directly comparable to values reported in downstream USGS literature
#     (e.g. Covert et al. 2020), which use this table as-is.
#
# "non-std" flagged rows (the "Toxicity value type/source" column, e.g.
# "non-std OPP", "non-std PPDB") are used AS PUBLISHED, not excluded: hundreds
# of compounds have no standard-duration ECOTOX data at all, so their
# published STC/MTC IS a non-std value -- excluding those would just delete
# the compound from coverage, not make the surviving data more standard.
# Nowell's table already represents the single best available value per
# compound; the flag documents a known bias direction, it isn't an
# instruction to discard the row. This also matches how downstream users
# (e.g. Covert et al. 2020) apply the table.
#
# Output  : Updates ../Data/taxotox_reference.fst in place.
#   stc_nowell_ng_L -- Nowell's published STC (5th percentile or minimum of
#                      individual toxicity values; see Nowell AppA.pdf sec 3-4)
#   n_stc_nowell    -- Nowell's "No. bioassays" count backing that STC
#   Compounds Nowell covers (via ECOTOX/OPP/PPDB) that TaxoTox's own ECOTOX
#   pull never surfaced are added as brand-new rows, the same pattern step I-2
#   (CompTox gap-fill) uses for compounds absent from ECOTOX.
#
# Prerequisites : taxotox_reference.fst must exist (run taxotox_install.R first)
# Authors : Yair Suari & Noam Gridish, 2025
# =============================================================================

library(dplyr)
library(readxl)
library(fst)

UG_TO_NG <- 1000  # source values are in ug/L; convert to ng/L

nowell_path <- "../Data/Nowell2014_AppB.xlsx"
if (!file.exists(nowell_path)) stop("Missing: ", nowell_path)

if (!file.exists("../Data/taxotox_reference.fst")) {
  stop("taxotox_reference.fst not found -- run taxotox_install.R first.")
}

# =============================================================================
# Parse one taxon sheet (Table B.1 or B.2). First 2 rows are title/notes; row 3
# is the header. Columns are addressed by position (not name) because the
# STC column's unit symbol (micro sign, U+00B5) doesn't round-trip reliably
# through Excel/R encodings -- position is more robust for this fixed,
# officially-published table layout. The sanity check below guards against
# silently misreading a future revision of the source file.
# =============================================================================
.read_nowell_sheet <- function(sheet, group) {
  raw <- read_excel(nowell_path, sheet = sheet, skip = 2)
  nm  <- names(raw)

  stopifnot(
    "Nowell AppB column layout changed: expected col 1 = Pesticide name" =
      grepl("Pesticide name", nm[1], ignore.case = TRUE),
    "Nowell AppB column layout changed: expected col 3 = CAS Registry Number" =
      grepl("CAS", nm[3], ignore.case = TRUE),
    "Nowell AppB column layout changed: expected col 4 = No. bioassays" =
      grepl("bioassays", nm[4], ignore.case = TRUE),
    "Nowell AppB column layout changed: expected col 9 = STC" =
      grepl("^STC", nm[9])
  )

  raw %>%
    transmute(
      chemical_name   = .data[[nm[1]]],
      cas_number_raw  = trimws(.data[[nm[3]]]),
      n_stc_nowell    = as.integer(.data[[nm[4]]]),
      stc_nowell_ng_L = as.numeric(.data[[nm[9]]]) * UG_TO_NG,
      ecotox_group    = group
    ) %>%
    # Drop rows with no assigned CAS number ("NAV" -- unregistered
    # degradates/mixtures in Nowell's table; cannot join without a CAS).
    filter(grepl("^[0-9]+(-[0-9]+)+$", cas_number_raw)) %>%
    mutate(cas_number = gsub("-", "", cas_number_raw)) %>%
    select(cas_number, chemical_name, ecotox_group, stc_nowell_ng_L, n_stc_nowell)
}

nowell_fish       <- .read_nowell_sheet("Table B.1 - Fish", "fish")
nowell_cladoceran <- .read_nowell_sheet("Table B.2 - Cladocerans", "crustacean")

message(sprintf("Nowell fish: %d compounds parsed", nrow(nowell_fish)))
message(sprintf("Nowell cladoceran: %d compounds parsed", nrow(nowell_cladoceran)))

nowell_stc <- bind_rows(nowell_fish, nowell_cladoceran)

# ---------------------------------------------------------------------------
# Sanity check -- guards against silently misparsing a future revision of
# Nowell2014_AppB.xlsx (e.g. reordered/renamed columns). Values below are
# taken directly from the published Table B.1/B.2 rows for these compounds.
# ---------------------------------------------------------------------------
.check_nowell <- function(cas, name, group, expected_stc_ug_L, tol = 0.01) {
  row <- nowell_stc[nowell_stc$cas_number == gsub("-", "", cas) &
                     nowell_stc$ecotox_group == group, ]
  if (nrow(row) == 0) {
    stop(sprintf("Nowell sanity check FAILED: %s (CAS %s, %s) not found in parsed table.",
                 name, cas, group))
  }
  parsed <- row$stc_nowell_ng_L[1] / UG_TO_NG
  if (abs(parsed - expected_stc_ug_L) / expected_stc_ug_L > tol) {
    stop(sprintf(paste(
      "Nowell sanity check FAILED for %s (CAS %s, %s): parsed STC = %s ug/L,",
      "expected %s ug/L. Nowell2014_AppB.xlsx structure may have changed --",
      "re-derive the column positions in .read_nowell_sheet() above."),
      name, cas, group, parsed, expected_stc_ug_L))
  }
}
.check_nowell("333-41-5",  "Diazinon",     "fish",       85)
.check_nowell("2921-88-2", "Chlorpyrifos", "crustacean", 0.05306)
message("Sanity check passed: Nowell AppB column mapping verified against Diazinon and Chlorpyrifos reference values.")

# =============================================================================
# Join into taxotox_reference.fst
# =============================================================================
ref <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
  mutate(cas_number = as.character(cas_number)) %>%
  select(-any_of(c("stc_nowell_ng_L", "n_stc_nowell")))  # safe re-run

# Add STC/n to existing (cas_number, ecotox_group) rows
ref <- ref %>%
  left_join(
    nowell_stc %>% select(cas_number, ecotox_group, stc_nowell_ng_L, n_stc_nowell),
    by = c("cas_number", "ecotox_group")
  )

# Compounds Nowell covers (via ECOTOX/OPP/PPDB) that TaxoTox's own ECOTOX pull
# never surfaced get added as new rows -- same pattern step I-2's CompTox
# gap-fill uses for compounds absent from ECOTOX.
existing_keys <- ref %>% select(cas_number, ecotox_group)

new_rows <- nowell_stc %>%
  anti_join(existing_keys, by = c("cas_number", "ecotox_group")) %>%
  mutate(n_ecotox = 0L)

ref <- bind_rows(ref, new_rows)

message("\nNowell STC coverage in reference table:")
message(sprintf("  fish       : %d / %d rows",
                sum(!is.na(ref$stc_nowell_ng_L[ref$ecotox_group == "fish"])),
                sum(ref$ecotox_group == "fish")))
message(sprintf("  crustacean : %d / %d rows",
                sum(!is.na(ref$stc_nowell_ng_L[ref$ecotox_group == "crustacean"])),
                sum(ref$ecotox_group == "crustacean")))
message(sprintf("  %d new rows added for compounds absent from TaxoTox's ECOTOX pull",
                nrow(new_rows)))

file.copy("../Data/taxotox_reference.fst",
         "../Data/taxotox_reference.fst.bak_pre_nowell_published", overwrite = TRUE)

write_fst(ref, "../Data/taxotox_reference.fst", compress = 50)
message(sprintf("\ntaxotox_reference.fst updated: %d rows, %d columns.", nrow(ref), ncol(ref)))
