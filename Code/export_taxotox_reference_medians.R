# Set working directory to this script's location.
.script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) tryCatch(
    dirname(normalizePath(sys.frames()[[1]]$ofile)),
    error = function(e) {
      f <- sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])
      if (length(f) && nzchar(f)) dirname(normalizePath(f))
      else getwd()
    }
  )
)
setwd(.script_dir)

# =============================================================================
# export_taxotox_reference_medians.R
# -----------------------------------------------------------------------------
# Purpose : Export cas_number / ecotox_group / median_lc50_ng_L (and a few
#           other useful columns) from ../Data/taxotox_reference.fst to CSV,
#           so the TaxoTox paper's ECOTOX median-validation check
#           (2MethodComparison/codes/add_ecotox_medians.py) can cross-
#           reference TaxoTox's own reference medians without needing R.
#
#           Replaces an earlier ad-hoc export that only ever existed in a
#           since-deleted Claude Code session scratchpad -- this version is
#           checked into the repo so it's reproducible.
#
# Input   : ../Data/taxotox_reference.fst
# Output  : ../../TaxoTox paper/data/external/taxotox_reference_medians.csv
# =============================================================================

suppressMessages(library(fst))

ref <- read_fst("../Data/taxotox_reference.fst")

out <- ref[, c("cas_number", "ecotox_group", "chemical_name", "n_ecotox",
               "median_lc50_ng_L", "hc5_ssd_ng_L")]

out_path <- "../../TaxoTox paper/data/external/taxotox_reference_medians.csv"
write.csv(out, out_path, row.names = FALSE)

cat("Wrote", nrow(out), "rows to", normalizePath(out_path), "\n")
