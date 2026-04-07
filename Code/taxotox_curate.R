# =============================================================================
# taxotox_curate.R
# -----------------------------------------------------------------------------
# Purpose : Interactive curation of temp_CAS — review compounds that were
#           encountered during app sessions but not yet in Known_CAS.fst.
#           Reads from Google Sheets when TAXOTOX_SHEET_ID env var is set
#           (deployed mode); otherwise reads local temp_CAS.fst (dev mode).
#           For each compound the researcher chooses:
#             [A]dd    — append to Known_CAS.fst (permanent)
#             [S]kip   — leave in temp_CAS.fst for the next curation run
#             [R]eject — remove from temp_CAS.fst permanently
#
# After curation, run taxotox_install.R (via update_ecotox.R) to rebuild
# taxotox_reference.fst and pick up the newly approved compounds.
#
# Authors : Yair Suari & Noam Gridish, 2025
# =============================================================================

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

library(dplyr)
library(fst)

known_path <- file.path(.script_dir, "..", "Data", "Known_CAS.fst")
temp_path  <- file.path(.script_dir, "..", "Data", "temp_CAS.fst")
sheet_id   <- Sys.getenv("TAXOTOX_SHEET_ID",
                          unset = "1pftfWQNfIStasPqH1CvpDSAH3yHxYmjhggvAwlBaxDA")
USE_SHEETS <- nzchar(sheet_id)

# ── Auth helper (mirrors app.R) ───────────────────────────────────────────────
.GS4_KEY_PATH <- file.path(.script_dir, "taxotox-service.json")
.gs4_auth_auto <- function() {
  if (file.exists(.GS4_KEY_PATH)) {
    googlesheets4::gs4_auth(path = .GS4_KEY_PATH)
  } else {
    googlesheets4::gs4_auth()   # interactive OAuth fallback
  }
}

# ── Load temp_CAS ─────────────────────────────────────────────────────────────

if (USE_SHEETS) {
  library(googlesheets4)
  .gs4_auth_auto()
  temp_cas <- tryCatch(
    googlesheets4::read_sheet(sheet_id, col_types = "cc") %>% as.data.frame(),
    error = function(e) stop("Cannot read Google Sheet: ", conditionMessage(e))
  )
  message(sprintf("Read %d row(s) from Google Sheet.", nrow(temp_cas)))
} else {
  if (!file.exists(temp_path)) {
    message("temp_CAS.fst not found — nothing to curate.")
    stop("temp_CAS.fst not found.", call. = FALSE)
  }
  temp_cas <- read.fst(temp_path, as.data.table = FALSE)
}

known_cas <- read.fst(known_path, as.data.table = FALSE)

if (nrow(temp_cas) == 0) {
  message("temp_CAS.fst is empty — nothing to curate.")
  stop("Nothing to curate.", call. = FALSE)
}

# ── Remove temp_CAS entries already present in Known_CAS ──────────────────────
already_known <- temp_cas %>%
  filter(CASRN %in% known_cas$CASRN | PREFERRED_NAME %in% known_cas$PREFERRED_NAME)

if (nrow(already_known) > 0) {
  cat(sprintf("Skipping %d compound(s) already in Known_CAS:\n", nrow(already_known)))
  cat(paste0("  ", already_known$PREFERRED_NAME, " (", already_known$CASRN, ")"),
      sep = "\n")
  cat("\n")
  temp_cas <- temp_cas %>%
    filter(!(CASRN %in% known_cas$CASRN | PREFERRED_NAME %in% known_cas$PREFERRED_NAME))
}

if (nrow(temp_cas) == 0) {
  message("All compounds already in Known_CAS — nothing left to review.")
  stop("Nothing left to review.", call. = FALSE)
}

cat(sprintf("=== TaxoTox Curation ===\n%d compound(s) to review.\n\n", nrow(temp_cas)))

# ── Interactive loop ───────────────────────────────────────────────────────────

added    <- character(0)
skipped  <- character(0)
rejected <- character(0)

for (i in seq_len(nrow(temp_cas))) {
  row <- temp_cas[i, ]

  cat(sprintf("[%d / %d]\n", i, nrow(temp_cas)))
  cat(sprintf("  Name  : %s\n", row$PREFERRED_NAME))
  cat(sprintf("  CASRN : %s\n", row$CASRN))
  if ("source" %in% names(row) && nzchar(as.character(row$source)))
    cat(sprintf("  Source: %s\n", row$source))
  if ("date" %in% names(row) && !is.na(row$date))
    cat(sprintf("  Date  : %s\n", as.character(row$date)))

  repeat {
    ans <- toupper(trimws(readline("  Action [A]dd / [S]kip / [R]eject: ")))
    if (ans %in% c("A", "S", "R")) break
    cat("  Please enter A, S, or R.\n")
  }

  if (ans == "A") {
    added <- c(added, row$PREFERRED_NAME)
    new_row <- data.frame(PREFERRED_NAME = row$PREFERRED_NAME,
                          CASRN          = gsub("-", "", trimws(row$CASRN)),
                          stringsAsFactors = FALSE)
    known_cas <- bind_rows(known_cas, new_row)
  } else if (ans == "R") {
    rejected <- c(rejected, row$PREFERRED_NAME)
  } else {
    skipped <- c(skipped, row$PREFERRED_NAME)
  }
  cat("\n")
}

# ── Write results ──────────────────────────────────────────────────────────────

if (length(added) > 0) {
  known_cas <- known_cas %>%
    distinct(CASRN,          .keep_all = TRUE) %>%
    distinct(PREFERRED_NAME, .keep_all = TRUE)
  write_fst(known_cas, known_path, compress = 50)
  cat(sprintf("Known_CAS.fst updated: %d compound(s) added.\n", length(added)))
  cat(paste0("  + ", added, collapse = "\n"), "\n\n")
}

# Keep only skipped compounds in temp_CAS
temp_remaining <- temp_cas %>% filter(PREFERRED_NAME %in% skipped)

if (USE_SHEETS) {
  # Overwrite the sheet with only the skipped rows (adds/rejects are removed)
  googlesheets4::write_sheet(data = temp_remaining, ss = sheet_id, sheet = 1)
  if (nrow(temp_remaining) == 0)
    message("Google Sheet cleared.")
  else
    message(sprintf("Google Sheet updated: %d skipped compound(s) retained.", nrow(temp_remaining)))
} else {
  write_fst(temp_remaining, temp_path, compress = 50)
}

cat(sprintf("Summary: %d added, %d skipped, %d rejected.\n\n",
            length(added), length(skipped), length(rejected)))

# ── Next step ─────────────────────────────────────────────────────────────────
# taxotox_install.R / update_ecotox.R will pick up newly added Known_CAS entries
# and populate taxotox_reference.fst (ECOTOX lookup + CompTox gap-fill) on the
# next full upgrade run. No reference-table update is needed here.

if (length(added) > 0) {
  message(sprintf(
    "%d compound(s) added to Known_CAS. Run taxotox_install.R (or update_ecotox.R) ",
    length(added)
  ), "to populate taxotox_reference.fst for the new entries.")
}

cat("Curation complete.\n")
