# =============================================================================
# taxotox_curate.R
# -----------------------------------------------------------------------------
# Purpose : Interactive curation of the pending-compounds Google Sheet —
#           review compounds that were encountered during app sessions but
#           not yet in Known_CAS.fst. Curation is Sheets-only (no local
#           temp_CAS.fst fallback): the Sheet is the single source of truth
#           so pending items persist across shinyapps.io restarts and
#           between dev machines, and there's exactly one queue to check.
#           For each compound the researcher chooses:
#             [A]dd    — append to Known_CAS.fst (permanent)
#             [S]kip   — leave in the Sheet for the next curation run
#             [R]eject — remove from the Sheet permanently
#
# This script is also invoked automatically (when running interactively)
# from taxotox_install.R's Step 0b if the Sheet has pending items -- run it
# standalone any time you want to curate without doing a full install.
#
# After curation, run taxotox_install.R to rebuild taxotox_reference.fst and
# pick up the newly approved compounds.
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
library(googlesheets4)

known_path <- file.path(.script_dir, "..", "Data", "Known_CAS.fst")
sheet_id   <- Sys.getenv("TAXOTOX_SHEET_ID",
                          unset = "1pftfWQNfIStasPqH1CvpDSAH3yHxYmjhggvAwlBaxDA")
if (!nzchar(sheet_id))
  stop("TAXOTOX_SHEET_ID is empty — curation is Sheets-only and needs a Sheet ID.", call. = FALSE)

# ── Auth helper (mirrors app.R) ───────────────────────────────────────────────
.GS4_KEY_PATH <- file.path(.script_dir, "taxotox-service.json")
.gs4_auth_auto <- function() {
  if (file.exists(.GS4_KEY_PATH)) {
    googlesheets4::gs4_auth(path = .GS4_KEY_PATH)
  } else {
    googlesheets4::gs4_auth()   # interactive OAuth fallback
  }
}

# ── Load pending compounds from the Google Sheet ──────────────────────────────

.gs4_auth_auto()
temp_cas <- tryCatch(
  googlesheets4::read_sheet(sheet_id, col_types = "ccc") %>% as.data.frame(),
  error = function(e) stop("Cannot read Google Sheet: ", conditionMessage(e))
)
message(sprintf("Read %d row(s) from Google Sheet.", nrow(temp_cas)))

known_cas <- read.fst(known_path, as.data.table = FALSE)

if (nrow(temp_cas) == 0) {
  message("Google Sheet is empty — nothing to curate.")
  stop("Nothing to curate.", call. = FALSE)
}

# ── Remove entries already present in Known_CAS ────────────────────────────────
# Compares normalised (digits-only CASRN, lowercased/trimmed name) values --
# Known_CAS.fst mixes dashed and dash-free CASRNs historically, and a raw
# string comparison silently misses matches across that formatting gap
# (confirmed: dozens of already-known compounds were slipping through the
# old exact-match filter for exactly this reason).
.norm_cas  <- function(x) gsub("[^0-9]", "", x)
.norm_name <- function(x) tolower(trimws(x))

# Standard dashed CASRN form (leading digits)-(2 digits)-(1 check digit),
# counting from the right -- mirrors app.R's .canon_casrn(). Keeping CASRN
# dashed matters beyond cosmetics: the CompTox API's DTXSID search only
# reliably resolves dashed CASRNs (~99% dashed vs. ~22% dash-free, confirmed
# against Data/dtxsid_map.csv), so a dash-free CASRN entering Known_CAS.fst
# quietly loses CompTox gap-fill coverage for that compound.
.canon_casrn <- function(x) {
  d <- .norm_cas(trimws(x))
  vapply(d, function(di) {
    if (nchar(di) < 4) return(di)
    n <- nchar(di)
    paste0(substr(di, 1, n - 3), "-", substr(di, n - 2, n - 1), "-", substr(di, n, n))
  }, character(1), USE.NAMES = FALSE)
}

already_known <- temp_cas %>%
  filter((nzchar(.norm_cas(CASRN)) & .norm_cas(CASRN) %in% .norm_cas(known_cas$CASRN)) |
         .norm_name(PREFERRED_NAME) %in% .norm_name(known_cas$PREFERRED_NAME))

if (nrow(already_known) > 0) {
  cat(sprintf("Skipping %d compound(s) already in Known_CAS:\n", nrow(already_known)))
  cat(paste0("  ", already_known$PREFERRED_NAME, " (", already_known$CASRN, ")"),
      sep = "\n")
  cat("\n")
  temp_cas <- temp_cas %>%
    filter(!((nzchar(.norm_cas(CASRN)) & .norm_cas(CASRN) %in% .norm_cas(known_cas$CASRN)) |
             .norm_name(PREFERRED_NAME) %in% .norm_name(known_cas$PREFERRED_NAME)))
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
                          CASRN          = .canon_casrn(row$CASRN),
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

# Keep only skipped compounds in the Sheet (adds/rejects are removed)
temp_remaining <- temp_cas %>% filter(PREFERRED_NAME %in% skipped)

googlesheets4::write_sheet(data = temp_remaining, ss = sheet_id, sheet = 1)
if (nrow(temp_remaining) == 0) {
  message("Google Sheet cleared.")
} else {
  message(sprintf("Google Sheet updated: %d skipped compound(s) retained.", nrow(temp_remaining)))
}

cat(sprintf("Summary: %d added, %d skipped, %d rejected.\n\n",
            length(added), length(skipped), length(rejected)))

# ── Next step ─────────────────────────────────────────────────────────────────
# taxotox_install.R will pick up newly added Known_CAS entries and populate
# taxotox_reference.fst (ECOTOX lookup + CompTox gap-fill) on the next full
# install run. No reference-table update is needed here.

if (length(added) > 0) {
  message(sprintf(
    "%d compound(s) added to Known_CAS. Run taxotox_install.R ",
    length(added)
  ), "to populate taxotox_reference.fst for the new entries.")
}

cat("Curation complete.\n")
