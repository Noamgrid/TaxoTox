# temp_v2_spotcheck.R — V-2 spot-check: HC5 values for well-studied compounds
# Prints hc5_ssd_ng_L, hc5_model_ng_L, median_lc50_ng_L for 5 benchmark compounds.
# Run once, inspect terminal output, then delete.

library(dplyr)
library(fst)

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

ref <- read.fst("../Data/taxotox_reference.fst", as.data.table = FALSE) %>%
  mutate(cas_number = gsub("-", "", as.character(cas_number)))

UG_TO_NG <- 1000  # published values are in µg/L; reference table is in ng/L

# pub_lc50_*  : typical median LC50 range from acute toxicity tests (µg/L)
#               used to validate median_lc50_ng_L
# pub_hc5_*   : HC5 / regulatory guideline value (µg/L)
#               used to validate hc5_ssd_ng_L
# All from ECOTOX published summaries, ANZG 2018, or EU EQS basis documents.
spotcheck <- tibble::tribble(
  ~cas_number,  ~compound,      ~pub_lc50_ug_L_low, ~pub_lc50_ug_L_high, ~pub_hc5_ug_L_low, ~pub_hc5_ug_L_high, ~pub_group,     ~pub_source,
  "1912-24-9",  "Atrazine",     100,                10000,               1,                  20,                 "fish/algae",   "ECOTOX; ANZG 2018",
  "2921-88-2",  "Chlorpyrifos", 0.01,               1,                   0.0001,             0.001,              "fish/crust",   "ECOTOX; ANZG 2018 marine",
  "330-54-1",   "Diuron",       1000,               100000,              0.5,                5,                  "fish/algae",   "ECOTOX; EU EQS basis",
  "52645-53-1", "Permethrin",   0.1,                10,                  0.01,               0.1,                "fish/crust",   "ECOTOX; literature",
  "7440-50-8",  "Copper",       10,                 1000,                0.1,                2,                  "fish/crust",   "ECOTOX; ANZG 2018 marine"
) %>%
  mutate(
    pub_lc50_ng_L_low  = pub_lc50_ug_L_low  * UG_TO_NG,
    pub_lc50_ng_L_high = pub_lc50_ug_L_high * UG_TO_NG,
    pub_hc5_ng_L_low   = pub_hc5_ug_L_low   * UG_TO_NG,
    pub_hc5_ng_L_high  = pub_hc5_ug_L_high  * UG_TO_NG
  )

cat("\n=== V-2 HC5 Spot-Check ===\n")
cat("All values in ng/L. Published ranges from ECOTOX literature / ANZG 2018 / EU EQS.\n")
cat("Flag <<< CHECK if value is >10x outside expected range.\n")

.fmt <- function(x) {
  if (is.na(x)) return("NA")
  if (x >= 1e5)  return(format(round(x), big.mark = ",", scientific = FALSE))
  if (x > 10)    return(as.character(round(x)))
  return(formatC(x, digits = 3, format = "g"))
}
.flag <- function(val, lo, hi) {
  if (is.na(val) || is.na(lo)) return("")
  if (val < lo / 10 || val > hi * 10) " <<< CHECK" else ""
}

for (i in seq_len(nrow(spotcheck))) {
  cas         <- spotcheck$cas_number[i]
  name        <- spotcheck$compound[i]
  pub         <- spotcheck$pub_source[i]
  grp         <- spotcheck$pub_group[i]
  lc50_lo     <- spotcheck$pub_lc50_ng_L_low[i]
  lc50_hi     <- spotcheck$pub_lc50_ng_L_high[i]
  hc5_lo      <- spotcheck$pub_hc5_ng_L_low[i]
  hc5_hi      <- spotcheck$pub_hc5_ng_L_high[i]

  rows <- ref %>% filter(cas_number == gsub("-", "", cas)) %>%
    select(ecotox_group, n_ecotox, median_lc50_ng_L, hc5_ssd_ng_L, hc5_model_ng_L, hc5_method)

  if (nrow(rows) == 0) {
    cat(sprintf("%-15s %-14s  *** NOT FOUND IN REFERENCE TABLE ***\n", cas, name))
    next
  }

  cat(sprintf("\n%s  (%s)  —  Source: %s\n", name, cas, pub))
  cat(sprintf("  Expected LC50 range : %s – %s ng/L  [%s]\n",
              .fmt(lc50_lo), .fmt(lc50_hi), grp))
  cat(sprintf("  Expected HC5  range : %s – %s ng/L  [%s]\n",
              .fmt(hc5_lo), .fmt(hc5_hi), grp))
  cat(sprintf("  %-12s %8s %14s %12s %14s %10s\n",
              "Group", "n_ecotox", "median_lc50", "hc5_ssd", "hc5_model", "method"))
  cat(sprintf("  %s\n", strrep("-", 75)))

  for (j in seq_len(nrow(rows))) {
    r <- rows[j, ]
    lc50_flag <- .flag(r$median_lc50_ng_L, lc50_lo, lc50_hi)
    hc5_flag  <- .flag(r$hc5_ssd_ng_L,    hc5_lo,  hc5_hi)
    cat(sprintf("  %-12s %8s %14s %12s %14s %10s%s%s\n",
                r$ecotox_group,
                ifelse(is.na(r$n_ecotox), "NA", as.character(r$n_ecotox)),
                .fmt(r$median_lc50_ng_L),
                .fmt(r$hc5_ssd_ng_L),
                .fmt(r$hc5_model_ng_L),
                ifelse(is.na(r$hc5_method), "NA", r$hc5_method),
                lc50_flag, hc5_flag))
  }
}
