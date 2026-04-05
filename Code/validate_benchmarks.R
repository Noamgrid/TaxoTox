# Set working directory to this script's location.
# Works in: RStudio (interactive/Source), VSCode source(), Rscript --file=
.script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),          # RStudio
  error = function(e) tryCatch(
    dirname(normalizePath(sys.frames()[[1]]$ofile)),             # source() in VSCode/R terminal
    error = function(e) {
      f <- sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])
      if (length(f) && nzchar(f)) dirname(normalizePath(f))     # Rscript --file=
      else getwd()
    }
  )
)
setwd(.script_dir)

# Purpose: Independent benchmark validation.
#   1. Query ECOTOX directly for LC50 / EC50 values (fish, algae, crustacean).
#   2. Convert units to ng/L, compute median per compound × taxonomic group.
#   3. Compare against medians in final_ecotox_data_ECOTOX_backup.fst
#      (produced by taxotox_install.R).
#   4. Write a self-contained HTML report.

library(RSQLite)
library(dplyr)
library(stringr)
library(tidyr)
library(fst)
library(htmltools)
library(DT)
library(plotly)

# ── 1. Query ECOTOX ───────────────────────────────────────────────────────────

conn <- dbConnect(RSQLite::SQLite(), ECOTOXr::get_ecotox_sqlite_file())

message("Querying ECOTOX...")
raw <- dbGetQuery(conn, "
SELECT
    c.cas_number,
    c.chemical_name,
    r.endpoint,
    r.conc1_mean,
    r.conc1_unit,
    CASE
        WHEN LOWER(s.ecotox_group) LIKE '%fish%'       THEN 'fish'
        WHEN LOWER(s.ecotox_group) LIKE '%algae%'      THEN 'algae'
        WHEN LOWER(s.ecotox_group) LIKE '%crustacean%' THEN 'crustacean'
    END AS ecotox_group
FROM results r
JOIN tests     t ON r.test_id        = t.test_id
JOIN species   s ON t.species_number = s.species_number
JOIN chemicals c ON t.test_cas       = c.cas_number
WHERE (LOWER(s.ecotox_group) LIKE '%fish%'
    OR LOWER(s.ecotox_group) LIKE '%algae%'
    OR LOWER(s.ecotox_group) LIKE '%crustacean%')
  AND r.conc1_mean > 0
  AND r.endpoint IN ('LC50','EC50','(log)LC50','(log)EC50',
                     'LC50*','EC50*','LC50/','EC50/')
")
dbDisconnect(conn)
message("  Rows returned: ", nrow(raw))

# ── 2. Unit conversion → ng/L ─────────────────────────────────────────────────

unit_factor <- c("ng/L" = 1, "ug/L" = 1e3, "mg/L" = 1e6, "g/L" = 1e9)

fresh <- raw %>%
  mutate(
    conc1_unit = str_remove(conc1_unit, "^AI\\s*"),
    conc1_unit = str_replace_all(conc1_unit, "\\bdm3\\b", "L"),
    conc1_unit = str_replace_all(conc1_unit, "\\bppt\\b", "ng/L"),
    conc1_unit = str_replace_all(conc1_unit, "\\bppm\\b", "mg/L"),
    conc1_unit = str_replace_all(conc1_unit, "\\bppb\\b", "ug/L"),
    conc1_unit = str_replace_all(conc1_unit, "0/00",      "g/L")
  ) %>%
  filter(conc1_unit %in% names(unit_factor)) %>%
  mutate(
    conc1_mean = as.numeric(conc1_mean),
    conc1_mean = if_else(str_detect(endpoint, "^\\(log\\)"), 10^conc1_mean, conc1_mean),
    conc_ng_L  = conc1_mean * unit_factor[conc1_unit],
    cas_number = as.character(cas_number)
  ) %>%
  filter(conc_ng_L > 0)

# ── 3. Median per compound × group ───────────────────────────────────────────

fresh_med <- fresh %>%
  group_by(cas_number, chemical_name, ecotox_group) %>%
  summarise(median_fresh_ngL = median(conc_ng_L, na.rm = TRUE),
            n_fresh          = n(),
            .groups = "drop")

message("Unique compound × group combinations (fresh): ", nrow(fresh_med))

# ── 4. Load FST and compute same medians ─────────────────────────────────────

fst_path <- if (file.exists("../Data/final_ecotox_data_ECOTOX_backup.fst"))
  "../Data/final_ecotox_data_ECOTOX_backup.fst" else "../Data/final_ecotox_data.fst"
message("FST source: ", fst_path)

fst_med <- read.fst(fst_path, as.data.table = FALSE) %>%
  mutate(cas_number = as.character(cas_number)) %>%
  group_by(cas_number, chemical_name, ecotox_group) %>%
  summarise(median_fst_ngL = median(conc_ng_L, na.rm = TRUE),
            n_fst          = n(),
            .groups = "drop")

message("Unique compound × group combinations (FST): ", nrow(fst_med))

# ── 5. Join and classify ──────────────────────────────────────────────────────

cmp <- full_join(fresh_med, fst_med,
                 by = c("cas_number", "ecotox_group"), suffix = c("", "_fst")) %>%
  mutate(
    chemical_name = coalesce(chemical_name, chemical_name_fst),
    ratio         = median_fresh_ngL / median_fst_ngL,
    log2_ratio    = log2(ratio),
    status = case_when(
      is.na(median_fresh_ngL) ~ "FST only",
      is.na(median_fst_ngL)   ~ "Fresh only",
      abs(log2_ratio) < 0.01  ~ "Match",
      TRUE                    ~ "Mismatch"
    )
  ) %>%
  select(cas_number, chemical_name, ecotox_group,
         n_fresh, median_fresh_ngL, n_fst, median_fst_ngL,
         ratio, log2_ratio, status) %>%
  arrange(desc(abs(coalesce(log2_ratio, 0))))

# ── 6. Build HTML report ──────────────────────────────────────────────────────

# 6a. Summary table
summary_df <- cmp %>%
  count(ecotox_group, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0)

# 6b. Scatter plot
plot_df <- cmp %>%
  filter(!is.na(median_fresh_ngL), !is.na(median_fst_ngL))

ax_range <- range(log10(c(plot_df$median_fresh_ngL, plot_df$median_fst_ngL)), na.rm = TRUE)

fig <- plot_ly(plot_df,
    x = ~median_fst_ngL, y = ~median_fresh_ngL,
    color = ~ecotox_group, colors = c("fish" = "#1f78b4", "algae" = "#33a02c", "crustacean" = "#e31a1c"),
    text = ~paste0("<b>", chemical_name, "</b><br>CAS: ", cas_number,
                   "<br>Group: ", ecotox_group,
                   "<br>Fresh median: ", signif(median_fresh_ngL, 3), " ng/L",
                   "<br>FST median: ",   signif(median_fst_ngL,   3), " ng/L",
                   "<br>Ratio (fresh/FST): ", round(ratio, 3)),
    hoverinfo = "text", type = "scatter", mode = "markers",
    marker = list(size = 5, opacity = 0.65)) %>%
  add_trace(
    x = 10^ax_range, y = 10^ax_range,
    type = "scatter", mode = "lines",
    line = list(color = "black", dash = "dash", width = 1),
    name = "1:1 line", inherit = FALSE, showlegend = TRUE) %>%
  layout(
    title = "Fresh extraction vs taxotox_install medians (ng/L, log scale)",
    xaxis = list(title = "taxotox_install FST median (ng/L)", type = "log"),
    yaxis = list(title = "Fresh extraction median (ng/L)",    type = "log"),
    legend = list(title = list(text = "Taxonomic group")))

# 6c. DT tables
fmt <- function(df) {
  df %>% mutate(across(c(median_fresh_ngL, median_fst_ngL, ratio), ~signif(.x, 4)),
                log2_ratio = round(log2_ratio, 3))
}

dt_opts <- list(pageLength = 30, scrollX = TRUE, autoWidth = TRUE)

dt_summary  <- datatable(summary_df,  rownames = FALSE, options = list(dom = "t"),
                          caption = "Compound × group counts by match status")
dt_mismatch <- datatable(fmt(filter(cmp, status == "Mismatch")),
                          rownames = FALSE, filter = "top", options = dt_opts,
                          caption = "Mismatches (|log2 ratio| ≥ 0.01)")
dt_fresh    <- datatable(fmt(filter(cmp, status == "Fresh only")),
                          rownames = FALSE, filter = "top", options = dt_opts,
                          caption = "In fresh extraction only (missing from FST)")
dt_fst      <- datatable(fmt(filter(cmp, status == "FST only")),
                          rownames = FALSE, filter = "top", options = dt_opts,
                          caption = "In FST only (missing from fresh extraction)")
dt_all      <- datatable(fmt(cmp), rownames = FALSE, filter = "top", options = dt_opts,
                          caption = "All compounds")

css <- "
  body { font-family: Arial, sans-serif; margin: 30px; color: #333; }
  h1   { color: #2c5f8a; border-bottom: 2px solid #2c5f8a; padding-bottom: 6px; }
  h2   { color: #4a4a4a; margin-top: 40px; }
  p.meta { color: #777; font-size: 0.9em; }
  .section { margin-bottom: 50px; }
"

report <- tagList(
  tags$html(lang = "en",
    tags$head(tags$meta(charset = "UTF-8"),
              tags$title("TaxoTox Benchmark Validation"),
              tags$style(HTML(css))),
    tags$body(
      tags$h1("TaxoTox Benchmark Validation"),
      tags$p(class = "meta",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"),
               " | FST: ", basename(fst_path),
               " | ECOTOX rows queried: ", nrow(raw))),

      tags$div(class = "section",
        tags$h2("1. Match Summary"),
        as.tags(dt_summary)),

      tags$div(class = "section",
        tags$h2("2. Fresh vs FST Scatter (log–log)"),
        as.tags(fig)),

      tags$div(class = "section",
        tags$h2("3. Mismatches"),
        as.tags(dt_mismatch)),

      tags$div(class = "section",
        tags$h2("4. Fresh Only (not in FST)"),
        as.tags(dt_fresh)),

      tags$div(class = "section",
        tags$h2("5. FST Only (not in fresh extraction)"),
        as.tags(dt_fst)),

      tags$div(class = "section",
        tags$h2("6. All Compounds"),
        as.tags(dt_all))
    )
  )
)

out <- "../../Testing_Taxotox/benchmark_validation.html"
save_html(report, out)
message("Report saved: ", normalizePath(out))
