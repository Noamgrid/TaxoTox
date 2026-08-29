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

# =============================================================================
# validate_reference.R
# -----------------------------------------------------------------------------
# Purpose : Diagnostic validation of the taxotox_reference.fst table.
#           Produces diagnostic plots and tables to verify the full install pipeline.
#
# Inputs  :
#   ../Data/final_ecotox_data.fst      – raw test results (LC50 + EC50 rows)
#   ../Data/taxotox_reference.fst      – pre-aggregated denominator table
#                                        (if not yet produced, built inline)
#   ../Data/hc5_scaling_factors.csv    – k per ecotox_group from step I-1b
#
# Output  : ../Docs/validate_reference.html   (self-contained, no server needed)
#
# Plots   :
#   1. HC5 (SSD) vs median LC50        — confirms proportionality; validates k
#   2. Scaled HC5 vs SSD HC5           — agreement plot; quantifies scaling error
#   3. LC50 vs EC50 medians            — endpoint interchangeability check
#   4. Distribution of k               — spread of HC5/median ratio per group
#   5. OPERA predicted vs ECOTOX median — CompTox gap-fill calibration (if I-2 ran)
#   6. Coverage by lc50_source          — how many compounds gained a denominator
#   7. Mode of Action group coverage    — compounds per group per ecotox_group (if I-3 ran)
#   8. Benchmark framework coverage     — compounds per national benchmark (if I-11 ran)
#   9. CAMA vs Standard PTI convergence — E_mix vs PTI on the repo's own sample workbooks
#                                         (if I-3 ran and Data/sample_*.xlsx are present)
#
# Authors : Yair Suari & Noam Gridish, 2025
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(fst)
library(htmltools)
library(DT)
library(plotly)

# ── 1. Load data ──────────────────────────────────────────────────────────────

message("Loading final_ecotox_data.fst...")
ecotox <- read.fst("../Data/final_ecotox_data.fst", as.data.table = FALSE) %>%
  mutate(cas_number = as.character(cas_number))

message("  Rows: ", nrow(ecotox),
        "  |  Endpoints: ",
        paste(names(table(ecotox$endpoint)), collapse = ", "))

# taxotox_reference.fst is written at step I-4.
# If it exists use it; otherwise reconstruct from final_ecotox_data + scaling factors
# so this script can be run immediately after steps I-1a/I-1b.
ref_path    <- "../Data/taxotox_reference.fst"
sf_path     <- "../Data/hc5_scaling_factors.csv"

if (file.exists(ref_path)) {
  message("Loading taxotox_reference.fst...")
  ref <- read.fst(ref_path, as.data.table = FALSE) %>%
    mutate(cas_number = as.character(cas_number))
} else {
  message("taxotox_reference.fst not found — reconstructing from final_ecotox_data...")

  library(ssdtools)

  # Mirrors Code/taxotox_install.R's fit_hc5_ssd_full() (step 7b) -- kept in
  # sync manually; if that function changes, update this one too. Also
  # extracts sigma_ssd_log10 (the fitted log-normal SSD's sdlog, converted
  # from ssdtools' native natural-log units to log10), needed for the
  # msPAF_EC50 diagnostic (de Zwart & Posthuma, 2005, Environ. Toxicol. Chem.
  # 24(11):2665-2679, Eq. 4-5) further down this script. sigma_ssd_lnorm
  # (natural-log, ssdtools' raw sdlog) is kept only for audit -- coalescing
  # it directly into sigma_ec50 mixed log bases with sigma_group_fallback
  # (log10 by construction) and badly overestimated msPAF_EC50.
  fit_hc5_ssd_full <- function(values) {
    tryCatch({
      df <- data.frame(Conc = values[values > 0 & is.finite(values)])
      if (nrow(df) < 5) {
        return(tibble(hc5_ssd_ng_L = NA_real_, sigma_ssd_lnorm = NA_real_, sigma_ssd_log10 = NA_real_))
      }
      fit <- ssd_fit_dists(df, dists = "lnorm")
      hc  <- ssd_hc(fit, proportion = 0.05, ci = FALSE)
      co  <- coef(fit)
      sdlog_val <- co$est[co$dist == "lnorm" & co$term == "sdlog"]
      sdlog_val <- if (length(sdlog_val) == 1) sdlog_val else NA_real_
      tibble(
        hc5_ssd_ng_L    = as.numeric(hc$est[[1]]),
        sigma_ssd_lnorm = sdlog_val,
        sigma_ssd_log10 = sdlog_val / log(10)
      )
    }, error = function(e) tibble(hc5_ssd_ng_L = NA_real_, sigma_ssd_lnorm = NA_real_, sigma_ssd_log10 = NA_real_))
  }

  ref <- ecotox %>%
    group_by(cas_number, ecotox_group) %>%
    summarise(
      chemical_name    = first(chemical_name),
      n_ecotox         = n(),
      median_lc50_ng_L = median(conc_ng_L, na.rm = TRUE),
      .groups = "drop"
    )

  eligible <- ecotox %>%
    filter(!is.na(conc_ng_L), conc_ng_L > 0) %>%
    group_by(cas_number, ecotox_group) %>%
    filter(n() >= 5) %>%
    ungroup()

  hc5_ssd <- eligible %>%
    group_by(cas_number, ecotox_group) %>%
    group_modify(~ fit_hc5_ssd_full(.x$conc_ng_L)) %>%
    ungroup()

  ref <- ref %>%
    left_join(hc5_ssd, by = c("cas_number", "ecotox_group")) %>%
    mutate(hc5_model_ng_L = NA_real_,
           hc5_method     = if_else(!is.na(hc5_ssd_ng_L), "SSD", NA_character_))

  if (file.exists(sf_path)) {
    scaling_factors <- read.csv(sf_path)
    # sigma_group_fallback (de Zwart & Posthuma, 2005, Eq. 4): present in
    # hc5_scaling_factors.csv once regenerated by taxotox_install.R step 7c
    # (this msPAF_EC50 change); tolerate older copies of the CSV that predate
    # it by deriving it inline from the same k = scaling_factor already there.
    if (!"sigma_group_fallback" %in% names(scaling_factors)) {
      scaling_factors$sigma_group_fallback <- log10(scaling_factors$scaling_factor) / qnorm(0.05)
    }
    ref <- ref %>%
      left_join(scaling_factors, by = "ecotox_group") %>%
      mutate(
        hc5_model_ng_L = if_else(
          !is.na(median_lc50_ng_L) & median_lc50_ng_L > 0,
          median_lc50_ng_L * scaling_factor,
          NA_real_
        ),
        hc5_method = case_when(
          !is.na(hc5_ssd_ng_L) & !is.na(hc5_model_ng_L) ~ "Both",
          !is.na(hc5_ssd_ng_L)                            ~ "SSD",
          !is.na(hc5_model_ng_L)                          ~ "Scaled",
          TRUE                                             ~ NA_character_
        ),
        # sigma_ec50: same fitted-over-fallback priority as
        # taxotox_install.R step 7c (Posthuma & de Zwart, 2012, Environ.
        # Toxicol. Chem. 31(10):2175-2181). Must use sigma_ssd_log10 (not
        # sigma_ssd_lnorm) -- see fit_hc5_ssd_full() comment above.
        sigma_ec50 = dplyr::coalesce(sigma_ssd_log10, sigma_group_fallback)
      ) %>%
      select(-scaling_factor, -n_calibration, -sigma_group_fallback)
  }
}

message("Reference rows: ", nrow(ref))

# Load scaling factors (for annotation lines on Plot 1)
scaling_factors <- if (file.exists(sf_path)) read.csv(sf_path) else NULL

# ── 2. Group colours (consistent across all plots) ────────────────────────────

group_colours <- c(fish = "#1f78b4", algae = "#33a02c", crustacean = "#e31a1c")

# ── 3. Plot 1 — HC5 (SSD) vs Median LC50 (log-log) ───────────────────────────
# Expected: tight log-linear relationship (slope ≈ 1 on log scale).
# Vertical offset = log10(k), i.e. the scaling factor.
# Points sized by n_ecotox; lines = OLS per group.

p1_df <- ref %>%
  filter(!is.na(hc5_ssd_ng_L), !is.na(median_lc50_ng_L),
         median_lc50_ng_L > 0, hc5_ssd_ng_L > 0)

ax_range1 <- range(log10(c(p1_df$median_lc50_ng_L, p1_df$hc5_ssd_ng_L)),
                   na.rm = TRUE)
ax_lim1   <- 10^(ax_range1 + c(-0.5, 0.5))

# Build one trace per group for the scatter
fig1 <- plot_ly()

for (grp in names(group_colours)) {
  d <- filter(p1_df, ecotox_group == grp)
  if (nrow(d) == 0) next

  fig1 <- fig1 %>%
    add_trace(
      data = d,
      x = ~median_lc50_ng_L, y = ~hc5_ssd_ng_L,
      type = "scatter", mode = "markers",
      name = grp,
      text = ~paste0("<b>", chemical_name, "</b><br>CAS: ", cas_number,
                     "<br>Group: ", ecotox_group,
                     "<br>n_ecotox: ", n_ecotox,
                     "<br>Median LC50: ", signif(median_lc50_ng_L, 3), " ng/L",
                     "<br>HC5 (SSD): ",  signif(hc5_ssd_ng_L, 3),    " ng/L"),
      hoverinfo = "text",
      marker = list(
        size   = ~pmin(6 + sqrt(n_ecotox), 20),
        color  = group_colours[grp],
        opacity = 0.7
      )
    )

  # OLS line per group (on log scale)
  if (nrow(d) >= 3) {
    lm_fit <- lm(log10(hc5_ssd_ng_L) ~ log10(median_lc50_ng_L), data = d)
    x_seq  <- 10^seq(ax_range1[1] - 0.5, ax_range1[2] + 0.5, length.out = 80)
    y_pred <- 10^predict(lm_fit, newdata = data.frame(median_lc50_ng_L = x_seq))

    fig1 <- fig1 %>%
      add_trace(
        x = x_seq, y = y_pred,
        type = "scatter", mode = "lines",
        name = paste(grp, "OLS"),
        line = list(color = group_colours[grp], width = 1.5, dash = "dot"),
        showlegend = TRUE, inherit = FALSE
      )
  }
}

# 1:1 reference line
fig1 <- fig1 %>%
  add_trace(
    x = ax_lim1, y = ax_lim1,
    type = "scatter", mode = "lines",
    name = "1:1 line",
    line = list(color = "black", dash = "dash", width = 1),
    inherit = FALSE, showlegend = TRUE
  )

# k-offset lines (y = x × k per group)
if (!is.null(scaling_factors)) {
  for (i in seq_len(nrow(scaling_factors))) {
    grp <- scaling_factors$ecotox_group[i]
    k   <- scaling_factors$scaling_factor[i]
    if (!grp %in% names(group_colours)) next

    fig1 <- fig1 %>%
      add_trace(
        x = ax_lim1, y = ax_lim1 * k,
        type = "scatter", mode = "lines",
        name = paste0(grp, " k=", round(k, 3)),
        line = list(color = group_colours[grp], dash = "solid", width = 1),
        inherit = FALSE, showlegend = TRUE
      )
  }
}

fig1 <- fig1 %>%
  layout(
    title  = "Plot 1 — HC5 (SSD-fitted) vs Median LC50 [log–log]",
    xaxis  = list(title = "Median LC50 (ng/L)", type = "log"),
    yaxis  = list(title = "HC5 SSD (ng/L)",     type = "log"),
    legend = list(title = list(text = "Group / line")),
    annotations = list(list(
      text = "Dashed lines = OLS per group  |  Solid coloured = y = x\u00d7k  |  Black dashed = 1:1",
      xref = "paper", yref = "paper", x = 0, y = -0.12,
      showarrow = FALSE, font = list(size = 11, color = "grey40")
    ))
  )

# ── 4. Plot 2 — SSD HC5 vs Scaled HC5 (agreement plot) ───────────────────────
# x: hc5_model_ng_L (scaling factor estimate), y: hc5_ssd_ng_L (full SSD)
# Diagonal = perfect agreement.  Spread reveals calibration uncertainty.

p2_df <- ref %>%
  filter(!is.na(hc5_ssd_ng_L), !is.na(hc5_model_ng_L),
         hc5_ssd_ng_L > 0, hc5_model_ng_L > 0)

ax_range2 <- range(log10(c(p2_df$hc5_model_ng_L, p2_df$hc5_ssd_ng_L)), na.rm = TRUE)
ax_lim2   <- 10^(ax_range2 + c(-0.3, 0.3))

fig2 <- plot_ly()
for (grp in names(group_colours)) {
  d <- filter(p2_df, ecotox_group == grp)
  if (nrow(d) == 0) next
  fig2 <- fig2 %>%
    add_trace(
      data = d,
      x = ~hc5_model_ng_L, y = ~hc5_ssd_ng_L,
      type = "scatter", mode = "markers",
      name = grp,
      text = ~paste0("<b>", chemical_name, "</b><br>CAS: ", cas_number,
                     "<br>HC5 (Scaled): ", signif(hc5_model_ng_L, 3), " ng/L",
                     "<br>HC5 (SSD): ",    signif(hc5_ssd_ng_L,   3), " ng/L",
                     "<br>Ratio (SSD/Scaled): ", round(hc5_ssd_ng_L / hc5_model_ng_L, 2)),
      hoverinfo = "text",
      marker = list(size = 6, color = group_colours[grp], opacity = 0.7)
    )
}
fig2 <- fig2 %>%
  add_trace(
    x = ax_lim2, y = ax_lim2,
    type = "scatter", mode = "lines",
    name = "Perfect agreement",
    line = list(color = "black", dash = "dash", width = 1.2),
    inherit = FALSE, showlegend = TRUE
  ) %>%
  layout(
    title = "Plot 2 — Scaled HC5 vs SSD-fitted HC5 (agreement)",
    xaxis = list(title = "HC5 Scaled (median \u00d7 k) [ng/L]", type = "log"),
    yaxis = list(title = "HC5 SSD-fitted [ng/L]",               type = "log"),
    legend = list(title = list(text = "Taxonomic group"))
  )

# ── 5. Plot 3 — LC50 vs EC50 medians per compound × group ────────────────────
# Computes separate medians for LC50 and EC50 from the raw ecotox file.
# Compounds with BOTH endpoint types are shown.
# Expected result: EC50 ≈ LC50 for acute endpoints (minor systematic offset
# is acceptable; large divergence would question combining the endpoints).

lc50_med <- ecotox %>%
  filter(endpoint == "LC50", !is.na(conc_ng_L), conc_ng_L > 0) %>%
  group_by(cas_number, chemical_name, ecotox_group) %>%
  summarise(median_lc50 = median(conc_ng_L, na.rm = TRUE),
            n_lc50      = n(),
            .groups     = "drop")

ec50_med <- ecotox %>%
  filter(endpoint == "EC50", !is.na(conc_ng_L), conc_ng_L > 0) %>%
  group_by(cas_number, ecotox_group) %>%
  summarise(median_ec50 = median(conc_ng_L, na.rm = TRUE),
            n_ec50      = n(),
            .groups     = "drop")

p3_df <- inner_join(lc50_med, ec50_med, by = c("cas_number", "ecotox_group")) %>%
  filter(median_lc50 > 0, median_ec50 > 0)

message("Plot 3: ", nrow(p3_df),
        " compound×group pairs with both LC50 and EC50 data")

ax_range3 <- range(log10(c(p3_df$median_lc50, p3_df$median_ec50)), na.rm = TRUE)
ax_lim3   <- 10^(ax_range3 + c(-0.5, 0.5))

fig3 <- plot_ly()
for (grp in names(group_colours)) {
  d <- filter(p3_df, ecotox_group == grp)
  if (nrow(d) == 0) next
  fig3 <- fig3 %>%
    add_trace(
      data = d,
      x = ~median_lc50, y = ~median_ec50,
      type = "scatter", mode = "markers",
      name = grp,
      text = ~paste0("<b>", chemical_name, "</b><br>CAS: ", cas_number,
                     "<br>Median LC50: ", signif(median_lc50, 3), " ng/L (n=", n_lc50, ")",
                     "<br>Median EC50: ", signif(median_ec50, 3), " ng/L (n=", n_ec50, ")",
                     "<br>EC50/LC50 ratio: ", round(median_ec50 / median_lc50, 2)),
      hoverinfo = "text",
      marker = list(size = 6, color = group_colours[grp], opacity = 0.7)
    )
}

# ── 5b. Plot 3b — EC50 vs IC50 medians per compound × group ─────────────────
# Same interchangeability check as Plot 3, applied to the EC50/IC50 pair added
# to the pooled median alongside LC50/EC50 (Docs §5.4). IC50 rows only exist
# in `ecotox` once taxotox_install.R's endpoint filter includes "IC50" (i.e.
# after a rebuild) -- fig3b stays NULL on old data so this section is simply
# omitted from the report rather than erroring.

fig3b <- NULL

ic50_med <- ecotox %>%
  filter(endpoint == "IC50", !is.na(conc_ng_L), conc_ng_L > 0) %>%
  group_by(cas_number, ecotox_group) %>%
  summarise(median_ic50 = median(conc_ng_L, na.rm = TRUE),
            n_ic50      = n(),
            .groups     = "drop")

p3b_df <- inner_join(ec50_med, ic50_med, by = c("cas_number", "ecotox_group")) %>%
  filter(median_ec50 > 0, median_ic50 > 0)

if (nrow(p3b_df) > 0) {
  # chemical_name isn't in ec50_med -- pull it back in for hover text
  p3b_df <- p3b_df %>%
    left_join(lc50_med %>% distinct(cas_number, chemical_name), by = "cas_number")

  message("Plot 3b: ", nrow(p3b_df), " compound×group pairs with both EC50 and IC50 data")

  ax_range3b <- range(log10(c(p3b_df$median_ec50, p3b_df$median_ic50)), na.rm = TRUE)
  ax_lim3b   <- 10^(ax_range3b + c(-0.5, 0.5))

  fig3b <- plot_ly()
  for (grp in names(group_colours)) {
    d <- filter(p3b_df, ecotox_group == grp)
    if (nrow(d) == 0) next
    fig3b <- fig3b %>%
      add_trace(
        data = d,
        x = ~median_ec50, y = ~median_ic50,
        type = "scatter", mode = "markers",
        name = grp,
        text = ~paste0("<b>", chemical_name, "</b><br>CAS: ", cas_number,
                       "<br>Median EC50: ", signif(median_ec50, 3), " ng/L (n=", n_ec50, ")",
                       "<br>Median IC50: ", signif(median_ic50, 3), " ng/L (n=", n_ic50, ")",
                       "<br>IC50/EC50 ratio: ", round(median_ic50 / median_ec50, 2)),
        hoverinfo = "text",
        marker = list(size = 6, color = group_colours[grp], opacity = 0.7)
      )
  }

  pair_counts3b <- p3b_df %>% count(ecotox_group) %>%
    mutate(label = paste0(ecotox_group, ": n=", n)) %>%
    pull(label) %>% paste(collapse = "  |  ")

  fig3b <- fig3b %>%
    add_trace(
      x = ax_lim3b, y = ax_lim3b,
      type = "scatter", mode = "lines",
      name = "EC50 = IC50",
      line = list(color = "black", dash = "dash", width = 1.2),
      inherit = FALSE, showlegend = TRUE
    ) %>%
    layout(
      title = "Plot 3b — Median EC50 vs Median IC50 per compound × group",
      xaxis = list(title = "Median EC50 (ng/L)", type = "log"),
      yaxis = list(title = "Median IC50 (ng/L)", type = "log"),
      legend = list(title = list(text = "Taxonomic group")),
      annotations = list(list(
        text = pair_counts3b,
        xref = "paper", yref = "paper", x = 0, y = -0.12,
        showarrow = FALSE, font = list(size = 11, color = "grey40")
      ))
    )
}

# Count pairs per group for annotation
pair_counts <- p3_df %>% count(ecotox_group) %>%
  mutate(label = paste0(ecotox_group, ": n=", n)) %>%
  pull(label) %>% paste(collapse = "  |  ")

fig3 <- fig3 %>%
  add_trace(
    x = ax_lim3, y = ax_lim3,
    type = "scatter", mode = "lines",
    name = "LC50 = EC50",
    line = list(color = "black", dash = "dash", width = 1.2),
    inherit = FALSE, showlegend = TRUE
  ) %>%
  layout(
    title = "Plot 3 — Median LC50 vs Median EC50 per compound \u00d7 group",
    xaxis = list(title = "Median LC50 (ng/L)", type = "log"),
    yaxis = list(title = "Median EC50 (ng/L)", type = "log"),
    legend = list(title = list(text = "Taxonomic group")),
    annotations = list(list(
      text = pair_counts,
      xref = "paper", yref = "paper", x = 0, y = -0.12,
      showarrow = FALSE, font = list(size = 11, color = "grey40")
    ))
  )

# ── 6. Plot 4 — Distribution of k (HC5/median LC50 ratio) ────────────────────
# Histogram of hc5_ssd_ng_L / median_lc50_ng_L, one panel per ecotox_group.
# Vertical line at group median (= k used in scaling factor).

ratio_df <- ref %>%
  filter(!is.na(hc5_ssd_ng_L), !is.na(median_lc50_ng_L),
         median_lc50_ng_L > 0, hc5_ssd_ng_L > 0) %>%
  mutate(ratio = hc5_ssd_ng_L / median_lc50_ng_L)

fig4 <- plot_ly()
for (grp in names(group_colours)) {
  d <- filter(ratio_df, ecotox_group == grp)
  if (nrow(d) == 0) next
  k_val <- median(d$ratio, na.rm = TRUE)
  fig4 <- fig4 %>%
    add_trace(
      x = ~ratio, data = d,
      type = "histogram",
      name = grp,
      marker = list(color = group_colours[grp], opacity = 0.65, line = list(width = 0.5)),
      nbinsx = 40,
      showlegend = TRUE
    ) %>%
    add_segments(
      x = k_val, xend = k_val, y = 0, yend = 50,
      type = "scatter", mode = "lines",
      line = list(color = group_colours[grp], dash = "dash", width = 2),
      name = paste0(grp, " k=", round(k_val, 3)),
      inherit = FALSE, showlegend = TRUE
    )
}
fig4 <- fig4 %>%
  layout(
    title      = "Plot 4 — Distribution of k = HC5 (SSD) / Median LC50 per group",
    xaxis      = list(title = "k ratio (HC5 / Median LC50)"),
    yaxis      = list(title = "Count"),
    barmode    = "overlay",
    legend     = list(title = list(text = "Group (dashed = k used)")),
    annotations = list(list(
      text = "Wide spread \u2192 higher uncertainty in Scaled HC5 estimates",
      xref = "paper", yref = "paper", x = 0, y = -0.12,
      showarrow = FALSE, font = list(size = 11, color = "grey40")
    ))
  )

# ── 7. Plot 5 — OPERA predicted vs ECOTOX median (calibration check) ─────────
# Only shown if I-2 has been run (predicted_lc50_ng_L column exists and is populated).
# Compounds where both values exist are plotted; the 1:1 line is the ideal.
# Systematic deviation indicates model bias; scatter quantifies prediction uncertainty.
# This plot is critical for deciding whether predicted values are reliable enough
# to use as TU denominators.

fig5 <- NULL
if ("predicted_lc50_ng_L" %in% names(ref)) {
  p5_df <- ref %>%
    filter(!is.na(predicted_lc50_ng_L), !is.na(median_lc50_ng_L),
           predicted_lc50_ng_L > 0, median_lc50_ng_L > 0,
           ecotox_group %in% c("fish", "crustacean"))  # algae has no prediction

  message("Plot 5: ", nrow(p5_df),
          " compound×group pairs with both OPERA prediction and ECOTOX median")

  if (nrow(p5_df) > 0) {
    ax_range5 <- range(log10(c(p5_df$predicted_lc50_ng_L, p5_df$median_lc50_ng_L)),
                       na.rm = TRUE)
    ax_lim5   <- 10^(ax_range5 + c(-0.5, 0.5))

    fig5 <- plot_ly()
    for (grp in c("fish", "crustacean")) {
      d <- filter(p5_df, ecotox_group == grp)
      if (nrow(d) == 0) next

      # RMSE on log10 scale as a performance metric
      rmse_log <- sqrt(mean((log10(d$predicted_lc50_ng_L) - log10(d$median_lc50_ng_L))^2,
                            na.rm = TRUE))

      fig5 <- fig5 %>%
        add_trace(
          data = d,
          x = ~median_lc50_ng_L, y = ~predicted_lc50_ng_L,
          type = "scatter", mode = "markers",
          name = paste0(grp, " (RMSE\u2081\u2080=", round(rmse_log, 2), ")"),
          text = ~paste0("<b>", chemical_name, "</b><br>CAS: ", cas_number,
                         "<br>Group: ", ecotox_group,
                         "<br>ECOTOX median: ",  signif(median_lc50_ng_L,    3), " ng/L",
                         "<br>OPERA predicted: ", signif(predicted_lc50_ng_L, 3), " ng/L",
                         "<br>Ratio (pred/obs): ", round(predicted_lc50_ng_L / median_lc50_ng_L, 2)),
          hoverinfo = "text",
          marker = list(size = 6, color = group_colours[grp], opacity = 0.7)
        )
    }
    fig5 <- fig5 %>%
      add_trace(
        x = ax_lim5, y = ax_lim5,
        type = "scatter", mode = "lines",
        name = "1:1 line (perfect prediction)",
        line = list(color = "black", dash = "dash", width = 1.2),
        inherit = FALSE, showlegend = TRUE
      ) %>%
      layout(
        title = "Plot 5 \u2014 OPERA Predicted vs ECOTOX Median LC50 [log\u2013log]",
        xaxis = list(title = "ECOTOX Median LC50 (ng/L) [observed]",  type = "log"),
        yaxis = list(title = "OPERA Predicted LC50 (ng/L) [predicted]", type = "log"),
        legend = list(title = list(text = "Group (RMSE on log\u2081\u2080 scale)")),
        annotations = list(list(
          text = paste0("RMSE\u2081\u2080 < 0.5 log-units = good agreement  |  ",
                        "RMSE\u2081\u2080 0.5\u20131.0 = acceptable  |  ",
                        "> 1.0 = use with caution"),
          xref = "paper", yref = "paper", x = 0, y = -0.12,
          showarrow = FALSE, font = list(size = 11, color = "grey40")
        ))
      )
  }
}

# ── 8. Coverage by lc50_source ────────────────────────────────────────────────
# Bar chart: how many compound×group rows fall into each lc50_source category.
# Shows the net gain from the CompTox gap-fill step.

fig6 <- NULL
if ("lc50_source" %in% names(ref)) {
  coverage_df <- ref %>%
    mutate(lc50_source = replace_na(as.character(lc50_source), "No denominator")) %>%
    count(ecotox_group, lc50_source) %>%
    mutate(lc50_source = factor(lc50_source,
             levels = c("ECOTOX_median", "Both", "CompTox_predicted", "No denominator")))

  source_colours <- c(
    "ECOTOX_median"     = "#1f78b4",
    "Both"              = "#33a02c",
    "CompTox_predicted" = "#ff7f00",
    "No denominator"    = "#cccccc"
  )

  fig6 <- plot_ly(coverage_df,
      x = ~ecotox_group, y = ~n, color = ~lc50_source,
      colors = source_colours,
      type = "bar",
      text = ~paste0(lc50_source, ": ", n, " compounds"),
      hoverinfo = "text") %>%
    layout(
      title  = "Plot 6 \u2014 Denominator Coverage by Source and Taxonomic Group",
      xaxis  = list(title = "Taxonomic group"),
      yaxis  = list(title = "Number of compound\u00d7group rows"),
      barmode = "stack",
      legend = list(title = list(text = "LC50 source")),
      annotations = list(list(
        text = paste0("Orange = CompTox-only rows added by gap-fill  |  ",
                      "Green = compounds present in both ECOTOX and CompTox  |  ",
                      "Grey = no denominator available"),
        xref = "paper", yref = "paper", x = 0, y = -0.12,
        showarrow = FALSE, font = list(size = 11, color = "grey40")
      ))
    )
}

# ── 9. Summary table ──────────────────────────────────────────────────────────

summary_tbl <- ref %>%
  group_by(ecotox_group) %>%
  summarise(
    n_rows         = n(),
    n_ecotox_data  = sum(!is.na(median_lc50_ng_L)),
    n_SSD          = sum(hc5_method == "SSD",    na.rm = TRUE),
    n_Both         = sum(hc5_method == "Both",   na.rm = TRUE),
    n_Scaled_only  = sum(hc5_method == "Scaled", na.rm = TRUE),
    pct_SSD        = round(100 * (n_SSD + n_Both) / n_rows, 1),
    .groups = "drop"
  )

# Add n_ecotox data-richness if available
if ("n_ecotox" %in% names(ref)) {
  richness <- ref %>%
    group_by(ecotox_group) %>%
    summarise(
      median_n_ecotox = round(median(n_ecotox, na.rm = TRUE), 1),
      n_ge5           = sum(n_ecotox >= 5, na.rm = TRUE),
      .groups = "drop"
    )
  summary_tbl <- left_join(summary_tbl, richness, by = "ecotox_group")
}

# Add CompTox coverage if I-2 has run
if ("predicted_lc50_ng_L" %in% names(ref)) {
  comptox_cov <- ref %>%
    group_by(ecotox_group) %>%
    summarise(n_predicted = sum(!is.na(predicted_lc50_ng_L)), .groups = "drop")
  summary_tbl <- left_join(summary_tbl, comptox_cov, by = "ecotox_group")
}

# Add MoA coverage if I-3 has run
if ("moa_group" %in% names(ref)) {
  moa_cov <- ref %>%
    group_by(ecotox_group) %>%
    summarise(
      # "Specific" = has a real Kramer et al. (2024) classification, i.e. not
      # the "unknown" catch-all (see Code/taxotox_install_moa.R).
      n_moa_specific    = sum(!is.na(moa_group) & moa_group != "unknown", na.rm = TRUE),
      # moa_group is guaranteed non-NA (always "unknown" rather than NA for
      # unclassified compounds -- see Code/taxotox_install_moa.R:144), so a
      # is.na() count would always read 0. Count the "unknown" label itself.
      n_moa_unknown     = sum(moa_group == "unknown", na.rm = TRUE),
      .groups = "drop"
    )
  summary_tbl <- left_join(summary_tbl, moa_cov, by = "ecotox_group")
}

dt_opts <- list(
  pageLength = 20, scrollX = TRUE, autoWidth = FALSE, dom = "t",
  initComplete = htmlwidgets::JS(
    "function(settings, json) { $(this.api().table().header()).css({'white-space':'nowrap'}); }"
  )
)

dt_summary <- datatable(
  summary_tbl, rownames = FALSE, options = dt_opts,
  caption = paste0(
    "Reference table coverage summary by taxonomic group.  ",
    "n_rows = total compound\u00d7group rows; ",
    "n_ecotox_data = rows with ECOTOX median LC50; ",
    "SSD = full Species Sensitivity Distribution fitted (n\u22655); ",
    "Both = SSD + scaling factor both computed; ",
    "Scaled = scaling factor only; ",
    "n_predicted = CompTox/OPERA predictions (fish + crustacean only); ",
    "n_moa_specific = rows with a real, curated Mode of Action classification; ",
    "n_moa_unknown = rows with no MoA classification available (Kramer et al. either ",
    "could not classify the compound or does not cover it).")
)

# Top mismatches between SSD and Scaled (for audit)
if (all(c("hc5_ssd_ng_L", "hc5_model_ng_L") %in% names(ref))) {
  mismatch_df <- ref %>%
    filter(!is.na(hc5_ssd_ng_L), !is.na(hc5_model_ng_L),
           hc5_ssd_ng_L > 0, hc5_model_ng_L > 0) %>%
    mutate(
      ratio_ssd_scaled = round(hc5_ssd_ng_L / hc5_model_ng_L, 3),
      log2_ratio       = round(log2(ratio_ssd_scaled), 3)
    ) %>%
    arrange(desc(abs(log2_ratio))) %>%
    select(cas_number, chemical_name, ecotox_group, n_ecotox,
           median_lc50_ng_L, hc5_ssd_ng_L, hc5_model_ng_L,
           ratio_ssd_scaled, log2_ratio) %>%
    mutate(across(c(median_lc50_ng_L, hc5_ssd_ng_L, hc5_model_ng_L), ~signif(.x, 4)))

  dt_mismatch <- datatable(
    mismatch_df, rownames = FALSE, filter = "top",
    options = list(pageLength = 25, scrollX = TRUE, autoWidth = TRUE),
    caption = "SSD vs Scaled HC5 — sorted by |log2(ratio)| (largest mismatch first)"
  )
} else {
  dt_mismatch <- NULL
}

# ── 8b. Plot 7 — Mode of Action coverage (if I-3 ran) ────────────────────────
# Stacked bar: number of unique compounds per moa_group × ecotox_group,
# filled by moa_source to show how each group was classified.

fig7 <- NULL
moa_summary_tbl <- NULL

if ("moa_group" %in% names(ref)) {

  # moa_group values come from Kramer et al. (2024) (Code/taxotox_install_moa.R)
  # and are data-driven, not a fixed list -- derive levels from what's actually
  # present, most-common first, so the stacked bar and legend stay correct as
  # the source data's category set changes. "unknown" (compounds with no MoA
  # classification, per Code/taxotox_install_moa.R) is pushed to the end.
  moa_group_levels <- ref %>%
    mutate(moa_group = replace_na(as.character(moa_group), "unknown")) %>%
    count(moa_group, sort = TRUE) %>%
    pull(moa_group)
  moa_group_levels <- c(setdiff(moa_group_levels, "unknown"),
                        intersect(moa_group_levels, "unknown"))

  moa_df <- ref %>%
    mutate(moa_group  = replace_na(as.character(moa_group),  "unknown"),
           moa_source = replace_na(as.character(moa_source), "not_in_Kramer2024")) %>%
    count(ecotox_group, moa_group, moa_source) %>%
    mutate(moa_group = factor(moa_group, levels = moa_group_levels))

  # moa_group is guaranteed non-NA (see Code/taxotox_install_moa.R:144); an
  # is.na() count here would always read 0, so count the "unknown" label.
  n_no_moa <- sum(ref$moa_group == "unknown", na.rm = TRUE)

  # No manual colour palette -- moa_group is now an open-ended, data-driven set
  # rather than a fixed 5-value list, so let plotly auto-cycle trace colours.
  fig7 <- plot_ly()
  for (grp in moa_group_levels) {
    d <- filter(moa_df, moa_group == grp)
    if (nrow(d) == 0) next
    fig7 <- fig7 %>%
      add_trace(
        data = d,
        x = ~ecotox_group, y = ~n,
        type = "bar", name = grp,
        text = ~paste0("<b>", moa_group, "</b><br>",
                       "Classified by: ", moa_source, "<br>",
                       "N rows: ", n),
        hoverinfo = "text"
      )
  }
  fig7 <- fig7 %>%
    layout(
      title   = "Plot 7 — Mode of Action Group Coverage by Taxonomic Group",
      xaxis   = list(title = "Taxonomic group"),
      yaxis   = list(title = "Number of compound×group rows"),
      barmode = "stack",
      legend  = list(title = list(text = "Mode of action group")),
      annotations = list(list(
        text = paste0(
          "Hover bars for classification source  |  ",
          "Rows without Mode of Action assignment (shown as unknown): ", n_no_moa),
        xref = "paper", yref = "paper", x = 0, y = -0.14,
        showarrow = FALSE, font = list(size = 11, color = "grey40")
      ))
    )

  # Summary table: compounds per moa_group × ecotox_group × source
  moa_summary_tbl <- ref %>%
    mutate(moa_group  = replace_na(as.character(moa_group),  "unknown"),
           moa_source = replace_na(as.character(moa_source), "not_in_Kramer2024")) %>%
    count(ecotox_group, moa_group, moa_source) %>%
    arrange(ecotox_group, desc(n)) %>%
    rename(
      `Taxonomic group`       = ecotox_group,
      `Mode of action group`  = moa_group,
      `Classification source` = moa_source,
      `N rows`                = n
    )
}

# ── 8c. Benchmark coverage ──────────────────────────────────────────────────
# Formerly Plot 8: a horizontal stacked bar comparing coverage across four
# national benchmark frameworks (US EPA, EU EQS, AU ANZG, CA CCME), sorted by
# how many frameworks covered each compound. EU EQS / AU ANZG / CA CCME were
# removed from Code/taxotox_install_benchmarks.R due to low compound coverage
# (~1-3% each vs. US EPA's ~9-10%) -- with only one framework left, a
# "coverage per framework" comparison plot no longer says anything a single
# coverage percentage doesn't already say more simply. taxotox_install_
# benchmarks.R prints that percentage directly in its own console output;
# fig8 stays NULL so the (still framework-agnostic) "Benchmarks: YES/NOT RUN"
# banner further down and the mismatch-table section numbering keep working
# unchanged.

fig8 <- NULL

bm_check_cols <- c("benchmark_usepa_fish_acute_ng_L",
                   "benchmark_usepa_crust_acute_ng_L",
                   "benchmark_usepa_algae_acute_ng_L")

# ── 8d. Plot 9 — CAMA vs Standard PTI convergence diagnostic ────────────────
# Checks a specific concern: CAMA's E_mix tends to correlate almost perfectly
# (rho ~ 1) with standard TU-sum PTI in practice, which can look suspicious --
# as if the Mode of Action grouping isn't doing anything. It isn't a bug: the
# closed form obtained by substituting Step 2 into Step 3 of the CAMA formula
# (Docs/TaxoTox_Technical_Methods.md §10.5.2) is
#   E_mix = 1 - exp(-Σ_g log(1 + group_TU_g))
# which reduces to 1 - exp(-PTI) ≈ PTI whenever every group_TU_g ≪ 1 --
# true for nearly every real water sample -- regardless of how many distinct
# MoA groups exist. This section measures that directly against the repo's
# own sample monitoring workbooks (Data/sample_*.xlsx) rather than relying on
# it being noticed by eye in the live app.
#
# Compound-name -> CASRN matching here replicates only the Known_CAS exact-
# match layer of app.R's CASRN pipeline (app.R:997-1006) -- sufficient for
# these curated sample files, which were built to resolve via that layer
# alone; it is not a substitute for the app's full PubChem/manual-entry flow.
#
# Two variants are computed per taxon, both with PTI and CAMA restricted to
# the SAME compound subset so the comparison isolates the CA-vs-CAMA formula
# effect from the "unknown"-coverage effect already covered by Plot 7/8b:
#   "all compounds"   -- what the app's own PTI and CAMA "all" sheets show
#   "known MoA only"  -- PTI and CAMA both restricted to compounds with a
#                        real Kramer et al. (2024) classification

fig9 <- NULL
cama_diag_tbl <- NULL

sample_files <- list.files("../Data", pattern = "^sample_.*\\.xlsx$", full.names = TRUE)

if ("moa_group" %in% names(ref) && length(sample_files) > 0 &&
    file.exists("../Data/Known_CAS.fst")) {

  known_cas <- read.fst("../Data/Known_CAS.fst", as.data.table = FALSE)

  # Mirrors app.R's .calc_tox() (app.R:1374-1386) and .calc_cama()
  # (app.R:1586-1621), computed together on the same compound subset so PTI
  # and E_mix are directly comparable row-for-row.
  .calc_pti_cama <- function(denom_tbl, user_data, exclude_unknown) {
    if (exclude_unknown) denom_tbl <- denom_tbl %>% filter(moa_group != "unknown")
    if (nrow(denom_tbl) == 0) return(NULL)

    conc_long <- denom_tbl %>%
      select(PREFERRED_NAME, moa_group, median_conc) %>%
      left_join(user_data, by = "PREFERRED_NAME") %>%
      pivot_longer(cols = -c(PREFERRED_NAME, moa_group, median_conc),
                   names_to = "Sample", values_to = "C") %>%
      mutate(TU = if_else(is.na(C) | is.na(median_conc), 0, C / median_conc))

    pti <- conc_long %>%
      group_by(Sample) %>%
      summarise(PTI = sum(TU, na.rm = TRUE), .groups = "drop")

    e_group <- conc_long %>%
      group_by(Sample, moa_group) %>%
      summarise(group_TU = sum(TU, na.rm = TRUE), .groups = "drop") %>%
      mutate(E_group = group_TU / (1 + group_TU))

    e_mix <- e_group %>%
      group_by(Sample) %>%
      summarise(E_mix        = 1 - prod(1 - E_group),
                max_group_TU = max(group_TU),
                n_groups     = n_distinct(moa_group[group_TU > 0]),
                .groups = "drop")

    left_join(pti, e_mix, by = "Sample")
  }

  cama_diag_rows <- list()

  for (f in sample_files) {
    d <- tryCatch(readxl::read_excel(f), error = function(e) NULL)
    if (is.null(d) || ncol(d) < 2) next
    names(d)[1] <- "PREFERRED_NAME"
    d <- d %>% mutate(across(-PREFERRED_NAME, ~suppressWarnings(as.numeric(.x))))

    matched <- known_cas[known_cas$PREFERRED_NAME %in% d$PREFERRED_NAME, ] %>%
      distinct(PREFERRED_NAME, .keep_all = TRUE) %>%
      mutate(cas_number = gsub("-", "", CASRN)) %>%
      select(PREFERRED_NAME, cas_number)
    if (nrow(matched) == 0) next

    ref_matched <- ref %>%
      filter(cas_number %in% matched$cas_number) %>%
      mutate(median_conc = dplyr::coalesce(median_lc50_ng_L, predicted_lc50_ng_L),
             moa_group   = if_else(is.na(moa_group) | moa_group == "", "unknown", moa_group)) %>%
      filter(!is.na(median_conc)) %>%
      # relationship = "many-to-many": a sample file can legitimately list two
      # synonym rows for the same CASRN (see Docs §11 "Synonym handling and
      # duplicate CASRN detection") -- app.R's own ref_matched join
      # (app.R:1363) has the same shape and accepts the same double-counting.
      left_join(matched, by = "cas_number", relationship = "many-to-many") %>%
      filter(!is.na(PREFERRED_NAME)) %>%
      select(PREFERRED_NAME, cas_number, ecotox_group, median_conc, moa_group)

    for (grp in c("algae", "crustacean", "fish")) {
      denom <- ref_matched %>% filter(ecotox_group == grp)
      if (nrow(denom) == 0) next

      res_all   <- .calc_pti_cama(denom, d, exclude_unknown = FALSE)
      res_known <- .calc_pti_cama(denom, d, exclude_unknown = TRUE)

      if (!is.null(res_all))
        cama_diag_rows[[length(cama_diag_rows) + 1]] <-
          res_all %>% mutate(dataset = basename(f), ecotox_group = grp, variant = "all compounds")
      if (!is.null(res_known))
        cama_diag_rows[[length(cama_diag_rows) + 1]] <-
          res_known %>% mutate(dataset = basename(f), ecotox_group = grp, variant = "known MoA only")
    }
  }

  if (length(cama_diag_rows) > 0) {
    # rho is undefined for all-zero (no-detect) samples -- exclude them.
    cama_diag_long <- bind_rows(cama_diag_rows) %>% filter(PTI > 0)

    if (nrow(cama_diag_long) > 0) {
      cama_diag_summary <- cama_diag_long %>%
        group_by(variant, ecotox_group) %>%
        summarise(
          n_samples       = n(),
          spearman_rho    = suppressWarnings(cor(PTI, E_mix, method = "spearman")),
          pearson_r       = suppressWarnings(cor(PTI, E_mix, method = "pearson")),
          max_abs_diff    = max(abs(E_mix - PTI / (1 + PTI))),
          n_high_group_TU = sum(max_group_TU >= 0.3),
          .groups = "drop"
        ) %>%
        arrange(variant, ecotox_group)

      cama_diag_tbl <- datatable(
        cama_diag_summary %>%
          mutate(across(c(spearman_rho, pearson_r, max_abs_diff), ~round(.x, 4))) %>%
          rename(Variant                              = variant,
                 `Taxonomic group`                     = ecotox_group,
                 `N samples`                           = n_samples,
                 `Spearman rho (E_mix vs PTI)`         = spearman_rho,
                 `Pearson r`                            = pearson_r,
                 `Max |E_mix - PTI/(1+PTI)|`           = max_abs_diff,
                 `N samples with a group TU ≥ 0.3` = n_high_group_TU),
        rownames = FALSE, options = list(pageLength = 10, dom = "t", scrollX = TRUE),
        caption = "CAMA E_mix vs standard PTI — pooled across all Data/sample_*.xlsx workbooks"
      )

      pti_range <- c(0, max(cama_diag_long$PTI, na.rm = TRUE))
      limit_line <- seq(pti_range[1], pti_range[2], length.out = 100)

      fig9 <- plot_ly(
        data = cama_diag_long %>% filter(variant == "known MoA only"),
        x = ~PTI, y = ~E_mix, color = ~ecotox_group, type = "scatter", mode = "markers",
        text = ~paste0(dataset, " | ", Sample,
                       "<br>max group TU: ", round(max_group_TU, 3),
                       "<br>n groups with TU>0: ", n_groups),
        hoverinfo = "text+x+y"
      ) %>%
        add_trace(
          data = data.frame(x = limit_line, y = limit_line / (1 + limit_line)),
          x = ~x, y = ~y, inherit = FALSE,
          type = "scatter", mode = "lines", name = "PTI/(1+PTI) (single-group limit)",
          line = list(dash = "dash", color = "black"), showlegend = TRUE, hoverinfo = "skip"
        ) %>%
        layout(
          title  = "Plot 9 — CAMA E_mix (known-MoA variant) vs Standard PTI (same subset)",
          xaxis  = list(title = "Standard PTI (Σ TU, known-MoA compound subset)"),
          yaxis  = list(title = "CAMA E_mix (known-MoA variant)"),
          legend = list(title = list(text = "Taxonomic group"))
        )
    }
  }
}

# ── 9. Build HTML report ──────────────────────────────────────────────────────

css <- "
  body  { font-family: Arial, sans-serif; margin: 30px; color: #333; max-width: 1400px; }
  h1    { color: #2c5f8a; border-bottom: 2px solid #2c5f8a; padding-bottom: 6px; }
  h2    { color: #4a4a4a; margin-top: 40px; }
  p.meta { color: #777; font-size: 0.9em; }
  .section { margin-bottom: 50px; }
  .note { background: #f5f5f5; border-left: 4px solid #2c5f8a;
          padding: 10px 14px; margin: 14px 0; font-size: 0.93em; }
"

make_section <- function(...) tags$div(class = "section", ...)

report <- tagList(
  tags$html(lang = "en",
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$title("TaxoTox Reference Validation"),
      tags$style(HTML(css))
    ),
    tags$body(
      tags$h1("TaxoTox Reference Validation — Diagnostic Report"),
      tags$p(class = "meta",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"),
               "  |  ECOTOX rows: ", nrow(ecotox),
               "  |  Reference rows: ", nrow(ref),
               "  |  Reference source: ",
               if (file.exists(ref_path)) "taxotox_reference.fst" else "reconstructed inline",
               if ("predicted_lc50_ng_L" %in% names(ref))
                 paste0("  |  CompTox gap-fill: YES (",
                        sum(!is.na(ref$predicted_lc50_ng_L)), " predictions)")
               else "  |  CompTox gap-fill: NOT RUN",
               if ("moa_group" %in% names(ref))
                 paste0("  |  Mode of Action: YES (",
                        sum(ref$moa_group != "unknown", na.rm = TRUE), " / ", nrow(ref),
                        " rows have a real classification)")
               else "  |  Mode of Action: NOT RUN",
               if (any(bm_check_cols %in% names(ref)))
                 paste0("  |  Benchmarks: YES")
               else "  |  Benchmarks: NOT RUN")),

      tags$p(
        "This report validates the ", tags$code("taxotox_reference.fst"), " table produced by ",
        tags$code("taxotox_install.R"), ". It covers four areas: ",
        tags$b("(1)"), " correctness of the HC5 computation pipeline (Species Sensitivity Distribution fitting and empirical scaling); ",
        tags$b("(2)"), " the assumption that LC50 and EC50 are interchangeable for Toxic Unit calculations; ",
        tags$b("(3)"), " the quality of CompTox/OPERA gap-fill predictions relative to measured ECOTOX data; and ",
        tags$b("(4)"), " coverage of the Mode of Action classification used by the CAMA method. ",
        "All plots are interactive — hover for compound details, click legend entries to isolate groups."
      ),

      # Summary table
      make_section(
        tags$h2("1. HC5 Coverage Summary"),
        tags$p(class = "note",
          tags$b("What this shows: "), "a single summary of coverage across all install steps. ",
          tags$b("n_ecotox_data"), " = rows with a measured ECOTOX median (the primary TU denominator). ",
          tags$b("SSD / Both / Scaled"), " = HC5 method coverage: SSD = full log-normal Species Sensitivity Distribution ",
          "fitted from ECOTOX data (n \u2265 5); Both = SSD + scaling factor computed; Scaled = scaling factor only. ",
          tags$b("n_predicted"), " = CompTox/OPERA predictions available (fish and crustacean only; algae has no suitable endpoint). ",
          tags$b("n_moa_assigned"), " = rows with a Mode of Action group assigned. ",
          tags$b("n_moa_specific"), " = rows with a real Kramer et al. (2024) Mode of Action classification ",
          "(e.g. AChE inhibition, Photosynthesis inhibition) rather than the “unknown” catch-all. ",
          "Columns not shown are absent because that install step has not been run yet."),
        as.tags(dt_summary)
      ),

      # Plot 1
      make_section(
        tags$h2("2. HC5 (SSD-fitted) vs Median LC50 [log\u2013log]"),
        tags$p(class = "note",
          tags$b("What this shows: "), "the relationship between the SSD-derived HC5 and the median LC50 ",
          "for all data-rich compound\u00d7group pairs (n \u2265 5). ",
          tags$b("Expected pattern: "), "a tight log-linear cloud with slope \u2248 1, ",
          "offset below the 1:1 line by log\u2081\u2080(k). ",
          "The solid coloured lines show the theoretical position y = x \u00d7 k for each group; ",
          "the dotted lines are OLS fits to the actual data. ",
          tags$b("What to look for: "),
          "OLS slope close to 1 (proportional scaling), ",
          "OLS intercept close to the k-line (scaling factor is unbiased), ",
          "low scatter around the fit (k is a good predictor). ",
          "Point size is proportional to n_ecotox."),
        as.tags(fig1)
      ),

      # Plot 2
      make_section(
        tags$h2("3. Scaled HC5 vs SSD-fitted HC5 (agreement)"),
        tags$p(class = "note",
          tags$b("What this shows: "), "direct comparison between the two HC5 estimates for all data-rich rows. ",
          "x-axis = hc5_model_ng_L (median \u00d7 k); y-axis = hc5_ssd_ng_L (full SSD). ",
          tags$b("Expected pattern: "), "points along the 1:1 diagonal. ",
          tags$b("What to look for: "),
          "systematic offset above or below the diagonal indicates that k over- or under-estimates HC5 ",
          "at the extremes of the LC50 range (heteroscedasticity). ",
          "Wide scatter quantifies the uncertainty in using the scaling factor for data-poor compounds. ",
          "This plot is empty if taxotox_reference.fst was not yet updated with the scaling factor (re-run install)."),
        as.tags(fig2)
      ),

      # Plot 3
      make_section(
        tags$h2("4. Median LC50 vs Median EC50 (endpoint interchangeability)"),
        tags$p(class = "note",
          tags$b("What this shows: "), "for compounds where both LC50 and EC50 data exist in ECOTOX, ",
          "this plots their separate medians against each other. ",
          tags$b("Why it matters: "), "TaxoTox combines LC50 and EC50 into a single median before computing TU. ",
          "This is standard practice (Backhaus & Faust, 2012) but rests on the assumption that the two ",
          "endpoint types produce similar concentration estimates for acute aquatic toxicity. ",
          tags$b("Expected pattern: "), "points near the 1:1 line. ",
          tags$b("What to look for: "),
          "systematic divergence (e.g. EC50 consistently lower than LC50 for a group) would indicate ",
          "the endpoints should be kept separate in the aggregation step."),
        as.tags(fig3)
      ),

      # Plot 3b — EC50 vs IC50 (conditional: only after a rebuild with IC50 in the filter)
      if (!is.null(fig3b)) {
        make_section(
          tags$h2("4b. Median EC50 vs Median IC50 (endpoint interchangeability)"),
          tags$p(class = "note",
            tags$b("What this shows: "), "for compounds where both EC50 and IC50 data exist in ECOTOX, ",
            "this plots their separate medians against each other — almost entirely an algae comparison, ",
            "since IC50 is rarely reported for fish/crustacean. ",
            tags$b("Why it matters: "), "TaxoTox pools IC50 into the same median as LC50/EC50 (Docs §5.4) ",
            "on the grounds that ECOTOX's IC50 rows share the same effect/measurement profile as its EC50 ",
            "rows for algae (population growth rate, biomass, chlorophyll, photosynthesis) — different ",
            "source studies label the identical growth-inhibition assay differently. ",
            tags$b("Expected pattern: "), "points near the 1:1 line. ",
            tags$b("What to look for: "),
            "systematic divergence would indicate IC50 should be kept separate from EC50 in the ",
            "aggregation step after all."),
          as.tags(fig3b)
        )
      } else tags$div(),

      # Plot 4
      make_section(
        tags$h2("5. Distribution of k = HC5 (SSD) / Median LC50"),
        tags$p(class = "note",
          tags$b("What this shows: "), "the distribution of the HC5/median ratio across all ",
          "data-rich compounds, separately per taxonomic group. ",
          "The dashed vertical line marks the group median (the k value used in the scaling factor). ",
          tags$b("Why it matters: "), "k is assumed constant within each group. ",
          "If the distribution is narrow and unimodal, this assumption holds well. ",
          tags$b("What to look for: "),
          "a narrow peak centred on the median line = low uncertainty in scaled HC5 estimates; ",
          "a wide or bimodal distribution = k is a poor predictor for some compound classes ",
          "and scaled HC5 estimates carry substantial uncertainty."),
        as.tags(fig4)
      ),

      # Plot 5 — CompTox validation (conditional)
      if (!is.null(fig5)) {
        make_section(
          tags$h2("6. OPERA Predicted vs ECOTOX Median LC50 (CompTox gap-fill calibration)"),
          tags$p(class = "note",
            tags$b("What this shows: "), "for compounds with both a measured ECOTOX median and an OPERA ",
            "model prediction, this plot compares the two on a log-log scale. ",
            tags$b("Why it matters: "), "the CompTox/OPERA predictions are used as TU denominators for ",
            "compounds absent from ECOTOX. This plot calibrates our confidence in those predictions ",
            "by testing them against measured data for compounds where we know the answer. ",
            tags$b("Expected pattern: "), "points near the 1:1 line. ",
            tags$b("Performance benchmark: "),
            "RMSE\u2081\u2080 < 0.5 log-units = good (within ~3\u00d7); ",
            "0.5\u20131.0 = acceptable for screening; > 1.0 = use with caution. ",
            "OPERA aquatic toxicity models typically achieve RMSE \u2248 0.6\u20130.8 log-units on external test sets ",
            "(Mansouri et al., 2018)."),
          as.tags(fig5)
        )
      } else tags$div(),

      # Plot 6 — coverage (conditional)
      if (!is.null(fig6)) {
        make_section(
          tags$h2("7. Denominator Coverage by Source"),
          tags$p(class = "note",
            tags$b("What this shows: "), "how many compound\u00d7group rows have a usable TU denominator, ",
            "broken down by data source. ",
            tags$b("Blue"), " = ECOTOX median only (most reliable). ",
            tags$b("Green"), " = compound present in both ECOTOX and CompTox (can be cross-validated). ",
            tags$b("Orange"), " = CompTox prediction only, no ECOTOX data (gap-fill rows, n_ecotox = 0). ",
            tags$b("Grey"), " = no denominator available (algae gap-fill not possible; some CAS not in CompTox). ",
            tags$b("Why it matters: "), "orange rows represent the net gain from the I-2 gap-fill step; ",
            "compounds in grey will still produce no TU contribution in the app."),
          as.tags(fig6)
        )
      } else tags$div(),

      # Plot 7 — Mode of Action coverage (conditional)
      if (!is.null(fig7)) {
        make_section(
          tags$h2("8. Mode of Action Group Coverage"),
          tags$p(class = "note",
            tags$b("What this shows: "), "how many compound×group rows in the reference table ",
            "have been assigned to a Mode of Action group, and which group. Groups come from Kramer ",
            "et al. (2024) (", tags$code("Code/taxotox_install_moa.R"), "), a curated external MoA ",
            "database keyed by CAS number — they are not a fixed list, so the legend reflects ",
            "whatever categories are actually present in the current reference table (e.g. ",
            tags$code("Neuromuscular system"), ", ", tags$code("Photosynthesis inhibition"), "). ",
            tags$code("unknown"), " marks compounds with no MoA classification available, either ",
            "because Kramer et al. couldn't classify them or because they're outside that dataset's ",
            "~3,400-chemical scope — hover a bar segment to see which case applies (",
            tags$code("moa_source"), " column). ",
            tags$b("Why it matters: "), "two separate effects push the Concentration Addition with Mode of ",
            "Action grouping (CAMA) method's output toward standard Concentration Addition (PTI). First, ",
            "a table dominated by 'unknown' means CAMA's 'all compounds' variant has few real groups to ",
            "partition by — see the 'known MoA' output variant for a view restricted to compounds with an ",
            "actual classification. Second, and independently of coverage, CAMA's formula mathematically ",
            "converges toward PTI whenever toxic-unit sums are low (Drescher & Boedeker, 1995) — true for ",
            "almost every real water sample, even with many real MoA groups present. Plot 9 below measures ",
            "both effects directly against the repo's own sample data."),
          as.tags(fig7),
          tags$br(),
          if (!is.null(moa_summary_tbl))
            as.tags(datatable(moa_summary_tbl, rownames = FALSE,
                              options = list(pageLength = 30, scrollX = TRUE, dom = "t",
                                             autoWidth = FALSE),
                              caption = "Mode of Action group counts per taxonomic group and classification source"))
          else tags$div()
        )
      } else tags$div(),

      # Plot 9 — CAMA vs standard PTI convergence diagnostic (conditional)
      if (!is.null(fig9)) {
        make_section(
          tags$h2("9. CAMA vs Standard PTI — Convergence Diagnostic"),
          tags$p(class = "note",
            tags$b("What this shows: "), "computed directly from the repo's own sample monitoring ",
            "workbooks (", tags$code("Data/sample_*.xlsx"), "), replicating app.R's PTI and CAMA logic: ",
            "how closely CAMA's mixture effect (", tags$code("E_mix"), ") tracks standard Concentration ",
            "Addition PTI, and why. ",
            tags$b("Why it matters: "), "at low toxic-unit levels, CA and IA mathematically converge ",
            "regardless of MoA grouping (Drescher & Boedeker, 1995) — substituting Step 2 into Step 3 of ",
            "the CAMA formula (§10.5.2) gives E_mix = 1 − exp(−Σ log(1+group_TU)) ≈ 1 − exp(−PTI) ≈ PTI ",
            "whenever every group's toxic-unit sum is ≪ 1 — true for almost every real water sample. The ",
            "dashed line is that single-group limit, PTI/(1+PTI); points hugging it indicate no MoA group ",
            "reached a toxic-unit sum large enough for the grouping to matter for that sample. ",
            tags$b("What to look for: "), "points pulling below the dashed line indicate samples where ",
            "CAMA is genuinely doing more than reproducing PTI — check the table's ",
            tags$code("N samples with a group TU ≥ 0.3"), " column for how often that regime is reached."),
          as.tags(fig9),
          tags$br(),
          if (!is.null(cama_diag_tbl)) as.tags(cama_diag_tbl) else tags$div()
        )
      } else tags$div(),

      # Mismatch table
      if (!is.null(dt_mismatch)) {
        make_section(
          tags$h2(if (!is.null(fig6) || !is.null(fig7) || !is.null(fig9))
                    "10. SSD vs Scaled HC5 — Largest Discrepancies"
                  else "6. SSD vs Scaled HC5 — Largest Discrepancies"),
          tags$p(class = "note",
            tags$b("What this shows: "), "the full table of compound\u00d7group rows where both HC5 methods ",
            "were computed, sorted by |log2(SSD/Scaled)| — largest disagreements first. ",
            "Use this table to identify specific compounds where the scaling factor performs poorly ",
            "and where caution is warranted when interpreting HC5-based TU outputs."),
          as.tags(dt_mismatch)
        )
      } else tags$div()
    )
  )
)

# ── 9. Write output ───────────────────────────────────────────────────────────

out_path <- "../Docs/validate_reference.html"
save_html(report, out_path)
message("Report saved: ", normalizePath(out_path))
