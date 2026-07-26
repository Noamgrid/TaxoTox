# =============================================================================
# TaxoTox — Taxonomic Toxicity Assessment Tool
# =============================================================================
#
# PURPOSE:
#   TaxoTox assesses the potential mixture toxicity of environmental pollutant
#   assemblages at monitoring stations using the Toxic Unit (TU) framework and
#   the Pollution Toxicity Index (PTI). For each sample and each taxonomic
#   group (algae, crustaceans, fish), it computes per-compound TUs by dividing
#   the measured concentration by the median acute effect concentration
#   (LC50/EC50) retrieved from the pre-processed US EPA ECOTOX database.
#   The PTI is the dimensionless sum of all TUs for a sample and taxonomic
#   group — a value ≥ 1 indicates that the mixture may cause acute effects,
#   while ≥ 0.1 indicates potential chronic risk.
#
# METHODOLOGY:
#   Toxic Unit (TU) for compound i in sample j, taxonomic group g:
#       TU_ijg = C_ij / median_LC50_ig
#   where C_ij  = measured concentration of compound i in sample j (ng/L)
#         median_LC50_ig = median of all LC50/EC50 values in ECOTOX for
#                          compound i and group g (ng/L).
#
#   Pollution Toxicity Index (PTI) for sample j and group g:
#       PTI_jg = Σ TU_ijg   (sum over all detected compounds i)
#
# CASRN RESOLUTION PIPELINE:
#   Compound names from the user's file are mapped to CAS Registry Numbers
#   (CASRNs) in two sequential layers:
#     Layer 1 — Known_CAS (exact, instant): curated internal table of
#               commonly monitored compounds. Confirmed automatically.
#     Layer 2 — PubChem REST API (optional, requires internet): queries
#               the PubChem compound database via the webchem R package.
#               Candidates are presented for user review in the
#               'CASRN Matching' tab. Enabled by the "Web search" checkbox.
#   Compounds still unresolved after all layers can be assigned CASRNs
#   manually in the 'CASRN Matching' tab → Manual CASRN Entry section.
#
# ECOTOX DATA PREPARATION:
#   final_ecotox_data.fst is built once by running taxotox_install.R.
#   That script queries the local ECOTOXr SQLite cache, retains only
#   LC50/EC50 endpoints for fish, algae, and crustaceans, applies unit
#   conversions to ng/L, and back-transforms (log)-reported values.
#   The resulting flat table is stored as an fst file for fast session
#   loading. taxotox_install.R's Step 0a offers to refresh the underlying
#   ECOTOX database (~3 months) before rebuilding, when run interactively.
#
# OUTPUT:
#   - Interactive per-sample and per-pollutant bar/box plots (ggplot2).
#   - Searchable DataTables with PTI colour-coded by risk threshold.
#   - Multi-sheet Excel workbook: original data, three toxicity sheets
#     (algae, crustacean, fish), CASRN report, toxicity coverage matrix.
#
# FILES:
#   Code/app.R               — this file (Shiny application)
#   Code/taxotox_install.R   — builds final_ecotox_data.fst from ECOTOX
#   Code/TaxoToxInstall_AFB.R— builds final_ecotox_data.fst from EPA
#                              Aquatic Life Benchmarks (validation only)
#   Code/install_dependencies.R — installs all required R packages
#   Data/Known_CAS.fst           — curated CASRN lookup table
#   Data/taxotox_reference.fst   — pre-aggregated denominators (median LC50, HC5,
#                                  predicted LC50, benchmarks, MoA) per compound × group
#   Data/final_ecotox_data.fst   — raw ECOTOX test results (install only, not deployed)
#
# AUTHORS: Noam Gridish, Yair Suari
#          School of Marine Sciences, Ruppin Academic Center, Israel
# LICENSE: MIT
# =============================================================================

# ── Dependencies ──────────────────────────────────────────────────────────────
# Core Shiny framework
library(shiny)
library(shinyjs)
# Excel import / export
library(openxlsx)
library(readxl)
# Tidy data manipulation
library(tidyr)
library(dplyr)
library(tidyverse)
library(purrr)
library(stringr)
# Fast columnar data storage (data.table and fst)
library(data.table)
library(bit64)
library(fst)
# Visualisation
library(ggplot2)
library(scales)
library(ggthemes)
library(ggpattern)
library(ggpubr)
library(ggpmisc)
library(hrbrthemes)
# Interactive tables in the UI
library(DT)
# webchem provides access to cheminformatics web databases (PubChem, CTS,
# ChemSpider, etc.) from R. Used here to query the PubChem REST API for
# CAS Registry Numbers when a compound is not found in Known_CAS.
# No API key is required for PubChem queries.
library(webchem)
# Miscellaneous
library(remotes)

# Allow uploads up to 500 MB (default is 5 MB)
options(shiny.maxRequestSize = 500 * 1024^2)

# =============================================================================
# USER INTERFACE
# =============================================================================
# The UI is a standard Shiny sidebarLayout:
#   Left panel  — three-step action buttons plus configuration controls.
#   Right panel — tabbed display: Instructions, CASRN Matching, Toxicity Plots,
#                 Toxicity Tables.
# =============================================================================

ui <- fluidPage(

    useShinyjs(),

    titlePanel("TaxoTox"),

    sidebarLayout(
        sidebarPanel(

            # ── Step 1: File upload ──────────────────────────────────────────
            fileInput("user_file", label = NULL,
                      buttonLabel = "1. Upload Samples Data",
                      placeholder = "No file selected",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv", ".xlsx", ".xls", ".txt", ".tsv"
                      ),
                      width = "100%"),
            uiOutput("step1_status"),

            # ── Data layout selector ─────────────────────────────────────────
            # Users commonly supply data in one of two orientations:
            #   "col" (default): rows = compounds, columns = stations.
            #   "row"          : rows = stations, columns = compounds.
            # When "row" is selected the matrix is transposed after loading so
            # that the rest of the pipeline always works on the canonical
            # (compounds × stations) layout.
            radioButtons("data_orientation", "Data layout:",
                         choices = c("Pollutant names in first column (default)" = "col",
                                     "Station names in first column (pollutants in header row)" = "row"),
                         selected = "col"),

            # ── Step 2: Run calculations ─────────────────────────────────────
            uiOutput("calc_button_ui"),
            uiOutput("step2_status"),

            # ── Step 3: Export results ───────────────────────────────────────
            uiOutput("download_button_ui"),
            hr(),

            checkboxInput("use_pubchem", "Web search for pollutant names on/off", value = TRUE),

            hr(),

            # ── Advanced Assessment panel (A-4) ──────────────────────────────
            checkboxInput("show_advanced", strong("Advanced Assessment"), value = FALSE),

            conditionalPanel(
                condition = "input.show_advanced == true",

                # Method checkboxes (A-5)
                tags$p(style = "margin-bottom:4px; font-weight:bold; font-size:13px;",
                       "Methods:"),
                checkboxInput("method_hc5",       "HC5-based TU/PTI",       value = FALSE),
                checkboxInput("method_benchmark", "Benchmark Hazard Index",  value = FALSE),
                checkboxInput("method_ia",        "Independent Action (IA)", value = FALSE),
                checkboxInput("method_cama",      "CAMA",                    value = FALSE),
                checkboxInput("method_nowell",    "Taxon-Sensitive PTI (Nowell et al. 2014)", value = FALSE),

                # Benchmark sub-selector (A-6) — below all methods, shown only when benchmark is checked
                # EU EQS / AU ANZG / CA CCME were removed due to low compound coverage
                # (<3% each in the current reference table) -- US EPA is the only
                # remaining framework, kept ~10% coverage.
                conditionalPanel(
                    condition = "input.method_benchmark == true",
                    hr(style = "margin:6px 0;"),
                    tags$p(style = "margin-bottom:4px; font-weight:bold; font-size:13px;",
                           "Benchmarks:"),
                    checkboxInput("bm_usepa", "US EPA",  value = TRUE),
                    hr(style = "margin:6px 0;")
                ),

                # Coverage warning (A-8) — rendered by server when coverage is low
                uiOutput("advanced_coverage_warning")
            ),

            # ── Reference Data download panel ─────────────────────────────────
            wellPanel(
                h4("Reference Data"),
                p(style = "font-size:0.85em; color:#555;",
                  "Download the underlying toxicity reference tables used by TaxoTox."),
                downloadButton("download_taxotox_reference",
                               "TaxoTox Reference Data (LC50/HC5/Benchmarks)",
                               class = "btn-default btn-block",
                               style = "white-space:normal; height:auto; padding:8px 12px;"),
                br(), br(),
                downloadButton("download_known_cas",
                               "Known CASRN Lookup Table",
                               class = "btn-default btn-block",
                               style = "white-space:normal; height:auto; padding:8px 12px;")
            )
        ),

        mainPanel(
           tabsetPanel(id = "main_tabs",

               # ── Instructions tab ─────────────────────────────────────────
               tabPanel("Instructions",
                        h3("Workflow"),
                        p("Follow the three steps in the sidebar to analyse your data:"),
                        tags$ol(
                            tags$li(
                                strong("Upload your data."), " Select the correct data layout using the radio buttons, then click ",
                                strong("'1. Upload Samples Data'"), ". The app immediately searches for CAS Registry Numbers (CASRNs) for all compounds and navigates to the CASRN Matching tab."
                            ),
                            tags$li(
                                strong("Review CASRN matches."), " In the ", strong("'CASRN Matching'"), " tab, verify any PubChem candidates using the checkboxes and click ", strong("'Submit Selected Lines'"), ". If any compounds remain unmatched, enter their CASRNs in the ", strong("'Manual CASRN Entry'"), " section below the table. Once all compounds are resolved, the '2. Calculate Toxicity' button becomes active."
                            ),
                            tags$li(
                                strong("Calculate and export."), " Click ", strong("'2. Calculate Toxicity'"), " to run the PTI calculation. Results appear in the ",
                                strong("'Toxicity Plots'"), " and ", strong("'Toxicity Tables'"), " tabs. Click ",
                                strong("'3. Download Results'"), " to export a multi-sheet Excel workbook."
                            )
                        ),

                        h4("Input File Format"),
                        p("Accepted formats: CSV, TSV, TXT, XLS, XLSX.",
                          strong(" All concentrations must be in ng/L."),
                          " The first column must be either compound names (Layout A) or sample/station identifiers (Layout B)."),

                        h5("Layout A — Pollutant names in first column (default)"),
                        p("Each row is a compound; each subsequent column is a monitoring station or sample identifier."),
                        tags$pre(paste(
                            "Compound    | Station1 | Station2 | ...",
                            "Caffeine    | 10.5     | 18.2     | ...",
                            "Atrazine    |  2.1     |  0.5     | ...",
                            "Bisphenol A |  2.1     |  0.3     | ...",
                            sep = "\n"
                        )),

                        h5("Layout B — Station names in first column"),
                        p("Each row is a monitoring station; compound names are in the header row. The app automatically transposes this matrix to Layout A before processing."),
                        tags$pre(paste(
                            "Station  | Caffeine | Atrazine | Bisphenol A | ...",
                            "Station1 | 10.5     |  2.1     |  2.1        | ...",
                            "Station2 | 18.2     |  0.5     |  0.3        | ...",
                            sep = "\n"
                        )),
                        p(em("Note: if the app detects that more than 75% of names in the first column are unrecognised, it will offer to transpose the data automatically.")),

                        h4("CASRN Identification"),
                        p("On upload, compound names are resolved to CAS Registry Numbers (CASRNs) via two sequential layers:"),
                        tags$ol(
                            tags$li(
                                strong("Known_CAS (instant, always on):"),
                                " a curated internal table of commonly monitored organic micropollutants. Exact-name matches are confirmed automatically with no user action required."
                            ),
                            tags$li(
                                strong("PubChem REST API (optional):"),
                                " live query of the PubChem compound database via the webchem R package. Enable with the 'Web search' checkbox. Requires an internet connection; no API key needed.",
                                " Candidate matches are shown in the CASRN Matching tab for review."
                            )
                        ),
                        p("Compounds still unresolved after both layers appear in the ",
                          strong("'Manual CASRN Entry'"), " section of the CASRN Matching tab. Enter the correct CASRN (e.g. 107-13-1) and click 'Add CASRN'. All newly confirmed CASRNs are logged for future integration into Known_CAS."),

                        h4("Toxicity Calculation"),
                        p("For each compound i, sample j, and taxonomic group g (algae, crustaceans, fish), the Toxic Unit is:"),
                        tags$pre("TU\u1d62\u2C7C\u1d4D = C\u1d62\u2C7C (ng/L)  /  median\u2080 LC50\u1d62\u1d4D (ng/L)"),
                        p("where the denominator is the median of all LC50 and EC50 values for compound i and group g in the pre-processed ECOTOX database (", em("final_ecotox_data.fst"), ")."),
                        p("The Pollution Toxicity Index (PTI) for sample j and group g is the sum of all TUs:"),
                        tags$pre("PTI\u2C7C\u1d4D = \u03a3\u1d62 TU\u1d62\u2C7C\u1d4D"),

                        h4("Risk Interpretation"),
                        tags$table(style = "border-collapse:collapse; margin-bottom:12px;",
                            tags$tr(
                                tags$th(style = "border:1px solid #ccc; padding:4px 10px; background:#f5f5f5;", "PTI"),
                                tags$th(style = "border:1px solid #ccc; padding:4px 10px; background:#f5f5f5;", "Interpretation")
                            ),
                            tags$tr(
                                tags$td(style = "border:1px solid #ccc; padding:4px 10px; color:#B71C1C; font-weight:bold;", "\u2265 1.0"),
                                tags$td(style = "border:1px solid #ccc; padding:4px 10px;", "Acute risk — mixture may cause lethal or immobilisation effects")
                            ),
                            tags$tr(
                                tags$td(style = "border:1px solid #ccc; padding:4px 10px; color:#E65100; font-weight:bold;", "0.1 \u2013 1.0"),
                                tags$td(style = "border:1px solid #ccc; padding:4px 10px;", "Chronic risk — sub-lethal effects on sensitive species possible")
                            ),
                            tags$tr(
                                tags$td(style = "border:1px solid #ccc; padding:4px 10px;", "< 0.1"),
                                tags$td(style = "border:1px solid #ccc; padding:4px 10px;", "Low risk — mixture unlikely to cause measurable toxicity")
                            )
                        ),
                        p(em("PTI values \u2265 0.1 are highlighted orange; \u2265 1.0 are highlighted red in the Toxicity Tables tab.")),

                        h4("Output"),
                        p("The downloaded Excel workbook contains six sheets:"),
                        tags$ul(
                            tags$li(strong("Original data"), " — the uploaded concentration matrix (ng/L)."),
                            tags$li(strong("Algae / Crustacean / Fish toxicity"), " — per-compound TUs and PTI for each sample."),
                            tags$li(strong("CASRN Report"), " — matched and unmatched compound list with CASRNs."),
                            tags$li(strong("Toxicity Coverage"), " — colour-coded matrix showing which taxonomic groups have ECOTOX data for each compound.")
                        )
                        ),
               tabPanel("CASRN Matching",
                        uiOutput("duplicate_casrn_warning"),
                        h4("Search Summary"),
                        verbatimTextOutput("cas_search_summary"),
                        hr(),
                        h4("Candidate Matches"),
                        p("PubChem matches (100%) are pre-checked. Review each suggestion and uncheck any you wish to reject, then click Submit."),
                        actionButton("submit_fuzzy_matches", "Submit Selected Lines", class = "btn-primary"),
                        br(),
                        DTOutput("fuzzy_match_table"),
                        hr(),
                        h4("Manual CASRN Entry"),
                        uiOutput("manual_entry_tab_ui")
                        ),

               tabPanel("Toxicity Plots",
                        div(style = "background:#fff8e1; border-left:4px solid #f9a825; padding:10px 14px; margin-bottom:14px; font-size:0.92em; color:#555;",
                            strong("Disclaimer: "),
                            "TaxoTox is a research screening tool. Outputs are indicative only and should not be used as the sole basis for regulatory decisions or risk management actions. ",
                            "Predicted LC50 values (marked \u007e) carry additional uncertainty. ",
                            "Benchmark values reflect the listed jurisdiction and may not apply to all study contexts; regulatory benchmarks are periodically revised and users should verify currency with the issuing authority. ",
                            "Results depend on the ECOTOX database version installed."
                        ),
                        h4("Top 10 Riskiest Samples (by PTI)"),
                        plotOutput("algae_sample_plot"),
                        plotOutput("crustacean_sample_plot"),
                        plotOutput("fish_sample_plot"),
                        hr(),
                        h4("Top 10 Most Toxic Pollutants per Taxonomic Group"),
                        plotOutput("algae_pollutant_plot"),
                        plotOutput("crustacean_pollutant_plot"),
                        plotOutput("fish_pollutant_plot")
                        ),
               tabPanel("Toxicity Tables",
                        h4("Algae"),
                        DTOutput("tox_table_algae"),
                        h4("Crustacean"),
                        DTOutput("tox_table_crustacean"),
                        h4("Fish"),
                        DTOutput("tox_table_fish")
                        )

           )
        )
    )
)

# =============================================================================
# SERVER
# =============================================================================

# ID of the Google Sheet used to log new compounds from app sessions.
# The Sheet must be shared with the service account in taxotox-service.json.
# Set TAXOTOX_SHEET_ID as an env var to override (e.g. for testing).
TAXOTOX_SHEET_ID <- "1pftfWQNfIStasPqH1CvpDSAH3yHxYmjhggvAwlBaxDA"
LOG_SHEET_NAME   <- "app_log"

# Hard cap on how many unresolved compound names a single upload will send to
# the PubChem API (two blocking HTTP calls each, ~1-2s/compound). Without this
# cap, a large file with many unresolved names (e.g. uploaded with the wrong
# data-orientation setting) triggers an unbounded serial network loop that
# freezes the whole session for hours. Names beyond the cap are routed
# straight to Manual CASRN Entry instead. Override via env var for tuning.
PUBCHEM_MAX_AUTO <- as.integer(Sys.getenv("TAXOTOX_PUBCHEM_MAX_AUTO", unset = 200))

server <- function(input, output, session) {

    # ── Reactive state ────────────────────────────────────────────────────────
    # All mutable application state is stored in a single reactiveValues object
    # to avoid scattered reactive dependencies and simplify reasoning about
    # data flow.
    v <- reactiveValues(
        user_data            = NULL,     # data.table: PREFERRED_NAME + one column per station (ng/L)
        p_vector             = NULL,     # character vector of compound names from the uploaded file
        summary_log          = c("Welcome to TaxoTox! Please upload a file to begin."),
        casrn_results        = NULL,     # (reserved for future use)
        fuzzy_to_review      = NULL,     # data.table of fuzzy matches awaiting user confirmation
        manual_to_fill       = NULL,     # data.table of compounds with no CASRN yet
        final_search_results = data.table(),  # accumulates confirmed PREFERRED_NAME / CASRN pairs
        tox_results          = NULL,     # named list of result data.tables (one per taxonomic group)
        hc5_results          = NULL,     # named list: HC5-based TU tables per group (advanced Excel)
        benchmark_results    = NULL,     # named list: one HI table per selected benchmark framework
        ia_results           = NULL,     # named list: IA effect-fraction tables per group
        cama_results         = NULL,     # data.frame: CAMA E_mix per sample
        plots                = NULL,     # named list of ggplot objects
        manual_additions     = data.table(PREFERRED_NAME = character(), CASRN = character()),
        pending_unfound      = NULL       # saved when orientation-check modal is shown
    )

    # ── Static reference data (loaded once at session start) ─────────────────
    # Known_CAS        : Curated table of compound name → CASRN mappings.
    # final_ecotox_data: Pre-processed ECOTOX data containing median acute
    #                    LC50/EC50 (conc_ng_L) per compound (cas_number) and
    #                    taxonomic group (ecotox_group).
    .data_dir          <- if (dir.exists("../Data")) "../Data" else "Data"
    Known_CAS          <- read.fst(file.path(.data_dir, "Known_CAS.fst"),         as.data.table = TRUE)
    taxotox_reference  <- read.fst(file.path(.data_dir, "taxotox_reference.fst"), as.data.table = FALSE) %>%
        mutate(cas_number = as.character(cas_number))

    # ── Reference Data download handlers ──────────────────────────────────────
    # Let users download the same reference tables the app uses internally,
    # as Excel workbooks, for transparency / offline reference.
    output$download_taxotox_reference <- downloadHandler(
        filename = function() "taxotox_reference_data.xlsx",
        content  = function(file) openxlsx::write.xlsx(taxotox_reference, file)
    )

    output$download_known_cas <- downloadHandler(
        filename = function() "taxotox_known_cas.xlsx",
        content  = function(file) openxlsx::write.xlsx(as.data.frame(Known_CAS), file)
    )

    # ── Button state management (D) ───────────────────────────────────────────
    # Disable buttons 2 and 3 at startup; enable reactively as prerequisites met.
    shinyjs::disable("start_toxicity_calc")
    shinyjs::disable("download_results")

    observe({
        casrn_ready <- nrow(v$final_search_results) > 0 &&
            (is.null(v$fuzzy_to_review) || nrow(v$fuzzy_to_review) == 0)
        if (casrn_ready) shinyjs::enable("start_toxicity_calc")
        else             shinyjs::disable("start_toxicity_calc")
    })

    observe({
        if (!is.null(v$tox_results)) shinyjs::enable("download_results")
        else                         shinyjs::disable("download_results")
    })

    # ── Advanced panel reactive helpers ──────────────────────────────────────────

    # Any advanced method selected?
    any_advanced <- reactive({
        input$show_advanced &&
            (input$method_hc5 || input$method_benchmark || input$method_ia || input$method_cama ||
             input$method_nowell)
    })

    # Dynamic button labels
    output$calc_button_ui <- renderUI({
        if (any_advanced())
            actionButton("start_toxicity_calc",
                         tagList("2. Calculate Basic", br(), "& Advanced Toxicity"),
                         class = "btn-default btn-block",
                         style = "margin-top:10px; white-space:normal; height:auto; padding:8px 12px;")
        else
            actionButton("start_toxicity_calc", "2. Calculate Toxicity",
                         class = "btn-default btn-block",
                         style = "margin-top:10px;")
    })

    output$download_button_ui <- renderUI({
        if (any_advanced())
            downloadButton("download_results",
                           tagList("3. Download Basic", br(), "& Advanced Results"),
                           class = "btn-default btn-block",
                           style = "margin-top:10px; white-space:normal; height:auto; padding:8px 12px;")
        else
            downloadButton("download_results", "3. Download Results",
                           class = "btn-default btn-block",
                           style = "margin-top:10px;")
    })

    # Coverage warning: fraction of matched compounds that have the required column
    output$advanced_coverage_warning <- renderUI({
        req(any_advanced(), v$tox_results)
        msgs <- character(0)

        matched_cas <- unique(gsub("-", "", v$final_search_results$CASRN))
        ref_matched <- taxotox_reference %>% filter(cas_number %in% matched_cas)
        n <- length(matched_cas)

        if (input$method_hc5) {
            n_covered <- ref_matched %>%
                filter(!is.na(hc5_ssd_ng_L) | !is.na(hc5_model_ng_L)) %>%
                distinct(cas_number) %>% nrow()
            pct <- round(100 * n_covered / max(n, 1))
            if (pct < 80) msgs <- c(msgs, sprintf("HC5: %d%% compound coverage", pct))
        }
        if (input$method_benchmark) {
            bm_cols <- c(
                if (input$bm_usepa) "benchmark_usepa_fish_acute_ng_L"  else NULL
            )
            for (col in bm_cols) {
                if (!col %in% names(ref_matched)) next
                n_covered <- ref_matched %>% filter(!is.na(.data[[col]])) %>%
                    distinct(cas_number) %>% nrow()
                pct <- round(100 * n_covered / max(n, 1))
                label <- switch(col,
                    benchmark_usepa_fish_acute_ng_L  = "US EPA")
                if (pct < 80) msgs <- c(msgs, sprintf("%s: %d%% compound coverage", label, pct))
            }
        }
        if (input$method_nowell) {
            n_covered <- ref_matched %>%
                filter(!is.na(stc_nowell_ng_L)) %>%
                distinct(cas_number) %>% nrow()
            pct <- round(100 * n_covered / max(n, 1))
            if (pct < 80) msgs <- c(msgs, sprintf("Nowell PTI: %d%% compound coverage", pct))
        }

        if (length(msgs) == 0) return(NULL)
        tags$div(
            style = "background:#fff3cd; border:1px solid #ffc107; border-radius:4px; padding:6px 8px; margin-top:6px; font-size:0.82em;",
            tags$strong("\u26a0 Low coverage:"),
            tags$ul(style = "margin:4px 0 0 0; padding-left:16px;",
                    lapply(msgs, tags$li))
        )
    })

    # ── Sidebar status badges (E) ─────────────────────────────────────────────
    output$step1_status <- renderUI({
        req(v$p_vector)
        n_total   <- length(v$p_vector)
        n_matched <- nrow(v$final_search_results)
        color <- if (n_matched == n_total) "#2a7a2a" else "#b36200"
        tags$p(
            style = paste0("color:", color, "; font-size:0.85em; margin:2px 0 6px 0;"),
            paste0("\u2713 ", n_matched, " / ", n_total, " compounds matched")
        )
    })

    output$step2_status <- renderUI({
        req(v$tox_results)
        tags$p(
            style = "color:#2a7a2a; font-size:0.85em; margin:2px 0 6px 0;",
            "\u2713 Toxicity calculated"
        )
    })

    # =========================================================================
    # HELPER FUNCTIONS
    # =========================================================================

    # ── load_user_file ────────────────────────────────────────────────────────
    # Reads the user-uploaded concentration file into a data.table, dispatching
    # on file extension. Supported: CSV, TXT (tab-delimited), TSV, XLS, XLSX.
    load_user_file <- function(file_path) {
        ext <- tools::file_ext(file_path)
        df <- switch(tolower(ext),
                     "csv"  = read.csv(file_path, stringsAsFactors = FALSE),
                     "txt"  = read.delim(file_path, stringsAsFactors = FALSE),
                     "tsv"  = read.delim(file_path, sep = "\t", stringsAsFactors = FALSE),
                     "xls"  = readxl::read_excel(file_path),
                     "xlsx" = readxl::read_excel(file_path),
                     stop("Unsupported file type: ", ext)
        )
        return(as.data.table(df))
    }

    # ── .gs4_auth_auto ────────────────────────────────────────────────────────
    # Authenticates with Google Sheets using the service account JSON bundled
    # alongside app.R (gitignored, uploaded to shinyapps.io with the app).
    # Falls back to interactive OAuth if the file is absent (dev machines that
    # haven't downloaded the key yet).
    .GS4_KEY_PATH <- "taxotox-service.json"
    .gs4_auth_auto <- function() {
        if (file.exists(.GS4_KEY_PATH)) {
            googlesheets4::gs4_auth(path = .GS4_KEY_PATH)
        } else {
            googlesheets4::gs4_auth()   # interactive OAuth — browser prompt
        }
    }

    # ── .ensure_log_tab / log_event ──────────────────────────────────────────
    # Creates the app_log tab in the Google Sheet if it does not already exist,
    # then appends one structured row per event. Both functions are fire-and-
    # forget: they wrap every Google Sheets call in tryCatch and never throw,
    # so a logging failure can never crash the app.
    .ensure_log_tab <- function(sheet_id) {
        tryCatch({
            .gs4_auth_auto()
            existing <- googlesheets4::sheet_names(sheet_id)
            if (!LOG_SHEET_NAME %in% existing) {
                googlesheets4::sheet_add(ss = sheet_id, sheet = LOG_SHEET_NAME)
                header <- data.frame(
                    timestamp   = character(), level       = character(),
                    event       = character(), message     = character(),
                    session_id  = character(), file_name   = character(),
                    n_compounds = integer(),   n_samples   = integer(),
                    stringsAsFactors = FALSE
                )
                googlesheets4::sheet_write(data = header, ss = sheet_id,
                                           sheet = LOG_SHEET_NAME)
            }
        }, error = function(e) invisible(NULL))
    }

    log_event <- function(level, event, message,
                          session_id  = NA_character_,
                          file_name   = NA_character_,
                          n_compounds = NA_integer_,
                          n_samples   = NA_integer_) {
        sheet_id <- Sys.getenv("TAXOTOX_SHEET_ID", unset = TAXOTOX_SHEET_ID)
        if (!nzchar(sheet_id)) return(invisible(NULL))
        row <- data.frame(
            timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            level       = level,
            event       = event,
            message     = as.character(message)[1],
            session_id  = as.character(session_id),
            file_name   = as.character(file_name),
            n_compounds = as.integer(n_compounds),
            n_samples   = as.integer(n_samples),
            stringsAsFactors = FALSE
        )
        tryCatch({
            .gs4_auth_auto()
            googlesheets4::sheet_append(ss = sheet_id, data = row,
                                        sheet = LOG_SHEET_NAME)
        }, error = function(e) invisible(NULL))
    }

    # ── Session logging init ──────────────────────────────────────────────────
    .session_id <- substr(session$token, 1, 8)
    local({
        sid <- Sys.getenv("TAXOTOX_SHEET_ID", unset = TAXOTOX_SHEET_ID)
        if (nzchar(sid)) .ensure_log_tab(sid)
    })
    log_event("INFO", "session_start", "New session", session_id = .session_id)

    # ── append_to_temp_cas ────────────────────────────────────────────────────
    # Appends newly confirmed PREFERRED_NAME / CASRN pairs to the curation
    # Google Sheet, so data persists across shinyapps.io server restarts and
    # curation always has one authoritative source. Curation is Sheets-only
    # now (no local temp_CAS.fst fallback) -- TAXOTOX_SHEET_ID has a real
    # default (see top of file), so a genuinely empty sheet_id only happens
    # if that default is deliberately overridden to "".
    append_to_temp_cas <- function(new_data) {
        if (!all(c("PREFERRED_NAME", "CASRN") %in% names(new_data)))
            stop("New data must contain PREFERRED_NAME and CASRN columns.")
        new_data <- new_data %>%
            mutate(CASRN = gsub("-", "", trimws(CASRN))) %>%
            select(PREFERRED_NAME, CASRN)

        sheet_id <- Sys.getenv("TAXOTOX_SHEET_ID", unset = TAXOTOX_SHEET_ID)

        if (!nzchar(sheet_id)) {
            warning("TAXOTOX_SHEET_ID is empty — compound not logged (curation is Sheets-only).")
            return(invisible(NULL))
        }

        tryCatch({
            library(googlesheets4)
            .gs4_auth_auto()
            googlesheets4::sheet_append(ss = sheet_id, data = new_data)
        }, error = function(e) {
            warning("Google Sheets write failed — compound not logged: ", conditionMessage(e))
        })
    }

    # ── pubchem_lookup ────────────────────────────────────────────────────────
    # Queries the PubChem REST API (via the webchem package) for CAS Registry
    # Numbers for compound names that could not be resolved by Known_CAS.
    #
    # For each name the function:
    #   1. Calls webchem::get_cid() to search PubChem by compound name and
    #      retrieve the best-matching PubChem Compound ID (CID).
    #   2. Calls webchem::pc_synonyms() to fetch all synonyms registered for
    #      that CID in PubChem. PubChem stores CAS Registry Numbers as synonyms
    #      in the canonical format: digits-digits-digit (e.g. "107-13-1").
    #   3. Extracts the first synonym that matches the CASRN regex pattern.
    #
    # Design notes:
    #   - No API key is required; PubChem's PUG REST endpoint is openly
    #     accessible.
    #   - The function wraps each compound query in tryCatch so a network
    #     failure or missing entry for one compound does not abort the loop.
    #   - 'match = "first"' in get_cid() retrieves only the highest-confidence
    #     PubChem hit, minimising false positives.
    #   - verbose = FALSE suppresses per-compound progress messages.
    #
    # Returns a data.table with columns:
    #   source_name  — the original compound name as supplied by the user.
    #   matched_name — same as source_name (PubChem was queried by exact name).
    #   CASRN        — the CAS Registry Number extracted from PubChem synonyms.
    #   source       — always "PubChem" (used to label rows in the review table).
    #
    # Returns an empty (0-row) data.table with the same column schema if no
    # names are resolved (preserves rbindlist compatibility).
    pubchem_lookup <- function(compound_names) {
        cas_pattern <- "^\\d{2,7}-\\d{2}-\\d$"   # canonical CASRN regex
        results <- list()
        for (name in compound_names) {
            tryCatch({
                # Step 1: name → PubChem CID
                cid_res <- webchem::get_cid(name, from = "name",
                                            domain = "compound",
                                            match  = "first",
                                            verbose = FALSE)
                if (is.null(cid_res) || all(is.na(cid_res$cid))) next
                cid <- cid_res$cid[!is.na(cid_res$cid)][1]

                # Step 2: CID → synonym list
                syns <- webchem::pc_synonyms(cid, from = "cid", verbose = FALSE)
                if (is.null(syns) || is.null(syns[[as.character(cid)]])) next
                syn_vec <- syns[[as.character(cid)]]

                # Step 3: extract first CASRN-pattern synonym
                cas_hits <- grep(cas_pattern, syn_vec, value = TRUE)
                if (length(cas_hits) == 0) next

                results[[length(results) + 1]] <- data.table(
                    source_name  = name,
                    matched_name = name,
                    CASRN        = cas_hits[1],
                    source       = "PubChem"
                )
            }, error = function(e) NULL)  # silently skip on network/parse errors
        }
        if (length(results) > 0) return(rbindlist(results))
        # Return empty table with correct schema for safe rbindlist() downstream
        data.table(source_name  = character(),
                   matched_name = character(),
                   CASRN        = character(),
                   source       = character())
    }

    # ── run_casrn_layers2to4 ──────────────────────────────────────────────────
    # Runs layer 2 (PubChem) of the CASRN pipeline on `unfound` (compounds not
    # matched in Layer 1 / Known_CAS). Builds the review table and updates
    # v$summary_log. Closes any open modal before starting.
    # Accesses v, input, session via closure.
    run_casrn_layers2to4 <- function(unfound) {
        removeModal()

        # ── Cap: bound the PubChem loop regardless of how it was reached ────
        # A large `unfound` count (e.g. from a file uploaded with the wrong
        # data-orientation setting) would otherwise trigger an unbounded
        # serial network loop in pubchem_lookup() below, freezing the app for
        # hours. Anything beyond the cap is simply skipped (not analysed) —
        # the user is notified clearly and the rest of the pipeline proceeds
        # normally with the capped subset.
        if (length(unfound) > PUBCHEM_MAX_AUTO) {
            n_skipped <- length(unfound) - PUBCHEM_MAX_AUTO
            unfound   <- unfound[seq_len(PUBCHEM_MAX_AUTO)]
            v$summary_log <- c(v$summary_log,
                paste0("Note: ", n_skipped, " unresolved compound name(s) exceeded the automatic ",
                       "lookup limit (", PUBCHEM_MAX_AUTO, " per upload) and were skipped — not ",
                       "matched or included in the toxicity calculation. If this count looks ",
                       "unexpectedly large, check that the data-layout radio button is set correctly."))
            showNotification(
                paste0("⚠ ", n_skipped, " compound(s) exceeded the automatic lookup limit (",
                       PUBCHEM_MAX_AUTO, ") and were skipped. Check your data-layout setting if this ",
                       "seems too high."),
                type = "warning", duration = 15)
        }

        withProgress(message = "Searching for CASRNs (layer 2)", value = 0, {
        tryCatch({

            # ── Layer 2: PubChem via webchem ─────────────────────────────────
            pubchem_rows <- data.table()
            if (length(unfound) > 0) {
                if (!isTRUE(input$use_pubchem)) {
                    v$summary_log <- c(v$summary_log, "Layer 2 (PubChem): skipped (web search is off).")
                } else {
                    setProgress(0.40, detail = paste("Layer 2: PubChem query (",
                                                     length(unfound), "compounds)..."))
                    v$summary_log <- c(v$summary_log,
                                       paste("Layer 2 (PubChem): querying",
                                             length(unfound), "compounds (may take a moment)..."))
                    wc <- pubchem_lookup(unfound)
                    if (nrow(wc) > 0) {
                        pubchem_rows <- data.table(
                            source_name  = wc$source_name,
                            matched_name = wc$matched_name,
                            distance     = NA_real_,
                            CASRN        = wc$CASRN,
                            source       = "PubChem"
                        )
                        unfound <- unfound[!unfound %in% wc$source_name]
                        v$summary_log <- c(v$summary_log,
                                           paste("Layer 2 (PubChem):", nrow(pubchem_rows),
                                                 "matches;", length(unfound), "remaining."))
                    } else {
                        v$summary_log <- c(v$summary_log,
                                           paste("Layer 2 (PubChem): no matches;",
                                                 length(unfound), "remaining."))
                    }
                }
            }

            # ── Build review table ───────────────────────────────────────────
            setProgress(0.90, detail = "Building review table...")

            known_cas_rows <- if (nrow(v$final_search_results) > 0) {
                data.table(
                    source_name  = v$final_search_results$PREFERRED_NAME,
                    matched_name = v$final_search_results$PREFERRED_NAME,
                    distance     = 0,
                    CASRN        = v$final_search_results$CASRN,
                    source       = "Known_CAS"
                )
            } else { data.table() }

            manual_rows <- if (length(unfound) > 0) {
                data.table(
                    source_name  = unfound,
                    matched_name = NA_character_,
                    distance     = NA_real_,
                    CASRN        = NA_character_,
                    source       = "Manual"
                )
            } else { data.table() }

            if (length(unfound) > 0)
                v$summary_log <- c(v$summary_log,
                                   paste(length(unfound),
                                         "compound(s) need manual CASRN entry."))

            review_rows <- rbindlist(
                list(manual_rows, pubchem_rows, known_cas_rows),
                use.names = TRUE, fill = TRUE
            )

            v$manual_to_fill <- data.table(PREFERRED_NAME = character())

            if (nrow(review_rows) > 0) {
                v$fuzzy_to_review <- review_rows
                n_need_action <- nrow(manual_rows)
                v$summary_log <- c(v$summary_log,
                                   paste("Done. Review table:", nrow(review_rows),
                                         "rows total;", n_need_action,
                                         "need attention (Manual)."),
                                   if (n_need_action == 0)
                                       ">> Ready! Click '2. Calculate Toxicity' to proceed."
                                   else
                                       ">> Review the table below, then click '2. Calculate Toxicity'.")
            } else {
                v$fuzzy_to_review <- NULL
                v$summary_log <- c(v$summary_log,
                                   "Done. All compounds auto-confirmed.",
                                   ">> Ready! Click '2. Calculate Toxicity' to proceed.")
            }

            setProgress(1.0, detail = "Complete.")
            log_event("INFO", "casrn_match",
                      paste0(nrow(v$final_search_results), " matched, ",
                             length(unfound), " unresolved"),
                      session_id  = .session_id,
                      file_name   = if (!is.null(input$user_file)) input$user_file$name else NA_character_,
                      n_compounds = length(v$p_vector))

            n_action <- nrow(manual_rows)
            if (n_action == 0) {
                showNotification("Matching complete \u2014 click '2. Calculate Toxicity'",
                                 type = "message", duration = 6)
            } else {
                showNotification(
                    paste0("Matching complete \u2014 ", n_action, " compound(s) need review"),
                    type = "warning", duration = 6)
            }

        }, error = function(e) {
            v$summary_log <- c(v$summary_log,
                               paste("ERROR:", e$message),
                               "Check the CASRN Matching tab for details.")
            showNotification(paste("Error:", e$message), type = "error", duration = 10)
            log_event("ERROR", "casrn_match", e$message,
                      session_id = .session_id,
                      file_name  = if (!is.null(input$user_file)) input$user_file$name else NA_character_)
        })
        }) # end withProgress
    }

    # ── make_sample_risk_plot ─────────────────────────────────────────────────
    # Produces a horizontal bar chart showing the PTI (Pollution Toxicity Index)
    # for the top-N riskiest monitoring stations for a given taxonomic group.
    # The subtitle names the three stations with the highest PTI.
    # Returns an empty plot with an informative title if no data are available.
    make_sample_risk_plot <- function(tox_cal, group_label) {
        if (is.null(tox_cal) || nrow(tox_cal) == 0)
            return(ggplot() + labs(title = paste(group_label, ": No data available")))
        n_show      <- min(10, nrow(tox_cal))
        top_samples <- tox_cal %>% arrange(desc(RQtg)) %>% slice(1:n_show)
        y_max       <- max(max(top_samples$RQtg, na.rm = TRUE) * 1.05, 1.15)
        ggplot(top_samples, aes(x = reorder(Sample, RQtg), y = RQtg)) +
            geom_col(fill = "steelblue") + theme_minimal() +
            geom_hline(yintercept = 0.1, linetype = "dashed", colour = "black", linewidth = 0.8) +
            geom_hline(yintercept = 1.0, linetype = "dashed", colour = "red",   linewidth = 0.8) +
            annotate("label", x = 0.55, y = 0.1,
                     label = "Chronic\neffect (0.1)",
                     angle = 90, hjust = 0, vjust = 0.5,
                     colour = "black", fill = "white", label.size = 0.3,
                     size = 4.2, fontface = "bold") +
            annotate("label", x = 0.55, y = 1.0,
                     label = "Acute\neffect (1.0)",
                     angle = 90, hjust = 0, vjust = 0.5,
                     colour = "red", fill = "white", label.size = 0.3,
                     size = 4.2, fontface = "bold") +
            coord_flip(clip = "off", ylim = c(0, y_max)) +
            labs(title    = paste(group_label, "\u2014 Top", n_show, "Riskiest Samples"),
                 x = "Sample", y = "Pollution Toxicity Index (PTI)")
    }

    # ── make_pollutant_plot ───────────────────────────────────────────────────
    # Produces a horizontal boxplot of per-compound TU distributions across all
    # monitoring stations, limited to the top-10 compounds by median TU.
    # This identifies the individual pollutants contributing most to mixture
    # toxicity across the dataset.
    make_pollutant_plot <- function(tox_cal, group_label) {
        if (is.null(tox_cal) || nrow(tox_cal) == 0)
            return(ggplot() + labs(title = paste(group_label, ": No data available")))
        long_data <- tox_cal %>%
            select(-RQtg) %>%
            pivot_longer(cols = -Sample, names_to = "Compound", values_to = "TU") %>%
            filter(!is.na(TU), TU > 0)
        top_compounds <- long_data %>%
            group_by(Compound) %>%
            summarise(med = median(TU, na.rm = TRUE), .groups = "drop") %>%
            arrange(desc(med)) %>% slice(1:min(10, n())) %>% pull(Compound)
        ggplot(long_data %>% filter(Compound %in% top_compounds),
               aes(x = reorder(Compound, TU, FUN = median), y = TU)) +
            geom_boxplot() + coord_flip() + theme_minimal() +
            labs(title = paste(group_label, "\u2014 Top", length(top_compounds), "Pollutants by TU"),
                 x = "Compound", y = "Toxic Unit (TU)")
    }

    # =========================================================================
    # EVENT: File upload
    # =========================================================================
    # Triggered when the user selects a file. Performs:
    #   1. File parsing and optional matrix transposition (layout B → layout A).
    #   2. Two-layer CASRN resolution pipeline:
    #        Layer 1 — exact match in Known_CAS (auto-confirmed, no review).
    #        Layer 2 — PubChem live API query via webchem (pre-checked).
    #   3. Navigation to the 'CASRN Matching' tab if any candidates exist,
    #      or directly to 'Manual Entry' if no matches were found at all.
    observeEvent(input$user_file, {
        v$summary_log <- c(paste0("File selected: ", input$user_file$name,
                                  " — searching for CASRNs..."))
        # Switch to CASRN Matching immediately so the user can watch the log
        updateTabsetPanel(session, "main_tabs", selected = "CASRN Matching")

        withProgress(message = "Searching for CASRNs", value = 0, {

        tryCatch({
            # ── Load & normalise ─────────────────────────────────────────────
            setProgress(0.05, detail = "Loading file...")
            v$user_data <- load_user_file(input$user_file$datapath)

            setProgress(0.10, detail = "Normalising layout...")
            if (input$data_orientation == "row") {
                station_names  <- as.character(v$user_data[[1]])
                compound_names <- names(v$user_data)[-1]
                transposed     <- as.data.table(t(v$user_data[, -1, with = FALSE]))
                names(transposed) <- station_names
                v$user_data <- cbind(data.table(PREFERRED_NAME = compound_names), transposed)
            } else {
                names(v$user_data)[1] <- "PREFERRED_NAME"
            }

            v$p_vector <- v$user_data[[1]]
            v$summary_log <- c(v$summary_log,
                               paste("Loaded", length(v$p_vector), "compounds from file."))
            log_event("INFO", "file_upload",
                      paste0(input$user_file$name, " — ", length(v$p_vector),
                             " compounds, ", ncol(v$user_data) - 1, " samples"),
                      session_id  = .session_id,
                      file_name   = input$user_file$name,
                      n_compounds = length(v$p_vector),
                      n_samples   = ncol(v$user_data) - 1L)

            # ── Unit plausibility check ──────────────────────────────────────
            numeric_vals <- unlist(v$user_data[, -1, with = FALSE])
            numeric_vals <- suppressWarnings(as.numeric(numeric_vals))
            numeric_vals <- numeric_vals[!is.na(numeric_vals) & is.finite(numeric_vals)]
            if (length(numeric_vals) > 0) {
                med_val <- median(numeric_vals)
                if (med_val > 1e6) {
                    showNotification(
                        paste0("Warning: median concentration is ", signif(med_val, 3),
                               " ng/L \u2014 this is unusually high. ",
                               "Check that units are ng/L, not \u00b5g/L or mg/L."),
                        type = "warning", duration = 15
                    )
                } else if (med_val > 0 && med_val < 0.001) {
                    showNotification(
                        paste0("Warning: median concentration is ", signif(med_val, 3),
                               " ng/L \u2014 this is unusually low. ",
                               "Check that units are ng/L, not pg/L."),
                        type = "warning", duration = 15
                    )
                }
            }

            # ── Layer 1: Known_CAS exact ─────────────────────────────────────
            setProgress(0.20, detail = paste("Layer 1: Known_CAS lookup (",
                                             length(v$p_vector), "compounds)..."))
            internal_list  <- Known_CAS[Known_CAS$PREFERRED_NAME %in% v$p_vector, ]
            # Deduplicate: Known_CAS may have multiple CASRNs per name; keep first.
            internal_list  <- unique(internal_list, by = "PREFERRED_NAME")
            p_vector_found <- internal_list$PREFERRED_NAME
            unfound        <- v$p_vector[!v$p_vector %in% p_vector_found]
            v$pending_unfound      <- unfound   # persisted so the modal buttons can use it
            v$final_search_results <- internal_list[, .(PREFERRED_NAME, CASRN)]
            v$summary_log <- c(v$summary_log,
                               paste("Layer 1 (Known_CAS):", length(p_vector_found),
                                     "exact matches;", length(unfound), "remaining."))

            # ── Orientation / structure check ────────────────────────────────
            # If >50 % of compound names were not found in Known_CAS, ask the
            # user to confirm the file is structured correctly before proceeding.
            pct_unfound <- length(unfound) / max(length(v$p_vector), 1)
            if (pct_unfound > 0.50) {
                v$pending_unfound <- unfound
                pct_label <- paste0(round(pct_unfound * 100), "%")
                showModal(modalDialog(
                    title = "\u26a0\ufe0f Low compound name recognition",
                    p(strong(paste0(length(unfound), " of ", length(v$p_vector),
                                   " (", pct_label, ") compound names were not recognised",
                                   " in the internal database."))),
                    p("Please confirm the file is structured correctly:"),
                    tags$ul(
                        tags$li("The ", strong("first column"), " should contain ",
                                strong("compound / pollutant names")),
                        tags$li("Remaining columns should each be a ",
                                strong("sample or station")),
                        tags$li("All concentration values should be in ", strong("ng/L"))
                    ),
                    p("If the layout is transposed (stations in first column), click ",
                      strong("'Transpose & Retry'")),
                    footer = tagList(
                        actionButton("orientation_transpose_btn",
                                     "Transpose & Retry", class = "btn-warning"),
                        actionButton("orientation_continue_btn",
                                     "Continue Anyway",   class = "btn-default")
                    ),
                    easyClose = FALSE
                ))
                setProgress(1.0)   # close the progress bar
                return()
            }

        }, error = function(e) {
            v$summary_log <- c(v$summary_log,
                               paste("ERROR:", e$message),
                               "Check the CASRN Matching tab for details.")
            showNotification(paste("Error:", e$message), type = "error", duration = 10)
            log_event("ERROR", "file_upload", e$message,
                      session_id = .session_id,
                      file_name  = input$user_file$name)
        })

        }) # end withProgress

        # Continue to layers 2–4, or auto-start if all compounds were found in Known_CAS
        if (length(v$pending_unfound) == 0) {
            v$summary_log <- c(v$summary_log,
                               "All compounds found in Known_CAS — starting calculation automatically.")
            shinyjs::enable("start_toxicity_calc")
            shinyjs::click("start_toxicity_calc")
        } else {
            run_casrn_layers2to4(v$pending_unfound)
        }
    })


    # ── Orientation modal: Transpose & Retry ─────────────────────────────────
    observeEvent(input$orientation_transpose_btn, {
        req(v$user_data)
        # Transpose: swap rows ↔ columns (compound names become header, stations become first col)
        station_names  <- as.character(v$user_data[[1]])
        compound_names <- names(v$user_data)[-1]
        transposed     <- as.data.table(t(v$user_data[, -1, with = FALSE]))
        names(transposed) <- station_names
        v$user_data <- cbind(data.table(PREFERRED_NAME = compound_names), transposed)
        v$p_vector  <- v$user_data[[1]]

        # Re-run Layer 1 on the transposed data
        internal_list  <- Known_CAS[Known_CAS$PREFERRED_NAME %in% v$p_vector, ]
        internal_list  <- unique(internal_list, by = "PREFERRED_NAME")
        p_vector_found <- internal_list$PREFERRED_NAME
        unfound        <- v$p_vector[!v$p_vector %in% p_vector_found]
        v$final_search_results <- internal_list[, .(PREFERRED_NAME, CASRN)]
        v$summary_log <- c(v$summary_log,
                           "Transposed data. Re-running CASRN search...",
                           paste("Layer 1 (Known_CAS):", length(p_vector_found),
                                 "exact matches;", length(unfound), "remaining."))
        run_casrn_layers2to4(unfound)
    })

    # ── Orientation modal: Continue Anyway ───────────────────────────────────
    observeEvent(input$orientation_continue_btn, {
        req(v$pending_unfound)
        v$summary_log <- c(v$summary_log, "Continuing with original orientation...")
        run_casrn_layers2to4(v$pending_unfound)
    })


    # ── CASRN Matching tab: search log ────────────────────────────────────────
    output$cas_search_summary <- renderPrint({ cat(v$summary_log, sep = "\n") })

    # ── Duplicate CASRN warning ───────────────────────────────────────────────
    # If two input compound names resolved to the same CASRN, their TUs will be
    # summed independently — effectively double-counting. Warn the user.
    output$duplicate_casrn_warning <- renderUI({
        req(nrow(v$final_search_results) > 0)
        dupes <- v$final_search_results %>%
            as.data.frame() %>%
            mutate(cas_nodash = gsub("-", "", CASRN)) %>%
            group_by(cas_nodash) %>%
            filter(n() > 1, nzchar(cas_nodash)) %>%
            arrange(cas_nodash)
        if (nrow(dupes) == 0) return(NULL)
        rows <- dupes %>%
            group_by(cas_nodash) %>%
            summarise(names = paste(PREFERRED_NAME, collapse = " / "),
                      casrn = first(CASRN), .groups = "drop")
        items <- lapply(seq_len(nrow(rows)), function(i)
            tags$li(sprintf('"%s" — both map to CASRN %s. Concentrations will be counted separately.',
                            rows$names[i], rows$casrn[i]))
        )
        tags$div(
            style = "background:#fff3cd; border:1px solid #ffc107; border-radius:4px; padding:8px 12px; margin-bottom:10px;",
            tags$strong("\u26a0 Duplicate CASRNs detected:"),
            tags$ul(style = "margin:6px 0 0 0;", items),
            tags$p(style = "margin:6px 0 0 0; font-size:0.88em;",
                   "If these are synonyms for the same compound, remove one from your input file before calculating.")
        )
    })

    # ── CASRN Matching tab: unified compound review table ────────────────────
    # Renders ALL compounds in a single DataTable, ordered so the ones needing
    # the most user attention appear first:
    #
    #   Row section  | Source    | Accept col       | CASRN col
    #   -------------|-----------|------------------|------------------
    #   1 (top)      | Manual    | checkbox (off)   | editable textInput
    #   2            | PubChem   | checkbox (on)    | static text
    #   3 (bottom)   | Known_CAS | "✓ Auto"         | static text
    #
    # The Accept column is rendered narrow (≈50 px).
    # Shiny input bindings are re-attached after every DataTables redraw so
    # that checkboxes and text inputs inside the table remain reactive.
    output$fuzzy_match_table <- renderDT({
        req(v$fuzzy_to_review)
        df <- as.data.frame(v$fuzzy_to_review)

        # ── Sort: Manual → PubChem → Known_CAS ───────────────────────────────
        df$source_order <- dplyr::case_when(
            df$source == "Manual"    ~ 1L,
            df$source == "PubChem"   ~ 2L,
            df$source == "Known_CAS" ~ 3L,
            TRUE                     ~ 4L
        )
        df <- df[order(df$source_order, df$distance, na.last = TRUE), ]

        # ── Match % display (character; "—" for Manual, "✓" for Known_CAS) ───
        df$match_disp <- dplyr::case_when(
            df$source == "Manual"                          ~ "\u2014",
            df$source == "Known_CAS"                       ~ "\u2713",
            is.na(df$distance) | df$distance == 0          ~ "100%",
            TRUE ~ paste0(round((1 - df$distance) * 100, 1), "%")
        )

        # ── Accept column ─────────────────────────────────────────────────────
        # Known_CAS: show a static "✓ Auto" label (already confirmed, no action).
        # PubChem: checkbox pre-checked.
        # Manual: checkbox unchecked.
        accept_col <- sapply(seq_len(nrow(df)), function(i) {
            src <- df$source[i]
            if (src == "Known_CAS") {
                "<span style='color:#2a7a2a; font-weight:bold;'>\u2713 Auto</span>"
            } else {
                pre <- src == "PubChem"
                as.character(checkboxInput(paste0("fuzzy_cb_", i), label = NULL, value = pre))
            }
        })

        # ── CASRN column ──────────────────────────────────────────────────────
        # Manual: editable textInput. All others: static text.
        casrn_col <- sapply(seq_len(nrow(df)), function(i) {
            src <- df$source[i]
            if (src == "Manual") {
                val <- if (is.na(df$CASRN[i])) "" else df$CASRN[i]
                as.character(
                    textInput(paste0("casrn_input_", i), label = NULL,
                              value = val, placeholder = "e.g. 107-13-1",
                              width = "150px")
                )
            } else {
                if (is.na(df$CASRN[i])) "\u2014" else df$CASRN[i]
            }
        })

        display_df <- data.frame(
            "\u2713"    = accept_col,
            "Compound"  = df$source_name,
            "Suggested" = ifelse(is.na(df$matched_name), "\u2014", df$matched_name),
            "Match %"   = df$match_disp,
            "CASRN"     = casrn_col,
            "Source"    = df$source,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        datatable(display_df, escape = FALSE, selection = "none",
                  options = list(
                      pageLength = 25,
                      autoWidth = FALSE,
                      columnDefs = list(
                          list(width = "45px",  targets = 0),    # Accept — narrow
                          list(width = "160px", targets = 4)     # CASRN — room for input
                      ),
                      preDrawCallback = JS("function() { Shiny.unbindAll(this.api().table().node()); }"),
                      drawCallback    = JS("function() { Shiny.bindAll(this.api().table().node()); }")
                  ))
    }, server = FALSE)

    # =========================================================================
    # EVENT: User confirms fuzzy matches
    # =========================================================================
    # Reads checkbox states, separates accepted from rejected candidates, adds
    # accepted pairs to final_search_results and temp_CAS, and moves rejected
    # names to the manual-entry queue.
    observeEvent(input$submit_fuzzy_matches, {
        req(v$fuzzy_to_review)

        # Apply the same sort order as the renderer so that checkbox/input
        # indices match the rows the user sees in the table.
        df <- as.data.frame(v$fuzzy_to_review)
        df$source_order <- dplyr::case_when(
            df$source == "Manual"    ~ 1L,
            df$source == "PubChem"   ~ 2L,
            df$source == "Known_CAS" ~ 3L,
            TRUE                     ~ 4L)
        df <- df[order(df$source_order,
                       ifelse(is.na(df$distance), Inf, df$distance)), ]

        confirmed  <- data.table()
        unresolved <- character()

        for (i in seq_len(nrow(df))) {
            src <- df$source[i]

            # Known_CAS rows are already auto-confirmed — skip them here.
            if (src == "Known_CAS") next

            checked <- isTRUE(input[[paste0("fuzzy_cb_", i)]])

            if (src == "Manual") {
                # CASRN comes from the editable text input.
                raw_val  <- input[[paste0("casrn_input_", i)]]
                casrn_val <- gsub("-", "", trimws(if (is.null(raw_val)) "" else raw_val))
                if (checked && nchar(casrn_val) > 0) {
                    confirmed <- rbindlist(list(confirmed,
                        data.table(PREFERRED_NAME = df$source_name[i],
                                   CASRN          = casrn_val)),
                        fill = TRUE)
                } else {
                    unresolved <- c(unresolved, df$source_name[i])
                }
            } else {
                # PubChem: CASRN is pre-filled in the data.
                casrn_val <- if (is.na(df$CASRN[i])) "" else df$CASRN[i]
                if (checked && nchar(casrn_val) > 0) {
                    confirmed <- rbindlist(list(confirmed,
                        data.table(PREFERRED_NAME = df$source_name[i],
                                   CASRN          = casrn_val)),
                        fill = TRUE)
                } else {
                    unresolved <- c(unresolved, df$source_name[i])
                }
            }
        }

        if (nrow(confirmed) > 0) {
            v$final_search_results <- rbindlist(
                list(v$final_search_results, confirmed),
                use.names = TRUE, fill = TRUE)
            append_to_temp_cas(confirmed)
            v$summary_log <- c(v$summary_log,
                paste("User confirmed", nrow(confirmed), "match(es)."))
        } else {
            v$summary_log <- c(v$summary_log, "No matches were confirmed.")
        }

        if (length(unresolved) > 0) {
            v$manual_to_fill <- rbind(v$manual_to_fill,
                                      data.table(PREFERRED_NAME = unresolved))
            v$summary_log <- c(v$summary_log,
                paste(length(unresolved), "compound(s) moved to manual entry."))
        }

        v$fuzzy_to_review <- NULL  # clear table so it disappears from the UI

        # If unresolved compounds remain, navigate to CASRN Matching (manual entry is there)
        if (!is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0)
            updateTabsetPanel(session, "main_tabs", selected = "CASRN Matching")
    })

    # =========================================================================
    # EVENT: User adds a CASRN manually
    # =========================================================================
    # Adds the entered PREFERRED_NAME / CASRN pair to final_search_results,
    # logs it for display, writes it to temp_CAS, and removes the compound from
    # the manual-entry queue.
    observeEvent(input$add_manual_casrn, {
        req(input$manual_name, input$manual_casrn)
        name  <- trimws(input$manual_name)
        casrn <- gsub("-", "", trimws(input$manual_casrn))
        if (name != "" && casrn != "") {
            new_entry              <- data.table(PREFERRED_NAME = name, CASRN = casrn)
            v$final_search_results <- rbind(v$final_search_results, new_entry)
            v$manual_additions     <- rbind(v$manual_additions,     new_entry)
            append_to_temp_cas(new_entry)
            v$manual_to_fill <- v$manual_to_fill[PREFERRED_NAME != name]
            v$summary_log    <- c(v$summary_log, paste("Manually added CASRN for", name))
            updateTextInput(session, "manual_casrn", value = "")
        }
    })

    output$manual_added_list <- renderDT({ v$manual_additions })

    # =========================================================================
    # EVENT: Calculate Toxicity
    # =========================================================================
    # Core computation. For each taxonomic group (algae, crustacean, fish):
    #   1. Filter the ECOTOX dataset to the matched compounds and group.
    #   2. Compute the median chronic effect concentration across all ECOTOX
    #      records for each compound (median_conc, ng/L).
    #   3. Join with user concentrations and divide each measured value by
    #      median_conc to obtain per-compound Toxic Units (TU).
    #   4. Sum TUs across compounds per station and multiply by 100 to obtain
    #      the Pollution Toxicity Index (PTI).
    observeEvent(input$start_toxicity_calc, {
        req(nrow(v$final_search_results) > 0)
        v$summary_log <- c(v$summary_log, "Starting toxicity calculations...")

        withProgress(message = "Calculating toxicity", value = 0, {
        tryCatch({
            # ── Build the matched compound lookup table ───────────────────────
            setProgress(0.05, detail = "Building compound lookup table...")
            p_final <- as.data.table(v$final_search_results) %>%
                na.omit() %>%
                distinct(PREFERRED_NAME, .keep_all = TRUE) %>%
                mutate(cas_number = gsub("-", "", CASRN)) %>%
                select(PREFERRED_NAME, cas_number)

            cas_search_list <- as.vector(p_final$cas_number)

            # ── Lookup denominators from reference table ──────────────────────
            # taxotox_reference has one row per (cas_number × ecotox_group).
            # Denominator priority: median_lc50_ng_L (ECOTOX) → predicted_lc50_ng_L (CompTox).
            setProgress(0.10, detail = "Looking up toxicity denominators...")
            ref_matched <- taxotox_reference %>%
                filter(cas_number %in% cas_search_list) %>%
                mutate(
                    median_conc  = dplyr::coalesce(median_lc50_ng_L, predicted_lc50_ng_L),
                    gap_fill     = is.na(median_lc50_ng_L) & !is.na(predicted_lc50_ng_L)
                ) %>%
                filter(!is.na(median_conc)) %>%
                left_join(p_final, by = "cas_number") %>%
                filter(!is.na(PREFERRED_NAME)) %>%
                select(PREFERRED_NAME, cas_number, ecotox_group, median_conc, gap_fill)

            # Names that used CompTox predicted LC50 (any group) — flagged with "~" in tables
            gap_fill_names <- ref_matched %>%
                filter(gap_fill) %>%
                distinct(PREFERRED_NAME) %>%
                pull(PREFERRED_NAME)

            # ── Per-group TU / PTI calculation ───────────────────────────────
            .calc_tox <- function(grp, progress_val, progress_label) {
                setProgress(progress_val, detail = progress_label)
                denom <- ref_matched %>%
                    filter(ecotox_group == grp) %>%
                    select(PREFERRED_NAME, median_conc)
                denom %>%
                    left_join(v$user_data, by = "PREFERRED_NAME") %>%
                    mutate(across(where(is.numeric) & !median_conc, ~ .x / median_conc)) %>%
                    select(-median_conc) %>%
                    pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>%
                    pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>%
                    mutate(RQtg = rowSums(across(where(is.numeric)), na.rm = TRUE))
            }

            tox_cal_algae      <- .calc_tox("algae",      0.25, "Algae: computing TUs...")
            tox_cal_crustacean <- .calc_tox("crustacean", 0.45, "Crustacean: computing TUs...")
            tox_cal_fish       <- .calc_tox("fish",       0.65, "Fish: computing TUs...")

            # ── Rename RQtg for display and reorder columns ──────────────────
            setProgress(0.82, detail = "Formatting results...")
            format_tox_output <- function(df) {
                df %>% rename("PTI" = RQtg) %>%
                    select(Sample, PTI, everything())
            }

            v$tox_results <- list(
                "Original data"       = v$user_data,
                "Algae toxicity"      = format_tox_output(tox_cal_algae),
                "Crustacean toxicity" = format_tox_output(tox_cal_crustacean),
                "Fish toxicity"       = format_tox_output(tox_cal_fish)
            )

            setProgress(0.92, detail = "Building plots...")
            v$plots <- list(
                algae_samples         = make_sample_risk_plot(tox_cal_algae,      "Algae"),
                crustacean_samples    = make_sample_risk_plot(tox_cal_crustacean, "Crustacean"),
                fish_samples          = make_sample_risk_plot(tox_cal_fish,       "Fish"),
                algae_pollutants      = make_pollutant_plot(tox_cal_algae,        "Algae"),
                crustacean_pollutants = make_pollutant_plot(tox_cal_crustacean,   "Crustacean"),
                fish_pollutants       = make_pollutant_plot(tox_cal_fish,         "Fish")
            )

            # ── A-9: HC5-based TU/PTI ────────────────────────────────────────
            # Uses hc5_ssd_ng_L preferentially; falls back to hc5_model_ng_L.
            # Runs only when the HC5 method is selected in the Advanced panel.
            if (input$method_hc5) {
                setProgress(0.86, detail = "HC5: computing TUs...")

                ref_hc5 <- taxotox_reference %>%
                    filter(cas_number %in% cas_search_list) %>%
                    mutate(
                        hc5_denom = dplyr::coalesce(hc5_ssd_ng_L, hc5_model_ng_L),
                        hc5_method_used = dplyr::case_when(
                            !is.na(hc5_ssd_ng_L)   ~ paste0(hc5_method, " (SSD)"),
                            !is.na(hc5_model_ng_L)  ~ paste0(hc5_method, " (Scaled)"),
                            TRUE                    ~ NA_character_
                        )
                    ) %>%
                    filter(!is.na(hc5_denom)) %>%
                    left_join(p_final, by = "cas_number") %>%
                    filter(!is.na(PREFERRED_NAME)) %>%
                    select(PREFERRED_NAME, ecotox_group, hc5_denom)

                .calc_hc5 <- function(grp) {
                    denom <- ref_hc5 %>%
                        filter(ecotox_group == grp) %>%
                        select(PREFERRED_NAME, hc5_denom)
                    denom %>%
                        left_join(v$user_data, by = "PREFERRED_NAME") %>%
                        mutate(across(where(is.numeric) & !hc5_denom, ~ .x / hc5_denom)) %>%
                        select(-hc5_denom) %>%
                        pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>%
                        pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>%
                        mutate(PTI = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
                        select(Sample, PTI, everything())
                }

                v$hc5_results <- list(
                    "HC5 Algae"      = .calc_hc5("algae"),
                    "HC5 Crustacean" = .calc_hc5("crustacean"),
                    "HC5 Fish"       = .calc_hc5("fish")
                )
            } else {
                v$hc5_results <- NULL
            }

            # ── A-10: Benchmark Hazard Index ──────────────────────────────────
            # HQ_i = C_i / benchmark_i;  HI = sum(HQ_i)
            # One sheet per selected framework; no taxonomic grouping.
            if (input$method_benchmark) {
                setProgress(0.88, detail = "Benchmark HI: computing...")

                # Helper: compute one HI sheet given a named benchmark column
                # filtered to a specific ecotox_group (or all groups if grp = NULL)
                .calc_hi <- function(col, grp = NULL) {
                    ref_bm <- taxotox_reference %>%
                        filter(cas_number %in% cas_search_list,
                               !is.na(.data[[col]]))
                    if (!is.null(grp))
                        ref_bm <- ref_bm %>% filter(ecotox_group == grp)
                    denom <- ref_bm %>%
                        group_by(cas_number) %>% slice(1) %>% ungroup() %>%
                        left_join(p_final, by = "cas_number") %>%
                        filter(!is.na(PREFERRED_NAME)) %>%
                        select(PREFERRED_NAME, bm_val = !!sym(col))
                    if (nrow(denom) == 0) return(NULL)
                    denom %>%
                        left_join(v$user_data, by = "PREFERRED_NAME") %>%
                        mutate(across(where(is.numeric) & !bm_val, ~ .x / bm_val)) %>%
                        select(-bm_val) %>%
                        pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "HQ") %>%
                        pivot_wider(names_from = PREFERRED_NAME, values_from = HQ) %>%
                        mutate(HI = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
                        select(Sample, HI, everything())
                }

                bm_results <- list()

                # US EPA — three group-specific sheets using acute columns
                if (input$bm_usepa) {
                    usepa_groups <- list(
                        list(grp = "fish",       col = "benchmark_usepa_fish_acute_ng_L",  label = "Bench. US EPA (Fish)"),
                        list(grp = "crustacean", col = "benchmark_usepa_crust_acute_ng_L", label = "Bench. US EPA (Crust)"),
                        list(grp = "algae",      col = "benchmark_usepa_algae_acute_ng_L", label = "Bench. US EPA (Algae)")
                    )
                    for (g in usepa_groups) {
                        if (!g$col %in% names(taxotox_reference)) next
                        df <- .calc_hi(g$col, g$grp)
                        if (!is.null(df)) bm_results[[g$label]] <- df
                    }
                }

                v$benchmark_results <- if (length(bm_results) > 0) bm_results else NULL
            } else {
                v$benchmark_results <- NULL
            }

            # ── A-11: Independent Action (IA, Hill slope n=1) ────────────────
            # E_i = C_i / (LC50_i + C_i)
            # E_mix = 1 - prod(1 - E_i)   [per sample, per taxonomic group]
            # Denominator: same coalesce priority as standard TU
            if (input$method_ia) {
                setProgress(0.91, detail = "IA: computing...")

                .calc_ia <- function(grp) {
                    denom <- ref_matched %>%
                        filter(ecotox_group == grp) %>%
                        select(PREFERRED_NAME, median_conc)
                    if (nrow(denom) == 0) return(NULL)

                    # Join concentrations, compute E_i per compound per sample
                    conc_long <- denom %>%
                        left_join(v$user_data, by = "PREFERRED_NAME") %>%
                        pivot_longer(cols      = -c(PREFERRED_NAME, median_conc),
                                     names_to  = "Sample",
                                     values_to = "C") %>%
                        mutate(E_i = C / (median_conc + C),
                               E_i = if_else(is.na(E_i) | is.nan(E_i), 0, E_i))

                    # E_mix per sample = 1 - prod(1 - E_i)
                    e_mix <- conc_long %>%
                        group_by(Sample) %>%
                        summarise(E_mix = 1 - prod(1 - E_i), .groups = "drop")

                    # Per-compound E_i wide for inspection
                    e_wide <- conc_long %>%
                        select(PREFERRED_NAME, Sample, E_i) %>%
                        pivot_wider(names_from = PREFERRED_NAME, values_from = E_i)

                    left_join(e_mix, e_wide, by = "Sample") %>%
                        select(Sample, E_mix, everything())
                }

                ia_algae      <- .calc_ia("algae")
                ia_crustacean <- .calc_ia("crustacean")
                ia_fish       <- .calc_ia("fish")

                v$ia_results <- Filter(Negate(is.null), list(
                    "IA Algae"      = ia_algae,
                    "IA Crustacean" = ia_crustacean,
                    "IA Fish"       = ia_fish
                ))
                if (length(v$ia_results) == 0) v$ia_results <- NULL
            } else {
                v$ia_results <- NULL
            }

            # ── A-12: CAMA ───────────────────────────────────────────────────
            # Step 1: within each MoA group → CA → group_TU = Σ(C_i / LC50_i)
            # Step 2: group fractional effect: E_group = group_TU / (1 + group_TU)
            # Step 3: between groups → IA: E_mix = 1 − ∏(1 − E_group)
            # One sheet per group; rows = samples; columns: E_mix, then one E_group
            # per MoA group. Computed twice per taxon: once with compounds of
            # unknown MoA included as their own group (honest default), and once
            # with them excluded entirely, so the effect of the unknown bucket on
            # the mixture result is visible rather than hidden inside E_mix.
            if (input$method_cama) {
                setProgress(0.94, detail = "CAMA: computing...")

                # Join moa_group onto ref_matched (already filtered to matched compounds).
                # moa_group/moa_source are populated by Code/taxotox_install_moa.R and
                # should already be "unknown" (never NA) for unclassified compounds;
                # the coalesce here is a defensive fallback only.
                ref_cama <- ref_matched %>%
                    left_join(
                        taxotox_reference %>%
                            select(cas_number, ecotox_group, moa_group) %>%
                            mutate(moa_group = if_else(is.na(moa_group) | moa_group == "",
                                                       "unknown", moa_group)),
                        by = c("cas_number", "ecotox_group")
                    )

                .calc_cama <- function(grp, exclude_unknown = FALSE) {
                    d <- ref_cama %>% filter(ecotox_group == grp)
                    if (exclude_unknown) d <- d %>% filter(moa_group != "unknown")
                    if (nrow(d) == 0) return(NULL)

                    # Long format: one row per compound × sample
                    conc_long <- d %>%
                        select(PREFERRED_NAME, moa_group, median_conc) %>%
                        left_join(v$user_data, by = "PREFERRED_NAME") %>%
                        pivot_longer(cols      = -c(PREFERRED_NAME, moa_group, median_conc),
                                     names_to  = "Sample",
                                     values_to = "C") %>%
                        mutate(TU = if_else(is.na(C) | is.na(median_conc), 0,
                                            C / median_conc))

                    # Step 1+2: group_TU and E_group per MoA × sample
                    e_group <- conc_long %>%
                        group_by(Sample, moa_group) %>%
                        summarise(group_TU = sum(TU, na.rm = TRUE), .groups = "drop") %>%
                        mutate(E_group = group_TU / (1 + group_TU))

                    # Step 3: E_mix per sample
                    e_mix <- e_group %>%
                        group_by(Sample) %>%
                        summarise(E_mix = 1 - prod(1 - E_group), .groups = "drop")

                    # Spread E_group per MoA group as additional columns
                    e_wide <- e_group %>%
                        select(Sample, moa_group, E_group) %>%
                        pivot_wider(names_from  = moa_group,
                                    values_from = E_group,
                                    names_prefix = "E_")

                    left_join(e_mix, e_wide, by = "Sample") %>%
                        select(Sample, E_mix, everything())
                }

                cama_algae      <- .calc_cama("algae")
                cama_crustacean <- .calc_cama("crustacean")
                cama_fish       <- .calc_cama("fish")

                cama_algae_known      <- .calc_cama("algae",      exclude_unknown = TRUE)
                cama_crustacean_known <- .calc_cama("crustacean", exclude_unknown = TRUE)
                cama_fish_known       <- .calc_cama("fish",       exclude_unknown = TRUE)

                v$cama_results <- Filter(Negate(is.null), list(
                    "CAMA Algae"      = cama_algae,
                    "CAMA Crustacean" = cama_crustacean,
                    "CAMA Fish"       = cama_fish,
                    # Excel worksheet names are capped at 31 characters --
                    # "CAMA Crustacean (MoA-known only)" is 32, so all three
                    # use the shorter "(known MoA)" suffix instead.
                    "CAMA Algae (known MoA)"      = cama_algae_known,
                    "CAMA Crustacean (known MoA)" = cama_crustacean_known,
                    "CAMA Fish (known MoA)"       = cama_fish_known
                ))
                if (length(v$cama_results) == 0) v$cama_results <- NULL
            } else {
                v$cama_results <- NULL
            }

            # ── A-13: Taxon-Sensitive PTI (Nowell et al. 2014) ────────────────
            # Denominator: sensitive toxicity concentration (STC) as PUBLISHED
            # by Nowell et al. (2014) Appendix B (see stc_nowell_ng_L in
            # taxotox_reference.fst / Code/taxotox_install_nowell.R), not
            # recomputed from TaxoTox's own ECOTOX pull -- this makes results
            # directly comparable to the literature (e.g. Covert et al. 2020).
            # Covers fish and cladocerans only (crustacean-group row is
            # restricted to Nowell's 17 cladoceran genera); no algae -- Nowell's
            # method does not cover algae at all.
            if (input$method_nowell) {
                setProgress(0.97, detail = "Nowell PTI: computing...")

                .calc_nowell <- function(grp) {
                    denom <- taxotox_reference %>%
                        filter(cas_number %in% cas_search_list,
                               ecotox_group == grp,
                               !is.na(stc_nowell_ng_L)) %>%
                        left_join(p_final, by = "cas_number") %>%
                        filter(!is.na(PREFERRED_NAME)) %>%
                        select(PREFERRED_NAME, stc_nowell_ng_L)
                    if (nrow(denom) == 0) return(NULL)

                    denom %>%
                        left_join(v$user_data, by = "PREFERRED_NAME") %>%
                        mutate(across(where(is.numeric) & !stc_nowell_ng_L,
                                     ~ .x / stc_nowell_ng_L)) %>%
                        select(-stc_nowell_ng_L) %>%
                        pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>%
                        pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>%
                        mutate(PTI = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
                        select(Sample, PTI, everything())
                }

                nowell_fish       <- .calc_nowell("fish")
                nowell_cladoceran <- .calc_nowell("crustacean")

                v$nowell_results <- Filter(Negate(is.null), list(
                    "Nowell PTI (Fish)"       = nowell_fish,
                    "Nowell PTI (Cladoceran)" = nowell_cladoceran
                ))
                if (length(v$nowell_results) == 0) v$nowell_results <- NULL
            } else {
                v$nowell_results <- NULL
            }

            setProgress(1.0, detail = "Done.")
            v$summary_log <- c(v$summary_log, "Toxicity calculations complete. See tabs for results.")
            log_event("INFO", "tox_calc",
                      paste0("Calculation complete — ", nrow(v$final_search_results), " compounds"),
                      session_id  = .session_id,
                      file_name   = input$user_file$name,
                      n_compounds = nrow(v$final_search_results),
                      n_samples   = ncol(v$user_data) - 1L)
            updateTabsetPanel(session, "main_tabs", selected = "Toxicity Plots")

        }, error = function(e) {
            v$summary_log <- c(v$summary_log, "An error occurred during toxicity calculation:", e$message)
            log_event("ERROR", "tox_calc", e$message,
                      session_id = .session_id,
                      file_name  = input$user_file$name)
        })
        }) # end withProgress
    })

    # ── Toxicity table outputs ────────────────────────────────────────────────
    # Format: 3 significant figures plain decimal (no scientific notation)
    sig3_render <- JS(
        "function(data, type) {",
        "  if (type === 'display' && data !== null && data !== '' &&",
        "      !isNaN(parseFloat(data)) && isFinite(data)) {",
        "    var v = parseFloat(data);",
        "    return v.toPrecision(3);",
        "  }",
        "  return data;",
        "}"
    )
    tox_dt_opts <- list(
        pageLength = 25,
        columnDefs = list(list(targets = "_all", render = sig3_render))
    )
    pti_style <- function(dt) {
        dt %>% formatStyle(
            "PTI",
            backgroundColor = styleInterval(c(0.1, 1.0),
                c("transparent", "#FFE0B2", "#FFCDD2")),
            color = styleInterval(c(0.1, 1.0),
                c("inherit", "#E65100", "#B71C1C")),
            fontWeight = styleInterval(c(0.1, 1.0),
                c("normal", "bold", "bold"))
        )
    }
    output$tox_table_algae      <- renderDT({
        req(v$tox_results)
        datatable(v$tox_results[["Algae toxicity"]],      rownames = FALSE, options = tox_dt_opts) %>% pti_style()
    })
    output$tox_table_crustacean <- renderDT({
        req(v$tox_results)
        datatable(v$tox_results[["Crustacean toxicity"]], rownames = FALSE, options = tox_dt_opts) %>% pti_style()
    })
    output$tox_table_fish       <- renderDT({
        req(v$tox_results)
        datatable(v$tox_results[["Fish toxicity"]],       rownames = FALSE, options = tox_dt_opts) %>% pti_style()
    })

    # ── Plot outputs ──────────────────────────────────────────────────────────
    output$algae_sample_plot         <- renderPlot({ req(v$plots); v$plots$algae_samples })
    output$crustacean_sample_plot    <- renderPlot({ req(v$plots); v$plots$crustacean_samples })
    output$fish_sample_plot          <- renderPlot({ req(v$plots); v$plots$fish_samples })
    output$algae_pollutant_plot      <- renderPlot({ req(v$plots); v$plots$algae_pollutants })
    output$crustacean_pollutant_plot <- renderPlot({ req(v$plots); v$plots$crustacean_pollutants })
    output$fish_pollutant_plot       <- renderPlot({ req(v$plots); v$plots$fish_pollutants })

    # =========================================================================
    # Download handler — multi-sheet Excel workbook
    # =========================================================================
    # Sheets written:
    #   "Original data"       — the user's uploaded concentration table.
    #   "Algae toxicity"      — TU per compound + PTI per station (algae).
    #   "Crustacean toxicity" — TU per compound + PTI per station (crustacean).
    #   "Fish toxicity"       — TU per compound + PTI per station (fish).
    #   "CASRN Report"        — matched and unmatched compound list.
    #   "Toxicity Coverage"   — which taxonomic groups have ECOTOX data for
    #                           each compound (colour-coded: green = data
    #                           available, red = no data).
    # ── Helper: build basic Excel workbook ───────────────────────────────────────
    .build_basic_wb <- function() {
        req(v$tox_results)

            # ── CASRN report ─────────────────────────────────────────────────
            casrn_report <- as.data.frame(v$final_search_results)[, c("PREFERRED_NAME", "CASRN")]
            names(casrn_report) <- c("Compound", "CASRN")
            casrn_report$Status <- "Matched"
            if (!is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0) {
                casrn_report <- rbind(casrn_report,
                    data.frame(Compound = v$manual_to_fill$PREFERRED_NAME,
                               CASRN    = NA_character_,
                               Status   = "Not found",
                               stringsAsFactors = FALSE))
            }
            casrn_report <- casrn_report[!duplicated(casrn_report$Compound), ]

            # ── Toxicity coverage matrix ──────────────────────────────────────
            p_final_rep <- as.data.frame(v$final_search_results) %>%
                mutate(cas_number = gsub("-", "", CASRN)) %>%
                distinct(PREFERRED_NAME, .keep_all = TRUE)

            # 3-state coverage: "ecotox" / "predicted" / "none"
            .cov_state <- function(grp) {
                has_ecotox    <- taxotox_reference$cas_number[taxotox_reference$ecotox_group == grp & !is.na(taxotox_reference$median_lc50_ng_L)]
                has_predicted <- taxotox_reference$cas_number[taxotox_reference$ecotox_group == grp & is.na(taxotox_reference$median_lc50_ng_L) & !is.na(taxotox_reference$predicted_lc50_ng_L)]
                function(cas) ifelse(cas %in% has_ecotox, "ecotox",
                                     ifelse(cas %in% has_predicted, "predicted", "none"))
            }
            cov_algae      <- .cov_state("algae")
            cov_crustacean <- .cov_state("crustacean")
            cov_fish       <- .cov_state("fish")

            coverage <- p_final_rep %>%
                mutate(Algae      = cov_algae(cas_number),
                       Crustacean = cov_crustacean(cas_number),
                       Fish       = cov_fish(cas_number)) %>%
                select(Compound = PREFERRED_NAME, Algae, Crustacean, Fish)

            if (!is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0) {
                coverage <- rbind(coverage,
                    data.frame(Compound   = v$manual_to_fill$PREFERRED_NAME,
                               Algae      = "none", Crustacean = "none", Fish = "none",
                               stringsAsFactors = FALSE))
            }

            # Replace state codes with display symbols
            coverage_display <- coverage %>%
                mutate(across(c(Algae, Crustacean, Fish), ~ case_when(
                    .x == "ecotox"    ~ "\u2713",
                    .x == "predicted" ~ "~",
                    TRUE              ~ "\u2717"
                )))

            # ── Build workbook ────────────────────────────────────────────────
            wb <- createWorkbook()

            # Legend sheet (first tab) — formulas and risk thresholds
            addWorksheet(wb, "Legend")
            legend_df <- data.frame(
                Item = c(
                    "Toxic Unit (TU)",
                    "Pollution Toxicity Index (PTI)",
                    "",
                    "PTI \u2265 1.0",
                    "PTI 0.1 \u2013 1.0",
                    "PTI < 0.1",
                    "",
                    "DISCLAIMER"
                ),
                Description = c(
                    "TU_ij = C_ij (ng/L) / median LC50_i (ng/L)   [compound i, sample j, per taxonomic group]",
                    "PTI_j = sum of all TU_ij   [sum over all detected compounds in sample j]",
                    "",
                    "Acute risk \u2014 mixture may cause lethal or immobilisation effects on sensitive species",
                    "Chronic risk \u2014 sub-lethal effects on sensitive species possible",
                    "Low risk \u2014 mixture unlikely to cause measurable toxicity",
                    "",
                    "TaxoTox is a research screening tool. Outputs are indicative only and should not be used as the sole basis for regulatory decisions or risk management actions. Predicted LC50 values (marked ~) carry additional uncertainty. Benchmark values reflect the listed jurisdiction and may not apply to all study contexts; regulatory benchmarks are periodically revised and users should verify currency with the issuing authority. Results depend on the ECOTOX database version installed."
                ),
                stringsAsFactors = FALSE
            )
            writeData(wb, "Legend", legend_df)
            setColWidths(wb, "Legend", cols = 1, widths = 25)
            setColWidths(wb, "Legend", cols = 2, widths = 80)

            for (sname in names(v$tox_results)) {
                addWorksheet(wb, sname)
                df <- v$tox_results[[sname]]
                # Transpose "Original data" to match the toxicity sheets:
                # rows = samples, columns = pollutants.
                if (sname == "Original data") {
                    df <- df %>%
                        pivot_longer(-PREFERRED_NAME,
                                     names_to  = "Sample",
                                     values_to = "Concentration_ng_L") %>%
                        pivot_wider(names_from  = PREFERRED_NAME,
                                    values_from = "Concentration_ng_L")
                }
                writeData(wb, sname, df)
            }
            addWorksheet(wb, "CASRN Report")
            writeData(wb, "CASRN Report", casrn_report)
            # ── Coverage summary rows ─────────────────────────────────────────
            cov_cols  <- coverage[, c("Algae", "Crustacean", "Fish")]
            n_total   <- nrow(coverage)
            n_ecotox  <- sapply(cov_cols, function(x) sum(x == "ecotox"))
            n_pred    <- sapply(cov_cols, function(x) sum(x == "predicted"))
            n_none    <- sapply(cov_cols, function(x) sum(x == "none"))
            pct_any   <- round(100 * (n_ecotox + n_pred) / n_total, 1)

            coverage_summary <- data.frame(
                Metric     = c("ECOTOX data (\u2713)",
                               "Predicted LC50 (~)",
                               "No data (\u2717)",
                               "% with any LC50"),
                Algae      = c(n_ecotox["Algae"],      n_pred["Algae"],      n_none["Algae"],      pct_any["Algae"]),
                Crustacean = c(n_ecotox["Crustacean"], n_pred["Crustacean"], n_none["Crustacean"], pct_any["Crustacean"]),
                Fish       = c(n_ecotox["Fish"],       n_pred["Fish"],       n_none["Fish"],       pct_any["Fish"]),
                stringsAsFactors = FALSE, row.names = NULL
            )

            # coverage_display is written starting at row 7
            # (4 summary rows + 1 header + 1 blank = rows 1-6 used above)
            COVERAGE_DATA_START_ROW <- 7L

            addWorksheet(wb, "Toxicity Coverage")
            writeData(wb, "Toxicity Coverage", coverage_summary, startRow = 1)
            writeData(wb, "Toxicity Coverage", coverage_display, startRow = COVERAGE_DATA_START_ROW)

            # Conditional formatting: green = ecotox or predicted, red = none
            green_style <- createStyle(fgFill = "#C6EFCE", halign = "CENTER", fontColour = "#276221")
            red_style   <- createStyle(fgFill = "#FFC7CE", halign = "CENTER", fontColour = "#9C0006")
            for (col_idx in 2:4) {
                vals       <- coverage[[names(coverage)[col_idx]]]
                green_rows <- which(vals %in% c("ecotox", "predicted")) + COVERAGE_DATA_START_ROW
                red_rows   <- which(vals == "none")                     + COVERAGE_DATA_START_ROW
                if (length(green_rows) > 0)
                    addStyle(wb, "Toxicity Coverage", green_style,
                             rows = green_rows, cols = col_idx, stack = FALSE)
                if (length(red_rows) > 0)
                    addStyle(wb, "Toxicity Coverage", red_style,
                             rows = red_rows,  cols = col_idx, stack = FALSE)
            }

            # ── Color legend below coverage table ─────────────────────────────
            legend_start_row <- COVERAGE_DATA_START_ROW + 1L + nrow(coverage_display) + 2L

            color_legend <- data.frame(
                Symbol  = c("\u2713", "~", "\u2717"),
                Meaning = c(
                    "ECOTOX data — median LC50/EC50 from the ECOTOX Knowledgebase used as TU denominator",
                    "Predicted LC50 — CompTox/OPERA model prediction used as gap-fill; treat results with lower confidence",
                    "No data — compound excluded from toxicity calculation for this taxonomic group"
                ),
                stringsAsFactors = FALSE
            )
            writeData(wb, "Toxicity Coverage", color_legend, startRow = legend_start_row)
            addStyle(wb, "Toxicity Coverage", green_style,
                     rows = legend_start_row + 1L, cols = 1L, stack = FALSE)
            addStyle(wb, "Toxicity Coverage", green_style,
                     rows = legend_start_row + 2L, cols = 1L, stack = FALSE)
            addStyle(wb, "Toxicity Coverage", red_style,
                     rows = legend_start_row + 3L, cols = 1L, stack = FALSE)
            setColWidths(wb, "Toxicity Coverage", cols = 1, widths = 20)
            setColWidths(wb, "Toxicity Coverage", cols = 2:4, widths = 15)
            setColWidths(wb, "Toxicity Coverage", cols = 2, widths = 90)

        wb
    }

    # ── Helper: build advanced Excel workbook ────────────────────────────────────
    .build_advanced_wb <- function() {
        wb <- createWorkbook()

        # ── HC5 sheets ───────────────────────────────────────────────────────────
        if (input$method_hc5) {
            if (!is.null(v$hc5_results) && length(v$hc5_results) > 0) {
                for (sname in names(v$hc5_results)) {
                    addWorksheet(wb, sname)
                    writeData(wb, sname, v$hc5_results[[sname]])
                }
            } else {
                addWorksheet(wb, "HC5 - No Data")
                writeData(wb, "HC5 - No Data", data.frame(
                    Note = "No HC5 values were available for any of the uploaded compounds. HC5 requires either ECOTOX SSD data (n \u2265 5) or a predicted LC50 for scaling.",
                    stringsAsFactors = FALSE
                ))
                setColWidths(wb, "HC5 - No Data", cols = 1, widths = 100)
            }
        }

        # ── Benchmark HI sheets ──────────────────────────────────────────────────
        if (input$method_benchmark) {
            if (!is.null(v$benchmark_results) && length(v$benchmark_results) > 0) {
                for (sname in names(v$benchmark_results)) {
                    addWorksheet(wb, sname)
                    writeData(wb, sname, v$benchmark_results[[sname]])
                }
            } else {
                addWorksheet(wb, "Bench. No Data")
                writeData(wb, "Bench. No Data", data.frame(
                    Note = paste0(
                        "No benchmark values were found for any of the uploaded compounds. ",
                        "The selected benchmark framework(s) do not cover the compounds in this dataset. ",
                        "Benchmark coverage is typically highest for legacy pesticides and priority substances; ",
                        "emerging contaminants are often not listed."
                    ),
                    stringsAsFactors = FALSE
                ))
                setColWidths(wb, "Bench. No Data", cols = 1, widths = 100)
            }
        }

        # ── IA sheets ────────────────────────────────────────────────────────────
        if (input$method_ia) {
            if (!is.null(v$ia_results) && length(v$ia_results) > 0) {
                for (sname in names(v$ia_results)) {
                    addWorksheet(wb, sname)
                    writeData(wb, sname, v$ia_results[[sname]])
                }
            } else {
                addWorksheet(wb, "IA - No Data")
                writeData(wb, "IA - No Data", data.frame(
                    Note = "No compounds with LC50 data were found for Independent Action calculation. IA requires at least one compound with a valid LC50 denominator.",
                    stringsAsFactors = FALSE
                ))
                setColWidths(wb, "IA - No Data", cols = 1, widths = 100)
            }
        }

        # ── CAMA sheets ──────────────────────────────────────────────────────────
        if (input$method_cama) {
            if (!is.null(v$cama_results) && length(v$cama_results) > 0) {
                for (sname in names(v$cama_results)) {
                    addWorksheet(wb, sname)
                    writeData(wb, sname, v$cama_results[[sname]])
                }
            } else {
                addWorksheet(wb, "CAMA - No Data")
                writeData(wb, "CAMA - No Data", data.frame(
                    Note = "No compounds with LC50 data were found for CAMA calculation. CAMA requires a valid LC50 denominator for at least one uploaded compound.",
                    stringsAsFactors = FALSE
                ))
                setColWidths(wb, "CAMA - No Data", cols = 1, widths = 100)
            }
        }

        # ── Nowell Taxon-Sensitive PTI sheets ───────────────────────────────────
        if (input$method_nowell) {
            if (!is.null(v$nowell_results) && length(v$nowell_results) > 0) {
                for (sname in names(v$nowell_results)) {
                    addWorksheet(wb, sname)
                    writeData(wb, sname, v$nowell_results[[sname]])
                }
            } else {
                addWorksheet(wb, "Nowell - No Data")
                writeData(wb, "Nowell - No Data", data.frame(
                    Note = paste0(
                        "No Taxon-Sensitive PTI (Nowell et al. 2014) values were found for the ",
                        "uploaded compounds. This method only covers the ~480 fish and ~460 ",
                        "cladoceran-genus pesticides/degradates Nowell et al. published STC values ",
                        "for; coverage is narrower than the Standard PTI or Benchmark HI methods."
                    ),
                    stringsAsFactors = FALSE
                ))
                setColWidths(wb, "Nowell - No Data", cols = 1, widths = 100)
            }
        }

        wb
    }

    # ── Download handler ─────────────────────────────────────────────────────────
    output$download_results <- downloadHandler(
        filename = function() {
            base <- tools::file_path_sans_ext(input$user_file$name)
            if (any_advanced()) paste0("taxotox_", base, ".zip")
            else                paste0("taxotox_", base, ".xlsx")
        },
        content = function(file) {
            req(v$tox_results)
            tryCatch({
            base_wb <- .build_basic_wb()

            if (any_advanced()) {
                # Write both workbooks to a temp dir and ZIP them
                tmp_dir   <- tempdir()
                base_name <- tools::file_path_sans_ext(input$user_file$name)
                basic_path    <- file.path(tmp_dir, paste0("taxotox_basic_",    base_name, ".xlsx"))
                advanced_path <- file.path(tmp_dir, paste0("taxotox_advanced_", base_name, ".xlsx"))

                saveWorkbook(base_wb, basic_path, overwrite = TRUE)

                adv_wb <- .build_advanced_wb()
                if (length(adv_wb$sheet_names) > 0) {
                    saveWorkbook(adv_wb, advanced_path, overwrite = TRUE)
                    zip(file, files = c(basic_path, advanced_path), flags = "-j")
                    showNotification(
                        paste0("This download contains TWO files zipped together: a ",
                               "BASIC workbook ('Algae/Crustacean/Fish toxicity' sheets, ",
                               "column 'PTI') and an ADVANCED workbook (e.g. 'Bench. US EPA ",
                               "(Fish)'). The Basic workbook's PTI column ALWAYS reports ",
                               "the Standard Concentration-Addition (median-LC50) result, ",
                               "regardless of which Advanced method(s) you selected. If you ",
                               "selected an Advanced method, your result is in the separate ",
                               "ADVANCED file - keep both and check the correct one."),
                        type     = "warning",
                        duration = 20
                    )
                } else {
                    # No advanced sheets produced — none of the selected methods had
                    # any compound coverage. Alert the user and fall back to basic only.
                    showNotification(
                        paste0("\u26a0 No advanced results produced: none of the selected ",
                               "methods (HC5 / Benchmark HI / IA / CAMA) found matching ",
                               "denominators for your compounds. Only the basic results ",
                               "file has been downloaded."),
                        type     = "warning",
                        duration = 12
                    )
                    saveWorkbook(base_wb, file, overwrite = TRUE)
                }
            } else {
                saveWorkbook(base_wb, file, overwrite = TRUE)
            }
            }, error = function(e) {
                log_event("ERROR", "export", e$message,
                          session_id = .session_id,
                          file_name  = input$user_file$name)
                stop(e)
            })
        }
    )

    # ── Manual Entry tab content ──────────────────────────────────────────────
    # Renders a self-contained entry form when unresolved compounds remain.
    # Shows a progress counter and a "nothing to do" confirmation when all
    # compounds are resolved. The form lets the user select a compound from a
    # dropdown, type its CASRN, and confirm — the compound is then removed from
    # the queue and added to final_search_results.
    output$manual_entry_tab_ui <- renderUI({

        # Case 1: fuzzy review still pending — ask user to finish that first
        if (!is.null(v$fuzzy_to_review) && nrow(v$fuzzy_to_review) > 0) {
            return(tagList(
                br(),
                div(class = "alert alert-info",
                    strong("Please submit the fuzzy matches in the 'CASRN Matching' tab first."))
            ))
        }

        n_remaining <- if (!is.null(v$manual_to_fill)) nrow(v$manual_to_fill) else 0
        n_added     <- nrow(v$manual_additions)

        # Case 2: nothing left to resolve
        if (n_remaining == 0) {
            return(tagList(
                br(),
                div(class = "alert alert-success",
                    icon("check-circle"),
                    strong(" All compounds resolved."),
                    if (n_added > 0)
                        paste0(" You manually added ", n_added, " CASRN(s).")
                    else
                        " No manual entry was needed."
                ),
                if (n_added > 0) tagList(h4("Manually added entries:"), DTOutput("manual_added_list"))
            ))
        }

        # Case 3: compounds still need a CASRN
        tagList(
            br(),
            div(class = "alert alert-warning",
                strong(paste0(n_remaining, " compound(s) still need a CASRN.")),
                " Enter each one below, then click 'Add CASRN'."
            ),
            fluidRow(
                column(5,
                    selectInput("manual_name", "Compound", choices = v$manual_to_fill$PREFERRED_NAME,
                                width = "100%")
                ),
                column(4,
                    textInput("manual_casrn", "CASRN (e.g., 123-45-6)", value = "", width = "100%")
                ),
                column(3,
                    br(),
                    actionButton("add_manual_casrn", "Add CASRN", class = "btn-primary",
                                 style = "margin-top:6px; width:100%;")
                )
            ),
            p(em("Tip: CASRNs can be found at ", tags$a("commonchemistry.cas.org", href = "https://commonchemistry.cas.org", target = "_blank"))),
            hr(),
            if (n_added > 0) tagList(h4("Added so far:"), DTOutput("manual_added_list"))
        )
    })
}

# =============================================================================
# Launch
# =============================================================================
shinyApp(ui = ui, server = server)
