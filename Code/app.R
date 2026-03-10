# =============================================================================
# TaxoTox — Taxonomic Toxicity Assessment Tool
# =============================================================================
#
# PURPOSE:
#   TaxoTox calculates the potential mixture toxicity of environmental pollutant
#   assemblages using the Toxic Unit (TU) framework and the Pollution Toxicity
#   Index (PTI). For each monitoring station and each taxonomic group (algae,
#   crustaceans, fish), it derives per-compound TUs by dividing the measured
#   concentration by the median chronic effect concentration retrieved from the
#   US EPA ECOTOX database. The PTI for a sample is then the sum of all TUs
#   multiplied by 100, providing a dimensionless index of chronic mixture risk.
#
# METHODOLOGY:
#   Toxic Unit (TU) for compound i in sample j:
#       TU_ij = C_ij / EC50_i
#   where C_ij is the measured concentration (ng/L) and EC50_i is the median
#   chronic effect concentration for the relevant taxonomic group (ng/L).
#
#   Pollution Toxicity Index (PTI) for sample j:
#       PTI_j = sum(TU_ij, for all i) * 100
#
# WORKFLOW:
#   1. User uploads a concentration table (CSV, TSV, XLS, XLSX).
#   2. The app automatically resolves compound names to CAS Registry Numbers
#      (CASRNs) via a four-layer lookup pipeline:
#        a. Exact match in the curated Known_CAS table (fastest, most reliable).
#        b. Exact name match in the full EPA DSSTox substance registry (covers
#           compounds in DSSTox but not yet in Known_CAS).
#        c. Live PubChem query via the webchem package — queries the PubChem
#           REST API by compound name and extracts the CAS number from the
#           synonym list (no API key required).
#        d. Jaro-Winkler fuzzy match against DSSTox (fallback for spelling
#           variants and abbreviations).
#      Matches from layers b–d are shown in the 'CASRN Matching' tab for user
#      review. Compounds unresolved by all layers go to Manual Entry.
#   3. The user clicks "Calculate Toxicity" to compute TUs and PTI values.
#   4. Results are displayed as interactive tables and bar/box plots, and
#      exported as a multi-sheet Excel workbook.
#
# DATA SOURCES:
#   - Known_CAS.fst   : Curated internal CASRN lookup table.
#   - DSSTox.fst      : EPA DSSTox substance registry (for fuzzy fallback).
#   - final_ecotox_data.fst : Pre-processed ECOTOX chronic endpoint data
#                             (median effect concentrations per compound and
#                             taxonomic group, in ng/L).
#
# AUTHORS: Noam [surname], Yair Suari
# CONTACT: [institution contact]
# LICENSE: [license]
# =============================================================================

# ── Dependencies ──────────────────────────────────────────────────────────────
# Core Shiny framework
library(shiny)
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
library(fst)
# String distance for fuzzy compound-name matching
library(stringdist)
# ECOTOX database interface (used during data pre-processing, not at runtime)
library(ECOTOXr)
library(RSQLite)
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
# webchem is optional — loaded softly so the app starts even when not installed.
# If absent, Layer 3 (PubChem lookup) is skipped and a notice is logged in the
# search summary. Install with: install.packages("webchem")
webchem_available <- requireNamespace("webchem", quietly = TRUE)
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
            actionButton("start_toxicity_calc", "2. Calculate Toxicity",
                         class = "btn-default btn-block",
                         style = "margin-top:10px;"),

            # ── Step 3: Export results ───────────────────────────────────────
            downloadButton("download_results", "3. Download Results",
                           class = "btn-default btn-block",
                           style = "margin-top:10px;"),
            hr(),

            # ── Advanced: fuzzy-match threshold ─────────────────────────────
            # Controls the Jaro-Winkler distance cutoff used when searching
            # DSSTox for compound names not found in the Known_CAS table.
            # Lower values require a closer match (stricter); higher values
            # accept more distant matches (more lenient, higher false-positive
            # rate). The default of 0.1 is conservative.
            numericInput("fuzzy_threshold",
                         "Fuzzy match sensitivity (0 = identical only, 0.3 = lenient)",
                         value = 0.1, min = 0.0, max = 0.5, step = 0.05),
            checkboxInput("use_pubchem", "Web search for pollutant names on/off", value = FALSE)
        ),

        mainPanel(
           tabsetPanel(id = "main_tabs",

               # ── Instructions tab ─────────────────────────────────────────
               tabPanel("Instructions",
                        h3("Workflow"),
                        p("Follow the three steps in the sidebar to analyse your data:"),
                        tags$ol(
                            tags$li(strong("Select your data layout"), " using the radio buttons, then click ",
                                    strong("'1. Upload Samples Data'"), " to load your file. The app will immediately attempt to identify all chemical compounds automatically."),
                            tags$li(strong("Resolve any unmatched compounds"), " in the 'CASRN Matching' tab if it appears, then click ",
                                    strong("'2. Calculate Toxicity'")),
                            tags$li("View results in the ", strong("'Toxicity Plots'"), " and ", strong("'Toxicity Tables'"),
                                    " tabs, then click ", strong("'3. Download Results'"), " to export an Excel workbook.")
                        ),

                        h4("Input File Format"),
                        p("Accepted formats: CSV, TSV, TXT, XLS, XLSX. Concentrations must be in ", strong("ng/L.")),

                        h5("Layout A — Pollutant names in first column (default)"),
                        p("Each row is a compound; each subsequent column is a monitoring station or sample."),
                        tags$pre(paste(
                            "Compound    | Station1 | Station2 | ...",
                            "Caffeine    | 10.5     | 18.2     | ...",
                            "Atrazine    |  2.1     |  0.5     | ...",
                            "Bisphenol A |  2.1     |  0.3     | ...",
                            sep = "\n"
                        )),

                        h5("Layout B — Station names in first column"),
                        p("Each row is a monitoring station; compound names are in the header row. The app transposes this automatically."),
                        tags$pre(paste(
                            "Station  | Caffeine | Atrazine | Bisphenol A | ...",
                            "Station1 | 10.5     |  2.1     |  2.1        | ...",
                            "Station2 | 18.2     |  0.5     |  0.3        | ...",
                            sep = "\n"
                        )),

                        h4("CASRN Identification"),
                        p("On upload, the app automatically resolves compound names to CAS Registry Numbers (CASRNs) via a four-layer pipeline:"),
                        tags$ol(
                            tags$li(strong("Known_CAS (exact):"), " a curated internal table of commonly monitored compounds. Matches are confirmed automatically — no review needed."),
                            tags$li(strong("DSSTox (exact):"), " the full EPA DSSTox substance registry. Compounds found here by exact name are added to the review table pre-checked."),
                            tags$li(strong("PubChem (API):"), " live query of the PubChem REST API via the webchem package. Requires an internet connection; no API key needed. Matches are pre-checked."),
                            tags$li(strong("DSSTox (fuzzy):"), " Jaro-Winkler string-distance search against DSSTox for spelling variants and abbreviations. These matches start unchecked — review each suggestion carefully.")
                        ),
                        p("All layer 2–4 candidates are shown together in the 'CASRN Matching' tab, grouped by source. Accept or reject each one using the checkboxes, then click Submit. Compounds rejected or not found by any layer go to the 'Manual Entry' tab."),

                        h4("Toxicity Calculation"),
                        p("For each compound and taxonomic group (algae, crustaceans, fish), the app retrieves chronic effect concentrations from the pre-processed ECOTOX database and computes the Toxic Unit:"),
                        tags$pre("TU = measured concentration (ng/L) / median ECOTOX effect concentration (ng/L)"),
                        p("The Pollution Toxicity Index (PTI) for each sample is the sum of all individual TUs multiplied by 100:"),
                        tags$pre("PTI = \u03a3 TU\u1d62 \u00d7 100"),
                        p("A PTI > 1 suggests the mixture may pose a chronic risk to the relevant taxonomic group."),

                        p(em("Tip: Lower the 'Fuzzy match sensitivity' value in the sidebar for stricter compound-name matching; raise it if legitimate compounds are being missed."))
                        ),
               tabPanel("CASRN Matching",
                        h4("Search Summary"),
                        verbatimTextOutput("cas_search_summary"),
                        hr(),
                        h4("Candidate Matches"),
                        p("DSSTox exact and PubChem matches (100%) are pre-checked. DSSTox fuzzy (partial) matches start unchecked — review and check any you wish to accept, then click Submit."),
                        DTOutput("fuzzy_match_table"),
                        br(),
                        actionButton("submit_fuzzy_matches", "Submit Selected Lines", class = "btn-primary"),
                        hr()
                        ),

               # ── Manual Entry tab ─────────────────────────────────────────
               # Shown when one or more compounds could not be matched
               # automatically. The tab is always present but shows a
               # "nothing to do" message when all compounds are resolved.
               tabPanel("Manual Entry",
                        uiOutput("manual_entry_tab_ui")
                        ),

               tabPanel("Toxicity Plots",
                        h4("Top 10 Riskiest Samples (by Risk Assessment score)"),
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
                        ),
           )
        )
    )
)

# =============================================================================
# SERVER
# =============================================================================

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
        plots                = NULL,     # named list of ggplot objects
        manual_additions     = data.table(PREFERRED_NAME = character(), CASRN = character()),
        pending_unfound      = NULL       # saved when orientation-check modal is shown
    )

    # ── Static reference data (loaded once at session start) ─────────────────
    # Known_CAS       : Curated table of compound name → CASRN mappings built
    #                   from prior runs and administrator-approved user entries.
    # DSSTox          : EPA DSSTox substance registry; used as fallback when a
    #                   compound name is absent from Known_CAS.
    # final_ecotox_data: Pre-processed ECOTOX data containing median chronic
    #                   effect concentrations (conc_ng_L) per compound
    #                   (cas_number) and taxonomic group (ecotox_group).
    Known_CAS         <- read.fst("../Data/Known_CAS.fst",         as.data.table = TRUE)
    DSSTox            <- read.fst("../Data/DSSTox.fst",            as.data.table = TRUE)
    final_ecotox_data <- read.fst("../Data/final_ecotox_data.fst", as.data.table = TRUE)

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

    # ── append_to_temp_cas ────────────────────────────────────────────────────
    # Persists newly confirmed PREFERRED_NAME / CASRN pairs to a session-level
    # temporary FST file (temp_CAS.fst). This file accumulates user-confirmed
    # fuzzy matches and manual entries across a session and can be reviewed by
    # administrators to expand the Known_CAS database over time.
    append_to_temp_cas <- function(new_data) {
        temp_cas_path <- "../Data/temp_CAS.fst"
        if (!all(c("PREFERRED_NAME", "CASRN") %in% names(new_data)))
            stop("New data must contain PREFERRED_NAME and CASRN columns.")
        temp_cas_dt <- if (file.exists(temp_cas_path)) {
            read.fst(temp_cas_path, as.data.table = TRUE)
        } else {
            data.table(PREFERRED_NAME = character(), CASRN = character())
        }
        updated_dt <- rbindlist(list(temp_cas_dt, new_data), use.names = TRUE, fill = TRUE)
        updated_dt <- unique(updated_dt, by = c("PREFERRED_NAME", "CASRN"))
        write.fst(updated_dt, temp_cas_path)
    }

    # ── fuzzy_match_non_interactive ───────────────────────────────────────────
    # For each name in source_names that was not found by exact lookup, computes
    # the Jaro-Winkler string distance to all entries in target_dt[[match_col]]
    # and returns the single closest match if it falls within 'threshold'.
    #
    # The Jaro-Winkler metric (range 0–1) gives extra weight to prefix
    # similarity, making it well-suited to chemical nomenclature where names
    # commonly share long common prefixes (e.g. "benzo[a]pyrene" vs
    # "benzo[b]pyrene"). A threshold of 0.1 corresponds to very close matches;
    # 0.3 accepts moderately different strings.
    #
    # Returns a data.table with columns: source_name, matched_name, distance,
    # CASRN. Returns an empty data.table if no matches pass the threshold.
    fuzzy_match_non_interactive <- function(source_names, target_dt, match_col, threshold = 0.2) {
        results <- list()
        for (i in seq_along(source_names)) {
            name <- source_names[i]
            if (is.na(name)) next
            valid_targets <- target_dt[[match_col]][!is.na(target_dt[[match_col]])]
            if (length(valid_targets) == 0) next
            distances <- stringdist(name, valid_targets, method = "jw")
            if (length(distances) > 0 && !all(is.na(distances))) {
                min_dist <- min(distances, na.rm = TRUE)
                if (!is.na(min_dist) && min_dist <= threshold) {
                    best_idx  <- which.min(distances)
                    best_match   <- valid_targets[best_idx]
                    casrn_number <- target_dt[get(match_col) == best_match, CASRN][1]
                    results[[length(results) + 1]] <- data.table(
                        source_name  = name,
                        matched_name = best_match,
                        distance     = min_dist,
                        CASRN        = casrn_number
                    )
                }
            }
        }
        if (length(results) > 0) return(rbindlist(results))
        return(data.table())
    }

    # ── pubchem_lookup ────────────────────────────────────────────────────────
    # Queries the PubChem REST API (via the webchem package) for CAS Registry
    # Numbers for compound names that could not be resolved by Known_CAS or
    # DSSTox exact match.
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
        # Return an empty table immediately if webchem is not installed
        if (!webchem_available) {
            return(data.table(source_name  = character(),
                              matched_name = character(),
                              CASRN        = character(),
                              source       = character()))
        }
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
    # Runs layers 2–4 of the CASRN pipeline on `unfound` (compounds not matched
    # in Layer 1 / Known_CAS). Builds the review table and updates v$summary_log.
    # Closes any open modal before starting. Accesses v, input, session via closure.
    run_casrn_layers2to4 <- function(unfound) {
        removeModal()
        withProgress(message = "Searching for CASRNs (layers 2\u20134)", value = 0, {
        tryCatch({

            # ── Layer 2: DSSTox exact ────────────────────────────────────────
            dsstox_exact_rows <- data.table()
            if (length(unfound) > 0) {
                setProgress(0.15, detail = paste("Layer 2: DSSTox exact (",
                                                 length(unfound), "remaining)..."))
                dsstox_exact_hits <- unique(
                    DSSTox[DSSTox$PREFERRED_NAME %in% unfound, .(PREFERRED_NAME, CASRN)],
                    by = "PREFERRED_NAME"
                )
                if (nrow(dsstox_exact_hits) > 0) {
                    dsstox_exact_rows <- data.table(
                        source_name  = dsstox_exact_hits$PREFERRED_NAME,
                        matched_name = dsstox_exact_hits$PREFERRED_NAME,
                        distance     = 0,
                        CASRN        = dsstox_exact_hits$CASRN,
                        source       = "DSSTox (exact)"
                    )
                    unfound <- unfound[!unfound %in% dsstox_exact_hits$PREFERRED_NAME]
                    v$summary_log <- c(v$summary_log,
                                       paste("Layer 2 (DSSTox exact):", nrow(dsstox_exact_rows),
                                             "matches;", length(unfound), "remaining."))
                } else {
                    v$summary_log <- c(v$summary_log,
                                       paste("Layer 2 (DSSTox exact): no matches;",
                                             length(unfound), "remaining."))
                }
            }

            # ── Layer 3: PubChem via webchem ─────────────────────────────────
            pubchem_rows <- data.table()
            if (length(unfound) > 0) {
                if (!webchem_available || !isTRUE(input$use_pubchem)) {
                    v$summary_log <- c(v$summary_log, if (!webchem_available)
                        "Layer 3 (PubChem): skipped \u2014 webchem not installed." else
                        "Layer 3 (PubChem): skipped (web search is off).")
                } else {
                    setProgress(0.45, detail = paste("Layer 3: PubChem query (",
                                                     length(unfound), "compounds)..."))
                    v$summary_log <- c(v$summary_log,
                                       paste("Layer 3 (PubChem): querying",
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
                                           paste("Layer 3 (PubChem):", nrow(pubchem_rows),
                                                 "matches;", length(unfound), "remaining."))
                    } else {
                        v$summary_log <- c(v$summary_log,
                                           paste("Layer 3 (PubChem): no matches;",
                                                 length(unfound), "remaining."))
                    }
                }
            }

            # ── Layer 4: DSSTox fuzzy ────────────────────────────────────────
            fuzzy_rows <- data.table()
            if (length(unfound) > 0) {
                setProgress(0.70, detail = paste("Layer 4: fuzzy matching (",
                                                 length(unfound), "remaining)..."))
                v$summary_log <- c(v$summary_log,
                                   paste("Layer 4 (DSSTox fuzzy): searching",
                                         length(unfound), "compounds..."))
                fm <- fuzzy_match_non_interactive(
                    unfound, DSSTox, "PREFERRED_NAME",
                    threshold = input$fuzzy_threshold
                )
                fm <- fm[fm$distance > 0, ]
                if (nrow(fm) > 0) {
                    fuzzy_rows <- data.table(
                        source_name  = fm$source_name,
                        matched_name = fm$matched_name,
                        distance     = fm$distance,
                        CASRN        = fm$CASRN,
                        source       = "DSSTox (fuzzy)"
                    )
                    v$summary_log <- c(v$summary_log,
                                       paste("Layer 4 (DSSTox fuzzy):", nrow(fuzzy_rows),
                                             "partial matches for review."))
                } else {
                    v$summary_log <- c(v$summary_log, "Layer 4 (DSSTox fuzzy): no matches.")
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

            truly_unfound <- unfound[!unfound %in% fuzzy_rows$source_name]
            manual_rows <- if (length(truly_unfound) > 0) {
                data.table(
                    source_name  = truly_unfound,
                    matched_name = NA_character_,
                    distance     = NA_real_,
                    CASRN        = NA_character_,
                    source       = "Manual"
                )
            } else { data.table() }

            if (length(truly_unfound) > 0)
                v$summary_log <- c(v$summary_log,
                                   paste(length(truly_unfound),
                                         "compound(s) need manual CASRN entry."))

            review_rows <- rbindlist(
                list(manual_rows, fuzzy_rows, pubchem_rows, dsstox_exact_rows, known_cas_rows),
                use.names = TRUE, fill = TRUE
            )

            v$manual_to_fill <- data.table(PREFERRED_NAME = character())

            if (nrow(review_rows) > 0) {
                v$fuzzy_to_review <- review_rows
                n_need_action <- nrow(manual_rows) + nrow(fuzzy_rows)
                v$summary_log <- c(v$summary_log,
                                   paste("Done. Review table:", nrow(review_rows),
                                         "rows total;", n_need_action,
                                         "need attention (Manual / fuzzy)."),
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

            n_action <- nrow(manual_rows) + nrow(fuzzy_rows)
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
        n_show     <- min(10, nrow(tox_cal))
        top_samples <- tox_cal %>% arrange(desc(RQtg)) %>% slice(1:n_show)
        top3_label  <- paste("Top 3:", paste(head(top_samples$Sample, min(3, nrow(top_samples))), collapse = ", "))
        ggplot(top_samples, aes(x = reorder(Sample, RQtg), y = RQtg)) +
            geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
            labs(title    = paste(group_label, "\u2014 Top", n_show, "Riskiest Samples"),
                 subtitle = top3_label,
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
    #   2. Four-layer CASRN resolution pipeline:
    #        Layer 1 — exact match in Known_CAS (auto-confirmed, no review).
    #        Layer 2 — exact match in DSSTox (pre-checked in review table).
    #        Layer 3 — PubChem live API query via webchem (pre-checked).
    #        Layer 4 — Jaro-Winkler fuzzy match in DSSTox (unchecked).
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

            # ── Orientation check ────────────────────────────────────────────
            # If >75 % of compound names were not found in Known_CAS the first
            # column is likely station names rather than compound names — the
            # user may have the wrong layout selected.
            pct_unfound <- length(unfound) / max(length(v$p_vector), 1)
            if (pct_unfound > 0.75) {
                v$pending_unfound <- unfound
                pct_label <- paste0(round(pct_unfound * 100), "%")
                showModal(modalDialog(
                    title = "\u26a0\ufe0f Possible orientation mismatch",
                    p(strong(paste0(length(unfound), " of ", length(v$p_vector),
                                   " (", pct_label, ") compound names were not recognised.")),
                      " This often means the file is in the opposite layout to the one selected:"),
                    tags$ul(
                        tags$li("If your ", strong("first column contains station names"),
                                " (Layout B), select ",
                                strong("\u201cStation names in first column\u201d"), " in the sidebar."),
                        tags$li("If your ", strong("first column contains compound names"),
                                " (Layout A), keep the current selection.")
                    ),
                    p("Would you like to ", strong("transpose the table and retry"), "?"),
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
        })

        }) # end withProgress

        # Continue to layers 2–4 (only reached when orientation check passes)
        run_casrn_layers2to4(v$pending_unfound)
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

    # ── CASRN Matching tab: unified compound review table ────────────────────
    # Renders ALL compounds in a single DataTable, ordered so the ones needing
    # the most user attention appear first:
    #
    #   Row section  | Source          | Accept col       | CASRN col
    #   -------------|-----------------|------------------|------------------
    #   1 (top)      | Manual          | checkbox (off)   | editable textInput
    #   2            | DSSTox (fuzzy)  | checkbox (off)   | editable textInput
    #                |                 |                  |   (pre-filled with
    #                |                 |                  |    fuzzy suggestion)
    #   3            | PubChem         | checkbox (on)    | static text
    #   4            | DSSTox (exact)  | checkbox (on)    | static text
    #   5 (bottom)   | Known_CAS       | "✓ Auto"         | static text
    #
    # The Accept column is rendered narrow (≈50 px).
    # Shiny input bindings are re-attached after every DataTables redraw so
    # that checkboxes and text inputs inside the table remain reactive.
    output$fuzzy_match_table <- renderDT({
        req(v$fuzzy_to_review)
        df <- as.data.frame(v$fuzzy_to_review)

        # ── Sort: Manual → fuzzy → PubChem → DSSTox exact → Known_CAS ────────
        df$source_order <- dplyr::case_when(
            df$source == "Manual"         ~ 1L,
            df$source == "DSSTox (fuzzy)" ~ 2L,
            df$source == "PubChem"        ~ 3L,
            df$source == "DSSTox (exact)" ~ 4L,
            df$source == "Known_CAS"      ~ 5L,
            TRUE                          ~ 6L
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
        # PubChem / DSSTox exact: checkbox pre-checked.
        # Manual / DSSTox fuzzy: checkbox unchecked.
        accept_col <- sapply(seq_len(nrow(df)), function(i) {
            src <- df$source[i]
            if (src == "Known_CAS") {
                "<span style='color:#2a7a2a; font-weight:bold;'>\u2713 Auto</span>"
            } else {
                pre <- src %in% c("PubChem", "DSSTox (exact)")
                as.character(checkboxInput(paste0("fuzzy_cb_", i), label = NULL, value = pre))
            }
        })

        # ── CASRN column ──────────────────────────────────────────────────────
        # Manual and DSSTox fuzzy: editable textInput (fuzzy pre-filled).
        # All others: static text.
        casrn_col <- sapply(seq_len(nrow(df)), function(i) {
            src <- df$source[i]
            if (src %in% c("Manual", "DSSTox (fuzzy)")) {
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
            df$source == "Manual"         ~ 1L,
            df$source == "DSSTox (fuzzy)" ~ 2L,
            df$source == "PubChem"        ~ 3L,
            df$source == "DSSTox (exact)" ~ 4L,
            df$source == "Known_CAS"      ~ 5L,
            TRUE                          ~ 6L)
        df <- df[order(df$source_order,
                       ifelse(is.na(df$distance), Inf, df$distance)), ]

        confirmed  <- data.table()
        unresolved <- character()

        for (i in seq_len(nrow(df))) {
            src <- df$source[i]

            # Known_CAS rows are already auto-confirmed — skip them here.
            if (src == "Known_CAS") next

            checked <- isTRUE(input[[paste0("fuzzy_cb_", i)]])

            if (src %in% c("Manual", "DSSTox (fuzzy)")) {
                # CASRN comes from the editable text input.
                raw_val  <- input[[paste0("casrn_input_", i)]]
                casrn_val <- trimws(if (is.null(raw_val)) "" else raw_val)
                if (checked && nchar(casrn_val) > 0) {
                    confirmed <- rbindlist(list(confirmed,
                        data.table(PREFERRED_NAME = df$source_name[i],
                                   CASRN          = casrn_val)),
                        fill = TRUE)
                } else {
                    unresolved <- c(unresolved, df$source_name[i])
                }
            } else {
                # PubChem / DSSTox exact: CASRN is pre-filled in the data.
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

        # If unresolved compounds remain, navigate to the Manual Entry tab
        if (!is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0)
            updateTabsetPanel(session, "main_tabs", selected = "Manual Entry")
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
        casrn <- trimws(input$manual_casrn)
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

        tryCatch({
            # ── Build the matched compound lookup table ───────────────────────
            # Remove hyphens from CASRNs for consistent matching with the
            # ECOTOX dataset (which stores CAS numbers without hyphens).
            p_final <- as.data.table(v$final_search_results) %>%
                na.omit() %>%
                distinct(PREFERRED_NAME, .keep_all = TRUE) %>%
                mutate(cas_number = gsub("-", "", CASRN)) %>%
                select(PREFERRED_NAME, cas_number)

            cas_search_list <- as.vector(p_final$cas_number)

            # ── Subset ECOTOX to matched compounds ───────────────────────────
            endpoint_data <- final_ecotox_data %>%
                mutate(cas_number = as.character(cas_number)) %>%
                filter(cas_number %in% cas_search_list)

            # ── Per-group TU / PTI calculation (identical for all three groups)
            # Steps (illustrated for algae; crustacean and fish are analogous):
            #   a. Filter endpoint_data to the group.
            #   b. Join compound names and compute median_conc per compound.
            #   c. Join user concentrations, divide each station column by
            #      median_conc → TU per compound per station.
            #   d. Pivot to (Sample × Compound) and sum TUs × 100 → PTI.

            # Algae ────────────────────────────────────────────────────────────
            endpoint_data_algae <- endpoint_data %>%
                filter(ecotox_group == "algae") %>%
                left_join(p_final, by = "cas_number", relationship = "many-to-many") %>%
                group_by(PREFERRED_NAME, cas_number) %>%
                summarise(median_conc = median(conc_ng_L, na.rm = TRUE), .groups = "drop")

            tox_cal_algae <- endpoint_data_algae %>%
                left_join(v$user_data, by = "PREFERRED_NAME") %>%
                mutate(across(where(is.numeric) & !c(cas_number, median_conc), ~ .x / median_conc)) %>%
                select(-c(cas_number, median_conc)) %>%
                pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>%
                pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>%
                mutate(RQtg = rowSums(across(where(is.numeric)), na.rm = TRUE))

            # Crustacean ───────────────────────────────────────────────────────
            endpoint_data_crustacean <- endpoint_data %>%
                filter(ecotox_group == "crustacean") %>%
                left_join(p_final, by = "cas_number", relationship = "many-to-many") %>%
                group_by(PREFERRED_NAME, cas_number) %>%
                summarise(median_conc = median(conc_ng_L, na.rm = TRUE), .groups = "drop")

            tox_cal_crustacean <- endpoint_data_crustacean %>%
                left_join(v$user_data, by = "PREFERRED_NAME") %>%
                mutate(across(where(is.numeric) & !c(cas_number, median_conc), ~ .x / median_conc)) %>%
                select(-c(cas_number, median_conc)) %>%
                pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>%
                pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>%
                mutate(RQtg = rowSums(across(where(is.numeric)), na.rm = TRUE))

            # Fish ─────────────────────────────────────────────────────────────
            endpoint_data_fish <- endpoint_data %>%
                filter(ecotox_group == "fish") %>%
                left_join(p_final, by = "cas_number", relationship = "many-to-many") %>%
                group_by(PREFERRED_NAME, cas_number) %>%
                summarise(median_conc = median(conc_ng_L, na.rm = TRUE), .groups = "drop")

            tox_cal_fish <- endpoint_data_fish %>%
                left_join(v$user_data, by = "PREFERRED_NAME") %>%
                mutate(across(where(is.numeric) & !c(cas_number, median_conc), ~ .x / median_conc)) %>%
                select(-c(cas_number, median_conc)) %>%
                pivot_longer(cols = 2:length(.), names_to = "Sample", values_to = "TU") %>%
                pivot_wider(names_from = PREFERRED_NAME, values_from = TU) %>%
                mutate(RQtg = rowSums(across(where(is.numeric)), na.rm = TRUE))

            # ── Rename RQtg for display and reorder columns ──────────────────
            # "Risk assessment" is used as the column header in tables and the
            # Excel export; RQtg is retained internally for plotting functions.
            format_tox_output <- function(df) {
                df %>% rename("Risk assessment" = RQtg) %>%
                    select(Sample, `Risk assessment`, everything())
            }

            v$tox_results <- list(
                "Original data"       = v$user_data,
                "Algae toxicity"      = format_tox_output(tox_cal_algae),
                "Crustacean toxicity" = format_tox_output(tox_cal_crustacean),
                "Fish toxicity"       = format_tox_output(tox_cal_fish)
            )

            v$plots <- list(
                algae_samples         = make_sample_risk_plot(tox_cal_algae,      "Algae"),
                crustacean_samples    = make_sample_risk_plot(tox_cal_crustacean, "Crustacean"),
                fish_samples          = make_sample_risk_plot(tox_cal_fish,       "Fish"),
                algae_pollutants      = make_pollutant_plot(tox_cal_algae,        "Algae"),
                crustacean_pollutants = make_pollutant_plot(tox_cal_crustacean,   "Crustacean"),
                fish_pollutants       = make_pollutant_plot(tox_cal_fish,         "Fish")
            )

            v$summary_log <- c(v$summary_log, "Toxicity calculations complete. See tabs for results.")
            updateTabsetPanel(session, "main_tabs", selected = "Toxicity Plots")

        }, error = function(e) {
            v$summary_log <- c(v$summary_log, "An error occurred during toxicity calculation:", e$message)
        })
    })

    # ── Toxicity table outputs ────────────────────────────────────────────────
    sci_render <- JS(
        "function(data, type) {",
        "  if (type === 'display' && data !== null && data !== '' &&",
        "      !isNaN(parseFloat(data)) && isFinite(data)) {",
        "    return parseFloat(data).toExponential(4);",
        "  }",
        "  return data;",
        "}"
    )
    tox_dt_opts <- list(
        pageLength = 25,
        columnDefs = list(list(targets = "_all", render = sci_render))
    )
    output$tox_table_algae      <- renderDT({
        req(v$tox_results)
        datatable(v$tox_results[["Algae toxicity"]],      rownames = FALSE, options = tox_dt_opts)
    })
    output$tox_table_crustacean <- renderDT({
        req(v$tox_results)
        datatable(v$tox_results[["Crustacean toxicity"]], rownames = FALSE, options = tox_dt_opts)
    })
    output$tox_table_fish       <- renderDT({
        req(v$tox_results)
        datatable(v$tox_results[["Fish toxicity"]],       rownames = FALSE, options = tox_dt_opts)
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
    output$download_results <- downloadHandler(
        filename = function() {
            paste0("taxotox_", tools::file_path_sans_ext(input$user_file$name), ".xlsx")
        },
        content = function(file) {
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

            ecotox_algae      <- unique(final_ecotox_data$cas_number[final_ecotox_data$ecotox_group == "algae"])
            ecotox_crustacean <- unique(final_ecotox_data$cas_number[final_ecotox_data$ecotox_group == "crustacean"])
            ecotox_fish       <- unique(final_ecotox_data$cas_number[final_ecotox_data$ecotox_group == "fish"])

            coverage <- p_final_rep %>%
                mutate(Algae      = cas_number %in% ecotox_algae,
                       Crustacean = cas_number %in% ecotox_crustacean,
                       Fish       = cas_number %in% ecotox_fish) %>%
                select(Compound = PREFERRED_NAME, Algae, Crustacean, Fish)

            if (!is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0) {
                coverage <- rbind(coverage,
                    data.frame(Compound   = v$manual_to_fill$PREFERRED_NAME,
                               Algae      = FALSE, Crustacean = FALSE, Fish = FALSE,
                               stringsAsFactors = FALSE))
            }

            # Replace logical with tick/cross symbols for readability
            coverage_display <- coverage %>%
                mutate(Algae      = ifelse(Algae,      "\u2713", "\u2717"),
                       Crustacean = ifelse(Crustacean, "\u2713", "\u2717"),
                       Fish       = ifelse(Fish,       "\u2713", "\u2717"))

            # ── Build workbook ────────────────────────────────────────────────
            wb <- createWorkbook()
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
            addWorksheet(wb, "Toxicity Coverage")
            writeData(wb, "Toxicity Coverage", coverage_display)

            # Conditional formatting: green = ECOTOX data present, red = absent
            green_style <- createStyle(fgFill = "#C6EFCE", halign = "CENTER", fontColour = "#276221")
            red_style   <- createStyle(fgFill = "#FFC7CE", halign = "CENTER", fontColour = "#9C0006")
            for (col_idx in 2:4) {
                col_name <- names(coverage)[col_idx]
                for (row_idx in seq_len(nrow(coverage))) {
                    addStyle(wb, "Toxicity Coverage",
                             if (isTRUE(coverage[[col_name]][row_idx])) green_style else red_style,
                             rows = row_idx + 1, cols = col_idx)
                }
            }
            saveWorkbook(wb, file, overwrite = TRUE)
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
