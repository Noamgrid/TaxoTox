#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(openxlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(tidyverse)
library(ggpattern)
library(ggpubr)
library(ggpmisc)
library(hrbrthemes)
library(remotes)
library(purrr)  
library(readxl)
library(stringr)
library(RSQLite)
library(data.table)
library(stringdist)
library(fst)
library(ECOTOXr)
library(DT)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("TaxoTox"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("user_file", "Choose your data file",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv",
                          ".xlsx",
                          ".xls",
                          ".txt",
                          ".tsv"
                      )),
            numericInput("fuzzy_threshold",
                         "Fuzzy match sensitivity (0 = identical only, 0.3 = lenient)",
                         value = 0.1, min = 0.0, max = 0.5, step = 0.05),
            actionButton("start_processing", "2. Find CASRNs",
                         title = "CASRN = Chemical Abstracts Service Registry Number: a unique identifier for each chemical compound"),
            hr(),
            uiOutput("interactive_casrn_ui"),
            hr(),
            actionButton("start_toxicity_calc", "3. Calculate Toxicity"),
            hr(),
            downloadButton("download_results", "4. Download Results")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           tabsetPanel(id = "main_tabs",
               tabPanel("Instructions",
                        h3("Workflow"),
                        p("1. Upload your data file using the 'Choose your data file'."),
                        h4("Expected File Format"),
                        p("The input file should be a table where the first column contains the chemical compound names. Subsequent columns should contain concentration data for each sample or site (in ng/L)."),
                        p(strong("Example Structure:")),
                        tags$pre(paste(
                            "Compound    |Sample1  |Sample2  |......",
                            "Caffeine    |10.5     |18.2     |",
                            "Atrazine    |2.1      |0.5      |1.1",
                            "Bisphenol A |2.1      |0.3      |0.5",
                            sep = "\n"
                        )),
                        p("2. Click 'Find CASRNs' to start the search. The app will find exact matches from the Known CASRNs database, then run a fuzzy search on anything unmatched."),
                        p("3. 'Interactive CASRN Matching' tab will enabble you to approve identification of polutant names similar to the ones in the input file. This match is resolved using fuzzy matches to the dsstox database of pollutants names and CASRNs."), 
                        p("4. For any remaining compounds, use the 'Manual Entry' section in the sidebar to add CASRNs one by one."),
                        p("5. Once all CASRNs are resolved, click 'Calculate Toxicity'."),
                        p("6. View the results in the 'Toxicity Plots' and 'Toxicity Tables' tabs."),
                        p("7. Click 'Download Results' to get the final Excel file.")
                        ),
               tabPanel("CASRN Search Summary", verbatimTextOutput("cas_search_summary")),
               tabPanel("Interactive CASRN Matching",
                        h4("Fuzzy Matches"),
                        p("All matches are pre-checked (accepted). Uncheck any row you want to reject, then click the button."),
                        DTOutput("fuzzy_match_table"),
                        br(),
                        actionButton("submit_fuzzy_matches", "Submit Selected Lines", class = "btn-primary"),
                        hr()
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
               tabPanel("Admin",
                        passwordInput("admin_password", "Admin Password"),
                        actionButton("admin_login", "Login"),
                        uiOutput("admin_panel_ui")
                        )
           )
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    
    v <- reactiveValues(
        user_data = NULL,
        p_vector = NULL,
        summary_log = c("Welcome to TaxoTox! Please upload a file to begin."),
        casrn_results = NULL,
        fuzzy_to_review = NULL,
        manual_to_fill = NULL,
        final_search_results = data.table(),
        tox_results = NULL,
        plots = NULL,
        manual_additions = data.table(PREFERRED_NAME=character(), CASRN=character()),
        admin_logged_in = FALSE
    )

    # Load static data once
    Known_CAS <- read.fst("../Data/Known_CAS.fst", as.data.table = TRUE)
    DSSTox <- read.fst("../Data/DSSTox.fst", as.data.table = TRUE)
    final_ecotox_data <- read.fst("../Data/final_ecotox_data.fst", as.data.table = TRUE)


    # Functions adapted from original script
    
    load_user_file <- function(file_path) {
        ext <- tools::file_ext(file_path)
        df <- switch(tolower(ext),
                     "csv" = read.csv(file_path, stringsAsFactors = FALSE),
                     "txt" = read.delim(file_path, stringsAsFactors = FALSE),
                     "tsv" = read.delim(file_path, sep = "\t", stringsAsFactors = FALSE),
                     "xls" = readxl::read_excel(file_path),
                     "xlsx" = readxl::read_excel(file_path),
                     stop("Unsupported file type: ", ext)
        )
        return(as.data.table(df))
    }

    
    # Helper function to append to the temporary CAS file
    append_to_temp_cas <- function(new_data) {
        temp_cas_path <- "../Data/temp_CAS.fst"
        
        # Ensure new_data has the correct columns
        if (!all(c("PREFERRED_NAME", "CASRN") %in% names(new_data))) {
            stop("New data must contain PREFERRED_NAME and CASRN columns.")
        }
        
        temp_cas_dt <- if (file.exists(temp_cas_path)) {
            read.fst(temp_cas_path, as.data.table = TRUE)
        } else {
            data.table(PREFERRED_NAME = character(), CASRN = character())
        }
        
        # Combine, ensure uniqueness, and write back
        updated_dt <- rbindlist(list(temp_cas_dt, new_data), use.names = TRUE, fill = TRUE)
        updated_dt <- unique(updated_dt, by = c("PREFERRED_NAME", "CASRN"))
        
        write.fst(updated_dt, temp_cas_path)
    }

    fuzzy_match_non_interactive <- function(source_names, target_dt, match_col, threshold = 0.2) {
        results <- list()
        for (i in 1:length(source_names)) {
            name <- source_names[i]
            if (is.na(name)) next

            valid_targets <- target_dt[[match_col]][!is.na(target_dt[[match_col]])]
            if (length(valid_targets) == 0) next

            distances <- stringdist(name, valid_targets, method = "jw")
            if (length(distances) > 0 && !all(is.na(distances))) {
                min_dist <- min(distances, na.rm = TRUE)
                if (!is.na(min_dist) && min_dist <= threshold) {
                    best_idx <- which.min(distances)
                    best_match <- valid_targets[best_idx]
                    casrn_number <- target_dt[get(match_col) == best_match, CASRN][1]
                    
                    results[[length(results) + 1]] <- data.table(
                        source_name = name,
                        matched_name = best_match,
                        distance = min_dist,
                        CASRN = casrn_number
                    )
                }
            }
        }
        if (length(results) > 0) return(rbindlist(results))
        return(data.table())
    }

    # Plot helper: bar chart of total RQtg per sample, top-10 riskiest, top-3 names as subtitle
    make_sample_risk_plot <- function(tox_cal, group_label) {
        if (is.null(tox_cal) || nrow(tox_cal) == 0)
            return(ggplot() + labs(title = paste(group_label, ": No data available")))
        n_show <- min(10, nrow(tox_cal))
        top_samples <- tox_cal %>% arrange(desc(RQtg)) %>% slice(1:n_show)
        top3_label  <- paste("Top 3:", paste(head(top_samples$Sample, min(3, nrow(top_samples))), collapse = ", "))
        ggplot(top_samples, aes(x = reorder(Sample, RQtg), y = RQtg)) +
            geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
            labs(title    = paste(group_label, "\u2014 Top", n_show, "Riskiest Samples"),
                 subtitle = top3_label,
                 x = "Sample", y = "Risk Assessment (RQtg)")
    }

    # Plot helper: top-10 pollutants by median TU, boxplot across all samples
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

    observeEvent(input$start_processing, {
        req(input$user_file)
        
        v$summary_log <- c("Processing started...")
        
        tryCatch({
            v$user_data <- load_user_file(input$user_file$datapath)
            names(v$user_data)[1] <- "PREFERRED_NAME"  # normalise: first column is always the compound name
            v$p_vector <- v$user_data[[1]]
            
            v$summary_log <- c(v$summary_log, paste("Loaded", length(v$p_vector), "compounds from user file."))
            
            # Step 2 - CASRN search
            internal_list <- Known_CAS[Known_CAS$PREFERRED_NAME %in% v$p_vector, ]
            p_vector_found <- internal_list$PREFERRED_NAME
            unfound <- v$p_vector[!v$p_vector %in% p_vector_found]
            
            v$final_search_results <- internal_list[, .(PREFERRED_NAME, CASRN)]
            
            v$summary_log <- c(v$summary_log, paste("Found", length(p_vector_found), "exact matches in Known_CAS database."))
            
            if (length(unfound) > 0) {
                v$summary_log <- c(v$summary_log, paste("Searching for", length(unfound), "compounds in DSSTox (fuzzy match)..."))
                fuzzy_matches <- fuzzy_match_non_interactive(
                    unfound,
                    DSSTox,
                    "PREFERRED_NAME",
                    threshold = input$fuzzy_threshold
                )
                
                if (nrow(fuzzy_matches) > 0) {
                    v$fuzzy_to_review <- fuzzy_matches
                    v$summary_log <- c(v$summary_log, paste("Found", nrow(v$fuzzy_to_review), "potential fuzzy matches for review."))
                } else {
                    v$summary_log <- c(v$summary_log, "No fuzzy matches found.")
                    v$fuzzy_to_review <- NULL
                }
                
                names_in_fuzzy_review <- v$fuzzy_to_review$source_name
                v$manual_to_fill <- data.table(PREFERRED_NAME = unfound[!unfound %in% names_in_fuzzy_review])
                
            } else {
                v$summary_log <- c(v$summary_log, "All compounds found in Known_CAS database. No fuzzy search needed.")
                v$fuzzy_to_review <- NULL
                v$manual_to_fill <- data.table(PREFERRED_NAME=character())
            }

            updateTabsetPanel(session, "main_tabs", selected = "Interactive CASRN Matching")
            
        }, error = function(e) {
            v$summary_log <- c(v$summary_log, "An error occurred:", e$message)
        })
    })

    output$cas_search_summary <- renderPrint({
        cat(v$summary_log, sep = "\n")
    })

    output$fuzzy_match_table <- renderDT({
        req(v$fuzzy_to_review)
        df <- as.data.frame(v$fuzzy_to_review)
        df <- cbind(
            Accept = sapply(seq_len(nrow(df)), function(i) {
                as.character(checkboxInput(paste0("fuzzy_cb_", i), label = NULL, value = TRUE))
            }),
            df
        )
        datatable(df, escape = FALSE, selection = 'none',
                  options = list(
                      pageLength = 10,
                      preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
                      drawCallback   = JS('function() { Shiny.bindAll(this.api().table().node()); }')
                  ))
    }, server = FALSE)

    observeEvent(input$submit_fuzzy_matches, {
        req(v$fuzzy_to_review)

        n <- nrow(v$fuzzy_to_review)
        accepted_idx <- which(sapply(seq_len(n), function(i) isTRUE(input[[paste0("fuzzy_cb_", i)]])))

        if (length(accepted_idx) > 0) {
            confirmed_matches <- v$fuzzy_to_review[accepted_idx, ]
            newly_confirmed <- confirmed_matches[, .(PREFERRED_NAME = source_name, CASRN)]
            v$final_search_results <- rbindlist(list(v$final_search_results, newly_confirmed), use.names = TRUE, fill = TRUE)
            append_to_temp_cas(newly_confirmed)
            v$summary_log <- c(v$summary_log, paste("User confirmed", nrow(newly_confirmed), "fuzzy matches."))

            rejected_names <- v$fuzzy_to_review[-accepted_idx, .(PREFERRED_NAME = source_name)]
            v$manual_to_fill <- rbind(v$manual_to_fill, rejected_names)
        } else {
            v$summary_log <- c(v$summary_log, "No fuzzy matches were confirmed.")
            v$manual_to_fill <- rbind(v$manual_to_fill, v$fuzzy_to_review[, .(PREFERRED_NAME = source_name)])
        }

        v$fuzzy_to_review <- NULL # Clear the review table
    })

    observeEvent(input$add_manual_casrn, {
        req(input$manual_name, input$manual_casrn)
        name <- trimws(input$manual_name)
        casrn <- trimws(input$manual_casrn)
        
        if (name != "" && casrn != "") {
            new_entry <- data.table(PREFERRED_NAME = name, CASRN = casrn)
            # Add directly to the final results
            v$final_search_results <- rbind(v$final_search_results, new_entry)
            # Track for display
            v$manual_additions <- rbind(v$manual_additions, new_entry)
            # Store in temp_CAS for admin review
            append_to_temp_cas(new_entry)
            # Remove from the list of those to fill
            v$manual_to_fill <- v$manual_to_fill[PREFERRED_NAME != name]

            v$summary_log <- c(v$summary_log, paste("Manually added CASRN for", name))

            # Clear input fields for next entry
            updateTextInput(session, "manual_casrn", value = "")
        }
    })

    output$manual_added_list <- renderDT({
        v$manual_additions
    })
    
    observeEvent(input$start_toxicity_calc, {
        req(nrow(v$final_search_results) > 0)
        
        v$summary_log <- c(v$summary_log, "Starting toxicity calculations...")

        tryCatch({
            
            p_final <- as.data.table(v$final_search_results) %>%
              na.omit() %>%
              distinct(PREFERRED_NAME, .keep_all = TRUE) %>% # Ensure unique compounds
              mutate(cas_number = gsub("-", "", CASRN)) %>%
              select(PREFERRED_NAME, cas_number)
            
            cas_search_list <- as.vector(p_final$cas_number)

            endpoint_data <- final_ecotox_data %>%
              mutate(cas_number = as.character(cas_number)) %>%
              filter(cas_number %in% cas_search_list)
            
            # Algae
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
              rowwise() %>%
              mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
              ungroup()
            
            # Crustacean
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
              rowwise() %>%
              mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
              ungroup()

            # Fish
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
              rowwise() %>%
              mutate(RQtg = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
              ungroup()
            
            # Rename RQtg -> "Risk assessment" and move to second column for output
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
                algae_samples       = make_sample_risk_plot(tox_cal_algae,       "Algae"),
                crustacean_samples  = make_sample_risk_plot(tox_cal_crustacean,  "Crustacean"),
                fish_samples        = make_sample_risk_plot(tox_cal_fish,        "Fish"),
                algae_pollutants    = make_pollutant_plot(tox_cal_algae,         "Algae"),
                crustacean_pollutants = make_pollutant_plot(tox_cal_crustacean,  "Crustacean"),
                fish_pollutants     = make_pollutant_plot(tox_cal_fish,          "Fish")
            )
            
            v$summary_log <- c(v$summary_log, "Toxicity calculations complete. See tabs for results.")
            updateTabsetPanel(session, "main_tabs", selected = "Toxicity Plots")

        }, error = function(e) {
            v$summary_log <- c(v$summary_log, "An error occurred during toxicity calculation:", e$message)
        })
    })
    
    output$tox_table_algae <- renderDT({ req(v$tox_results); v$tox_results[["Algae toxicity"]] })
    output$tox_table_crustacean <- renderDT({ req(v$tox_results); v$tox_results[["Crustacean toxicity"]] })
    output$tox_table_fish <- renderDT({ req(v$tox_results); v$tox_results[["Fish toxicity"]] })
    
    output$algae_sample_plot       <- renderPlot({ req(v$plots); v$plots$algae_samples })
    output$crustacean_sample_plot  <- renderPlot({ req(v$plots); v$plots$crustacean_samples })
    output$fish_sample_plot        <- renderPlot({ req(v$plots); v$plots$fish_samples })
    output$algae_pollutant_plot    <- renderPlot({ req(v$plots); v$plots$algae_pollutants })
    output$crustacean_pollutant_plot <- renderPlot({ req(v$plots); v$plots$crustacean_pollutants })
    output$fish_pollutant_plot     <- renderPlot({ req(v$plots); v$plots$fish_pollutants })

    output$download_results <- downloadHandler(
        filename = function() {
            paste0("taxotox_", tools::file_path_sans_ext(input$user_file$name), ".xlsx")
        },
        content = function(file) {
            req(v$tox_results)
            write.xlsx(v$tox_results, file)
        }
    )
    
    # ── Admin panel ───────────────────────────────────────────────────────────
    # Change ADMIN_PASSWORD before deployment!
    ADMIN_PASSWORD <- "TaxoTox2025!"

    observeEvent(input$admin_login, {
        if (isTRUE(input$admin_password == ADMIN_PASSWORD)) {
            v$admin_logged_in <- TRUE
        } else {
            v$admin_logged_in <- FALSE
            showNotification("Incorrect password.", type = "error")
        }
    })

    output$admin_panel_ui <- renderUI({
        if (!isTRUE(v$admin_logged_in)) return(NULL)
        temp_path <- "../Data/temp_CAS.fst"
        if (!file.exists(temp_path))
            return(tagList(br(), p("No pending entries in temp_CAS — nothing to promote.")))
        tagList(
            h4("Pending entries in temp_CAS (awaiting promotion to Known_CAS)"),
            p("Select rows to approve, then click the button. Unselected rows stay in temp_CAS."),
            DTOutput("admin_temp_table"),
            br(),
            actionButton("admin_approve", "Promote selected to Known_CAS", class = "btn-success")
        )
    })

    output$admin_temp_table <- renderDT({
        req(v$admin_logged_in)
        temp_path <- "../Data/temp_CAS.fst"
        if (!file.exists(temp_path)) return(data.table())
        datatable(read.fst(temp_path, as.data.table = TRUE),
                  selection = "multiple", options = list(pageLength = 20))
    })

    observeEvent(input$admin_approve, {
        req(v$admin_logged_in)
        temp_path <- "../Data/temp_CAS.fst"
        if (!file.exists(temp_path)) { showNotification("temp_CAS.fst not found.", type = "warning"); return() }
        temp_dt <- read.fst(temp_path, as.data.table = TRUE)
        sel     <- input$admin_temp_table_rows_selected
        if (is.null(sel) || length(sel) == 0) { showNotification("No rows selected.", type = "warning"); return() }
        to_promote  <- temp_dt[sel, ]
        known_dt    <- read.fst("../Data/Known_CAS.fst", as.data.table = TRUE)
        updated_known <- unique(rbindlist(list(known_dt, to_promote[, .(PREFERRED_NAME, CASRN)]),
                                          use.names = TRUE, fill = TRUE),
                                by = c("PREFERRED_NAME", "CASRN"))
        write.fst(updated_known, "../Data/Known_CAS.fst")
        remaining <- temp_dt[-sel, ]
        if (nrow(remaining) > 0) write.fst(remaining, temp_path) else file.remove(temp_path)
        showNotification(paste(nrow(to_promote), "entries promoted to Known_CAS.fst."), type = "message")
    })

    output$interactive_casrn_ui <- renderUI({
      if (is.null(v$fuzzy_to_review) && !is.null(v$manual_to_fill) && nrow(v$manual_to_fill) > 0) {
        tagList(
          hr(),
          h4("Manual Entry"),
          p("Add CASRNs for remaining compounds."),
          selectInput("manual_name", "Select Compound", choices = v$manual_to_fill$PREFERRED_NAME),
          textInput("manual_casrn", "Enter CASRN (e.g., 123-45-6)"),
          actionButton("add_manual_casrn", "Add CASRN"),
          h5("Manually Added:"),
          DTOutput("manual_added_list")
        )
      }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
